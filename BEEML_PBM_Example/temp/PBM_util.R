require(minpack.lm)
require(MASS)
require(seqLogo)



# extracts the top matrix from the final round from consensus output file
get.consensus.top.pwm <- function(consensus.output.file) {
  # read all lines of the file into an array
  out = readLines(consensus.output.file, n=-1)

  # find the line that starts the list of matrices from the final cycle
  start = grep("THE LIST OF MATRICES FROM FINAL CYCLE",out)[1]
  # first matrix starts 8 lines below ...
  mtx.start = start + 8

  if (mtx.start > length(out)) {
    mtx.consensus = NA
  } else {
    mtx.consensus = c()
    for (j in 0:3) {
      line.no = mtx.start + j
      mtx.consensus[j+1] = substr(out[line.no],gregexpr(" +",out[line.no])[[1]][2] + attr(gregexpr(" +",out[line.no])[[1]],"match.length")[2], nchar(out[line.no]))
    }
    mtx.consensus = sapply(mtx.consensus,function(x) strsplit(x,split=" +"))
    mtx.consensus = lapply(mtx.consensus, as.integer)
    # get rows into ACGT order (consensus has it as ATCG)
    mtx.consensus = mtx.consensus[c(1,3,4,2)]
    # mtx.consensus is now a count matrix
    mtx.consensus = matrix(unlist(mtx.consensus),nrow=4,byrow=T)
    # add pseudo count, turn into frequency matrix
    mtx.consensus = apply(mtx.consensus + 1, 2 ,function(x) x/sum(x))
    # turn frenquency matrix into energy matrix (ignore background)
    mtx.consensus = as.vector(apply(-log(mtx.consensus),2,function(x) x-min(x)))
  }
  mtx.consensus
}

pretty.print.pwm <- function(pwm, comment, file, sign_digit = 2) {
  pretty.pwm = round(matrix(pwm,4),sign_digit)
  rownames(pretty.pwm) = c("A","C","G","T")
  write(paste("#", comment) , file=file)
  write.table(pretty.pwm,file=file,append=T,col.names=F,sep="\t",quote=F)
}


get.high.information.positions <- function(mtx, k, energy.mtx = T) {
  if (energy.mtx) {
    freq.mtx = apply(matrix(exp(-mtx),4),2,function(x) x/sum(x))
  } else {
    freq.mtx = mtx
  }
  entropy = apply(freq.mtx,2,function(x) sum(x*-log(x)))
  start = which.min(sapply(1:(length(entropy)-k+1), function(x) sum(entropy[x:(x+k-1)])))
  as.vector(apply(-log(freq.mtx[,start:(start+k-1)]),2,function(x) x-min(x)))
}

get.kmer.intensity.ratios <- function(kmer,seqs,intensities, ratio=T, both.orientations =T) {
  pos = c()
  int = c()

  # palidromic kmers only occur on 16 sequences
  matches.idx = grep(kmer,seqs)
  for (idx in matches.idx) {
    location = as.integer(regexpr(kmer,seqs[idx]))
    pos = c(pos, location)
    int = c(int, intensities[idx])
  }

  if (both.orientations) {
    kmer.revcomp = decode01(rev(encode01(kmer)))
    if ( kmer != kmer.revcomp) {
      matches.idx = grep(kmer.revcomp, seqs)
      for (idx in matches.idx) {
        location = as.integer(regexpr(kmer.revcomp,seqs[idx]))
        pos = c(pos, location)
        int = c(int, intensities[idx])
      }
    }
  }
  if (ratio) int = int/mean(int, na.rm=T)
  cbind(pos,int)
}

est.orientation.weight <- function(pwm, seqs.01, seqs, values, top=25, smooth=F, kmers = NULL) {
   seqs = seqs[!is.na(values)]
   values = values[!is.na(values)]

   if (is.null(kmers)){
     energy = seqs.01 %*% pwm
     ## get the top i = 25 kmers according to the matrix,estimate position effect based on them
     kmers = apply(seqs.01[order(energy)[1:top],],1,decode01)
   } else {
     top = length(kmers)
   }


   l = unlist(sapply(1:25, function(q) {
               kmer = kmers[q]
               kmer.rc = decode01(rev(encode01(kmer)))
               idx.kmer = grep(kmer, seqs)
               idx.kmer.rc = grep(kmer.rc, seqs)
               pos.kmer = unlist(gregexpr(kmer, seqs[idx.kmer]))
               pos.kmer.rc = unlist(gregexpr(kmer.rc, seqs[idx.kmer.rc]))
               common = unique(intersect(pos.kmer,pos.kmer.rc))
               sapply(unique(intersect(pos.kmer,pos.kmer.rc)),function(z) mean(values[idx.kmer[pos.kmer==z]]) / mean(values[idx.kmer.rc[pos.kmer.rc==z]]) )
             }))
   median(l,na.rm=T)
 }

est.position.weights <- function(pwm, seqs.01, seqs, values, kmers = NULL, top=25, smooth=F,both.orientations =T) {
  seqs = seqs[!is.na(values)]
  values = values[!is.na(values)]

  if (is.null(kmers)) {
    energy = seqs.01 %*% pwm
    # get the top kmers according to the matrix,estimate position effect based on them
    kmers = apply(seqs.01[order(energy)[1:top],],1,decode01)
  } else {
    # use given kmers to estimate position effect
    top = length(kmers)
  }
  # weights are averaged over top 100 8mers, effect of position on intensity level
  d = get.kmer.intensity.ratios(kmers[1], seqs, values, both.orientations = both.orientations)
  for(j in 2:top) d = rbind(d,get.kmer.intensity.ratios(kmers[j], seqs, values))
    pos = sort(unique(d[,1]))
    for (j in 1:length(pos))
      if (j != pos[j]) {
        j = j - 1
        break
    }
  weights = sapply(1:j, function(idx) mean(d[d[,1]==idx,2],na.rm=T))

  # do some reasonably agreesive smoothing on the weights make f value smaller to do less smoothing
  if (smooth) weights = lowess(1:length(weights),weights,f=1/5)$y
  # since we are extending into the constant region, assume the last weight holds for the rest of them
  weights/sum(weights)
}

est.bg.weights <- function(values, num.bins=200, lowess.f=0.1) {

  ## estimate lower half of background distribution (upper ground is convoluted with binding signal)
  ## make a histogram of
  norm.values = (values - median(values))/mad(values)
##  bg.values = values[norm.values < 1 & values > 0]
  bg.values = values[norm.values < 2 ]

  h = hist(bg.values,nclass=num.bins,plot=F)

  ## lowess smoothing to get better density estimate
  total = lowess(h$mids,h$counts,f=lowess.f)
  x.lower = total$x[1:which.max(total$y)]
  y.lower = total$y[1:which.max(total$y)]

  # assume background distribution is symmetrical, so upper tail is just the mirror image of the lower tail
  x.upper = rev(sapply(x.lower, function(x) x + 2 * (max(x.lower - x))))
  y.upper = rev(y.lower)
  x.upper = x.upper[2:length(x.upper)]
  y.upper = y.upper[2:length(y.upper)]
  bg = list(x = c(x.lower, x.upper), y = c(y.lower, y.upper))
  if (length(total$x) > length(bg$x)) {
    total = list(x = total$x[1:length(bg$x)], y = total$y[1:length(bg$y)])
  } else {
    bg = list(x=bg$x[1:length(total$x)], y = bg$y[1:length(total$y)])
  }
  bin.width = (bg$x[2]-bg$x[1])
  bin.weight = rep(0,length(bg$x))
  start = which.max(bg$y)
  total = sapply(start:length(bg$x), function(x) sum(values > bg$x[x] - bin.width/2 & values < bg$x[x] + bin.width/2))
  bg.exp = bg$y[start:length(bg$x)]
  bin.weight[start:length(bg$x)] = (total-bg.exp)/total
  bin.weight[bin.weight<0] = 0
#  bin.weight = (total$y-bg$y)/total$y
  bin.weight = round(bin.weight,2)

  probe.weights = rep(0,length(values))
  bin.begin = bg$x[bin.weight>0] - bin.width/2
  bin.end = bg$x[bin.weight>0] + bin.width/2
  for (i in 1:length(bin.begin)) probe.weights[values > bin.begin[i] & values < bin.end[i]] = bin.weight[bin.weight>0][i]
  probe.weights[values > bin.end[length(bin.end)]] = 1
  probe.weights
}

predict.occupancy <- function(mtx,seqs.01,idx.f, idx.r, pos.weights, orientation.weight = 1, boltzmann=F) {
  if (boltzmann) {
    p = exp(-seqs.01 %*% mtx)
  } else {
    p = 1/(1+exp(seqs.01%*%mtx))
  }
  pf = matrix(p[idx.f],ncol=ncol(idx.f))
  pnf = 1-pf
  pr = matrix(p[idx.r],ncol=ncol(idx.r))
  (pf + pnf * pr*orientation.weight) %*% pos.weights
}

beeml.ob.f.sum <- function(para, free.pos, target, idx.f, idx.r, position.weights, seqs.01, probe.weights, orientation.weight = 1, boltzmann=F, lambda = 0.1){
  pred = predict.occupancy(para, seqs.01, idx.f, idx.r, position.weights, orientation.weight, boltzmann)
  # don't penalize mu
  if (boltzmann) {
    pen = lambda * para
  } else {
    pen = lambda * para[1:(length(para)-1)]
  }
  c(lm(target~pred)$residuals * probe.weights, pen)
}

beeml.ob.f.pali.sum <- function(par, halfsite.len, spacer.len, free.pos, target, idx.f, idx.r, position.weights, seqs.01, probe.weights, orientation.weight = 1, boltzmann=F, lambda=0.1){

  mtx = make.full.from.halfsite(par[1:(3*halfsite.len)],halfsite.len, spacer.len,free.pos)

  if (spacer.len > 0) {
    free.pos = c(free.pos[1:(4*halfsite.len)], free.pos[(4*halfsite.len+1):(4*(halfsite.len+spacer.len))], rev(free.pos[1:(4*halfsite.len)]))
  } else {
    free.pos = c(free.pos, rev(free.pos))
  }
  if (boltzmann) {
    pred = predict.occupancy(mtx[free.pos], seqs.01, idx.f, idx.r, position.weights, orientation.weight, boltzmann)
  } else {
    mu = par[length(par)]
    pred = predict.occupancy(c(mtx[free.pos], mu), cbind(seqs.01,-1), idx.f, idx.r, position.weights, orientation.weight, boltzmann)
  }

  c(lm(target~pred)$residuals * probe.weights, lambda * mtx[free.pos])
}


get.beeml.solution.sum <- function(values, idx.f, idx.r, position.weights, seed.mtx, seqs.01 = seqs.n.01, palidromic = F, nprint = 0, ignore.background=F, halfsite.len=5, spacer.len=0, orientation.weight = 1, boltzmann=F, lambda = 0, maxiter = 15, initial.mu = 2) {
   # get rid of na values.
   good.idx = !is.na(values)
  good.values = values[good.idx]
   good.idx.f = idx.f[good.idx,]
   good.idx.r = idx.r[good.idx,]

   if (ignore.background) {
     bg.weights = rep(1,length(good.values))
   } else {
     # have to take the square root, since nls.lm takes residuals rather than square of residuals, so weights have to be adjusted accordingly
     bg.weights = sqrt(est.bg.weights(good.values))
   }
   # make sure consensus of seed matrix is set to 0
   seed.mtx = as.vector(apply(matrix(seed.mtx,4),2,function(x) x-min(x)))
   f.pos = seed.mtx>0

   if (palidromic) {
     mtx = make.full.from.halfsite(seed.mtx[f.pos],halfsite.len, spacer.len,f.pos)
     mtx = as.vector(apply(matrix(mtx, 4),2,function(x) x-min(x)))
     seqs.f  = seqs.01[,mtx>0]

     ## sol = nls.lm(runif(sum(f.pos)+1),  fn = beeml.ob.f.pali, target = values[good.idx]
     if (boltzmann) {
       sol = nls.lm(seed.mtx[f.pos], fn = beeml.ob.f.pali.sum, target = good.values, halfsite.len = halfsite.len, spacer.len = spacer.len, free.pos = f.pos, idx.f = good.idx.f, idx.r = good.idx.r, position.weights = position.weights, seqs.01 = seqs.f, probe.weights = bg.weights, boltzmann = boltzmann,lambda = lambda, control = nls.lm.control(nprint=nprint,maxiter=maxiter))
     } else {
       sol = nls.lm(c(seed.mtx[f.pos], initial.mu), fn = beeml.ob.f.pali.sum, target = good.values, halfsite.len = halfsite.len, spacer.len = spacer.len, free.pos = f.pos, idx.f = good.idx.f, idx.r = good.idx.r, position.weights = position.weights, seqs.01 = seqs.f, probe.weights = bg.weights, boltzmann = boltzmann,lambda = lambda, control = nls.lm.control(nprint=nprint,maxiter=maxiter))
     }
   }
   else {
     seqs.f = seqs.01[,f.pos]
     if (boltzmann) {
        sol = nls.lm(seed.mtx[f.pos], fn = beeml.ob.f.sum, target = good.values, free.pos = f.pos, idx.f = good.idx.f, idx.r = good.idx.r, position.weights = position.weights, seqs.01 = seqs.f, probe.weights = bg.weights, boltzmann = boltzmann, lambda = lambda, orientation.weight  = orientation.weight , control = nls.lm.control(nprint=nprint,maxiter=maxiter))
     } else {
       seqs.f = cbind(seqs.f, -1)
       sol = nls.lm(c(seed.mtx[f.pos], initial.mu), fn = beeml.ob.f.sum, target = good.values, free.pos = f.pos, idx.f = good.idx.f, idx.r = good.idx.r, position.weights = position.weights, seqs.01 = seqs.f, probe.weights = bg.weights, boltzmann = boltzmann, lambda = lambda, orientation.weight  = orientation.weight , control = nls.lm.control(nprint=nprint,maxiter=maxiter))
     }
   }

   mtx = seed.mtx

   if (boltzmann) {
     mtx[f.pos] = sol$par
     mu = NA
   } else {
     mtx[f.pos] = sol$par[1:(length(sol$par)-1)]
   }

   if (palidromic) mtx = make.full.from.halfsite(mtx[f.pos], halfsite.len, spacer.len, f.pos)


   if (!boltzmann) {
     mu = sol$par[length(sol$par)] - sum(apply(matrix(mtx,4),2,min))
   }
   mtx = as.vector(apply(matrix(mtx,4),2,function(x)x-min(x)))
   list(sol = sol, mtx = mtx, mu = mu)
}

mu.ob.f <- function(mu, mtx, target, seqs.01, pos.weights, f.weight = 1) {
  pred = predict.occupancy(mtx = mtx, mu = mu, seqs.01 = seqs.01, weights = pos.weights, f.weight = f.weight, boltzmann = F)
  summary(lm(target~pred))$r.squared
}

make.full.from.halfsite <- function(par, halfsite.len, spacer.len, free.pos) {
  halfsite.idx = 1:(4*halfsite.len)
  halfsite.values = par[1:(3*halfsite.len)]

  if(spacer.len>0) {
    spacer.idx = (4*halfsite.len+1):(4*halfsite.len+4*spacer.len)
    spacer.values = par[(3*halfsite.len+1):(3*halfsite.len+3*spacer.len)]
    flip.halfsite.idx = (max(spacer.idx)+1):(length(halfsite.idx)+max(spacer.idx))
  } else {
    flip.halfsite.idx = (max(halfsite.idx)+1):(length(halfsite.idx)+max(halfsite.idx))
  }

  mtx = rep(0,4*(2*halfsite.len+spacer.len))
  mtx[halfsite.idx[free.pos[halfsite.idx]]] = halfsite.values
  if (spacer.len > 0) {mtx[spacer.idx[free.pos[spacer.idx]]] = spacer.values}
  mtx[flip.halfsite.idx[rev(free.pos[halfsite.idx])]] = rev(halfsite.values)
  mtx
}
