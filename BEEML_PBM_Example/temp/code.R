## different ways of encoding DNA sequences. Assume only ACGT are present, no ambigouity letters such as "N"

## matrix used to transform 01 encoded sequences to wyk basis (from 4 to 3 dimensions)
mk.wyk.tr.mtx <- function(){
  tr.mtx = matrix(c(1,-1,-1,-1,1,-1,-1,-1,1,1,1,1),nrow=4,byrow=T)
  rownames(tr.mtx) = c("A","C","G","T")
  colnames(tr.mtx) = c("W","Y","K")
  tr.mtx
}

## matrix used to transform 01 encoded sequences to wyk basis
mk.wyk.di.tr.mtx <- function(){
  tr.mtx = mk.wyk.tr.mtx()
  tr.mtx.di = matrix(0,nrow=16,ncol=15)
  for(x in 1:4) for(y in 1:4) tr.mtx.di[(x-1)*4+y,] = c(tr.mtx[x,],tr.mtx[y,], as.vector(tr.mtx[y,] %o% tr.mtx[x,]))
  colnames(tr.mtx.di) = c("W1","Y1","K1","W2","Y2","K2","WW","WY","WK","YW","YY","YK","KW","KY","KK")
  rownames(tr.mtx.di) = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
  tr.mtx.di
}

mk.01.tr.mtx <- function(){
  tr.mtx=matrix(0,nrow=4,ncol=4)
  diag(tr.mtx) = 1
  rownames(tr.mtx) = c("A","C","G","T")
  colnames(tr.mtx) = c("A","C","G","T")
  tr.mtx
}

mk.01.di.tr.mtx <- function(){
  tr.mtx=matrix(0,nrow=16,ncol=16)
  diag(tr.mtx) = 1
  rownames(tr.mtx) = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
  colnames(tr.mtx) = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
  tr.mtx
}

calc.p <- function(s.e,ns.e,mu) (exp(-s.e)+exp(-ns.e))/(exp(-s.e)+exp(-ns.e)+exp(-mu))

encode01 <- function(kmer,consensus = NULL){
    kmer <- toupper(kmer)
    n <- nchar(kmer)
    out <- rep(0,4*n)
    alph <- c("A","C","G","T")
    out[seq(0,n*4-1,4) + match(strsplit(kmer,"")[[1]], alph)] = 1
    # drop the columns corresponding to consensus base so that columns are linearly independent
    if (!is.null(consensus)) {
      if (nchar(consensus) != n) {
        cat("consensus sequence is not the same length as sequence")
      } else {
       out = out[!encode01(consensus)]
      }
    }
    as.integer(out)
}

decode01 <- function(seq.01,filter=NULL,alph=c("A","C","G","T")) {
    if (is.null(filter)) filter = rep(1:4,length(seq.01)/4)
    paste(alph[seq.01*filter],collapse="")
}

# Example: AAGC first identify column corresponsing to AA, then columns that end with G in the third position and C in the 4th position
find.consensus.col.idx.di <- function(consensus) {
  alph <- c("A","C","G","T")
  # the column corresponding to the first di-nucleotide of the consensus sequence
  d = (which(alph==substr(consensus,1,1))-1) * 4 + which(alph==substr(consensus,2,2))
  # all the columns containing rest of the sequences
  i = 3
  while (i <= nchar(consensus)){
    d = c(d, seq((i-2)*16 + which(alph==substr(consensus,i,i)), (i-2)*16+16, 4))
    i = i + 1
  }
  d
}

encode.01.di <- function(kmer,consensus=NULL) {
  out = c()
  len = nchar(kmer)
  if (len >= 2) {
    tr.mtx = mk.01.di.tr.mtx()
    out = tr.mtx[rownames(tr.mtx)==substr(kmer,1,2),]
    i = 2
    while (i < len) {
      out = c(out, tr.mtx[rownames(tr.mtx)==substr(kmer,i,i+1),])
      i = i + 1
    }

    # drop the columns corresponding to consensus base so that columns are linearly independent
    if (!is.null(consensus)) {
      if (nchar(consensus) != len) {
        cat("consensus sequence is not the same length as sequence")
      } else {
        out = out[-find.consensus.col.idx.di(consensus)]
      }
    }
  }
  as.integer(out)
}

make.kmers.01 <- function(k) {
    out = vector("integer",(4^k)*4*k)
    for (i in 1:k) {
        for (j in 0:3){
	    idx = (i-1)*4+(j+1)
            idx.start = (idx - 1) * 4^k + 1
            idx.end = idx.start + 4^k - 1
            x = vector("integer",4^i)
	    start = j*(4^(i-1))+1
	    end = start + 4^(i-1) - 1
	    x[start:end] = as.integer(1)
	    out[idx.start:idx.end] = rep(x, 4^(k-i))
        }
    }
    matrix(out,nrow=4^k)
}

encode.wyk <- function(kmer){
  as.vector(t(mk.wyk.tr.mtx()) %*% matrix(encode01(kmer),4)) 
}

decode.wyk <- function(seq.wyk, alph = c("A","C","G","T")) {
  idx = apply(matrix(seq.wyk,3),2,function(x) ifelse(sum(x) == 3,4, (1:3)[x==1]))
  paste(alph[idx],collapse="")
}

# di-nucleotide encoding 
encode.wyk.di <- function(kmer) {
  out = c()
  len = nchar(kmer)
  if (len >= 2) {
    tr.mtx = mk.wyk.di.tr.mtx()
    out = tr.mtx[rownames(tr.mtx)==substr(kmer,1,2),]
    i = 3
    while (i <= len) {
      out = c(out, tr.mtx[rownames(tr.mtx)==substr(kmer,i-1,i),][4:15])
      i = i + 1
    }
  }
  as.integer(out)
}

# very slow, but too lazy to figure out the repeating patterns like in make.kmers.01
make.kmers.wyk <- function(k) {
  t(apply(make.kmers.01(k),1,function(x) as.integer(t(mk.wyk.tr.mtx()) %*% matrix(x,4))))
}

make.kmers.wyk.di <- function(k) {
  ## first make wyk encoded pwm
  out = make.kmers.wyk(k)
  ## add all the interaction terms (outer product of pwm columns in wyk encoding)
  for( i in seq(1, by=3, length.out = k-1)) {
    out = cbind(out, t(apply(out[,i:(i+5)],1,function(x) as.integer(x[1:3] %o% x[4:6]))))
  }
  out
}

find.idx.01 <- function(seq.01,weights = NULL,filter=NULL) {
    if (is.null(filter)) filter = rep(1:4,length(seq.01)/4)
    bases = seq.01 * filter
    bases = bases[bases>0]-1
    if (is.null(weights)) weights = 4^(0:(length(bases)-1))
    as.integer(bases %*% weights + 1)
}

find.kmer.idx <- function(kmer){
    kmer.01 = encode01(kmer)
    find.idx.01(kmer.01)
}


count.unique <- function(s,counts=NULL) {
    s.idx = unique(s)
    s.idx.ch = as.character(s.idx)
    s.ch = as.character(s)
    s.idx.map.total.counts = new.env(hash=T,parent=emptyenv())
    s.idx.map.num.occ = new.env(hash=T,parent=emptyenv())
    for(i in s.idx.ch) assign(i,0,envir=s.idx.map.total.counts)
    for(i in s.idx.ch) assign(i,0,envir=s.idx.map.num.occ)
    if (is.null(counts)) counts = rep(1,length(s))
    for (i in 1:length(s.ch)) assign(s.ch[i], get(s.ch[i],envir=s.idx.map.total.counts)+counts[i],envir=s.idx.map.total.counts)
    for (i in 1:length(s.ch)) assign(s.ch[i], get(s.ch[i],envir=s.idx.map.num.occ)+1,envir=s.idx.map.num.occ)
    s.counts = sapply(s.idx.ch,function(x) get(x,envir=s.idx.map.total.counts))
    s.occ = sapply(s.idx.ch,function(x) get(x,envir=s.idx.map.num.occ))
    list(idx = s.idx, counts = s.counts, num.occ = s.occ)
}

#log odds matrix with pseudocounts (defaults to 1)
make.matrix.lo.01 <- function(seqs.01, s, pseudo=1){
  counts = apply(seqs.01[s,],2,sum)
  lo = apply(matrix(counts,4),2,function(x) ifelse(x==max(x),0,log( (max(x)+pseudo) / (x+pseudo) )))
  as.vector(lo)
}

# calculation partition function for full binding model, including mu and non-specific binding
# estimate energy distribution by FFT (convolve() function), then calculate binding probability
# of energy levels to save time. Energy levels are rounded to 2 decimal places by default
est.partition.fn <- function(mtx, ns.e, mu, bg.f = NULL,precision=2) {
   m.rounded = matrix(as.integer(round(mtx*10^precision)),4)
   m = apply(m.rounded,2,function(x) x-min(x))
   if(is.null(bg.f)) bg.f = matrix(0.25,nrow=nrow(m),ncol=ncol(m))
   col.idx = 1
   eng.dist = rep(0,max(m[,col.idx])+1)
   for (base.idx in 1:4) eng.dist[m[base.idx,col.idx]+1] = eng.dist[m[base.idx,col.idx]+1]+bg.f[base.idx,col.idx]
   col.idx = 2
   y = rep(0,max(m[,col.idx]+1))
   for (base.idx in 1:4) y[m[base.idx,col.idx]+1] = y[m[base.idx,col.idx]+1]+bg.f[base.idx,col.idx]
   eng.dist = convolve(eng.dist,rev(y),type="o")
     
   for (col.idx in 3:(ncol(m))){
        y = rep(0,max(m[,col.idx]+1))
        for (base.idx in 1:4) y[m[base.idx,col.idx]+1] = y[m[base.idx,col.idx]+1]+bg.f[base.idx,col.idx]
        eng.dist = convolve(eng.dist,rev(y),type="o")
 
   }
   eng.dist[eng.dist<1e-15] = 0
   rounded.max = sum(apply(m.rounded,2,max))
   rounded.min = sum(apply(m.rounded,2,min))
   eng = (rounded.min:rounded.max) / (10^precision)
   bound.weights = exp(-eng) + exp(-ns.e)
   unbound.weights = exp(-mu)
   sum((bound.weights)/(bound.weights+unbound.weights) * eng.dist) * (4^ncol(m))
}

sse.fn.lm.est <- function(free.para.val, free.para.pos, mtx, seqs, sample.size,target.counts,precision=2,prior=NULL) {
#    free.para.val = round(free.para.val,precision)
    n.par = length(free.para.val)
  #  mtx[free.para.pos] = free.para.val[1:(n.par-1)]
    mtx[free.para.pos] = free.para.val[1:(n.par-2)]
    ns.e = free.para.val[n.par-1]
    mu = free.para.val[n.par]
    obs.p = calc.p(seqs %*% mtx,ns.e,mu)
    # only able to handle prior in terms of biased base probabilities
    # not implmented right now
    z = est.partition.fn(mtx=mtx,ns.e=ns.e,mu=mu,precision=precision)
    model.counts = obs.p/z*sample.size
    (target.counts-model.counts)
}

sse.fn.lm <- function(free.para.val, free.para.pos, mtx, seqs, s.idx, sample.size,target.counts,prior=NULL,s.idx.prior=NULL,s.idx.weights=NULL) {
  #    free.para.val = round(free.para.val,precision)
      n.par = length(free.para.val)
    #  mtx[free.para.pos] = free.para.val[1:(n.par-1)]
      mtx[free.para.pos] = free.para.val[1:(n.par-2)]
      ns.e = free.para.val[n.par-1]
      mu = free.para.val[n.par]
      bind.p = calc.p(seqs %*% mtx,ns.e,mu)
      if (!is.null(prior)) bind.p = bind.p * prior
      z = sum(bind.p)
      model.counts = bind.p[s.idx] / z * sample.size
      if (!is.null(s.idx.prior)) model.counts = model.counts * (s.idx.prior/prior[s.idx])
      if (!is.null(s.idx.weights)) {
          model.counts = model.counts*s.idx.weights
          target.counts = target.counts * s.idx.weights
      }
      (target.counts-model.counts)
}

library(seqLogo)
print.logo <- function(mtx) seqLogo(apply(matrix(mtx,4),2,function(x) exp(-x)/sum(exp(-x))))
