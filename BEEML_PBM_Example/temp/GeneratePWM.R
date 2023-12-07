
#Usage: R CMD BATCH '--args deBruijn_file.txt SNW_PWM_file.txt output_prefix' GeneratePWM.R

## data files are from uniprobe database http://thebrain.bwh.harvard.edu/uniprobe/downloads.php
## first column is normalized probe intensity, second column is probe sequence

#data.file = "./Plagl1_0972.2_v1_deBruijn.txt"
#seed.pwm.file = "./Plagl1_pwm_primary.txt"


args <- commandArgs(trailingOnly = TRUE)

if (length(args)<3)
{
	stop("Missing arguments")
}

data.file = args[1]
seed.pwm.file = args[2]
output.prefix = args[3]

source("/data/bulyk/pipelines/universal_PBM_analysis/BEEML_PBM_Example/code.R")
source("/data/bulyk/pipelines/universal_PBM_analysis/BEEML_PBM_Example/PBM_util.R")

L = 10
all.seqs.01.l = make.kmers.01(L)
all.seqs.01.8 = make.kmers.01(8)

## list of 8mers that does not include reverse complement of itself
i=1
seen = new.env(hash=T)
for(i in 1:nrow(all.seqs.01.8)) {
  word.revcomp.seen = exists(decode01(rev(all.seqs.01.8[i,])), envir=seen)
  if ( !word.revcomp.seen ) assign(decode01(all.seqs.01.8[i,]),i,envir=seen)
}
x = ls(envir=seen)
good.8mer.idx = sapply(x, function(idx) get(idx,envir=seen))

print(good.8mer.idx[1:10])


rm(i)
rm(x)
rm(seen)
rm(word.revcomp.seen)



pbm.data = read.table(data.file, as.is=T)


## process sequence. this only needs to be done once for each array design, takes ~8 minutes on an imac (2006 model)
## probe sequences from Badis et. al. have 36 long variable region, I take the sequences of the first 40 bases into account
## which includes some of the costant region.
var.len = 40
num.lmers = var.len - L + 1

var.seqs = sapply(pbm.data[,2], function(seq) substr(seq, 1, var.len))
var.seqs.01 = t(sapply(var.seqs, encode01))

## for each 8mer, figure out which probe contains it
## pbm.8mer.probes.idx is a list, each element of which is the vector of probe indices that contain a particular 8mer, in either orientation
num.8mers = var.len - 8 + 1
pbm.8mer.probes.idx = list()
pbm.8mer.probes.idx[4^8+1]=NA
start.pos = seq(1, by=4, length.out=num.8mers)
for(probe.idx in 1:nrow(var.seqs.01)) {
  for (start in start.pos) {
    word = var.seqs.01[probe.idx, start:(start+(4*8-1))]
    word.idx = find.idx.01(word)
    word.revcomp.idx = find.idx.01(rev(word))
    pbm.8mer.probes.idx[[word.idx]] = c(pbm.8mer.probes.idx[[word.idx]], probe.idx)
    pbm.8mer.probes.idx[[word.revcomp.idx]] = c(pbm.8mer.probes.idx[[word.revcomp.idx]], probe.idx)
  }
}
pbm.8mer.probes.idx = lapply(pbm.8mer.probes.idx[1:(4^8)],unique)
pbm.probes.8mer.idx.f = t(apply(var.seqs.01, 1, function(x) sapply(start.pos, function(i) find.idx.01(x[i:(i+8*4-1)]))))
pbm.probes.8mer.idx.r = t(apply(var.seqs.01, 1, function(x) {y = rev(x);sapply(rev(start.pos), function(i) find.idx.01(y[i:(i+8*4-1)]))}))

## for each probe, figure out which L-mer are on it.
start.pos = seq(1, by=4, length.out=num.lmers)
pbm.probes.lmer.idx.f = t(apply(var.seqs.01, 1, function(x) sapply(start.pos, function(i) find.idx.01(x[i:(i+L*4-1)]))))
pbm.probes.lmer.idx.r = t(apply(var.seqs.01, 1, function(x) {y = rev(x);sapply(rev(start.pos), function(i) find.idx.01(y[i:(i+L*4-1)]))}))


## figure out the position weights, the same sequence seem to have more signal farther away from the glass
pbm.8mer.medians = unlist(lapply(pbm.8mer.probes.idx, function(x) median(pbm.data[x,1],na.rm=T)))
num.top.8mers = 50
top.8mers = apply(all.seqs.01.8[order(pbm.8mer.medians,decreasing=T)[1:num.top.8mers],],1,decode01)
pbm.position.weights = est.position.weights(0, 0, var.seqs, pbm.data[,1],kmers = top.8mers)
pbm.position.weights = c(pbm.position.weights, rep(pbm.position.weights[length(pbm.position.weights)], 40 - length(pbm.position.weights)))

## 8mer median intensities explain ~80% of variance of probe intensities?
pbm.8mers.pred = apply(cbind(pbm.probes.8mer.idx.f, pbm.probes.8mer.idx.r), 1, function(x) max(pbm.8mer.medians[x] %*% c(pbm.position.weights[1:ncol(pbm.probes.8mer.idx.f)], pbm.position.weights[1:ncol(pbm.probes.8mer.idx.r)])))
#pbm.8mers.pred.rsqr = round(cor(pbm.8mers.pred, pbm.data[,1])^2,2)


#pdf(paste(output.prefix, 'BEEML_analysis.pdf', sep="_"))

#plot(pbm.8mers.pred, pbm.data[,1], xlab="Predicted intensity", ylab="Observed intensity", main=paste('8-mer median R^2 =',pbm.8mers.pred.rsqr))

## take L contigous columns with highest information from uniprobe pwm to use as starting position for optimization
seed.mtx = get.high.information.positions(as.matrix(read.table(seed.pwm.file,skip=1)[,-1]), L, F)
pbm.sol = get.beeml.solution.sum(values = pbm.data[,1]/sd(pbm.data[,1],na.rm=T), idx.f = pbm.probes.lmer.idx.f, idx.r = pbm.probes.lmer.idx.r, position.weights = pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)], seed.mtx = seed.mtx, seqs.01 = all.seqs.01.l, palidromic = F, nprint = 1, lambda = 0.1)

#pbm.sol

## how well does the beeml regression model perform?
## fit at the probe level
pbm.beeml.pred = predict.occupancy(cbind(pbm.sol$mtx, pbm.sol$mu), all.seqs.01.l, pbm.probes.lmer.idx.f, pbm.probes.lmer.idx.r, pbm.position.weights[1:ncol(pbm.probes.lmer.idx.f)])
#pbm.beeml.pred.rsqr = round(cor(pbm.beeml.pred, pbm.data[,1])^2,2)
## fit at the 8mer median intensity level, excluding reverse complements
pbm.beeml.pred.8mers = unlist(lapply(pbm.8mer.probes.idx, function(x) median(pbm.beeml.pred[x],na.rm=T)))
#pbm.beeml.pred.8mers.rsqr = round(cor(pbm.beeml.pred.8mers[good.8mer.idx], pbm.8mer.medians[good.8mer.idx])^2,2)

#plot(pbm.beeml.pred, pbm.data[,1], xlab="Predicted probe intensity", ylab="Observed probe intensity", main=paste('BEEML-PBM R^2 =',pbm.beeml.pred.rsqr))
#plot(pbm.beeml.pred.8mers[good.8mer.idx], pbm.8mer.medians[good.8mer.idx], xlab="Predicted 8-mer intensity", ylab="Observed 8-mer intensity", main=paste('BEEML-PBM 8-mer R^2 =',pbm.beeml.pred.8mers.rsqr))


#dev.off()


# Create data frames

exp.name <- strsplit(output.prefix, '/', fixed = TRUE)[[1]]
exp.name <- exp.name[length(exp.name)] #Experiment name is the prefix basename

kmer.df <- data.frame(Experiment = exp.name, k.mer = names(good.8mer.idx), Observed.MI = pbm.8mer.medians[good.8mer.idx], Predicted.MI = pbm.beeml.pred.8mers[good.8mer.idx])

write.table(kmer.df, file = paste(output.prefix, 'kmer_df.txt', sep="_"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

probe.df <- data.frame(Experiment = exp.name, Probe.Seq = pbm.data[,2], Observed.MI  = pbm.data[,1], BEEML.Pred.MI = pbm.beeml.pred, Kmer.Pred.MI =  pbm.8mers.pred)

write.table(probe.df, file = paste(output.prefix, 'probe_df.txt', sep="_"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")




# Write

## make a plot of the results
#par(mfrow=c(2,2),pch=16,col=rgb(0,0,1,0.8),bty="n")
#hist(pbm.data[,1],breaks="Scott",main="Histogram of Probe Intensities")
#plot(pbm.8mers.pred, pbm.data[,1], main=paste("Probe Intensities from 8mer medians\n R^2 =", pbm.8mers.pred.rsqr), xlab="8mer predicted probe intensities", ylab="Probe Intensities")
#plot(pbm.beeml.pred.8mers, pbm.8mer.medians, main=paste("BEEML Fit on 8mer Median Intensities\n R^2 =", pbm.beeml.pred.8mers.rsqr), xlab="BEEML predicted 8mer medians", ylab="PBM 8mer Median Intensities")
#plot(pbm.beeml.pred, pbm.data[,1], main=paste("Probe Intensities from BEEML model\n R^2 =", pbm.beeml.pred.rsqr), xlab="BEEML predicted probe intensities", ylab="PBM Probe Intensities")

#png(paste(output.prefix,"-SnW.png",sep=""),width = 750, height = 600)
#seqLogo(as.matrix(read.table(seed.pwm.file,skip=1)[,-1]),2)
#dev.off()

#png(paste(output.prefix,"-BEEML.png",sep=""),width = 750, height = 600)
#seqLogo(apply(matrix(pbm.sol$mtx,4),2,function(x) exp(-x)/sum(exp(-x))))
#dev.off()

pretty.print.pwm(matrix(pbm.sol$mtx,4), 'bla', paste(output.prefix,"_BEEML_pwm.txt",sep=""))

freqmatrix=t(apply(exp(-matrix(pbm.sol$mtx,4)),1,"/",colSums(exp(-matrix(pbm.sol$mtx,4)))))
pretty.print.pwm(freqmatrix, 'bla', paste(output.prefix,"_BEEML_frequency_matrix.txt",sep=""), 6)
