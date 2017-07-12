# CALCULATE DXY FROM PHYLIP ALIGNMENTS FOR EACH LOCUS

library(ape)
source("/home/thom/research/scripts/dxy.R")
# dxy.R output (matrix format): return(c(pi,pi.x,pi.y,pi.xy,dxy))
DIR <- "/home/thom/research/seq/coal/akor_pstI/full/bbmap/stacks/rxstacks/phased_alaskanThreespines/fasta696/"
suffix <- ".fasta"

# pops1   <- c("bt","rs")
# pops2   <- c("bp","rs")
# pops3   <- c("bp","bt")
# pops4	<- c("bt|bp","rs")

pops1   <- c("rs","bt")
pops2   <- c("rs","bp")
pops3   <- c("bt","bp")
pops4	<- c("rs","bt|bp")

alns <- list.files(DIR)
alns <- alns[grepl(suffix,alns)]

n.alns <- length(alns)

locus  <- gsub(suffix,"",alns)

# pre-allocate memory for bigass matrix
# cols: pi, pi.rs, pi.bt, pi.bp, pi.fw, dxy.rsbt, dxy.rsbp, dxy.btbp ,dxy.rsfw
dstats <- matrix(nrow = n.alns, ncol = 9)

for (i in 1:n.alns) {
    aln	  <-	    alns[i]
    d	  <-	    c(2,2,2,2,2,2,2,2,2)
    # ADD IN PI, PI.RS, PI.BT, DXY.RSBT
    try(
	d[c(1,2,3,6)]  <-	    dxy(path.to.aln = paste0(DIR,aln,collapse=""), pops1, as.matrix = TRUE)[c(1,2,3,5)])
    # ADD IN PI.BP, DXY.RSBP
    try(
	d[c(4,7)]  <-	    dxy(path.to.aln = paste0(DIR,aln,collapse=""), pops2, as.matrix = TRUE)[c(3,5)])
    # ADD IN DXY.BTBP
    try(
	d[8]	<-	    dxy(path.to.aln = paste0(DIR,aln,collapse=""), pops3, as.matrix = TRUE)[5])
    # ADD IN PI.FW, DXY.RSFW
    try(
	d[c(5,9)]	<-	    dxy(path.to.aln = paste0(DIR,aln,collapse=""), pops4, as.matrix = TRUE)[c(3,5)])
    dstats[i,] <- d
    if (i %% 1000 == 0) {cat(paste0("Just got through ",as.character(i)," loci\r"))}
}
cat(paste0("Just got through ",as.character(i)," loci\n"))

dstats <- data.frame(locus,dstats)
names(dstats) <- c("locus","pi", "pi.rs", "pi.bt", "pi.bp", "pi.fw", "dxy.rsbt", "dxy.rsbp", "dxy.btbp" ,"dxy.rsfw")
write.table(dstats, file="/home/thom/research/seq/coal/akor_pstI/full/bbmap/stacks/rxstacks/phased_alaskanThreespines/dstats.tsv",
          sep = "\t", quote = FALSE, col.names = T, row.names=F)
