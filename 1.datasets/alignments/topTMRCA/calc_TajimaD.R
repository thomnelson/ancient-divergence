source("~/scripts/TajimaD.R")

loci	<-	read.table("/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/3.data/0.datasets/alignments/topTMRCA/loci.txt")[,1]
DIR		<-	"/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/3.data/0.datasets/alignments/topTMRCA/"

poplabels	<-	c("bt",'bp','rs')
n.loci	<-	length(loci)

results	<-	matrix(nrow = n.loci, ncol = 4)

pb		<-	txtProgressBar(1, length(loci), style = 3)
for (i in 1:length(loci)) {
	locus	<-	loci[i]
	path.to.aln	<-	paste0(DIR,locus,".phylip")
	res	<-	TajimaD(path.to.aln,poplabels)
	results[i,]	<-	res
	setTxtProgressBar(pb, i)

}