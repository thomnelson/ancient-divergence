#####################################
#
# THIS R SCRIPT CONTAINS ALL ANALYSES
#   AND PLOTTING CODE FOR THE MANUSCRIPT
#   ""
#
#####################################

# SPECIFY A HOME DIRECTORY WHERE THIS DIRECTORY EXISTS

HOME <- "/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/4.data/"

# FUNCTION FROM GLAZER et al. FOR CONVERTING BETWEEN GENOME ASSEMBLIES

source(paste0(HOME,"2.scripts/convertCoordinate.R"))
newScafs    <-    paste0(HOME,"NewScaffoldOrder.csv")
convertCoordinate(1, 21730178, 'old2new', newScafs)

# TABLE OF COORDINATES/LENGTHS FOR CHROMOSOMES

LGs.info    <-    read.table(paste0(HOME,"GacNewLGs.tsv"), sep = "\t", 
                                header = TRUE, stringsAsFactors = FALSE)
    group    <-    2    # ***column in which to find 'group' label (You specify)
    bp        <-    3    # ***column in which to find chromosome-specific bp position (You specify)

# DRAW LINES AT TO SEPARATE CHROMOSOMES ON WHOLE-GENOME PLOTS

chr.breaks    <-    function() {
    for (i in 2:length(LGs.info$groups)) {
        xright    <- LGs.info$LGs.breaks[i] / 1000000
        abline(v = xright, col = rgb(0.5,0.5,0.5,0.5))
    }
}

# genomic2Mb CONVERT 'GENOMIC' COORDINATE SYSTEM BACK INTO c(chr, Mb)

genomic2Mb <- function(coord, LGs = LGs.info) {
	if (coord > 20000) {cat("genomic coordinate should be in megabases! Exiting...\n"); break }
	chr.start.pos <- max(LGs.info$LGs.breaks[LGs.info$LGs.breaks <= coord * 1000000])
	chromosome    <- LGs.info$groups[LGs.info$LGs.breaks == chr.start.pos]
	new.coord     <- coord - (chr.start.pos / 1000000)
	return(list(chromosome, new.coord))
}

# PERMUTATION TEST P-VALUES

perm.test    <-    function(samp, source, n.perms = 1000, tail = 'lower', write.every = 100, graph.it = FALSE) {
    # takes two vectors, one from the sample you're testing, one representing all values
    # option 'tail' can be either 'lower' or 'upper'
    # first, get the mean and n of the sample
    samp.mean    <-    mean(samp)
    source.mean  <-    mean(source)
    samp.n       <-    length(samp)
    # next, set up an object to hold the resampled means
    rand.means    <-    NULL
    pval        <-    NULL
    # sample the source vector n.perms times, sampling samp.n entries
    for (i in 1:n.perms) {
        rand.mean    <-    mean(sample(source, samp.n))
        rand.means    <-    append(rand.means, rand.mean)
        if (i %% write.every == 0) {
            rand.means    <-    rand.means[order(rand.means)]
            ifelse(tail == 'lower', 
                   pval    <-    (length(rand.means[rand.means < samp.mean]) / length(rand.means)),
                   pval    <-    (length(rand.means[rand.means > samp.mean]) / length(rand.means)))
            cat(paste0("p-value after ",i," replicates: ",pval,"\n"))
        }
    }
    rand.means    <-    rand.means[order(rand.means)]
    ifelse(tail == 'lower', 
           pval    <-    (length(rand.means[rand.means < samp.mean]) / length(rand.means)),
           pval    <-    (length(rand.means[rand.means > samp.mean]) / length(rand.means)))
    result    <-    matrix(nrow = 1, ncol = 5)
    result[1,] <- c(source.mean, samp.mean, samp.n, i, pval)
    result    <-    as.data.frame(result)
    names(result) <- c("source.mean", "samp.mean", "samp.n", "n.perms", "p.value")
    result$p.val.string <- ifelse(tail == 'lower',
                                  paste0("≤ ", length(rand.means[rand.means < samp.mean])+1, "/", length(rand.means)),
                                  paste0("≤ ", length(rand.means[rand.means > samp.mean])+1, "/", length(rand.means))
                                  )
    if (graph.it == TRUE) {
    hist(c(rand.means,samp.mean), main = '', xlab = "sample means, permuted and estimated", breaks = 50)
        abline(v = samp.mean, col = 'blue', lwd = 2)
    }
    return(result)
}


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###    LOAD DATA
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

cat("Loading data from file...\n")
    coal <- read.table(
            paste0(HOME,"1.datasets/coalescence.tsv"), 
            sep = "\t", header = TRUE, stringsAsFactors=FALSE)
        coal$chr    <-    as.character(coal$chr)
        # coal        <-    coal[coal$chr != "19",]
    coal            <-    coal[order(coal$chr, coal$Mb),]

summary(coal$fst.rsbt[coal$fst.rsbt <= 1])
summary(coal$fst.rsbp[coal$fst.rsbp <= 1])
summary(coal$fst.bpbt[coal$fst.bpbt <= 1])
chromsNoD <- coal[coal$chr %in% c(2,6,15),]
summary(chromsNoD$fst.rsbt[chromsNoD$fst.rsbt <= 1])
summary(chromsNoD$fst.rsbp[chromsNoD$fst.rsbp <= 1])
summary(chromsNoD$fst.bpbt[chromsNoD$fst.bpbt <= 1])

cat("Loading segregating sites from file...\n")
    segsites    <-    read.table(
            paste0(HOME,"1.datasets/segsites.tsv"), 
            sep = "\t", header = TRUE, stringsAsFactors=FALSE)
     segsites    <-    segsites[segsites$locus %in% intersect(segsites$locus, coal$locus),]

inv1    <-    coal[coal$chr == 1 & coal$Mb > 26 & coal$Mb < 26.5 & coal$cons.sort == 'D',]
eda        <-    coal[coal$chr == 4 & coal$Mb > 12.79 & coal$Mb < 12.82 & coal$cons.sort == 'D',]

d    <-    coal[coal$cons.sort == 'D',]
# d[d$threespine.scaled > 7.468,c(2,4,14,40)]


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###    GENOME-WIDE SUMMARIES
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

hist(coal$threespine.scaled, breaks = seq(0, 15, by = 0.1), xlim = c(0,15))
hist(coal$boot.scaled, breaks = seq(0, 15, by = 0.1), xlim = c(0,15))
hist(coal$bepa.scaled, breaks = seq(0, 15, by = .1), xlim = c(0,15))
hist(coal$rs.scaled, breaks = seq(0, 15, by = 0.1), xlim = c(0,15))
# hist(coal$suqi, breaks = seq(0, .2, by = 0.001), xlim = c(0,0.05))

    summary(coal$threespine.scaled[coal$cons.sort != "N"])
    summary(coal$boot.scaled[coal$cons.sort != "N"])
    summary(coal$bepa.scaled[coal$cons.sort != "N"])
    summary(coal$fw.scaled[coal$cons.sort != "N"])
    summary(coal$rs.scaled[coal$cons.sort != "N"])

#    PLOT DENSITIES

bw = 0.5 # in MYA

wg    <-    density(coal$threespine.scaled[coal$cons.sort != "N"], 
                bw = bw, kernel = 'gaussian')
dl    <-    density(coal$threespine.scaled[coal$cons.sort == "D"],
                bw = bw, kernel = 'gaussian')
dl.not    <-    density(coal$threespine.scaled[coal$cons.sort != "D" & coal$cons.sort != "N"],
                bw = bw, kernel = 'gaussian')

plot(0,0, type = 'n', xlim = c(0,15), ylim = c(0,0.30), xlab = "tmrca", ylab = 'density')
    polygon(x = c(wg$x, rev(wg$x)), y = c(wg$y, rep(0,length(wg$y))),
            col = rgb(0.9,0.9,0.9))
    polygon(x = c(dl.not$x, rev(dl.not$x)), y = c(dl.not$y, rep(0,length(dl.not$y))),
            col = rgb(0.9,0.9,0.9, 0.25))
    polygon(x = c(dl$x, rev(dl$x)), y = c(dl$y, rep(0,length(dl$y))),
            col = rgb(0.0,0.0,0.0, 0.25))

# PRINT OUT GENOME-WIDE STATISTICS WITH PERMUTED P-VALUES

pi.permute    <-    perm.test(coal$pi[coal$cons.sort == "D"], source = coal$pi, 
                              n.perms = 10000, tail = 'upper', write.every = 100000)
pi.bt.permute    <-    perm.test(coal$pi.bt[coal$cons.sort == "D"], source = coal$pi.bt, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
pi.bp.permute    <-    perm.test(coal$pi.bp[coal$cons.sort == "D"], source = coal$pi.bp, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
pi.fw.permute    <-    perm.test(coal$pi.fw[coal$cons.sort == "D"], source = coal$pi.fw, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
pi.rs.permute    <-    perm.test(coal$pi.rs[coal$cons.sort == "D"], source = coal$pi.rs, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
theta.w.permute    <-    perm.test(coal$theta.w[coal$cons.sort == "D"], 
                              source = coal$theta.w.threespine, 
                              n.perms = 10000, tail = 'upper', write.every = 100000)
theta.w.bt.permute    <-    perm.test(coal$theta.w.bt[coal$cons.sort == "D"], 
                              source = coal$theta.w.bt, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
theta.w.bp.permute    <-    perm.test(coal$theta.w.bp[coal$cons.sort == "D"], 
                              source = coal$theta.w.bp, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
theta.w.fw.permute    <-    perm.test(coal$theta.w.fw[coal$cons.sort == "D"], 
                              source = coal$theta.w.fw, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
theta.w.rs.permute    <-    perm.test(coal$theta.w.rs[coal$cons.sort == "D"], 
                              source = coal$theta.w.rs, 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
dxy.rsfw.permute    <-    perm.test(coal$dxy.rsfw[coal$cons.sort == "D"], source = coal$dxy.rsfw, 
                              n.perms = 10000, tail = 'upper', write.every = 100000)
tmrca.permute    <-    perm.test(coal$threespine.scaled[coal$cons.sort == "D"], source = coal$threespine.scaled[coal$cons.sort != "N"], 
                              n.perms = 10000, tail = 'upper', write.every = 100000)
tmrca.fw.permute    <-    perm.test(coal$fw.scaled[coal$cons.sort == "D"], source = coal$fw.scaled[coal$cons.sort != "N"], 
                              n.perms = 10000, tail = 'lower', write.every = 100000)
tmrca.rs.permute    <-    perm.test(coal$rs.scaled[coal$cons.sort == "D"], source = coal$rs.scaled[coal$cons.sort != "N"], 
                              n.perms = 10000, tail = 'lower', write.every = 100000)

cat(paste0("\npi\n",
           "genome-wide mean:   ",mean(coal$pi),'\n',
           "genome-wide median: ",median(coal$pi),'\n',
           "recM mean:          ",mean(coal$pi[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$pi[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", pi.permute$p.val.string,"\n",
           "\npi Boot\n",
           "genome-wide mean:   ",mean(coal$pi.bt),'\n',
           "genome-wide median: ",median(coal$pi.bt),'\n',
           "recM mean:          ",mean(coal$pi.bt[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$pi.bt[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", pi.bt.permute$p.val.string,"\n",
           "\npi Bear Paw\n",
           "genome-wide mean:   ",mean(coal$pi.bp),'\n',
           "genome-wide median: ",median(coal$pi.bp),'\n',
           "recM mean:          ",mean(coal$pi.bp[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$pi.bp[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", pi.bp.permute$p.val.string,"\n",
           "\npi Freshwater combined\n",
           "genome-wide mean:   ",mean(coal$pi.fw),'\n',
           "genome-wide median: ",median(coal$pi.fw),'\n',
           "recM mean:          ",mean(coal$pi.fw[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$pi.fw[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", pi.fw.permute$p.val.string,"\n",
           "\npi Rabbit Slough\n",
           "genome-wide mean:   ",mean(coal$pi.rs),'\n',
           "genome-wide median: ",median(coal$pi.rs),'\n',
           "recM mean:          ",mean(coal$pi.rs[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$pi.rs[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", pi.rs.permute$p.val.string,"\n",
           "\ndxy RS-FW\n",
           "genome-wide mean:   ",mean(coal$dxy.rsfw),'\n',
           "genome-wide median: ",median(coal$dxy.rsfw),'\n',
           "recM mean:          ",mean(coal$dxy.rsfw[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$dxy.rsfw[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", dxy.rsfw.permute$p.val.string,"\n",
           "\ntmrca\n",
           "genome-wide mean:   ",mean(coal$threespine.scaled[coal$cons.sort != "N"]),'\n',
           "genome-wide median: ",median(coal$threespine.scaled[coal$cons.sort != "N"]),'\n',
           "recM mean:          ",mean(coal$threespine.scaled[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$threespine.scaled[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", tmrca.permute$p.val.string,"\n",
           "\ntmrca Freshwater\n",
           "genome-wide mean:   ",mean(coal$fw.scaled[coal$cons.sort != "N"]),'\n',
           "genome-wide median: ",median(coal$fw.scaled[coal$cons.sort != "N"]),'\n',
           "recM mean:          ",mean(coal$fw.scaled[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$fw.scaled[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", tmrca.fw.permute$p.val.string,"\n",
           "\ntmrca Rabbit Slough\n",
           "genome-wide mean:   ",mean(coal$rs.scaled[coal$cons.sort != "N"]),'\n',
           "genome-wide median: ",median(coal$rs.scaled[coal$cons.sort != "N"]),'\n',
           "recM mean:          ",mean(coal$rs.scaled[coal$cons.sort == "D"]),'\n',
           "recM median:        ",median(coal$rs.scaled[coal$cons.sort == "D"]),'\n',
           "permuted p-value of difference in means: ", tmrca.rs.permute$p.val.string,"\n"
    )
)


print(ttest <- t.test(coal$pi.fw, coal$pi.rs))
print(ttest <- t.test(coal$pi.bt, coal$pi.rs))
print(ttest <- t.test(coal$pi.bp, coal$pi.rs))
print(ttest <- t.test(coal$pi.bt, coal$pi.bp))

print(ttest <- t.test(coal$pi, coal$dxy.rsfw))
wilcox.test(coal$pi, coal$dxy.rsfw)

# T-TESTS USING CHROMOSOMES WITH LEAST EVIDENCE OF DIVERGENT SELECTION

neu <- coal[coal$cons.sort != 'N',]
neu <- neu[neu$chr %in% c(3,6,15),]

print(ttest <- t.test(neu$pi.fw[neu$cons.sort != "N"], neu$pi.rs[neu$cons.sort != "N"]))
print(ttest <- t.test(neu$pi.bt[neu$cons.sort != "N"], neu$pi.rs[neu$cons.sort != "N"]))
print(ttest <- t.test(neu$pi.bp[neu$cons.sort != "N"], neu$pi.rs[neu$cons.sort != "N"]))
print(ttest <- t.test(neu$pi.bt[neu$cons.sort != "N"], neu$pi.bp[neu$cons.sort != "N"]))
print(ttest <- t.test(neu$pi.bt[neu$cons.sort != "N"], neu$pi.bp[neu$cons.sort != "N"]))

# MANN-WHITNEY TESTS FOR DIFFERENCE IN PAIRWISE FST

wilcox.test(neu$fst.rsbt[neu$fst.rsbt <= 1], neu$fst.bpbt[neu$fst.bpbt <= 1])
wilcox.test(neu$fst.rsbt[neu$fst.rsbt <= 1], neu$fst.rsbp[neu$fst.rsbp <= 1])
wilcox.test(neu$fst.rsbp[neu$fst.rsbp <= 1], neu$fst.bpbt[neu$fst.bpbt <= 1])

# HISTOGRAMS OF N. SEGREGATING SITES AND DISTANCES

hist(segsites$segsites, xlim = c(0,40), breaks = 0:155, 
    xlab = 'segregating sites per locus', main = '')
    length(segsites$segsites[segsites$segsites <= 40])  / length(segsites$segsites)
    max(segsites$segsites, na.rm = T)
    abline(v = mean(segsites$segsites, na.rm = T), lty = 'dashed')
    abline(v = median(segsites$segsites, na.rm = T), lwd = 2)

chrs    <-    as.character(1:21)
dists    <-    NULL
for (chr in chrs) {
    cat(paste0("Scanning chr ",chr,"\n"))
    bps = coal$bp[coal$chr == chr]
    for (i in 2:length(bps)) {
        d    <-    bps[i] - bps[i-1]
        dists    <-    append(dists,d)
    }    
}
dists.kb    <-    dists / 1000

hist(dists.kb, xlim = c(0,50), breaks = seq(0,500,by=1),
    xlab = 'distance between adjacent RAD loci (cutsites)', main="")
    abline(v = mean(dists.kb), lty = 'dashed')
    abline(v = median(dists.kb), lwd = 2)
    length(dists.kb[dists.kb <= 50])  / length(dists.kb)
    max(dists.kb)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###                         SLIDING WINDOWS                             ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

options(scipen=999)

# MANIPULATE ONLY THE DATA I NEED TO

c        <-    coal[, c(2:44,46)]
c$chr    <-    as.numeric(c$chr)
c$fst.rsbp[c$fst.rsbp == 99] <- 0
c$fst.rsbt[c$fst.rsbt == 99] <- 0
c$fst.rsfw[c$fst.rsfw == 99] <- 0
c$fst.bpbt[c$fst.bpbt == 99] <- 0
w.size    <-    0.5 # size of windows, in Mb
w.step    <-    0.1 # step size, in Mb


# MOVING CHROMOSOME BY CHROMOSOME
n.cols    <-    dim(c)[2]

windowed    <-    NULL

for (i in 1:21) {
    chrom    <-    c[c$chr == i,]
    start    <-    0
    end        <-    max(chrom$Mb)
    windows    <-    seq(start, end, by = w.step)
    n.rows    <-    length(windows)
    cat(paste0("Steppin\' thru chromosome ",i,"\n"))
    # CREATE MATRIX TO ALLOCATE MEMORY
    result    <-    matrix(nrow= n.rows, ncol = n.cols)
    for (j in 1:n.rows) {
        w.start    <-    windows[j]
        w.end    <-    w.start + w.size
        x        <-    chrom[ chrom$Mb > w.start & chrom$Mb <= w.end , ]
        res        <-    apply(x, 2, mean, na.rm=T)
        result[j,]    <-    res
        }
    windowed    <-    rbind(windowed, result)
}

windowed    <-    as.data.frame(windowed)
names(windowed)    <-    names(c)

plot(windowed$genomic, windowed$threespine.scaled, type = 'l')
plot(windowed$genomic, windowed$fst.rsbt, type = 'l')


# USING GENOMIC COORDINATES
n.cols    <-    dim(c)[2]

start    <-    0
end        <-    max(chrom$genomic)
windows    <-    seq(start, end, by = w.step)
n.rows    <-    length(windows)
# CREATE MATRIX TO ALLOCATE MEMORY
windowed    <-    matrix(nrow= n.rows, ncol = n.cols)

pb    <-    txtProgressBar(1, n.rows, style = 3)
for (i in 1:n.rows) {
    w.start    <-    windows[i]
    w.end    <-    w.start + w.size
    x        <-    c[ c$genomic > w.start & c$genomic <= w.end , ]
    res        <-    apply(x, 2, mean, na.rm=T)
    windowed[i,]    <-    res
    setTxtProgressBar(pb, i)
}

windowed    <-    as.data.frame(windowed)
names(windowed)    <-    names(c)

plot(windowed$genomic, windowed$threespine.scaled, type = 'l') ;chr.breaks()
plot(windowed$genomic, windowed$dxy.rsfw, type = 'l') ;chr.breaks()
plot(windowed$genomic, windowed$rs.scaled, type = 'l') ; chr.breaks()
plot(windowed$genomic, windowed$fw.scaled, type = 'l') ; chr.breaks()
plot(windowed$genomic, windowed$fst.rsbt, type = 'l') ; chr.breaks()
plot(windowed$genomic, windowed$pi, type = 'l') ; chr.breaks()
plot(windowed$genomic, windowed$theta.w.threespine, type = 'l') ; chr.breaks()
plot(windowed$genomic, windowed$tajima.D.threespine, type = 'l') ; chr.breaks()
plot(windowed$genomic, windowed$tajima.D.rs, type = 'l', col = rgb(1,.5,.5,.75))
    chr.breaks()
    lines(windowed$genomic, windowed$tajima.D.fw, col = rgb(0.5,.5,1,.5))
plot(windowed$genomic, windowed$tajima.D.fw, type = 'l') ; chr.breaks()




###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###         RELATIVE AND ABSOLUTE SEQUENCE DIVERGENCE                   ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###



summary(coal$fst.rsbp[coal$fst.rsbp <=1], na.rm=T)
summary(coal$fst.rsbt[coal$fst.rsbt <=1], na.rm=T)
summary(coal$fst.bpbt[coal$fst.bpbt <=1], na.rm=T)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### FST -- WHOLE GENOME
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

coal    <-    coal[order(coal$genomic),]
bw = 0.5

c              <-    coal[coal$fst.rsbt <= 1,]
fst.sm.rsbt    <-    ksmooth(x = c$genomic, y = c$fst.rsbt, bandwidth = bw, kernel='normal')
c              <-    coal[coal$fst.rsbp <= 1,]
fst.sm.rsbp    <-    ksmooth(x = c$genomic, y = c$fst.rsbp, bandwidth = bw, kernel='normal')
c              <-    coal[coal$fst.bpbt <= 1,]
fst.sm.bpbt    <-    ksmooth(x = c$genomic, y = c$fst.bpbt, bandwidth = bw, kernel='normal')

ymin    <-    0
ymax    <-    1
plot(0,0, type = 'n', xlim = c(0,450), ylim = c(ymin,ymax))
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(c$genomic[c$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
        lines(fst.sm.rsbt, col = 'blue3')
        lines(fst.sm.rsbp, col = 'steelblue3')
        lines(fst.sm.bpbt, col = 'gray50')


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### FST -- SINGLE CHROMOSOMES
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

chr <- 1
bw = 0.5

coal    <-    coal[order(coal$genomic),]

c    <-    coal[coal$fst.rsbt <= 1 & coal$chr == chr,]
fst.sm.rsbt    <-    ksmooth(x = c$Mb, y = c$fst.rsbt, bandwidth = bw, kernel='normal')
c    <-    coal[coal$fst.rsbp <= 1 & coal$chr == chr,]
fst.sm.rsbp    <-    ksmooth(x = c$Mb, y = c$fst.rsbp, bandwidth = bw, kernel='normal')
c    <-    coal[coal$fst.bpbt <= 1 & coal$chr == chr,]
fst.sm.bpbt    <-    ksmooth(x = c$Mb, y = c$fst.bpbt, bandwidth = bw, kernel='normal')


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### DXY -- WHOLE GENOME
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

bw = 0.5 

c    <-    coal
#c    <-    coal[coal$fst.rsbt <= 1,]
dxy.sm.rsbt    <-    ksmooth(x = l$genomic, y = l$dxy.rsbt, bandwidth = bw, kernel='normal')
#c    <-    coal[coal$dxy.rsbp <= 1,]
dxy.sm.rsbp    <-    ksmooth(x = l$genomic, y = l$dxy.rsbp, bandwidth = bw, kernel='normal')
#c    <-    coal[coal$dxy.bpbt <= 1,]
dxy.sm.bpbt    <-    ksmooth(x = l$genomic, y = l$dxy.btbp, bandwidth = bw, kernel='normal')

ymin    <-    0
ymax    <-    0.03
plot(0,0, type = 'n', xlim = c(0,450), ylim = c(ymin,ymax))
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(l$genomic[l$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
        lines(dxy.sm.rsbt, col = 'blue3')
        lines(dxy.sm.rsbp, col = 'steelblue3')
        lines(dxy.sm.bpbt, col = 'gray50')


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###                     SEQUENCE DIVERSITY STATISTICS                   ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

bw    <-    0.5

pi.threespine    <-    ksmooth(coal$genomic, coal$pi, 
                            bandwidth = bw, kernel = 'normal')
pi.rs            <-    ksmooth(coal$genomic, coal$pi.rs, 
                            bandwidth = bw, kernel = 'normal')
pi.fw            <-    ksmooth(coal$genomic, coal$pi.fw, 
                            bandwidth = bw, kernel = 'normal')

theta.threespine    <-    ksmooth(coal$genomic, coal$pi, 
                            bandwidth = bw, kernel = 'normal')
theta.rs            <-    ksmooth(coal$genomic, coal$theta.w.rs, 
                            bandwidth = bw, kernel = 'normal')
theta.fw            <-    ksmooth(coal$genomic, coal$theta.w.fw, 
                            bandwidth = bw, kernel = 'normal')

c                    <-    coal[!(is.na(coal$tajima.D.threespine)),]
tajima.D.threespine    <-    ksmooth(c$genomic, c$tajima.D.threespine, 
                            bandwidth = bw, kernel = 'normal')
c                    <-    coal[!(is.na(coal$tajima.D.rs)),]
tajima.D.rs            <-    ksmooth(c$genomic, c$tajima.D.rs, 
                            bandwidth = bw, kernel = 'normal')
c                    <-    coal[!(is.na(coal$tajima.D.fw)),]
tajima.D.fw            <-    ksmooth(c$genomic, c$tajima.D.fw, 
                            bandwidth = bw, kernel = 'normal')


# PI 

ymin    <-    0
ymax    <-    0.02
xmin    <-    0
xmax    <-    max(pi.threespine$x)
plot(0,0, type = 'n', xlim = c(xmin,xmax), ylim = c(ymin,ymax), ylab = 'pi')
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(coal$genomic[coal$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
    lines(pi.threespine, lwd = 1)
    lines(pi.rs, lwd = 1, col = 'steelblue3')
    lines(pi.fw, lwd = 1, col = 'red3')

# WATTERSON'S THETA

ymin    <-    0
ymax    <-    0.02
xmin    <-    0
xmax    <-    max(theta.threespine$x)
plot(0,0, type = 'n', xlim = c(xmin,xmax), ylim = c(ymin,ymax), ylab = 'Watterson\'s theta')
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(coal$genomic[coal$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
    lines(theta.threespine, lwd = 1)
    lines(theta.rs, lwd = 1, col = 'steelblue3')
    lines(theta.fw, lwd = 1, col = 'red3')

# TAJIMA'S D

ymin    <-    -2
ymax    <-    3
xmin    <-    0
xmax    <-    max(tajima.D.threespine$x)
# xmin <- 330
# xmax <- 360
plot(0,0, type = 'n', xlim = c(xmin,xmax), ylim = c(ymin,ymax), ylab = "Tajima's D")
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(coal$genomic[coal$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
    lines(tajima.D.threespine, lwd = 1)
    lines(tajima.D.rs, lwd = 1, col = 'red3')
    lines(tajima.D.fw, lwd = 1, col = 'steelblue3')

# SEQUENCE DIVERSITY AT D VS ALL

d    <-    coal[coal$cons.sort == 'D',]

hist(d$pi, breaks = seq(0,0.2, by = 0.0005), xlim = c(0,0.02),
        xlab = 'pi', main = '', col = rgb(0.5,0.5,0.5))
    par(new=T)
    hist(coal$pi, breaks = seq(0,0.2, by = 0.0005), xlim = c(0,0.02),
        xaxt='n',yaxt='n', main = 'n', xlab = '', ylab = '',
        col = rgb(1,1,1,0.25))

hist(d$tajima.D.threespine, breaks = seq(-5,5, by = 0.1), xlim = c(-4,4),
        xlab = 'Tajima\'s D', main = '', col = rgb(.5,.5,.5))
    par(new=T)
    hist(coal$tajima.D.threespine, breaks = seq(-5,5, by = 0.1), xlim = c(0-4,4),
        xaxt='n',yaxt='n', main = '', xlab = '', ylab = '',
        col = rgb(1,1,1,0.25))

hist(d$tajima.D.rs, breaks = seq(-5,5, by = 0.1), xlim = c(-4,4),
        xlab = 'Tajima\'s D', main = '', col = rgb(.5,.5,.5))
    par(new=T)
    hist(coal$tajima.D.rs, breaks = seq(-5,5, by = 0.1), xlim = c(0-4,4),
        xaxt='n',yaxt='n', main = '', xlab = '', ylab = '',
        col = rgb(1,1,1,0.25))

hist(d$tajima.D.fw, breaks = seq(-5,5, by = 0.1), xlim = c(-4,4),
        xlab = 'Tajima\'s D', main = '', col = rgb(.5,.5,.5))
    par(new=T)
    hist(coal$tajima.D.fw, breaks = seq(-5,5, by = 0.1), xlim = c(0-4,4),
        xaxt='n',yaxt='n', main = '', xlab = '', ylab = '',
        col = rgb(1,1,1,0.25))




###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###         POPULATION/ECOTYPE-SPECIFIC LINEAGE SORTING                 ###
###                 USING MCC TREE TOPOLOGIES                           ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# sorting categories: BT only, BP only, FW combined, RS
# 3to9 lineage sorting must be observed in > 50% of sampled trees

c            <-    coal[coal$mono.3to9 == 1,]
ncats        <-    4
# btsort        <-    c$genomic[c$mono.bt == 1 & c$mono.fw == 0]
# bpsort        <-    c$genomic[c$mono.bp == 1 & c$mono.fw == 0]
# fwsort        <-    c$genomic[c$mono.fw == 1]
fwsort        <-    c$genomic[c$mono.fw == 1 & c$mono.oc == 0]
# fwunsort    <-    c$genomic[c$mono.bt == 1 & c$mono.bp == 1 & c$mono.fw == 1]
rssort        <-    c$genomic[c$mono.oc == 1 & c$mono.fw == 0]
recM        <-    c$genomic[c$mono.fw == 1 & c$mono.oc == 1]

xmin    <-    -50
xmax    <-    max(coal$genomic)
ymin    <-    -1
ymax    <-    ncats

plot(0,0,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax), yaxt = 'n',
        xlab="",ylab="",bty='n')
    # segments(x0 = btsort,x1=btsort, y0=4,y1=3, col = "blue3")
    segments(x0 = fwunsort,x1=fwunsort, y0=4,y1=3, col = "green3")
#    segments(x0 = bpsort,x1=bpsort, y0=3,y1=2, col = "steelblue3")
    segments(x0 = fwsort,x1=fwsort, y0=2,y1=1, col = "steelblue")
    segments(x0 = rssort,x1=rssort, y0=1,y1=0, col = "red3")
    segments(x0 = recM,x1=recM, y0=-1,y1=0, col = "gray50")
    text(x = 0, y = 3.5, labels = paste0("sep: ",length(fwunsort)), pos = 2)
    # text(x = 0, y = 3.5, labels = paste0("bt: ",length(btsort)), pos = 2)
    # text(x = 0, y = 2.5, labels = paste0("bp: ",length(bpsort)), pos = 2)
    text(x = 0, y = 1.5, labels = paste0("same: ",length(fwsort)), pos = 2)
    text(x = 0, y = 0.5, labels = paste0("rs: ",length(rssort)), pos = 2)
    text(x = 0, y = -0.5, labels = paste0("recM: ",length(recM)), pos = 2)

# plot cutoff vs n loci with that level of evidence

cutoffs    <-    seq(0,1,by=0.01)

nloci    <-    matrix(nrow = length(cutoffs)-1, ncol = ncats)

for (i in 2:length(cutoffs))    {
    cutoff    <-    cutoffs[i]
    bt    <-    length(btsort    <-    l$genomic[l$mono.bt >= cutoff])
    bp    <-    length(bpsort    <-    l$genomic[l$mono.bp >= cutoff])
    fw    <-    length(fwsort    <-    l$genomic[l$mono.fw >= cutoff])
    rs    <-    length(rssort    <-    l$genomic[l$mono.oc >= cutoff])
    nloci[i-1,]    <-    c(bt,bp,fw,rs)
}

ymax    <-    max(nloci)
xvals    <-    cutoffs[2:length(cutoffs)]
plot(0,0,type='n',xlim=range(cutoffs),ylim=c(0,ymax),
    xlab='% sampled trees with given clade',ylab='RAD loci (n)')
    lines(x = xvals, y = nloci[,1], col = 'blue3')
    lines(x = xvals, y = nloci[,2], col = 'steelblue3')
    lines(x = xvals, y = nloci[,3], col = 'steelblue', lwd = 2)
    lines(x = xvals, y = nloci[,4], col = 'red3')


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### LOGIT REGRESSION OF "IS.RECIPROCALLY.MONOPHYLETIC" ONTO FST ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# COMBINED FW
x        <-    coal[coal$fst.rsfw <= 1,]
    x$recM    <-    rep(0, length(x$bp))
    x$recM[x$mono.fw == 1 & x$mono.oc == 1]    <-    1
    fst.rsfw    <-    x$fst.rsfw[order(x$fst.rsfw)]

model     <- glm(recM ~ fst.rsfw, family=binomial(link="logit"), data = x)
    rec.fit    <-    predict(model, data.frame(len=fst.rsfw), type = 'response')

plot(x$fst.rsfw, x$recM, pch = 20, col = rgb(.5,.5,.5,.25))
    lines(fst.rsfw, rec.fit)


# RS-BL
x        <-    coal[coal$fst.rsbt <= 1,]
    x$recM    <-    rep(0, length(x$bp))
    x$recM[x$mono.fw == 1 & x$mono.oc == 1]    <-    1
    fst.rsbt    <-    x$fst.rsbt[order(x$fst.rsbt)]

model     <- glm(recM ~ fst.rsbt, family=binomial(link="logit"), data = x)
    rec.fit    <-    predict(model, data.frame(len=fst.rsbt), type = 'response')

plot(x$fst.rsbt, x$recM, pch = 20, col = rgb(.5,.5,.5,.25))
    lines(fst.rsbt, rec.fit)


# RS-BP
x        <-    coal[coal$fst.rsbp <= 1,]
    x$recM    <-    rep(0, length(x$bp))
    x$recM[x$mono.fw == 1 & x$mono.oc == 1]    <-    1
    fst.rsbp    <-    x$fst.rsbp[order(x$fst.rsbp)]

model     <- glm(recM ~ fst.rsbp, family=binomial(link="logit"), data = x)
    rec.fit    <-    predict(model, data.frame(len=fst.rsbp), type = 'response')

plot(x$fst.rsbp, x$recM, pch = 20, col = rgb(.5,.5,.5,.25))
    lines(fst.rsbp, rec.fit)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###                  NON-OVERLAPPING GENOMIC WINDOWS                    ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

library(lmodel2)

all.window    <-    NULL
window        <-    0.25 #size, in Mb

for (i in 1:21) {
#    i = 1
    d = coal[coal$chr == i,]
    cat(paste0("Scanning chromosome",as.character(i),"...\n"))
    start = 0
    wind.end = start + window
    end = max(d$Mb)
    
    Mb  = NULL
    genomic  = NULL
    pi    = NULL
    pi.bt     = NULL
    pi.bp     = NULL
    pi.rs     = NULL
    pi.fw     = NULL
    dxy.rsbt     = NULL
    dxy.rsbp     = NULL
    fst.rsbt     = NULL
    fst.rsbp     = NULL
    d.fract         = NULL
    threespine.scaled     = NULL
    threespine.noD        = NULL
    while (start < end) {
        loci = d[d$Mb >= start & d$Mb < wind.end,]
        nl    <-    length(loci[,1])
        Mb  <- append(Mb,  mean(loci$Mb))
        genomic  <- append(genomic,  mean(loci$genomic))
        pi    <- append(pi, mean(loci$pi))
        pi.bt <- append(pi.bt, mean(loci$pi.bt))
        pi.bp <- append(pi.bp, mean(loci$pi.bp))
        pi.rs <- append(pi.rs, mean(loci$pi.rs))
        dxy.rsbt <- append(dxy.rsbt, mean(loci$dxy.rsbt))
        dxy.rsbp <- append(dxy.rsbp, mean(loci$dxy.rsbp))
        fst.rsbt <- append(fst.rsbt, mean(loci$fst.rsbt[loci$fst.rsbt <=1]))
        fst.rsbp <- append(fst.rsbp, mean(loci$fst.rsbp[loci$fst.rsbp <=1]))
        threespine.scaled <- append(threespine.scaled, mean(loci$threespine.scaled))
        threespine.noD    <- append(threespine.noD, mean(loci$threespine.scaled[loci$cons.sort <= 'D' & loci$fst.rsfw <= 0.5]))
        sort    <-    loci$cons.sort
        nd    <-    length(sort[sort == "D"])
        df    <-    nd / nl
        d.fract    <-    append(d.fract, df)
        start <- start + window
        wind.end <- wind.end + window
    }
    chr = rep(i, length(Mb))
    x = data.frame(chr,Mb,genomic,pi,pi.bt,pi.bp,pi.rs,dxy.rsbt,dxy.rsbp,fst.rsbt,fst.rsbp,threespine.scaled,threespine.noD,d.fract)
    all.window    <-    rbind(all.window,x)

}


plot(all.window$fst.rsbt, all.window$dxy.rsbt, type = 'p', pch = 20, col = rgb(.5,.5,.5,.5),
     xlab = "Fst, RS vs BL", ylab = "dXY, RS vs BL", ylim = c(0,0.03), xlim = c(0,0.8))
plot(all.window$fst.rsbp, all.window$dxy.rsbp, type = 'p', pch = 20, col = rgb(.5,.5,.5,.5),
     xlab = "Fst, RS vs BP", ylab = "dXY, RS vs BP", ylim = c(0,0.03), xlim = c(0,0.8))
plot(all.window$genomic, all.window$d.fract, type = 'p', pch = 20, col = rgb(.5,.5,.5,.5))
plot(all.window$threespine.scaled[all.window$d.fract != 0], all.window$d.fract[all.window$d.fract != 0], 
     type = 'p', pch = 20, col = rgb(.5,.5,.5,.5))

plot(all.window$threespine.scaled, all.window$threespine.noD, type = 'p', pch = 20, col = rgb(.5,.5,.5,.5),
     xlab = 'TMRCA, all loci', ylab = "TMRCA, differentiation outliers removed")


dxy.rma.bt	<-	lmodel2(dxy.rsbt ~ fst.rsbt, data = all.window, 
					range.x = "interval", range.y = 'relative', nperm = 99)
dxy.rma.bp	<-	lmodel2(dxy.rsbp ~ fst.rsbp, data = all.window, 
					range.x = "interval", range.y = 'relative', nperm = 99)



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###                             CHEVRON PLOTS                           ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

        ###~~~~~~~~~~~~~~~~~~~~~###
        ###     FW COMBINED
        ###~~~~~~~~~~~~~~~~~~~~~###

d    <-    coal[coal$cons.sort == 'D',]
nd    <-    dim(d)[1]
eda    <-    coal[coal$locus == 30517,]

d    <-    d[order(d$pi.fw),]

div        <-    c(d$pi.rs,d$dxy.rsfw,d$pi.fw)
pop        <-    c(rep("Marine",nd), rep("Total",nd), rep("Freshwater",nd))
pop        <-    factor(pop, levels = c("Marine","Total","Freshwater"))
d.cat    <-    data.frame(pop,div)

ylim = c(0,0.05)
xlim = c(0.5,3.5)

plot(0,0, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab='',bty='n')
    abline(h = median(coal$dxy.rsfw, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'black')
    abline(h = median(coal$pi, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'gray50')
    abline(h = median(coal$pi.rs, na.rm=T), lwd = 3, lty = 'dashed', col = 'red3')
    abline(h = median(coal$pi.fw, na.rm=T), lwd = 3, lty = 'dashed', col = 'blue3')
    segments(x0 = 1, x1 = 2, y0 = d$pi.rs, y1 = d$dxy.rsfw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 2, x1 = 3, y0 = d$dxy.rsfw, y1 = d$pi.fw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 1, x1 = 2, y0 = eda$pi.rs, y1 = eda$dxy.rsfw,
            col = 'goldenrod3', lwd = 2)
    segments(x0 = 3, x1 = 2, y0 = eda$pi.fw, y1 = eda$dxy.rsfw,
            col = 'goldenrod3', lwd = 2)
    boxplot(div ~ pop, d.cat, border = c('red3','black','blue3'), add = TRUE,
            ylim = ylim, ylab = "pi.rs, dxy, pi.bt",
            pch = 20, cex = 0.5)


        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###     FW COMBINED - highlighting all in inv1
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

library(wesanderson)
inv.col <- wes_palette("Darjeeling")[2]

d     <-    coal[coal$cons.sort == 'D',]
nd    <-    dim(d)[1]
dim(inv1)

d    <-    d[order(d$pi.fw),]

div      <-    c(d$pi.rs,d$dxy.rsfw,d$pi.fw)
pop      <-    c(rep("Marine",nd), rep("Total",nd), rep("Freshwater",nd))
pop      <-    factor(pop, levels = c("Marine","Total","Freshwater"))
d.cat    <-    data.frame(pop,div)

ylim = c(0,0.05)
xlim = c(0.5,3.5)

plot(0,0, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab='',bty='n')
    abline(h = median(coal$dxy.rsfw, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'black')
    abline(h = median(coal$pi, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'gray50')
    abline(h = median(coal$pi.rs, na.rm=T), lwd = 3, lty = 'dashed', col = 'red3')
    abline(h = median(coal$pi.fw, na.rm=T), lwd = 3, lty = 'dashed', col = 'blue3')
    segments(x0 = 1, x1 = 2, y0 = d$pi.rs, y1 = d$dxy.rsfw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 2, x1 = 3, y0 = d$dxy.rsfw, y1 = d$pi.fw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 1, x1 = 2, y0 = inv1$pi.rs, y1 = inv1$dxy.rsfw,
            col = inv.col, lwd = 2)
    segments(x0 = 3, x1 = 2, y0 = inv1$pi.fw, y1 = inv1$dxy.rsfw,
            col = inv.col, lwd = 2)
    boxplot(div ~ pop, d.cat, border = c('red3','black','blue3'), add = TRUE,
            ylim = ylim, ylab = "pi.rs, dxy, pi.bt",
            pch = 20, cex = 0.5)


        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###     FW COMBINED - highlighting all in inv21
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

library(wesanderson)
inv.col <- wes_palette("Darjeeling")[5]

d     <-    coal[coal$cons.sort == 'D',]
nd    <-    dim(d)[1]
inv21    <-    coal[coal$chr == 21 & coal$Mb > 10 & coal$Mb < 11.7 & coal$cons.sort == 'D',]

d    <-    d[order(d$pi.fw),]

div      <-    c(d$pi.rs,d$dxy.rsfw,d$pi.fw)
pop      <-    c(rep("Marine",nd), rep("Total",nd), rep("Freshwater",nd))
pop      <-    factor(pop, levels = c("Marine","Total","Freshwater"))
d.cat    <-    data.frame(pop,div)

ylim = c(0,0.05)
xlim = c(0.5,3.5)

plot(0,0, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab='',bty='n')
    abline(h = median(coal$dxy.rsfw, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'black')
    abline(h = median(coal$pi, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'gray50')
    abline(h = median(coal$pi.rs, na.rm=T), lwd = 3, lty = 'dashed', col = 'red3')
    abline(h = median(coal$pi.fw, na.rm=T), lwd = 3, lty = 'dashed', col = 'blue3')
    segments(x0 = 1, x1 = 2, y0 = d$pi.rs, y1 = d$dxy.rsfw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 2, x1 = 3, y0 = d$dxy.rsfw, y1 = d$pi.fw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 1, x1 = 2, y0 = inv21$pi.rs, y1 = inv21$dxy.rsfw,
            col = inv.col, lwd = 2)
    segments(x0 = 3, x1 = 2, y0 = inv21$pi.fw, y1 = inv21$dxy.rsfw,
            col = inv.col, lwd = 2)
    boxplot(div ~ pop, d.cat, border = c('red3','black','blue3'), add = TRUE,
            ylim = ylim, ylab = "pi.rs, dxy, pi.bt",
            pch = 20, cex = 0.5)


        ###~~~~~~~~~~~~~~~~~~~~~###
        ###         TMRCA
        ###~~~~~~~~~~~~~~~~~~~~~###

d    <-    coal[coal$cons.sort == 'D',]
nd    <-    dim(d)[1]
eda    <-    coal[coal$locus == 30517,]

div        <-    c(d$rs.scaled,d$threespine.scaled,d$fw.scaled)
pop        <-    c(rep("Marine",nd), rep("Total",nd), rep("Freshwater",nd))
pop        <-    factor(pop, levels = c("Marine","Total","Freshwater"))
d.cat    <-    data.frame(pop,div)

ylim = c(0,15)
xlim = c(0.5,3.5)

plot(0,0, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab='',bty='n')
    abline(h = median(coal$threespine.scaled, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'gray50')
    abline(h = median(coal$rs.scaled, na.rm=T), lwd = 3, lty = 'dashed', col = 'red3')
    abline(h = median(coal$fw.scaled, na.rm=T), lwd = 3, lty = 'dashed', col = 'blue3')
    segments(x0 = 1, x1 = 2, y0 = d$rs.scaled, y1 = d$threespine.scaled,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 2, x1 = 3, y0 = d$threespine.scaled, y1 = d$fw.scaled,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 1, x1 = 2, y0 = eda$rs.scaled, y1 = eda$threespine.scaled,
            col = 'goldenrod3', lwd = 2)
    segments(x0 = 3, x1 = 2, y0 = eda$fw.scaled, y1 = eda$threespine.scaled,
            col = 'goldenrod3', lwd = 2)
    boxplot(div ~ pop, d.cat, border = c('red3','black','blue3'), add = TRUE,
            ylim = ylim, ylab = "pi.rs, dxy, pi.bt",
            pch = 20, cex = 0.5)


# TAKE A LOOK AT TOP TMRCA

d    <-    d[order(d$threespine.scaled),]

tail(d, n = 10)

l    <-    tail(d$locus, n = 10)
for (i in l) {cat(paste0(i,"\n"))}



        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###        HIGHEST 1% FST       ###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

d    <-    coal[coal$fst.rsfw <= 1 & coal$cons.sort != 'N',]
d    <-    d[order(d$fst.rsfw),]
d    <-    tail(d, n = (dim(d)[1] * 0.1))
    
    # OR

nd    <-    dim(d)[1]
eda    <-    coal[coal$locus == 30517,]

d    <-    d[order(d$pi.fw),]

div        <-    c(d$pi.rs,d$dxy.rsfw,d$pi.fw)
pop        <-    c(rep("Marine",nd), rep("Total",nd), rep("Freshwater",nd))
pop        <-    factor(pop, levels = c("Marine","Total","Freshwater"))
d.cat    <-    data.frame(pop,div)

ylim = c(0,0.05)
xlim = c(0.5,3.5)

plot(0,0, xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = '', ylab='',bty='n')
    abline(h = median(coal$dxy.rsfw, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'black')
    abline(h = median(coal$pi, na.rm=T)   , lwd = 3, lty = 'dashed', col = 'gray50')
    abline(h = median(coal$pi.rs, na.rm=T), lwd = 3, lty = 'dashed', col = 'red3')
    abline(h = median(coal$pi.fw, na.rm=T), lwd = 3, lty = 'dashed', col = 'blue3')
    segments(x0 = 1, x1 = 2, y0 = d$pi.rs, y1 = d$dxy.rsfw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 2, x1 = 3, y0 = d$dxy.rsfw, y1 = d$pi.fw,
            col = rgb(0.25,0.25,0.25,0.1))
    segments(x0 = 1, x1 = 2, y0 = eda$pi.rs, y1 = eda$dxy.rsfw,
            col = 'goldenrod3', lwd = 2)
    segments(x0 = 3, x1 = 2, y0 = eda$pi.fw, y1 = eda$dxy.rsfw,
            col = 'goldenrod3', lwd = 2)
    boxplot(div ~ pop, d.cat, border = c('red3','black','blue3'), add = TRUE,
            ylim = ylim, ylab = "pi.rs, dxy, pi.bt",
            pch = 20, cex = 0.5)









###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#                 GENOME-WIDE TMRCA SUMMARIES
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# GRAND MEAN, ETC
summary(coal$threespine.scaled[coal$cons.sort != "N"])

# MIDDLE 90%
t    <-    coal$threespine.scaled[coal$cons.sort != "N"]
t    <-    sort(t)
tl    <-    length(t)
t5    <-    round(tl * 0.05)
t[t5 + 1] ; t[tl - t5]

print("Boot Lake:"); summary(coal$boot.scaled[coal$cons.sort != "N"])
print("Bear Paw Lake:"); summary(coal$bepa.scaled[coal$cons.sort != "N"])
print("Freshwater:"); summary(coal$fw.scaled[coal$cons.sort != "N"])
print("Rabbit Slough:"); summary(coal$rs.scaled[coal$cons.sort != "N"])

# ANOVA, TUKEY TEST FOR TMRCA BY POPULATION

bt    <-    coal$boot.scaled[coal$cons.sort != "N"]
bp    <-    coal$bepa.scaled[coal$cons.sort != "N"]
rs    <-    coal$rs.scaled[coal$cons.sort != "N"]

tmrca    <-    c(log(bt),log(bp),log(rs))
pop        <-    c(rep("BT", length(bt)),rep("BP", length(bp)),rep("RS", length(rs)))

tbyp    <-    data.frame(pop,tmrca)
boxplot(tmrca ~ pop, tbyp)
summary(pop.aov    <-    aov(tmrca ~ pop, tbyp))
TukeyHSD(pop.aov)

# T-TEST FOR FW VS MARINE

fw    <-    coal$fw.scaled[coal$cons.sort != "N"]
t.test(rs,fw)






###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###                         HAPLOTYPE NETWORKS                          ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

library(ape)
library(pegas)

d    <-    d[order(d$threespine.scaled),]
l    <-    tail(d$locus, n = 10)


plotHaplo    <-    function(locus = test.locus) {
    require(ape)
    require(pegas)
    DIR    <-    paste0(HOME,"/1.datasets/alignments/all/")
#    locus <- "757087"
    aln    <-    paste0(DIR,locus,".phylip")
    aln    <-    read.dna(aln, format='sequential')    
    pops    <-    c('rs', 'bt', 'bp', 'suqi')
    label.ind    <-    rownames(aln)
    label.pop    <-    rep('n', length(label.ind))
        for (pop in pops)    {
            label.pop[grep(pop, label.ind)]    <-    pop
        }
    rownames(aln)    <-    label.pop
    
    
    # identify haplotypes present in alignment and compute network
    haplos    <-    haplotype(aln, what = 'labels', decreasing = TRUE)
    rownames(haplos)    <-    seq(1:length(rownames(haplos)))
    network    <-    haploNet(h = haplos)
    
    # Calculate conversion so that haplotype frequencies in network
    #  are described by the area of the circle, not the diameter
    
    areas        <-    attr(network, 'freq')
    diameters    <-    2 * sqrt(areas / pi)
    
    # Create contingency table so that haplotypes are colored by population
    #  (still not sure how this code works) 
    
    assigns    <-    stack(setNames(attr(haplos, "index"), rownames(haplos)))
    assigns <-    as.matrix(stack(setNames(attr(haplos, "index"), rownames(haplos))))
    assigns    <-    data.frame(assigns)
    assigns$ind    <-    as.character(assigns$ind); assigns$ind <- as.numeric(assigns$ind)
        assigns    <-    table(hap = assigns$ind, pop = rownames(aln)[assigns$values])
        assigns    <-    as.data.frame(assigns); assigns$hap <- as.numeric(as.character(assigns$hap))
        assigns1    <-    assigns[assigns$pop == pops[1], ]
        assigns2    <-    assigns[assigns$pop == pops[2], ]
        assigns3    <-    assigns[assigns$pop == pops[3], ]
        assigns4    <-    assigns[assigns$pop == pops[4], ]
        assigns    <-    data.frame(assigns1$Freq,assigns2$Freq,assigns3$Freq,assigns4$Freq)
        names(assigns)    <-    pops
    
    plot(network, main = locus,
        size = diameters, pie = as.matrix(assigns), scale.ratio = 0.5,
        bg = c('red', 'blue3','steelblue','gray50'), fast = TRUE, labels = FALSE
    )
    
}

plotHaplo(locus = l[10])

# NEAR ATP1A1
atp1    <-    "475382_inv1"
atp2    <-    "75297_inv1"

plotHaplo(atp1)
plotHaplo(atp2)

file    <-    paste0(HOME,"network_475382_chr1_26257782.ps")
postscript(file = file, width = 5, height = 5)
plotHaplo(atp1)
dev.off()

file    <-    paste0(HOME,"network_75297_chr1_26312886.ps")
postscript(file = file, width = 5, height = 5)
plotHaplo(atp2)
dev.off()


# NEAR EDA
file    <-    paste0(HOME,"network_30517_chr4_12808683.ps")
postscript(file = file, width = 5, height = 5)
eda1    <-    "30517"
plotHaplo(eda1)
dev.off()




###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###                                                                     ###
###                     PERMUTE KERNEL SMOOTHING TO                     ###
###                        GENERATE INTERVALS OUTSIDE                   ###
###                            CONFIDENCE BANDS                         ###
###                                                                     ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


bw      <-    0.5
conf    <-    0.99
perms   <-    100
xvals   <-    coal$genomic[coal$cons.sort != "N"]
yvals   <-    coal$threespine.scaled[coal$cons.sort != "N"]

conf.int    <-    function(xvals,yvals,conf = 0.95, perms = 1000, bw = 0.25){
    # FUNCTIONS TO SORT VALUES AND EXTRACT VALS AT, e.g., 2.5% AND 97.5%
    # NEED TO DEFINE SEPARATELY TO INSERT INTO 'APPLY' FUNCTION
    lwr    <-    (1-conf) / 2
    upr    <-    1 - lwr
    
    lim1    <-    function (x) {
        limit    <-    lwr
        test    <-    x[!(is.na(x))]
        sorted    <-    test[order(test)]
        n        <-    length(sorted)
        result    <-    round(n * limit)
        result    <-    sorted[result]
        return(result)
    }
    lim2    <-    function (x) {
        limit    <-    upr
        test    <-    x[!(is.na(x))]
        sorted    <-    test[order(test)]
        n        <-    length(sorted)
        result    <-    round(n * limit)
        result    <-    sorted[result]
        return(result)
    }
    
    # KSMOOTH REAL VALUES (& TO ID SIZE OF MATRIX NEEDED)
    est        <-    ksmooth(xvals, yvals, bandwidth = bw, kernel = 'normal')
    slots    <-    length(est$x)
    
    # RANDOMIZE YVALUES AND RE-SMOOTH WITH SAME X VALUES
    randos    <-    matrix(nrow = slots, ncol = perms)
    cat(paste0("Randomizing and re-smoothing y values ",perms," times...\n"))
    for (i in 1:perms) {
        rand    <-    sample(x = yvals, size = length(yvals), replace = FALSE)
        rand.sm    <-    ksmooth(xvals, rand, bandwidth = bw, kernel = 'normal')
        randos[,i]    <-    rand.sm$y
        if (i %% 100 == 0) {cat(paste0("Done with iteration ",i,"\r"))}
    }
    cat(paste0("\nComputing confidence intervals for each window\n"))
    lower        <-    as.vector(apply(X=randos, MARGIN=1,FUN=lim1))
    upper        <-    as.vector(apply(X=randos, MARGIN=1,FUN=lim2))
    limits        <-    cbind(lower,upper)
    estimate    <-    as.vector(est$y)
    Mb            <-    as.vector(est$x)
    band        <-    data.frame(Mb,estimate,limits)
    # names(band)    <-    c("Mb","estimate","lower","upper")
    band$outlier    <-    rep(0,dim(band)[1])
        band$outlier[band$estimate < band$lower]    <-    -1
        band$outlier[band$estimate > band$upper]    <-    1
    # plot(band$Mb, band$lower,ylim = c(0,10), type = 'l', lwd = 2, col = 'gray50')
        # lines(band$Mb, band$upper, lwd = 2, col = 'gray50')
        # lines(band$Mb, band$est, lwd = 1)
        # oob.lo    <-    band[band$outlier == -1,]
        # oob.hi    <-    band[band$outlier == 1,]
        # points(oob.lo$Mb, rep(1,length(oob.lo$Mb)), pch = 3)
        # points(oob.hi$Mb, rep(8,length(oob.hi$Mb)), pch = 3)
    
    # FIND INTERVALS - LOOP THROUGH 'BAND' DATAFRAME TO FIND
    #   CONSECUTIVE ROWS WITH SAME 'OUTLIER' STATUS
    
    intervals    <-    NULL
    start        <-    band$Mb[1]
    end            <-    band$Mb[1]
    current        <-    band$outlier[1]
    cat("Finding genomic intervals for significantly high and low values...\n")
    for (i in 2:dim(band)[1])    {
        stat    <-    band$outlier[i]
        if (stat != current)    {
            end        <-    band$Mb[i - 1]
            result    <-    c(start, end, current)
            intervals    <-    rbind(intervals,result)
            start    <-    band$Mb[i]
            end        <-    band$Mb[i]
            current    <-    stat
        }
    }
    rownames(intervals)    <-    NULL
    intervals    <-    data.frame(intervals)
    names(intervals)    <-    c("Mb.start","Mb.end","signif")
    intervals$signif[intervals$signif == -1]    <-    'lo'
    intervals$signif[intervals$signif ==  0]    <-    'ns'
    intervals$signif[intervals$signif ==  1]    <-    'hi'
    return(intervals)
}


int     <-    conf.int(xvals = coal$genomic[coal$cons.sort != "N"], yvals = coal$threespine.scaled[coal$cons.sort != "N"], perms = 1000, conf = 0.999)

int.noD <-    conf.int(xvals = coal$genomic[coal$cons.sort != "D" & coal$cons.sort != "N" & coal$fst.rsfw <= 0.5], yvals = coal$threespine.scaled[coal$cons.sort != "D" & coal$cons.sort != "N" & coal$fst.rsfw <= 0.5], perms = 1000, conf = 0.999)

write.table(int, paste0(HOME,"1.datasets/confidenceBands_threespineScaledTMRCA_99.9percent.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


int.hi <- int[int$signif == 'hi',]
int.noD.hi <- int.noD[int.noD$signif == 'hi',]

plot(0,0, xlim = c(0,430), ylim = c(0,1), type = 'n')
    rect(xleft = int.hi$Mb.start, xright = int.hi$Mb.end,
         ybottom = 0.6, ytop = 0.9, col = 'gray25')
    rect(xleft = int.noD.hi$Mb.start, xright = int.noD.hi$Mb.end,
         ybottom = 0.1, ytop = 0.4, col = 'gray25')
    chr.breaks()
    text(x = 0,y = 0.75, labels = "w/ D")
    text(x = 0,y = 0.25, labels = "no D")


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###
###         GENOME SCAN FIGURE
###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

# GENOME-WIDE TMRCA

bw = 0.5

    c    <-    coal[coal$cons.sort != "N",]  
    threespine.smooth    <-    ksmooth(c$genomic, c$threespine.scaled,
                                        kernel = 'normal', bandwidth = bw)
    # CONFIDENCE INTERVALS FOR INCREASED TMRCA
    conf.int  <-  read.table(paste0(HOME,"1.datasets/confidenceBands_threespineScaledTMRCA_99.9percent.tsv"), sep = '\t', header =T)
    ci    <-    conf.int[conf.int$signif == 'hi',]
    spans <- sum(ci$Mb.end - ci$Mb.start)
    percent.hi <- spans / max(coal$genomic)

ymin = 0 ; ymax = 11
xlim = c(0,max(c$genomic))
# xlim = c(60, 200)
plot(0,0, type = 'p', col = 'white', xlim = xlim, ylim = c(ymin,ymax),
        xlab = "Reference Genome, Mb", ylab = expression("T"["MRCA"]*", Mya"))
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(c$genomic[c$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
        for(i in 0:10) {
            abline(h = i, lwd = 1, col = 'gray75', lty = "dotted")
        }
    lines(threespine.smooth, lwd = 1, col = "black")


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### DESCRIPTIVE STATS OF OUTLIER REGIONS
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

int <- read.table(paste0(HOME,"1.datasets/confidenceBands_threespineScaledTMRCA_99.9percent.tsv"), header = T, sep = '\t')
    names(int) <- c("genomic.start","genomic.end","signif")
    int$Mb.start <- int$genomic.start
    int$Mb.end   <- int$genomic.end
    int$chr      <- rep(22, length(int$Mb.start))
    int$tmrca    <- rep(0, length(int$Mb.start))

c <- coal[coal$cons.sort != "N",]

for (i in 1:length(int$Mb.start)) {
	new.start <- genomic2Mb(int$genomic.start[i], LGs <- LGs.info)
	chr       <- new.start[[1]]
	new.coord <- new.start[[2]]
	span      <- int$genomic.end[i] - int$genomic.start[i]
	new.end   <- new.coord + span
	int$Mb.start[i]	<- new.coord
	int$Mb.end[i]   <- new.end
	int$chr[i]      <- chr
	int$tmrca[i]    <- mean(c$threespine.scaled[c$chr == chr & c$Mb > new.coord & c$Mb < new.end])
}


int.hi <- int[int$signif == 'hi',]
int.ns <- int[int$signif == 'ns',]

inv1.start  <- convertCoordinate(1, 21940000, 'old2new', newScafs)[[2]] / 1000000
inv1.end    <- convertCoordinate(1, 21490000, 'old2new', newScafs)[[2]] / 1000000

inv11.start <- convertCoordinate(11, 5400000, 'old2new', newScafs)[[2]] / 1000000
inv11.end   <- convertCoordinate(11, 5900000, 'old2new', newScafs)[[2]] / 1000000

inv21.start <- convertCoordinate(21, 5790000, 'old2new', newScafs)[[2]] / 1000000
inv21.end   <- convertCoordinate(21, 7480000, 'old2new', newScafs)[[2]] / 1000000

inv1     <- c[c$chr == 1  & c$Mb >= inv1.start  & c$Mb <= inv1.end,]
inv11    <- c[c$chr == 11 & c$Mb >= inv11.start & c$Mb <= inv11.end,]
inv21    <- c[c$chr == 21 & c$Mb >= inv21.start & c$Mb <= inv21.end,]

# number of blocks of divergent loci in each inversion 
# (stretches of divergent loci separated by topologies A, B, or C)
# inv1:  7 (separated by 'A' topologies: 2)
# inv11: 8 (separated by 'A' topologies: 2)
# inv21: 31 (separated by 'A' topologies: 10)

# PLOT A CHROMOSOME WITH SEGMENTS FOR INTERVALS AND A KSMOOTHED LINE

chromosome <- 4

chr     <- c[c$chr == chromosome,] ; chr <- chr[order(chr$Mb),]
chr.sm  <- ksmooth(chr$Mb, chr$threespine.scaled, bandwidth = 0.5, kernel = 'normal')
hi      <- int[int$signif == 'hi' & int$chr == chromosome,]
ns      <- int[int$signif == 'ns' & int$chr == chromosome,]

plot(0,0, xlim = range(chr.sm$x), ylim = c(0,10))
	segments(x0 = ns$Mb.start, x1 = ns$Mb.end, y0 = ns$tmrca, y1 = ns$tmrca, col = 'gray50', lwd = 2)
	segments(x0 = hi$Mb.start, x1 = hi$Mb.end, y0 = hi$tmrca, y1 = hi$tmrca, col = 'goldenrod', lwd = 2)
	lines(chr.sm)

# look at the RAD loci in the intervals surrounding Eda and the potentially higher ones elsewhere on chr4
# (They're not)

eda.int     <- c[c$chr == 4 & c$Mb >= 12.695363 & c$Mb <= 12.985418,]
higher1.int <- c[c$chr == 4 & c$Mb >= 20.825691 & c$Mb <= 21.985911,]
higher2.int <- c[c$chr == 4 & c$Mb >= 26.345525 & c$Mb <= 26.811371,]

t.test(eda.int$threespine.scaled, higher1.int$threespine.scaled)
t.test(eda.int$threespine.scaled, higher2.int$threespine.scaled)

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### TMRCA - DIVERGENT LOCI VS OTHER CHR
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

percentile    <- 99

top.fw        <-    coal[coal$fst.rsfw <= 1 & coal$cons.sort != "N",]
    top.fw        <-    top.fw[order(top.fw$fst.rsfw),]
    n            <-    dim(top.fw)[1]
    top.fw        <-    top.fw[round(n * (percentile/100)):n,]

rest    <-    coal[!(coal$locus %in% top.fw$locus) & coal$cons.sort != "N",]
    rest    <-    rest[rest$fst.rsbt <= 1 & rest$fst.rsbp <= 1,]
other    <-    coal[coal$dist2d == 50000 & coal$fst.rsbt <= 1 & coal$fst.rsbp <= 1 & coal$cons.sort != "N",]

par(mfrow = c(2,3))
hi.all        <-    top.fw$threespine.scaled
hi.all.log    <-    top.fw$log.total
    hi.all.cat    <-    rep("hi.fst", length(hi.all))
    lo.all        <-    rest$threespine.scaled
    lo.all.log    <-    rest$log.total
    lo.all.cat    <-    rep("lo.fst", length(lo.all))
    tmrca    <-    c(lo.all,hi.all)
    log.t    <-    c(lo.all.log, hi.all.log)
    cat        <-    c(lo.all.cat,hi.all.cat)
    threespine    <-    data.frame(cat,tmrca)
    three.log    <-    data.frame(cat,log.t)
    T.ALL    <-    t.test(tmrca ~ cat, data = three.log)
    boxplot(tmrca ~ cat , data = threespine, pch = 20, ylim = c(0,12),
        ylab = expression("T"["MRCA"]), xlab = "threespine")

hi.rs        <-    top.fw$rs.scaled
    hi.rs.cat    <-    rep("hi.fst", length(hi.rs))
    lo.rs        <-    rest$rs.scaled
    lo.rs.cat    <-    rep("lo.fst", length(lo.rs))
    tmrca    <-    c(lo.rs,hi.rs)
    cat        <-    c(lo.rs.cat,hi.rs.cat)
    threespine    <-    data.frame(cat,tmrca)
    T.RS    <-    t.test(tmrca ~ cat, data = threespine)
    boxplot(tmrca ~ cat , data = threespine, pch = 20, ylim = c(0,12),
        ylab = expression("T"["MRCA"]), xlab = "RS")

hi.bt        <-    top.fw$boot.scaled
    hi.bt.cat    <-    rep("hi.fst", length(hi.bt))
    lo.bt        <-    rest$boot.scaled
    lo.bt.cat    <-    rep("lo.fst", length(lo.bt))
    tmrca    <-    c(lo.bt,hi.bt)
    cat        <-    c(lo.bt.cat,hi.bt.cat)
    threespine    <-    data.frame(cat,tmrca)
    T.BT    <-    t.test(tmrca ~ cat, data = threespine)
    boxplot(tmrca ~ cat , data = threespine, pch = 20, ylim = c(0,12),
        ylab = expression("T"["MRCA"]), xlab = "BT")

hi.bp        <-    top.fw$bepa.scaled
    hi.bp.cat    <-    rep("hi.fst", length(hi.bp))
    lo.bp        <-    rest$bepa.scaled
    lo.bp.cat    <-    rep("lo.fst", length(lo.bp))
    tmrca    <-    c(lo.bp,hi.bp)
    cat        <-    c(lo.bp.cat,hi.bp.cat)
    threespine    <-    data.frame(cat,tmrca)
    T.BP    <-    t.test(tmrca ~ cat, data = threespine)
    boxplot(tmrca ~ cat , data = threespine, pch = 20, ylim = c(0,12),
        ylab = expression("T"["MRCA"]), xlab = "BP")

hi.fw        <-    top.fw$fw.scaled
    hi.fw.cat    <-    rep("hi.fst", length(hi.fw))
    lo.fw        <-    rest$fw.scaled
    lo.fw.cat    <-    rep("lo.fst", length(lo.fw))
    tmrca    <-    c(lo.fw,hi.fw)
    cat        <-    c(lo.fw.cat,hi.fw.cat)
    threespine    <-    data.frame(cat,tmrca)
    T.FW    <-    t.test(tmrca ~ cat, data = threespine)
    boxplot(tmrca ~ cat , data = threespine, pch = 20, ylim = c(0,12),
        ylab = expression("T"["MRCA"]), xlab = "FW")
T.ALL
T.BT
T.BP
T.RS
T.FW        



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### GENOME-WIDE PI
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


bw = 0.25
    
    boot.smooth    <-    ksmooth(coal$genomic, coal$pi.bt, kernel = 'normal', bandwidth = bw)
    bepa.smooth    <-    ksmooth(coal$genomic, coal$pi.bp, kernel = 'normal', bandwidth = bw)
    rs.smooth    <-    ksmooth(coal$genomic, coal$pi.rs, kernel = 'normal', bandwidth = bw)
    threespine.smooth    <-    ksmooth(coal$genomic, coal$pi,
                                        kernel = 'normal', bandwidth = bw)
    freshwater.smooth    <-    ksmooth(coal$genomic, coal$pi.fw, 
                                        kernel = 'normal', bandwidth = bw)

ymin = 0 ; ymax = 0.015
xlim = c(0,max(coal$genomic))
# xlim = c(60, 200)
plot(0,0, type = 'p', col = 'white', xlim = xlim, ylim = c(ymin,ymax),
        xlab = "Reference Genome, Mb", ylab = expression(pi))
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(coal$genomic[coal$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
    lines(threespine.smooth, lwd = 1, col = "black")
    lines(bepa.smooth, lwd = 0.5, col = rgb(0, (92/255), (142/255)))
    lines(boot.smooth, lwd = 0.5, col = rgb(0, (155/255) , (219/255)))
    lines(freshwater.smooth, lwd = 1, col = rgb((28/255), (117/255), (188/255)))
    lines(rs.smooth, lwd = 1, col = 'red')
    # legend(x = 'topright', 
        # col = c(rgb(0, (155/255) , (219/255)), rgb(0, (92/255), (142/255)), pch = 20,
                # col = rgb((28/255), (117/255), (188/255)), 'red', 'black'),
        # legend=c("Boot","Bear","Fresh","Rabbit","All"))

outfile="/Users/thom/Dropbox/Shared_with_Bill/0.Dissertation_Chapters/2.HistoryOfMolecularDivergence/02.FiguresAndTables/01.Drafts/genomescan_pi.ps"
postscript(file = outfile, width = 7.5, height = 3)
dev.off()


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### GENOME-WIDE FST
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


bw = 0.5
cfst    <-    coal[coal$fst.rsbt <= 1 & coal$fst.rsbp <= 1 & coal$fst.bpbt <= 1,]
    btrs.smooth    <-    ksmooth(cfst$genomic, cfst$fst.rsbt, kernel = 'normal', bandwidth = bw)
    bprs.smooth    <-    ksmooth(cfst$genomic, cfst$fst.rsbp, kernel = 'normal', bandwidth = bw)
    btbp.smooth    <-    ksmooth(cfst$genomic, cfst$fst.bpbt, kernel = 'normal', bandwidth = bw)

ymin = 0 ; ymax = 0.8
xlim = c(0,max(coal$genomic))
# xlim = c(60, 200)
plot(0,0, type = 'p', col = 'white', xlim = xlim, ylim = c(ymin,ymax),
        xlab = "Reference Genome, Mb", ylab = expression("F"["ST"]))
        # PLOT CHROMOSOME BANDS, IF GENOMIC
            for (i in 2:length(LGs.info$groups)) {
                xleft    <- LGs.info$LGs.breaks[i - 1] / 1000000
                xright    <- LGs.info$LGs.breaks[i] / 1000000
                xmean    <-    (xright + xleft) / 2
                ybottom    <- ymin-1
                ytop    <- ymax+1
                rect(xleft = xleft, ybottom = ybottom, ytop = ytop, xright = xright,
                    col = ifelse(i %% 2 == 0,"white", "gray90"), border = NA)
                text(x = xmean, y = ymax, labels = LGs.info$groups[i-1], pos = 1)
            }
            start21    <-    LGs.info$LGs.breaks[21] / 1000000
            end21    <-    max(coal$genomic[coal$chr == '21'])
            mid21    <-    (start21 + end21) / 2
            text(x = mid21, y = ymax, labels = "21", pos = 1)
    lines(btrs.smooth, lwd = 1, col = "black")
    lines(bprs.smooth, lwd = 0.5, col = "gray25")
    lines(btbp.smooth, lwd = 0.5, col = rgb(0, (155/255) , (219/255)))











