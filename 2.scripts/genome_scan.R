#####################################
#
# SOME USEFUL FUNCTIONS AND RESOURCES
#
#####################################

# FUNCTION FROM GLAZER et al. FOR CONVERTING BETWEEN GENOME ASSEMBLIES

source("/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/4.data/2.scripts/convertCoordinate.R")
newScafs    <-    "/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/4.data/NewScaffoldOrder.csv"
convertCoordinate(1, 21730178, 'old2new', newScafs)

# TABLE OF COORDINATES/LENGTHS FOR CHROMOSOMES

# ORIGINAL ASSEMBLY
LGs.info    <-    read.table("/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/4.data/GacLGs.tsv", sep = "\t", 
                                header = TRUE, stringsAsFactors = FALSE)
    group    <-    2    # ***column in which to find 'group' label (You specify)
    bp        <-    3    # ***column in which to find chromosome-specific bp position (You specify)

# GLAZER ET AL ASSEMBLY
LGs.info    <-    read.table("/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/4.data/GacNewLGs.tsv", sep = "\t", 
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

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###    LOAD DATA
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

cat("Loading data from file...\n")    
    coal <- read.table(
            "/Users/thom/Dropbox/Shared_with_Bill/9.Manuscripts/1.AncientAdaptation/4.data/1.datasets/coalescence.tsv", 
            sep = "\t", header = TRUE, stringsAsFactors=FALSE)
        coal$chr    <-    as.character(coal$chr)
        # coal        <-    coal[coal$chr != "19",]
    coal            <-    coal[order(coal$chr, coal$Mb),]


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
### PLOT WHOLE GENOME
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

xvals <- coal$genomic[coal$fst.rsbt <= 1]
yvals <- coal$fst.rsbt[coal$fst.rsbt <= 1]
bw    <- 0.5 #bandwidth for smoothing (Mb)

xlim <- c(0,450)
ylim <- c(0,1)

sm.vals    <-    ksmooth(x = xvals[order(xvals)], y = yvals[order(xvals)], bandwidth = bw, kernel='normal')

ymin    <-    0
ymax    <-    1
plot(0,0, type = 'n', xlim = xlim, ylim = ylim) # create blank plotting region
        # PLOT CHROMOSOME BOXES
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
    # PLOT DATA AS POINTS AND LINES
        points(xvals, yvals, col = rgb(0.5,0.5,0.5, 0.1), pch = 20)
        lines(sm.vals$x, sm.vals$y, col = 'blue3')






