
# functions for analysis.R
# source("./scripts/utils.R")

# could put libraries and includes here
library(methylKit)
# getData
# getMethylationStats
# getSampleID
# extract

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

# get the list of files in the data directory and names
# input:
#   datadir: path to directory containing methylcall files
# returns:
#   named character vector of file paths
getDataFiles <- function(datadir) {
  file.list <- list.files(datadir, full=T)
  names(file.list) <- sub(".mincov0.txt.gz","", unlist(lapply(strsplit(file.list, "/"), "[[", 3)))
  file.list <- as.list(file.list)
}

# reads in methylcall and returns a granges object
# input:
#   single list item with file path
# returns:
#   granges object of methylcall file
getCpgGr <- function(mcfile) {
  # library(GenomicRanges)
  # library(BSgenome.Hsapiens.UCSC.hg19)
  # see if there is a header
  cpg.frame <- read.table(mcfile[[1]], nrows=1, as.is=T)
  header <- FALSE
  if ('chrBase' %in% cpg.frame) { header <- TRUE }
  # read file
  # cpg.frame <- read.table(mcfile[[1]], header=header, nrows=sample(1000:10000,1))
  cpg.frame <- read.table(mcfile[[1]], header=header)
  tmp <- rep("+", nrow(cpg.frame))
  tmp[cpg.frame[,4] == "R"] <- "-"
  cpg.frame[,4] <- as.factor(tmp)
  cpg.frame[,6] <- round(cpg.frame[,5]*cpg.frame[,6]/100) # numCs
  cpg.frame[,7] <- round(cpg.frame[,5]*cpg.frame[,7]/100) # numTs
  mcgr   <- GRanges(seqnames=cpg.frame[,2],
                      ranges=IRanges(cpg.frame[,3],end=cpg.frame[,3]),
                      strand=cpg.frame[,4],
                      coverage=cpg.frame[,5],
                      numC=cpg.frame[,6],
                      numT=cpg.frame[,7])
  chr.len       <- seqlengths(Hsapiens)
  seqinfo(mcgr) <- seqinfo(Hsapiens)[names(seqlengths(mcgr))]
  mcgr
}

# plot the number of CpGs per sample as barplot
# input:
#   named numeric vector
# returns:
#   null
# output:
#   barplot
plotNumberOfCpGs <- function(x) {
  barplot(x)
}


if(F) {
# some stuff for dealing with myobj
coverage_vals <- lapply(myobj, function(i){ x <- getData(i)$coverage })
names(coverage_vals) <- getSampleID(myobj)

# maybe want to filter for coverage values that make sense
filt_coverage_vals <- lapply(coverage_vals, function(i){
  x <- i[i >= 10 & i <= 500]
})
boxplot(filt_coverage_vals)

library(vioplot)
orig_margins <- par()$mar
par(mar=c(7.1,4.1,4.1,2.1))
plot(0, type="n", xlim=c(0,7), ylim=c(0,500), axes=F, main="all graphs must have labels", xlab="",ylab="")
for (i in 1:length(filt_coverage_vals)){ vioplot(filt_coverage_vals[[i]], h=10, at=i, add=T) }
axis(2)
axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2)

dev.new()
plot(0, type="n", xlim=c(0,7), ylim=c(1,3), axes=F, main="all graphs must have labels", xlab="",ylab="")
for (i in 1:length(filt_coverage_vals)){ vioplot(log10(filt_coverage_vals[[i]]), h=0.1, at=i, add=T) }

par(mfrow=c(3,2))
for (i in 1:length(filt_coverage_vals)){ plot(density(filt_coverage_vals[[i]], bw=0.1)) }

par(mfrow=c(1,1))
plot(0, type="n", xlim=c(0,500), ylim=c(0,0.21), axes=F, main="all graphs must have labels", xlab="",ylab="")
for (i in 1:length(filt_coverage_vals)){ lines(density(filt_coverage_vals[[i]], bw=0.1)) }
axis(1)
axis(2)


par(mfrow=c(3,2))
for (i in 1:length(coverage_vals)){ plot(density(coverage_vals[[i]], bw=0.1)) }
axis(1)
axis(2)

par(mfrow=c(1,1))
plot(0, type="n", xlim=c(0,500), ylim=c(0,3.5), axes=F, main="all graphs must have labels", xlab="",ylab="")
for (i in 1:length(coverage_vals)){ lines(density(coverage_vals[[i]], bw=0.1)) }
axis(2)
axis(1)
}

