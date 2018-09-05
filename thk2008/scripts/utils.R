
# functions for analysis.R
# source("./scripts/utils.R")

# could put libraries and includes here
library(methylKit)
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
  cpg.frame <- read.table(mcfile[[1]], header=header, nrows=sample(1000:10000,1))
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
