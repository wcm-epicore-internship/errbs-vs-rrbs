#!Rscript

# purpose
#   an example script for a brief analysis of methylseq datasets
#
# for now, while we are working on it we'll source the code
#   source("./scripts/analysis.R")
#
# input:
#   /path/to/methylcall/data
#
# eventually we may want to use it from the command line
# usage:
#   Rscript --vanilla analysis.R /path/to/methylcall/data

# libraries and includes
library(methylKit)
library(GenomicRanges)

source("scripts/utils.R")

## eventually we may do this
## check command line
# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input directory).n", call.=FALSE)
# }
# dataDir <- args[1]

# for now i'll hardcode the data directory
dataDir <- "../data"
# check that it exists
if (dir.exists(dataDir)) {
  cat("given data directory [", dataDir, "]\n")
} else {
  stop("cannot find data directory [", dataDir, "]\n")
}

# while i'm working on a script i like to save some of the data structures
# that take a long time to create so i dont have to create them every time.
rdataDir <- "rdata"
if (! dir.exists(rdataDir)) {
  dir.create(rdataDir)
  cat("created rdata directory [", rdataDir, "]\n")
}

# maybe the result will just be a pdf,
# but think about saving results to a directory
resDir <- "results"
# or even better, with a date,
# if you want to compare how your tests evolve
# these will probably need to be removed
# resDir <- paste("results", format(Sys.time(), "%Y%m%d"), sep="-")
if (! dir.exists(resDir)) {
  dir.create(resDir)
  cat("created results directory [", resDir, "]\n")
}

# get the data files
file.list <- getDataFiles(dataDir)

cat("collecting data, generating mcgrlist\n")
mcgrlist <- lapply(file.list, getCpgGr)
# while testing you may wish to save myobj so you dont have to recreate it every time
# save(mcgrlist, file="mcgrlist.rda")
rda <- paste(rdataDir, "mcgrlist.rda", sep="/")
cat("saving mcgrlist to [", rda, "]\n")
save(mcgrlist, file=rda)

# now that we have some data
# we can start looking at it
# ...
x <- unlist(lapply(mcgrlist,length))
plotNumberOfCpGs(x)
dev.new()

hist(mcgrlist[[1]]$coverage, col="blue")
dev.new()

plot(density(mcgrlist[[1]]$coverage))
dev.new()

hist(mcgrlist[[3]]$numC / mcgrlist[[3]]$coverage)
dev.new()

pdf(paste(resDir,"plot1.pdf", sep="/"))
qqplot(mcgrlist[[4]]$numC,mcgrlist[[4]]$numT)
dev.off()