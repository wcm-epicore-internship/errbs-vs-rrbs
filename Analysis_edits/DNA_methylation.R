#R packages
library(methylKit)
library(GenomicRanges)
library(rmarkdown)
library(knitr)
library(vioplot)

#set up some variables
analysis.directory  <- "/Users/shazedaomar/Desktop/Internship_notes"
working.directory   <- paste(analysis.directory, "MethylCall_Data_Plots", sep="/")
data.file.directory <- paste(analysis.directory, "MethylCall_Data_Epicore", sep="/")

# (thk2008) for me on my machine
if (Sys.getenv()[["USER"]] == "thk2008") {
  analysis.directory  <- "/scratch001/thk2008/wcm-epicore-internship/errbs-vs-rrbs/Analysis_edits"
  working.directory   <- paste(analysis.directory, "MethylCall_Data_Plots", sep="/")
  data.file.directory <- paste(analysis.directory, "MethylCall_Data_Epicore", sep="/")

  # (thk2008) ive already run this so for testing
  # i want to turn it off
  run_methylation <- FALSE
}

# (thk2008) should check that these exist of throw an error if not
if (file.exists(working.directory)) {
  cat("set working directory to [",working.directory,"]\n")
  setwd(working.directory)
} else {
  stop("Cannot find working directory [",working.directory,"]")
}

#methread for Multiple Files - load myobj.rda that is already created
#if not myobj will be created

# (thk2008) if the data object myobj exists in the R environment, yay, it's already loaded
# otherwise if the saved data object r data file exists in the current working directory, load it
# otherwise create the data object myobj from the raw methylcall data files and save it
if (!exists("myobj")) {
  cat("Did not find myobj\n")
  if (file.exists("myobj.rda")) {
    cat("myobj.rda loading\n")
    load("myobj.rda")
  } else {
    if (file.exists(data.file.directory)) {
      #read in files - path file stored in variable data.file.path
      # (thk2008) path to methylcall files from alignment pipeline are stored in variable data.file.directory
      file.names.list <- list.files(path = data.file.directory, pattern = ".mincov0.txt.gz", full.names = TRUE)
    } else {
      stop("Cannot find data file directory [",data.file.directory,"]")
    }
    cat("creating and saving myobj\n")
    myobj <- methRead(as.list(file.names.list),
      sample.id=as.list(gsub(".mincov0.txt.gz", "", basename(file.names.list))), #gsub finds all the files that end in ("*.mincov0.txt.gz"),
      # Basename removes the ending of the file. -- are you sure? see ?basename
      assembly="hg19",
      # treatment=rep(0,6), #rep function - 0 = same -- (thk2008) what does 0 mean?.  6 = number of items in the file.names.list -- (thk2008) then use a variable or the length of the vector
      treatment=rep(0,length(file.names.list)),
      context="CpG",
    )
    save(myobj,file = "myobj.rda")
    cat ("myobj saved \n")
  }
}


#Methylation function:Counts the number of rows, descriptive statistics, percent methylation distribution, histogram with the percentage of CpG for each dataset.
methylation <- function(myobj) {
  for (i in 1:length(myobj)) {
    sample.output <- myobj[[i]]@sample.id
    cat("sample", sample.output, "\n", sep=": ")
    cat("number of sites", nrow(myobj[[i]]), "\n", sep=": ")
    pdf(paste(sample.output, "methylation_output.pdf", sep="_"))
    getMethylationStats(myobj[[i]])
    getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE)
    getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE)
    dev.off() #pdf saved in current directory
  }
}

#to execute function
if (run_methylation) {
  methylation(myobj)
}

#Output number of rows/ CpGs for each sample
samples <- sapply(myobj, nrow)
names(samples) <- getSampleID(myobj)
# (thk2008) this appears to be unnecessary. where is this used?
# file.list <- as.list(samples)
# file.list

#filtering the data by read coverage to avoid bias (obtain coverage that is not too high or too low)
if (!exists("filter_myobj")) {
  cat("Did not find filter_myobj\n")
  if (file.exists("filter_myobj.rda")) {
    cat("filter_myobj.rda loading\n")
    load("filter_myobj.rda")
  } else {
    cat("creating and saving ")
    filter_myobj = filterByCoverage(myobj, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
    save(filter_myobj, file = "filter_myobj.rda")
    cat ("filter_myobj.rda saved \n")
  }
}

#Merging sample to find bases covered in samples
if (!exists("meth")) {
  cat("Did not find meth obj\n")
  if(file.exists("meth.rda")) {
    cat("meth.rda loading\n")
    load("meth.rda")
  } else {
    cat("creating and saving ")
    meth=unite(myobj, destrand=FALSE)
    save(meth, file = "meth.rda")
    cat ("meth.rda saved \n")
  }
}

# (thk2008) this is already done above
#load meth.rda and filter_myobj.rda that is already created
# load("/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Plots/meth.rda")
# load("/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Plots/filter_myobj.rda")

#Using perc.meth from the Methylation package to compare data
perc.meth=percMethylation(meth)

# (thk2008) this appears to be extra
# hist(perc.meth)
dens <- apply(perc.meth, 2, density)

# (thk2008) this plot needs labeling. title, xaxis, yaxis
# and eventually a discussion of what it is and means
plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
mapply(lines, dens, col=1:length(dens))
legend("topright", legend=names(dens), fill=1:length(dens))

#for loop and function to obtain correlation and Cluster Samples
# (thk2008) What are you trying to do here?
# where is 'i' used to iterate?
# each of these functions plots all of the data
# it looks like each function is run for the length of meth (so 22 times)
# so it's creating 22 correlation graphs,
# and 22 cluster plots, and 22 PCA variance plots and 22 PCA plots
# what is the purpose of break here?
Analysis <- function(meth) {

  for (i in 1:length(meth)) {
    getCorrelation(meth,plot=TRUE)
    clusterSamples(meth, dist="correlation", plot=TRUE) #find similarity in method
    #Principal component analysis - brings out variation and strong pattern in sample (helps to visualize/explore)
    PCASamples(meth, screeplot=TRUE)
    PCASamples(meth)

    break
  }
}


#To execute the function
# (thk2008) do not run this until it is fixed
if (F) {
  Analysis(meth)
}

#plots and graphs - 1) barplot of all samples
barplot(samples,
  main="Analysis",
  xlab="Methylation Method",
  ylab="Frequency",
  col= "blue")
dev.new()

if(F) {
# (thk2008) does this work? see ?pdf
pdf(paste(barplot(samples), "barplot_samples.pdf", sep = "_"))
#dev.off()
}

# 2) plot the samples
# (thk2008) does this work? fix the labels.
# when i run it i get an error.
# the x-axis uses the names from a data object that is different from
# the one being plotted. this could lead to trouble if the names are different.
# and it doesn't exist in the environment yet
plot(samples,
  main="Analysis",
  main = "",
  xlab = "" ,
  ylab = "frequency",
  type = "o",  #plot type
  col = "red") #plot color
axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2) #changes the labels 1,2,3,4,5 with sample.id names (works with this line of code)
dev.new()

if(F) {
# (thk2008) save as above, does this work?
pdf(paste(plot(samples), "plot_samples.pdf", sep = "_"))
#dev.off()
}


# !!!
# (thk2008) there is a similar data object above: filter_myobj
# actually, from the very beginning, myobj <- methRead(...)
# filters for >= 10x coverage by default so we dont need this
# !!!


# plot the data in that object instead of this below
# if (F) {
# # (thk2008) this comment does match what the code is doing.
# # what could it say?
# #Store and output the sample.id
# coverage_vals <- lapply(myobj, function(i){
#   x <- getData(i)$coverage
# })
# names(coverage_vals) <- getSampleID(myobj)
# names(coverage_vals)

# # (thk2008) this comment does match what the code is doing. the boxplot is below it
# # what could it say?
# #boxplot of samples
# filt_coverage_vals <- lapply(coverage_vals, function(i){
#   x <- i[i >= 10 & i <= 500]
# })
# boxplot(filt_coverage_vals)
# }

filt_coverage_vals <- lapply(filter_myobj, function(i){
  x <- getData(i)$coverage
})
names(filt_coverage_vals) <- getSampleID(filter_myobj)
boxplot(filt_coverage_vals)


#function for graphs
# (thk2008) the plots get overwritten
# maybe consider dev.new() between plots or write them to a file
plots <- function(filt_coverage_vals) {


  orig_margins <- par()$mar  # (thk2008) this needs to be done only once
  cl <- rainbow(6)           # (thk2008) save for this, it can be used for all the plots

  for (i in 1:length(filt_coverage_vals))  {

    par(mar=c(7.1,4.1,4.1,2.1)) # (thk2008) changing the margins to see sample names
    # (thk2008) need to correct the labels
    plot(0, type="n", xlim=c(0,7), ylim=c(0,500), axes=F, main="Methylation analysis", xlab="Method", ylab="Coverage")
    for (i in 1:length(filt_coverage_vals)) {
      vioplot(filt_coverage_vals[[i]], h=10, at=i, add=T, col=cl[i])
    }
    axis(2)
    axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2)

    plot(0, type="n", xlim=c(0,7), ylim=c(1,3), axes=F, main="all graphs must have labels", xlab="",ylab="")
    for (i in 1:length(filt_coverage_vals)) {
      vioplot(log10(filt_coverage_vals[[i]]), h=0.1, at=i, add=T, col=cl[i])
    }


    par(mfrow=c(3,2))
    for (i in 1:length(filt_coverage_vals)) {
      plot(density(filt_coverage_vals[[i]], bw=0.1), col=cl[i])
    }

    par(mfrow=c(1,1))
    # (thk2008) careful, plotting filt_coverage_vals but legend has names from a different data object
    plot(0, type="n", xlim=c(0,500), ylim=c(0,0.21), axes=F, main="", xlab="",ylab="")
    title(main = "Methylation Coverage", sub = "Coverage (0-500)", xlab = NULL, ylab = "C methylated",
    line = NA, outer = FALSE, legend("topright", legend=names(coverage_vals), fill=cl))
    for (i in 1:length(filt_coverage_vals)) {
      lines(density(filt_coverage_vals[[i]],bw=0.1), col=cl[i])
    }
    axis(1)
    axis(2)

    # (thk2008) why the switch to coverage_vals from filt_coverage_vals?
    par(mfrow=c(3,2))
    for (i in 1:length(coverage_vals)){
      plot(density(coverage_vals[[i]], bw=0.1), col=cl[i])
    }
    axis(1) # (thk2008) does this do something?
    axis(2) # (thk2008) does this do something?


    par(mfrow=c(1,1))
    # (thk2008) careful, plotting filt_coverage_vals but legend has names from a different data object
    # also be careful cutting and pasting code. it will perpetuate the discrepancy or worse an error.
    plot(0, type="n", xlim=c(0,500), ylim=c(0,3.5), axes=F, main="", xlab="",ylab="")
    title(main = "Methylation Coverage", sub = "Coverage (0-500)", xlab = NULL, ylab = "C methylated",
    line = NA, outer = FALSE, legend("topright", legend=names(coverage_vals), fill=cl))

    for (i in 1:length(coverage_vals)) {
      lines(density(coverage_vals[[i]],bw=0.1), col=cl[i])
    }
    axis(2)
    axis(1)

    break # (thk2008) what is the purpose break? does it do something?
  }
}

#execute the function
# (thk2008) does this work? i get an error when i try to run it.
plots(filt_coverage_vals)


#find the bandwidth for each sample - highest bw = 3.63, lowest bw = 0.1786
plot(density(filt_coverage_vals[[i]])) #N= 12646218, Bandwidth=0.1786
plot(density(filt_coverage_vals[[1]])) #N= 6693861,  Bandwidth=2.488
plot(density(filt_coverage_vals[[2]])) #N= 3819823,  Bandwidth=3.63
plot(density(filt_coverage_vals[[3]])) #N= 6001227,  Bandwidth=1.806
plot(density(filt_coverage_vals[[4]])) #N= 4028562,  Bandwidth=1.379
plot(density(filt_coverage_vals[[5]])) #N= 5336906,  Bandwidth=0.5154
plot(density(filt_coverage_vals[[6]])) #N= 12646218, Bandwidth=0.1786


#loop for plot(density(filter_coverage_values)) showing bw measures

for (i in 1:length(filt_coverage_vals))  {
# (thk2008) this was already done above
#   filt_coverage_vals <- lapply(coverage_vals, function(i){
#     x <- i[i >= 10 & i <= 500]
#   })
  plot(density(filt_coverage_vals[[i]]))
}
# (thk2008)
# also you can extract the bw from the return value from density()
# for example:
# bandwidths <- sapply(filt_coverage_vals, function(i) { density(i)$bw })

x <- myobj[[1]]$numCs / myobj[[1]]$coverage * 100
hist(x,
  main="Number of Cs",
  xlab="% methylation per base",
  ylab="Frequency",
  col= "red")

# (thk2008) do we really want to plot the methylation value for each CpG site?
# a histogram is sufficient: hist(x)
# also see ?smoothScatter
# smoothScatter(x) # and tinker with the parameters

x <- myobj[[i]]$numCs / myobj[[i]]$coverage *100
plot(x,
main="Number of Cs",
xlab="% methylation per base",
ylab="Frequency",
col= "yellow")

#library(BSgenome.Hsapiens.UCSC.hg19)
#Using genomicRanges
gr=GRanges(seqnames=c("Chr1"),
ranges=IRanges(start=c(50,150,200), end = c(100,200,300)), #ranges contain the start and end position of the genome
strand=c("+", "-","-"))

gr[1,]
seqnames(gr)
mcols(gr)

#need help
cgr <- GRanges(seqnames=myobj$chr[,2],
  ranges=IRanges(myobj$start[,3],end=myobj$end[,3]),
  strand=myobj$strand[,4],
  coverage=myobj$coverage[,5],
  numC=myobj$numCs[,6],
  numT=myobj$numTs[,7])

chr.len    <- seqlengths(Hsapiens)
seqinfo(cgr) <- seqinfo(Hsapiens)[names(seqlengths(cgr))]
cgr
