#R packages
library(methylKit)
library(GenomicRanges)
library(knitr)
library(vioplot)

#read in and list files
analysis.directory <- "/Users/shazedaomar/Desktop/Internship_notes"
working.directory <- paste(analysis.directory, "MethylCall_Data_Plots", sep="/")
data.file.directory <- paste(analysis.directory, "MethylCall_Data_Epicore", sep="/")
list_of_files <- list.files(path = data.file.directory , pattern = "*.mincov0.txt.gz", full.names = TRUE)
list_of_files


#check working directory
if (file.exists(working.directory)) {
    cat("set working directory to [", working.directory, "]\n")
    setwd(working.directory)
} else{
    stop("cannot find working directory [", working.directory,"]")
}

#methread for Multiple Files - load myobj.rda that is already created.
#Otherwise, create the data object 'myobj' from the raw methylcall data files and save it.
load("myobj.rda")


#If myobj is not created, create and save the data object myobj from the raw methylcall data files
if  (!exists("myobj")) {
    cat("Did not find myobj\n")
    if(file.exists("myobj.rda")) {
        cat("myobj.rda loading\n")
        load("myobj.rda")
    } else {
        cat("creating and saving myobj\n")
        myobj=methRead(as.list(file.names.list),
        sample.id=as.list(gsub(".mincov0.txt.gz", "", basename(file.names.list))), #gsub finds all the files that end in ("*.mincov0.txt.gz"), Basename removes the ending of the file.
        assembly="hg19",
        treatment=rep(0,6), #rep function - 0 = same , 6 = number of items in the file.names.list
        context="CpG",
        dbdir = "methylDB"
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
        pdf(paste(sample.output,"methylation_output.pdf",sep="_"))
        getMethylationStats(myobj[[i]])
        getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
        getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
        dev.off() #pdf saved in current directory
    }
}

#to execute function
methylation(myobj)


#filtering the data by read coverage to avoid bias (obtain coverage that is not too high or too low)
if  (!exists("filter_myobj")) {
    cat("Did not find filter_myobj\n")
    if(file.exists("filter_myobj.rda")) {
        cat("filter_myobj.rda loading\n")
        load("filter_myobj.rda")
    } else {
        cat("creating and saving ")
        filter_myobj = filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
        hi.count=NULL,hi.perc=99.9)
        save(filter_myobj,file = "filter_myobj.rda")
        cat ("filter_myobj.rda saved \n")
    }
}

#Merging sample to find bases covered in samples
if  (!exists("meth")) {
    cat("Did not find meth.myobj\n")
    if(file.exists("meth.rda")) {
        cat("meth.rda loading\n")
        load("meth.rda")
    } else {
        cat("creating and saving ")
        meth=unite(myobj, destrand=FALSE)
        save(meth,file = "meth.rda")
        cat ("meth.rda saved \n")
    }
}


#load meth.rda and filter_myobj.rda that is already created
load("/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Plots/meth.rda")
load("/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Plots/filter_myobj.rda")


#Using perc.meth from the Methylation package to compare data
perc.meth=percMethylation(meth)
dens <- apply(perc.meth, 2, density)
plot(NA, xlab = "Number of CpG Unit ", ylab= "Percent (%)", main= "Distribution of Methylation ",
xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
mapply(lines, dens, col=1:length(dens))
legend("topright", legend=names(dens), fill=1:length(dens))


#for loop and function to obtain correlation and Cluster Samples
#creates one correlation chart, cluster of sample chart and PCA sample chart.
#break ends the for loop
Analysis <- function(meth) {
    
    getCorrelation(meth,plot=TRUE)
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE) #find similarity in method
    #Principal component analysis - brings out variation and strong pattern in sample (helps to visualize/explore)
    PCASamples(meth, screeplot=TRUE)
    PCASamples(meth) #is this important
    
}

#To execute the function
Analysis(meth)


#Store and output the sample.id and nrow for each sample (CpGs)
coverage_vals <- lapply(myobj, function(i){
    x <- getData(i)$coverage #?
})
names(coverage_vals) <- getSampleID(myobj) #extra names from sample myobj
file.list <- (sapply(myobj,nrow)) #number of samples
names(file.list) <- names(coverage_vals)
file.list


#plots and graphs - barplot of all samples
barplot(file.list,
        main="Analysis",
        xlab="Methylation Method",
        ylab="Frequency",
        col= "blue")


#filter the coverage 
filt_coverage_vals <- lapply(coverage_vals, function(i){
    x <- i[i >= 10 & i <= 500]
})


#graphs
plots <- function(filt_coverage_vals) {
    
    orig_margins <- par()$mar  # (thk2008) this needs to be done only once
    cl <- rainbow(6)           # (thk2008) save for this, it can be used for all the plots
    
    
    #vioplot graph
    par(mar=c(7.1,4.1,4.1,2.1)) # changing the margins to see sample names
    # need to correct the labels
    plot(0, type="n", xlim=c(0,7), ylim=c(0,500), axes=F, main="Methylation analysis", xlab="", ylab="Coverage")
    for (i in 1:length(filt_coverage_vals)) {
        vioplot(filt_coverage_vals[[i]], h=10, at=i, add=T, col=cl[i])
    }
    axis(2)
    axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2)
    
    
    #Bandwith graphs
    par(mfrow=c(3,2))
    for (i in 1:length(filt_coverage_vals)) {
        plot(density
        (filt_coverage_vals[[i]], bw=0.1), col=cl[i])
    }
    
    
    #methylation coverage graph
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
    
    #graph
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
    
}


#execute the function
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
    filt_coverage_vals <- lapply(coverage_vals, function(i){
        x <- i[i >= 10 & i <= 500]
    })
    plot(density(filt_coverage_vals[[i]]))
    
}


x <- myobj[[1]]$numCs / myobj[[1]]$coverage *100
hist(x,
main="Number of Cs",
xlab="% methylation per base",
ylab="Frequency",
col= "red")

x <- myobj[[i]]$numCs / myobj[[i]]$coverage *100
plot(x,
main="Number of Cs",
xlab="% methylation per base",
ylab="Frequency",
col= "yellow")

library(BSgenome.Hsapiens.UCSC.hg19)
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

chr.len      <- seqlengths(Hsapiens)
seqinfo(cgr) <- seqinfo(Hsapiens)[names(seqlengths(cgr))]
cgr
