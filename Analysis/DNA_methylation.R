#R packages
library(methylKit)
library(GenomicRanges)
library(knitr)
library(vioplot)

#read in files - path file stored in variable data.file.path
data.file.path <- "/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Epicore"
data <- data.file.path
file.names.list <- list.files(path = data, pattern = "*.mincov0.txt.gz", full.names = TRUE)
file.names.list

#set working directory
setwd("~/Desktop/Internship_notes/MethylCall_Data_Plots")
getwd()

#methread for Multiple Files - load myobj.rda that is already created
load("/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Plots/myobj.rda")

#if not myobj will be created
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


#Output number of rows for each sample
samples <- sapply(myobj,nrow)
names(samples) <- getSampleID(myobj)
file.list <- as.list(samples)
file.list

#filtering the data by read coverage to avoid bias (obtain coverage that is not too high or too low)
filter_myobj = filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
hi.count=NULL,hi.perc=99.9)
head(filter_myobj, n=20)

#Merging sample to find bases covered in samples
meth=unite(myobj, destrand=FALSE)
head(meth) #view the file meth

#save meth file so that it doesn't have to be created
if  (!exists("meth.rda")) {
    cat("Did not find meth.myobj\n")
    if(file.exists("meth.rda")) {
        cat("meth.rda loading\n")
        load("meth.rda")
    } else {
        cat("creating and saving meth.rda\n")
        meth.rda = filter_myobj = filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
        hi.count=NULL,hi.perc=99.9)
        meth=unite(myobj, destrand=FALSE)
        save(myobj,file = "meth.rda")
        cat ("meth.rda saved \n")
    }
}

#load meth.rdathat is already created
load("/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Plots/meth.rda")

#creates a methylBase object where only CpGs covered with at least 1 sample per group will be returned
meth.min=unite(myobj,min.per.group=1L)
head(meth.min)


#getMethylationStats doesn't give enough information to plot nor extract sample
#meth and percentMethylationStat extra data
x <- meth #rda save
perc.meth=percMethylation(x)
hist(perc.meth)

dens <- apply(perc.meth, 2, density)
plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
mapply(lines, dens, col=1:length(dens))
legend("topright", legend=names(dens), fill=1:length(dens))

#for loop and function to obtain correlation and Cluster Samples
Analysis <- function(meth) {
    
    for (i in 1:length(meth)) {
        getCorrelation(meth,plot=TRUE)
        clusterSamples(meth, dist="correlation", method="ward", plot=TRUE) #find similarity in method
        #Principal component analysis - brings out variation and strong pattern in sample (helps to visualize/explore)
        PCASamples(meth, screeplot=TRUE)
        PCASamples(meth) #is this important
        
    break
    }
}


#To execute the function
Analysis(meth)


#plots and graphs - 1) barplot of all samples
barplot(samples,
main="Analysis",
xlab="Methylation Method",
ylab="Frequency",
col= "blue")
dev.new()
pdf(paste(barplot(samples), "barplot_samples.pdf", sep = "_"))
#dev.off()

# 2) plot the samples
plot(samples,
main="Analysis",
main = "",
xlab = "" ,
ylab = "frequency",
type = "o", #plot type
col = "red") #plot color
axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2) #changes the labels 1,2,3,4,5 with sample.id names (works with this line of code)
dev.new()
pdf(paste(plot(samples), "plot_samples.pdf", sep = "_"))
#dev.off()

#Store and output the sample.id
coverage_vals <- lapply(myobj, function(i){
  x <- getData(i)$coverage
})
names(coverage_vals) <- getSampleID(myobj)
names(coverage_vals)

#boxplot of samples
filt_coverage_vals <- lapply(coverage_vals, function(i){
    x <- i[i >= 10 & i <= 500]
})
boxplot(filt_coverage_vals)

#function for graphs
plots <- function(filt_coverage_vals) {
    
    for (i in 1:length(filt_coverage_vals))  {
        
        orig_margins <- par()$mar
        cl <- rainbow(6)
        par(mar=c(7.1,4.1,4.1,2.1))
        plot(0, type="n", xlim=c(0,7), ylim=c(0,500), axes=F, main="Methylation analysis", xlab="Method",ylab="Coverage")
        for (i in 1:length(filt_coverage_vals)){ vioplot(filt_coverage_vals[[i]], h=10, at=i, add=T, col=cl[i]) }
        axis(2)
        axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2)
        
        plot(0, type="n", xlim=c(0,7), ylim=c(1,3), axes=F, main="all graphs must have labels", xlab="",ylab="")
        for (i in 1:length(filt_coverage_vals)){ vioplot(log10(filt_coverage_vals[[i]]), h=0.1, at=i, add=T, col=cl[i]) }
        
        
        par(mfrow=c(3,2))
        cl <- rainbow(6)
        for (i in 1:length(filt_coverage_vals)){ plot(density(filt_coverage_vals[[i]], bw=0.1), col=cl[i]) }
        
        par(mfrow=c(1,1))
        cl <- rainbow(6)
        plot(0, type="n", xlim=c(0,500), ylim=c(0,0.21), axes=F, main="", xlab="",ylab="")
        title(main = "Methylation Coverage", sub = "Coverage (0-500)", xlab = NULL, ylab = "C methylated",
        line = NA, outer = FALSE, legend("topright", legend=names(coverage_vals), fill=cl))
        for (i in 1:length(filt_coverage_vals)){ lines(density(filt_coverage_vals[[i]],bw=0.1), col=cl[i]) }
        axis(1)
        axis(2)
        
        
        par(mfrow=c(3,2))
        cl <- rainbow(6)
        for (i in 1:length(coverage_vals)){ plot(density(coverage_vals[[i]], bw=0.1), col=cl[i]) }
        axis(1)
        axis(2)
        
        
        par(mfrow=c(1,1))
        cl <- rainbow(6)
        plot(0, type="n", xlim=c(0,500), ylim=c(0,3.5), axes=F, main="", xlab="",ylab="")
        title(main = "Methylation Coverage", sub = "Coverage (0-500)", xlab = NULL, ylab = "C methylated",
        line = NA, outer = FALSE, legend("topright", legend=names(coverage_vals), fill=cl))
        for (i in 1:length(coverage_vals)) {lines(density(coverage_vals[[i]],bw=0.1), col=cl[i]) }
        axis(2)
        axis(1)
        
        
        break
    }}

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

chr.len      <- seqlengths(Hsapiens)
seqinfo(cgr) <- seqinfo(Hsapiens)[names(seqlengths(cgr))]
cgr
