#R packages
library(methylKit)
library(GenomicRanges)
install.packages("vioplot")
library(vioplot)

#read in files
#path file in variable data.file.path
data.file.path <- "/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Epicore"
data <- data.file.path
file.names.list <- list.files(path = data, pattern = "*.mincov0.txt.gz", full.names = TRUE)
file.names.list

#set working directory
setwd("~/Desktop/Internship_notes/MethylCall_Data_Plots")
getwd()


#Explain
#rda <- paste(rdataDir, "mcgrlist.rda", sep="/")
#cat("saving mcgrlist to [", rda, "]\n")
#save(mcgrlist, file=rda)
#but then you will have to load it

#if myobj.rda is saved it gives this error
#Did not find myobj Error in readChar(con, 5L, useBytes = TRUE) : cannot open the connection
#In addition: Warning message:
#In readChar(con, 5L, useBytes = TRUE) : cannot open compressed file 'myobj.rda ',probable reason 'No such file or directory'
#does not let me keep on working

if  (!exists("myobj")) {
    cat("Did not find myobj\n")
    if(file.exists("myobj.rda")) {
        cat("myobj.rda loaded")
        load("myobj.rda\n")
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

#filter the data by read coverage to avoid bias (obtain coverage that is not too high or too low)
filter_myobj = filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
hi.count=NULL,hi.perc=99.9)
head(filter_myobj, n=20)

#Merging sample to find bases covered in samples 
meth=unite(myobj, destrand=FALSE)
head(meth) #view the file meth

#creates a methylBase object where only CpGs covered with at least 1 sample per group will be returned
meth.min=unite(myobj,min.per.group=1L)
head(meth.min)

#find the correlation of all samples
getCorrelation(meth,plot=TRUE)

#Cluster the samples to find similarity in methods
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

#Principal component analysis - brings out variation and strong pattern in sample (helps to visualize/explore) 
PCASamples(meth, screeplot=TRUE)
PCASamples(meth)

#Using genomicRanges - still working on this
getCpG <- function() 
gr=GRanges(seqnames=c("Chr1")),
ranges=IRanges(start=c(50,150,200), end = c(100,200,300)),
strand=c("+", "-")

#barplot of all samples - #? pdf isn't opening
barplot(samples)
dev.new()
pdf(paste(barplot(samples), "barplot_samples.pdf", sep = "_"))
#dev.off()

#plot the samples
plot(samples,
main = "DNA methylation Samples",
xlab = "" , #trying to label 1,2,3,4,5, with sample.id names
ylab = "frequency",
type = "o", #plot type
col = "red") #plot color
dev.new()
pdf(paste(plot(samples), "plot_samples.pdf", sep = "_"))
#dev.off()

#Store and output the sample.id
coverage_vals <- lapply(myobj, function(i){
  x <- getData(i)$coverage
})
names(coverage_vals) <- getSampleID(myobj)
names(coverage_vals) 

#using plots to analyze the data
#filter coverage values ?
filt_coverage_vals <- lapply(coverage_vals, function(i){
       x <- i[i >= 10 & i <= 500]
   })
boxplot(filt_coverage_vals)
 

install.packages("vioplot")
library(vioplot)
orig_margins <- par()$mar
par(mar=c(7.1,4.1,4.1,2.1))
plot(0, type="n", xlim=c(0,7), ylim=c(0,500), axes=F, main="DNA Methylation Samples", xlab="",ylab="")
for (i in 1:length(filt_coverage_vals)){ vioplot(filt_coverage_vals[[i]], h=10, at=i, add=T) }
axis(2)
axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2)

par(mfrow=c(3,2))
for (i in 1:length(filt_coverage_vals)){ plot(density(filt_coverage_vals[[i]], bw=0.1)) }

par(mfrow=c(3,2))
for (i in 1:length(filt_coverage_vals)){ plot(density(filt_coverage_vals[[i]], bw=0.1)) }
pdf(paste(names(coverage_vals) ,"plot1_output.pdf",sep="_"))
dev.off()

par(mfrow=c(1,1))
plot(0, type="n", xlim=c(0,500), ylim=c(0,0.21), axes=F, main="DNA Methylatin analysis", xlab="",ylab="")
for (i in 1:length(filt_coverage_vals)){ lines(density(filt_coverage_vals[[i]], bw=0.1)) }
axis(1)
axis(2)

par(mfrow=c(3,2))
for (i in 1:length(coverage_vals)){ plot(density(coverage_vals[[i]], bw=0.1)) }
axis(1)
axis(2)

par(mfrow=c(1,1))
plot(0, type="n", xlim=c(0,500), ylim=c(0,3.5), axes=F, main="", xlab="",ylab="")
for (i in 1:length(coverage_vals)){ lines(density(coverage_vals[[i]], bw=0.1)) }
axis(2)
axis(1)


#function for graphs 
install.packages("vioplot")
library(vioplot)
filt_coverage_vals <- lapply(coverage_vals, function(i){
  x <- i[i >= 10 & i <= 500]
})
boxplot(filt_coverage_vals)
library(vioplot) #package 

plots <- function(filt_coverage_vals, myobj) { 
 
  for (i in 1:length(filt_coverage_vals, bw=0.1))  {
  
  par(mar=c(7.1,4.1,4.1,2.1))
  plot(0, type="n", xlim=c(0,7), ylim=c(0,500), axes=F, main="DNA Methylation Samples", xlab="",ylab="")
  #for (i in 1:length(filt_coverage_vals)){ vioplot(filt_coverage_vals[[i]], h=10, at=i, add=T) }
  axis(2)
  axis(1, labels=names(filt_coverage_vals), at=seq(1,6), las=2)
    
  par(mfrow=c(3,2)) 
  plot(density(filt_coverage_vals[[i]], bw=0.1)) 
  pdf(paste(sample.output,"plot.pdf",sep="_")) #send the location setwd()
  plot(density(coverage_vals[[i]], bw=0.1))
  axis(1)
  axis(2)
  pdf(paste(sample.output,"lines.pdf",sep="_"))

  par(mfrow=c(1,1))
  plot(0, type="n", xlim=c(0,500), ylim=c(0,3.5), axes=F, main="insert title", xlab="",ylab="")
  lines(density(coverage_vals[[i]], bw=0.1)) 
  axis(2)
  axis(1)
  
  plot(0, type="n", xlim=c(0,500), ylim=c(0,0.21), axes=F, main="insert title", xlab="",ylab="")
  lines(density(filt_coverage_vals[[i]], bw=0.1)) 
  axis(1)
  axis(2)

  }
} 

