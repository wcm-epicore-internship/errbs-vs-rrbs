#R packages
library(methylKit)
library(GenomicRanges)

compressed.files <- FALSE

#read in files
#store path into variable data.file.path
#data.file.path <- "/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Epicore"
data <- data.file.path
file.names.list <- list.files(path = data, pattern = "*.mincov0.txt.gz", full.names = TRUE)
file.names.list

#methread for Multiple Files
#save it so that it doesn't have to be run everytime - saves more time
rda <- paste(rdataDir, "mcgrlist.rda", sep="/")
cat("saving mcgrlist to [", rda, "]\n")
save(mcgrlist, file=rda)
# but then you will have to load it

if  (!exists("myobj")) {
    cat("Did not find myobj\n")
    elseif(file.exists("myobj.rda")) {
        load("myobj.rda\n")
        cat("myobj.rda loaded")
   else { 
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
} 


#creates a path to the folder methylDB -- compressed and indexed files .tbi and .bgz
#print(myobj[[1]]@dbpath) #need to run the dbtype = "tabix", command to get result

#Output the file names and nrows for each sample
samples <- sapply(myobj,nrow)
names(samples) <- getSampleID(myobj)
file.list <- as.list(samples)
file.list

#myobj[[i]]@sample.id

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

#view where the pdf is located
getwd()

#barplot of all samples 
barplot(samples)
#pdf(paste(barplot(samples), "barplot_samples.pdf", sep = "_"))
#dev.off()

    
#plot the samples
plot(samples,
main = "DNA methylation Samples ",
xlab = "samples",
ylab = "frequency",
#dev.off()
)


#filtering the data by coverage to avoid bias (obtain coverage that is not too high or too low)
filter_myobj = filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
hi.count=NULL,hi.perc=99.9)
head(filter_myobj, n=20)


#Merging sample to find common sites
meth=unite(myobj, destrand=FALSE)
head(meth) #view the file meth

#creates a methylBase object where only CpGs covered with at least 1 sample per group will be returned
meth.min=unite(myobj,min.per.group=1L)
head(meth.min)

#correlation
getCorrelation(meth,plot=TRUE)
#get numbers to compare
