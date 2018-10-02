#R packages
library(methylKit)
library(GenomicRanges)

#read in file 
file.names.list <- list.files(path = "/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Epicore", pattern = "*.mincov0.txt.gz", full.names = TRUE)
head(file.names.list)

#methread for MULTIPLE FILES
#read the files to a methylRawListDB object: myobjDB
# save the obj so that it doesn't have to be run everytime??  
myobj=methRead(as.list(file.names.list),
sample.id=as.list(gsub(".mincov0.txt.gz", "", basename(file.names.list))), #gsub finds all the files that end in ("*.mincov0.txt.gz"), Basename removes the ending of the file.
assembly="hg19",
treatment=rep(0,6), #rep function - 0 = same , 6 = number of items in the file.names.list
context="CpG",
dbdir = "methylDB"
)

#creates a path to the folder methylDB -- compressed and indexed files
print(myobj[[1]]@dbpath) #need to run the dbtype = "tabix", command to get result


#Output the file names and nrows for each sample  
samples <- sapply(myobj,nrow) 
names(samples) <- c("ERRBS", "methylcaptureA", "methylcaptureB", "methylcaptureC", "RRBS", "WGBS")
file.list <- as.list(samples)
file.list


#using a function to find nrow, statistics and histogram (condense code)
methylation <- function(myobj, files.names.lists) {
  
  for (i in 1:length(file.names.list)) {
    print(basename(file.names.list[[i]]))
    print(nrow(myobj[[i]]))
    getMethylationStats(myobj[[i]])
    getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
    
  }
   
}

#obtaining nrow using loop
for (i in 1:length(file.names.list)) {
    print(nrow(myobj[[i]]))
}

#descriptive statistics on the ERRBS, RRBS, WGBS, methylcaptureA, methylcaptureB, methylcaptureC
#Option2
for (i in 1:length(file.names.list)) {
    print(basename(file.names.list[[i]]))
    getMethylationStats(myobj[[i]])
}

#percent methylation distribution

for (i in 1:length(file.names.list)) {
    print(basename(file.names.list[[i]]))
    getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
    
}

#Create a histogram with the percentage of CpG for each dataset
for (i in 1:length(file.names.list)) {
    print(basename(file.names.list[[i]]))
    getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
}


#plot the samples
plot(samples,
          main = "DNA methylation Samples ",
          xlab = "",
          ylab = ""
          
)

#filtering the data by coverage
filter_myobj = filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
hi.count=NULL,hi.perc=99.9)
head(filter_myobj, n=20)


#Merging sample to find common sites
meth=unite(myobj, destrand=FALSE)
head(meth) #view the file meth

# creates a methylBase object where only CpGs covered with at least 1 sample per group will be returned
meth.min=unite(myobj,min.per.group=1L)
head(meth.min)

#correlation
getCorrelation(meth,plot=TRUE)
