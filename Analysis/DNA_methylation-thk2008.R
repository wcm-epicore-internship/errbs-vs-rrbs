#R packages
library(methylKit)
library(GenomicRanges)

use.compressed.files <- FALSE

#read in file
# data_path <- ""/Users/shazedaomar/Desktop/Internship_notes/MethylCall_Data_Epicore""
data_path <- "./data"
file.names.list <- list.files(path = data_path, pattern = "*.mincov0.txt.gz", full.names = TRUE)
head(file.names.list)

#methread for MULTIPLE FILES
#read the files to a methylRawListDB object: myobjDB
# save the obj so that it doesn't have to be run everytime??
# (thk2008) indeed, while testing you may wish to save myobj so you dont have to recreate it every time

# rda <- paste(rdataDir, "mcgrlist.rda", sep="/")
# cat("saving mcgrlist to [", rda, "]\n")
# save(mcgrlist, file=rda)
# but then you will have to load it

if (!exists("myobj")) {
    cat("cannot find myobj\n")
  if (file.exists("myobj.rda")) {
    cat("loading myobj.rda\n")
    load("myobj.rda")
  } else {
    cat("generating myobj\n")
    # create myobj and save it
    myobj=methRead(as.list(file.names.list),
      sample.id=as.list(gsub(".mincov0.txt.gz", "", basename(file.names.list))), #gsub finds all the files that end in ("*.mincov0.txt.gz"), Basename removes the ending of the file.
      assembly="hg19",
      treatment=rep(0,6), #rep function - 0 = same , 6 = number of items in the file.names.list
      context="CpG",
      dbdir = "methylDB"
    )
    cat("saving myobj\n")
    save(myobj, file="myobj.rda")
  }
}

#creates a path to the folder methylDB -- compressed and indexed files
if (use.compressed.files) {
  print(myobj[[1]]@dbpath) #need to run the dbtype = "tabix", command to get result
}

#for loop to count the number of rows in the files - function for nrow, statistics and coverage.
loop <- function() {
  for (i in 1:length(file.names.list))
  print(nrow(myobj[[i]]))
}

head(loop)

for (i in 1:length(file.names.list)) {
    print(nrow(myobj[[i]]))
}

#Output the file names and nrows for each sample
samples <- sapply(myobj,nrow)
# names(samples) <- getSampleID(myobj) # c("ERRBS", "methylcaptureA", "methylcaptureB", "methylcaptureC", "RRBS", "WGBS")
names(samples) <- c("ERRBS", "methylcaptureA", "methylcaptureB", "methylcaptureC", "RRBS", "WGBS")
file.list <- as.list(samples)

#using a function to find nrow, statistics and histogram (condense code)
methylation <- function(myobj) {
    for (i in 1:length(myobj)) {
        #print(basename(file.names.list[[i]]))
        sample.name <- myobj[[i]]@sample.id
        cat("sample", sample.name, "\n", sep=": ")
        cat("number of sites", nrow(myobj[[i]]), "\n", sep=": ")
        pdf(paste(sample.name,"methylation.pdf",sep="_"))
         getMethylationStats(myobj[[i]])
         getMethylationStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
         getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
        dev.off()
    }
}



#descriptive statistics on the ERRBS, RRBS, WGBS, methylcaptureA, methylcaptureB, methylcaptureC
#Option2
for (i in 1:length(file.names.list)) {
    print(basename(file.names.list[[i]]))
    getMethylationStats(myobj[[i]])
}

#Create a histogram with the percentage of CpG for each dataset
for (i in 1:length(file.names.list)) {
    print(basename(file.names.list[[i]]))
    getCoverageStats(myobj[[i]],plot=TRUE,both.strands=FALSE)
}

#barplot to graph the files
#plot the data that's in my obj or samples
barplot(table(samples),  #how to add names of files ERRBS,etc
        main="DNA Methylation files",
        xlab="",
        ylab="",
        col="red"
)

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
