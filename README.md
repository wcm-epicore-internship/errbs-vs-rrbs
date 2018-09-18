## Purpose 
The purpose of this project is to develop a software for the analysis of multiple DNA methylation data obtained from Next Generation Sequencing Technology by inputting, annotating and generating data. The software will be able to evaluate DNA methylation data for visualization, location of methylation sites and coverage.

## Background 
DNA methylation at cytosine is an epigenetic mark critical in mammalian cells for a variety of biological processes, including but not limited to imprinting, X chromosome inactivation, development, and regulation of gene expression (Garrett-Bakelman, 2015). DNA methylation can occur at different regions in the genome in which the cytosine nucleotide base follows a guanine base or CpG dinucleotides. CpG islands are regions of the DNA where the frequency of cytosine- guanine dinucleotide bases is higher compared to other regions in the genome. DNA methylation changes in these locations occur at the promoter region where transcription binding factor results in gene silencing and altered gene expressions. Levels of methylation in the DNA can be identified using methods such as Enhanced Reduced Representation Bisulfite Sequencing (ERRBS), Reduced Representation Bisulfite Sequencing (RRBS), Whole-Genome Bisulfite sequencing (WGBS) and Methylome capture. The ERRBS protocol has been used to profile cytosine methylation in DNA from human, mouse and other organisms, producing an increase in the biologically relevant genomic loci covered (Garrett-Bakelman, 2015).

In order to create a software and scripts to distinguish the variation in the DNA methylation protocols, a concrete understanding of DNA methylation is required. A comparison between the protocols can demonstrate variation in the methods of analyzing DNA methylation which impacts the resulting DNA methylation file. Currently in the field of Epigenomics, methods of analyzing DNA methylation are limited in aiding researchers to choose the best protocol for a given sample. The software will be beneficial in finding the differences and commonalities in the protocols from DNA methylation files. Upon completion of this project, skills gained include knowledge of R programming language, the use of Bioconductor and R packages and downloading data from one or more servers.

## Description of Project
The methylome of an organism’s genome can be obtained using different protocols which can then be sequenced and aligned to produce methylation files. The file contains information at base resolution where DNA methylation has occurred including chromosome base, strand, coverage and the frequency of methylation or unmethylation. To build a software that compares DNA methylation protocols, I will use the instructions available on the GitHub repository which contains 

a schedule and a Readme file. I will update the Readme file to include information about the different protocols, the software and codes used to analyze the data. The software will use data from the Epicore Server to perform analysis using R packages such as MethylKitand GenomicRanges. MethylKit is a multi-threaded R package used to rapidly analyze and characterize
I will run the scripts in R Studio to validate the accuracy and consistency of the ERRBS method coverage.
My final presentation will have the following format: (1) Introduction to DNA methylation, CpG islands and methods used for DNA methylation analysis. (2) Explain the difference in protocols and perform analysis using R packages such as methylKit and GenomicRanges. (3) Describe the software or scripts that I have created. (4) Use the scripts in R Studio and various platforms, servers, etc. (5) Explain the advantages or disadvantages of the script and ERRBS method. (6) Skills or knowledge gained from the project/internship.


## ERRBS and RRBS libraries
__Reduced representation bisulfite sequencing (RRBS):__ an efficient method for quantitative, base-pair resolution of cytosine methylation across the genome.

__Protocol:__
1. Enzyme Digestion: Digest the genomic DNA with Msp1 enzyme targeting 5’CCGG3’ sequences. It is a Methylation insensitive digestive enzyme, meaning whether or not the CG is methylated, the Msp1 enzyme will cleave it. 
1. End repair, where Msp1 digested fragments have sticky ends, the 3’ terminal of the ends of the strands must be filled in. A tailing, the addition of an adenosine on both strands, necessary for adapter ligation. Then Ligate (close off) with methylated adaptors, so the flow cell can recognize the sequences. If they were not methylated, once treated with the bisulfate the C -> U and will no longer be recognized (and attach) to the flow cell.
1. Size Selection/Reduced representation: Cut out a small section of the genome, 40-200 bp fragments (2.5 of the genome, but CpG rich). Avoids a lot of sequencing, while focusing on CpG rich regions that are more likely to cause changes in gene expression (CpG islands).
1. Bisulfate conversion: distinguish between methylated and un-methylated cytosine. 
1. PCR amplification and sequencing. 

__Enhanced Reduced representation bisulfite sequencing (ERRBS):__ provides biochemical and bioinformatic methodological improvements that generate more extensive and balanced coverage of methylation (improving the coverage of regions outside CpG islands).

The main difference between RRBS and ERRBS ibrary preparation protocols is the size selection is no longer done in ERRBS. 
RRBS only interrogates CpGs within short MspI delimited fragments between 40 to 220 bp, making it biased towards representing CpG islands, which typically contain more densely clustered MspI sites. In order to enhance the capture of regions beyond CpG islands, MspI fragments ranging from 70–320 bp are selected instead. This enhanced RRBS (ERRBS) method yielded a 75% increase in coverage of CpG sites with a 54% increase in coverage of CpG shores. 


## Analysis 
__Read in methylcall files:__

``` errbs <- read.table("/Methyl/File/Path.txt", header = TRUE) ```

__Plot methylation values:__

``` Hist(errbs$freqC, main="Methylation-RRBS Data", xlab="FreqC", col="red", border="black") ```

__Comparing site calls__

__Calculate the amount of sites:__

``` nrows(errbs) ``` 

__Calculate the amount of common sites:__ 

``` 
common <- intersect(errbs$chrBase,rrbs$chrBase)
length(common) 
```

