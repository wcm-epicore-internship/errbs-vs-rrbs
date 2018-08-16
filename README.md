## Purpose 
The purpose of this project is to develop a software for the analysis of multiple DNA methylation files obtained from Next Generation Sequencing Technology by inputting, annotating and generating data. The software will be able to specify DNA methylation data for visualization, location of methylation sites and coverage. 

## Background 
DNA methylation is an epigenetic mechanism that is characterized by the addition of a methyl (CH3) group to the cytosine base of the DNA molecule. Methylation occurs at different regions in the genome in which the cytosine nucleotide base follows a guanine base or CpG dinucleotides.  CpG islands are regions of the DNA where the frequency of cytosine-guanine dinucleotide bases is higher compared to other regions. Methylation in these locations occur at the promoter region where transcription binding factor results in gene silencing and altered gene expression. Levels of methylation in the DNA are identified by methods such as Enhanced Reduced Representation Bisulfite Sequencing (ERRBS), Reduced Representation Bisulfite Sequencing (RRBS), Whole-Genome Bisulfite sequencing (WGBS) and Methylome capture. 
In order to create software for DNA methylation comparisons, the following steps are required: 
1.	Establish a concrete understanding of DNA methylation.  
2.	Learn and compare the different methods used to analyze DNA methylation. 
3.	Using R, analyze the DNA methylation files obtained from Pubshare database using Bioconductor and R packages such as methylKit. 
4.	 Create scripts necessary that can read DNA methylation files, filter and reformat the data, compute statistics, plot the coverage of CpG and determine methylation levels for multiple datasets. 

## Description of Project
The methylome of an organism’s genome can be obtained using different protocols which can then be sequenced and aligned to produce methylation files. The file contains information at base resolution where DNA methylation has occurred including chromosome base, strand, coverage and the frequency of methylation or unmethylation. To analyze DNA methylation files, the current methods do not allow annotation of multiple files simultaneously. In order to build the software, I will use the instructions available on the GitHub repository which contains a schedule and a Readme file. I will update the Readme file to include information about the current  methods, the software and codes used to analyze the files.  The software will use data from Pubshare database to perform analysis using R packages such as methylKit by gathering more advanced knowledge in programming and understanding the software that are already available. I will run script in R to validate its accuracy and consistency.
My final presentation will have the following format: (1) Introduction to DNA methylation, 
CpG islands and methods used for DNA methylation analysis. (2) Explain the difference in current methods and perform analysis on DNA methylation files using methylKit. (3) Describe the software that I have created, i.e. bash scripts, R. (4) Use of the scripts in programs such as R on various platforms, servers, etc. (5) Advantages or disadvantages of the script(s). (6) Skills or knowledge gained from the project/internship.  


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

