## Purpose 
The purpose of this project is to develop a software for the analysis of multiple DNA methylation data obtained from Next Generation Sequencing Technology by inputting, annotating and generating data. The software will be able to evaluate DNA methylation data for visualization, location of methylation sites and coverage.

## Background 
The human genome contains essential information that is used to maintain the expression of genes and the regulation of proteins in different organ systems. Underlying mechanisms can directly or indirectly affect this process resulting in gene alterations. Epigenetics is a process that creates alterations in the gene expressions without modifying the DNA sequence1. The DNA methylation mechanism that is crucial for development with an important role in many vital processes including X-chromosome inactivation, and suppression of repetitive element transcription and, when dysregulated, contributes to diseases like cancer2. DNA methylation influences gene expression by affecting the nitrogenous cytosine base and areas where the cytosine follows a guanine base, commonly called CpG islands (CGIs) or short interspersed sequences.

Many tools developed can significantly differ the common genomic pattern by being GC-rich, CpG-rich, and predominantly nonmethylated3. Methylation is regulated by DNA methyltransferases (Dnmts) which shifts the methyl group from the S-adenyl methionine (SAM) to the fifth carbon cytosine residue forming a 5-methylcytosine4. The detection of DNA methylation uses bisulfite sequencing in protocols such as Enhanced Reduced Representation Bisulfite Sequencing (ERRBS)5, Reduced Representation Bisulfite Sequencing (RRBS)6, Whole-Genome Bisulfite sequencing (WGBS) and Methylome Capture Sequencing (targeted methylome sequencing)7. The hypothesis conducted will evaluate the Enhanced Reduced Representation Bisulfite Sequencing (ERRBS) protocol which provides improved coverage of DNA methylation sites across specific regions of the genome compared to other protocols.

Bisulfite sequencing detection is used to identify segments of methylation or nonmethylation in the DNA, it uses sodium bisulfite to find the amination reactions of cytosine and 5-methylcytosine(5mC)5. Cytosines in single-stranded DNA are converted into uracil and recognized as thymine in PCR amplification, but 5mCs are immune to this conversion and remain as cytosines allowing it to be analyzed from unmethylated cytosines8. Both the RRBS and ERRBS protocol are similar in that they use bisulfite sequencing and restriction enzymes. The RRBS protocol uses methylation-dependent restriction enzymes (MSREs) in fragments of 500-600bps7 whereas the ERRBS is an advanced protocol derived from RRBS which produces more in-depth methylation coverage of CpG islands (CGIs) or cytosines. The ERRBS protocol uses a restriction enzyme that generates a low molecular fragment weight in library preparation5 and yields approximately ten percent of genomic CpG sites9. This provides enrichment in CpG islands and CpG shores, promoters, exons, introns and intergenic regions9.

The Whole genome bisulfite sequencing (WGBS) is produced by Next Generation Sequencing technology (NGS) that can capture all the cytosine in the genome at single-nucleotide resolution but has several drawbacks amplified with increasing sample numbers7. While the WGBS is increasingly accessible for both primary and clinical research10, its cost has remained substantial and limits the widespread use for multi-sample comparison of large methylomes such as mammals11. An alternative to WGBS is Methylome Capture Sequencing or targeted methylome sequencing (TMS). Methylome Capture Sequencing is based on the capture of methylated DNA using the methyl-binding domain of methyl CpG binding protein and subsequent NGS of eluted DNA, allowing the DNA to be separated by salt gradients according to its CpG methylation density12. An aligner tool, Bismark uses a combination of bisulfite treatment of DNA and high throughput sequencing (BS-Seq) for efficient time analysis, read mapping and methylation calling13.

## Description of Project
The methylome of an organism’s genome can be obtained using different protocols which can then be sequenced and aligned to produce methylation files. The file contains information at base resolution where DNA methylation has occurred including chromosome base, strand, coverage and the frequency of methylation or unmethylation. To build a software that compares DNA methylation protocols, I will use the instructions available on the GitHub repository which contains 

a schedule and a Readme file. I will update the Readme file to include information about the different protocols, the software and codes used to analyze the data. The software will use data from the Epicore Server to perform analysis using R packages such as MethylKitand GenomicRanges. MethylKit is a multi-threaded R package used to rapidly analyze and characterize
data from many methylation experiments at once (Akalin, 2012). I will run the scripts in R Studio to validate the accuracy and consistency of the ERRBS method coverage.
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

