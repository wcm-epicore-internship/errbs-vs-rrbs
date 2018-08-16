# Schedule

**Week 1 - Week 3:** Read background articles on DNA methylation, CpG islands and methylation library protocols such as ERRBS, RRBS, WGBS and Methylome capture. 
**Week 4 – Week 6:** Compare methylation library protocols. What's the difference in the library preparation protocols?
- The [ERRBS paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002781)
**Week 7 – Week 9:** Read about the MethylKit and install into R studio. Run analysis on data from datasets such as ERRBS, RRBS, WGBS, etc. Include reading methylcall files, compute statistics, plot methylation values and compare methylation site calls. 
**Methylcall data**
- We use [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) to align the sequencing reads (fastq) and genreate methylation calls
  - It would be useful to look at [Bismark paper](https://academic.oup.com/bioinformatics/article/27/11/1571/216956) and their [RRBS guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf)
  - The authors of the errbs paper developed, [methylKit](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3491415/), software for analyzing methyllation data

**Week 10 - Week 15:** Create bash scripts that can perform numerous analysis on the DNA methylation files including coverage, filtering, identifying the methylated value of cytosine or CpG islands and compute statistics. 
Week 16 - Week 18: Update the GitHub repository, Readme file. Prepare the final report and presentation. 

**Find the datasets for each**
- in R or Python
  - read in methylcall files
  - plot methylation values
  - compare site calls
    - how man sites
    - how many common
  
