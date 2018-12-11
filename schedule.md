# Schedule

**Week 1 - Week 3:** Read background articles on DNA methylation, CpG islands and methylation library protocols such as ERRBS, RRBS, WGBS and Methylome capture.

**Week 4 – Week 6:** Compare methylation library protocols. What's the difference in the library preparation protocols?
- The [ERRBS paper](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002781)

**Week 7 – Week 9:** Read about the MethylKit and install into R studio. Run analysis on data using library protocol data and include reading in methylcall files, computing statistics, plot methylation values and compare methylation site calls.

**Week 10 - Week 15:** Run analysis on the DNA Methylation data such as ERRBS, RRBS, WGBS and methylome capture. Create a script that can perform numerous analysis on DNA methylation data including reading DNA methylation files, filtering and reformatting the data, computing and graphing statistics, finding and plotting the coverage of CpG and determining methylation levels for multiple datasets.

**Week 16 - Week 18:** Create a .RMD file using Rscript, update the GitHub repository and prepare the final report. 

**Methylcall data**
- We use [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/) to align the sequencing reads (fastq) and genreate methylation calls
  - It would be useful to look at [Bismark paper](https://academic.oup.com/bioinformatics/article/27/11/1571/216956) and their [RRBS guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/RRBS_Guide.pdf)
  - The authors of the errbs paper developed, [methylKit](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3491415/), software for analyzing methyllation data

**Find the datasets for each**
- I have palaced some data here (you may wish to copy it somewhere local):
  - /athena/epicore/users/scratch/thk2008/methyl_seq_comparison_data/wcm-epicore-internship
- in R or Python
  - read in methylcall files
  - plot methylation values
  - compare site calls
    - how many sites
    - how many common

