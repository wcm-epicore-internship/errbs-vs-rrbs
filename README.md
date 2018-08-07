# ERRBS-vs-RRBS
A comparison of ERRBS methylation data with RRBS methylation data

### Background 
__DNA methylation__ is an epigenetic mechanism in which methyl groups are added to the DNA molecule. Methylation can change the activity of a DNA segment without changing the sequence, playing a role in gene expression regulation. 
__Bisulfite sequencing__ (also known as bisulphite sequencing) is the use of bisulfite treatment of DNA to determine its pattern of methylation. In animals it predominantly involves the addition of a methyl group to the carbon-5 position of cytosine residues of the dinucleotide CpG. 
__Bisulfite sequencing__ applies routine sequencing methods on bisulfite-treated genomic DNA to determine methylation status at CpG dinucleotides. Treatment of DNA with bisulfite converts cytosine residues to uracil, but leaves 5-methylcytosine residues unaffected. Therefore DNA that has been treated with bisulfate retains only methylated cytosines thus yielding single-nucleotide resolution information about the methylation status of a segment of DNA. 

### ERRBS and RRBS libraries
__Reduced representation bisulfite sequencing (RRBS):__ an efficient method for quantitative, base-pair resolution of cytosine methylation across the genome.

__Protocol:__
1. Enzyme Digestion: Digest the genomic DNA with Msp1 enzyme targeting 5’CCGG3’ sequences. It is a Methylation insensitive digestive enzyme, meaning whether or not the CG is methylated, the Msp1 enzyme will cleave it. 
1. End repair, where Msp1 digested fragments have sticky ends, the 3’ terminal of the ends of the strands must be filled in. A tailing, the addition of an adenosine on both strands, necessary for adapter ligation. Then Ligate (close off) with methylated adaptors, so the flow cell can recognize the sequences. If they were not methylated, once treated with the bisulfate the C -> U and will no longer be recognized (and attach) to the flow cell.
1. Size Selection/Reduced representation: Cut out a small section of the genome, 40-200 bp fragments (2.5 of the genome, but CpG rich). Avoids a lot of sequencing, while focusing on CpG rich regions that are more likely to cause changes in gene expression (CpG islands).
1. Bisulfate conversion: distinguish between methylated and un-methylated cytosine. 
1. PCR amplification and sequencing. 
