
# Simulate recombinant genomes for breakpoint enrichment analysis

# Purpose
The ranges that we can infer where recombination breakpoints occured in recombinant genomes are large (4 kb long on average, but vary from 67 nucleotidesbases to 17 kb long), so we quanified the breakpoint frequency as the number of overlapping breakpoint ranges at any one point in the genome when across all recombinant genomes. To determine if the frequencies we observed across the genome were enriched relative to a null model that assumes no bias in the location of breakpoints, we simulated recombinant genomes with different numbers of breakpoints per genome and compared the distribution of breakpoint frequencies.

This script simulates recombinant genomes by picking two random parent clades and picking random locations throughout the genome to place breakpoints. Nucleotides for the 37 clade-defining SNPs of these simulated genomes are then assigned by first randomly picked which parent clade the very first cdSNP originated from, then assigned the remaining cdSNPs of the recombinant according to the locations of the breakpoints. This script is also capable of simulating an infinite number of recombination breakpoints by assigning nucleotides for the 37 clade-defining SNPs randomly and independently. After simulating the desired number of recombinants that passed the same logic that was used to identify  real recombinants in GISAID data (≥2 cdSNPs supporting recombination and at least one pair of predicted parent clades with no conflicting cdSNPs), the fold change in the number of overlapping breakpoint ranges in the real recombinants relative to the simulated recombinants was calculated in 10 nucleotide long bins across the genome. 95% confidence intervals are calculated by including data from all batches of simulated genomes. 

# Usage
Edit the "USER DEFINED VARIABLES" section at the top of the script to select the number of breakpoints, number of genomes to simulate per batch, and number of batches to sumulate. To simulate an infinite number of breakpoints, adjust the breakpoint count to any value greater than 100. 

