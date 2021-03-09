
# Summary of accession numbers included in analysis, and of putative recombinants identified

`genome_info.txt`
* This file has the GISAID accession number of every genome included in the analyses of our manuscript, along with basic summary statistics about those genomes after they have been processed and aligned to the reference genome (including trimmed genome length, N count, nearest clade, and distance to the nearest clade).

`high_quality_reference_genomes.txt`
* This file contains the accession numbers of the high quality reference genomes that were selected for minimizing redundancy while maximizing diversity to subcluster the NextStrain clades and identify clade-defining SNPs

`unique_genomes.min_cdSNP_distance.txt`
* This file contains the cdSNP profiles and the nearest clade and distance to that clade of all unique genomes included in the analysis.

`recombinants.min_diffSNP_distance.txt`
* This file lists all of the putatively recombinant genomes identified by our analysis, along with the clade-defining SNP profile of that genome, and how that profile is similar and different relative to all 14 clades. The file also has a summary of the clade-defining SNP profiles of the 14 clades. This file also summarizes detected amino acid substitutions in each recombinant sequence.

`recombinant_parent_clade_info.txt`
* This file lists all of the possible parent clade pairs for the putatively recombinant genomes identified by our analysis, along with the distance to predicted parent clades, the number of SNPs supporting recombination, and the minimum number of break points required to explain that genotype as recombination between the two parent clades.
