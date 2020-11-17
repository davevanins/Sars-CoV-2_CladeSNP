
# SARS-CoV-2 cladeSNP_blast

# Purpose
cladeSNP_blast is a local alignment based version of the full cladeSNP pipeline that is designed to screen SARS-CoV-2 genomes for recombinants without needing to make a whole genome, global alignment. It works by trimming and quality screening the input genomes before using BLASTn to find 37 clade defining SNPs, then screening those SNP profiles for genomes with atypical combinations of SNPs that could be indicative of recombination.

Running cladeSNP_blast on a 3.5 GHz single-core CPU with 10,000 GISAID genomes as input took approximately 15 minutes to complete. 

# Dependencies

* NCBI Blastn

# Usage
All template sequences are included in the blast_template_files directory, so the only input is the your input genome fasta formatted file.  To use this script, simply run "python cladeSNP_blast.py your_input_sequence_file.fasta" from the cladeSNP_blast directory.

# Interpreting the output

There are four main output files that summarize any information supporting recombination

* putative_recombinant_info.txt - is the output for information on all possible parent clade pairs for each putative recombinant, including the number of predicted recombination breakpoints and SNPs supporting recombination 
* recombinants.min_diffSNP_distance.txt - another output for genomes that pass recombinant screen, including the clade defining SNP profile of the putative recombinant and how that profile differs from each of the 14 clades' cdSNP profiles
* all_genomes.min_cdSNP_distance.txt - reports the minimum distance of all input genomes to the nearest clade based on clade-defining SNPs. In current implementation, N's are not counted towards distance
* zero_cdSNP_distance.clade_assignment.txt - outputs only those genomes with clade-defining SNP profiles that perfectly match one of the 14 clades' profiles
