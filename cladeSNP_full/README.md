
# SARS-CoV-2 cladeSNP

# Purpose
This pipeline is designed to take the raw SARS-CoV-2 genomes from GISAID, align them to a reference alignment, then identifies putative recombinant genomes by highlighting genomes with unusual combinations of clade-defining SNPs and automatically determining if recombination between any two parent clades could explain the query genome.

# Dependencies

* MAFFT (initially designed an operated using v7.464)
* NCBI Blastn

# Usage
At the beginning of each script is a "user defined variables" section which needs to be updated to reflect the directory of the project, name of the input file downloaded from GISAID, and name of the reference alignment to which the genomes should be aligned.
To replicate the results from our manuscript, we recommend downloading the NCBI SARS-CoV-2 reference genome (accession NC_045512) to align all GISAID sequences against, but you can use any fasta formatted alignment. We included the trimmed version of this genome in this section to replicate our analysis exactly ("SARSCoV2_trimmed.fasta").

Scripts should be run in the following order:

`01_genome_rename.py`
* Genomes uploaded to GISAID often have different standards for how to handle ambiguous nucleotides, as DNA vs RNA sequences, etc... so this script reformats all sequences in the input file and 

`02_trim_genomes.py`
* This script takes the blast output file and uses that information to consistently cut each genome and generate a set of unique sequences to align to the reference alignment. It finishes by submitting a MAFFT alignment job.

`03_alignment_parse_cladeSNPs.py`
* This script generates the final output for manually identifying putative recombinant genomes. After searching each genome for the clade-defining SNPs and comparing them to the SNP profiles of each clade.


# Interpreting the output

There are four main output files that summarize any information supporting recombination

* putative_recombinant_info.txt - is the output for information on all possible parent clade pairs for each putative recombinant, including the number of predicted recombination breakpoints and SNPs supporting recombination 
* recombinants.min_diffSNP_distance.txt - another output for genomes that pass recombinant screen, including the clade defining SNP profile of the putative recombinant and how that profile differs from each of the 14 clades' cdSNP profiles
* all_genomes.min_cdSNP_distance.txt - reports the minimum distance of all input genomes to the nearest clade based on clade-defining SNPs. In current implementation, N's are not counted towards distance
* zero_cdSNP_distance.clade_assignment.txt - outputs only those genomes with clade-defining SNP profiles that perfectly match one of the 14 clades' profiles

For examples of what validated recombinant sequences look like in the output file, we recommend include the following sequences in your analysis for reference:
* EPI_ISL_468407
* EPI_ISL_475584
* EPI_ISL_479572
* EPI_ISL_444583

# Generating clade-defining SNPs for clades you chose
`alignment_parse_groupSNPs.py` 
* The script takes a reference alignment and a text file that contains a key for what clade each genome in the reference alignment is part of in the format: "genome \t cladeID \n". The output is a file containing the SNPs that differentiate the clades, and a separate file indicating the prevalence of those SNPs in each clade

# Notes
* If there are any sequences which you wish to exclude from your analysis, make a file called "genomes_to_exclude.txt" and put it in the working directory with the GISAID 'ESP_ISL_####' accession numbers, each listed on a new line separately.
* Likewise, if there are genomes you want to be sure are included, make a file called "genomes_to_include.txt" and put it in the working directory. This will not over-rule any quality filtering steps, but when the genome sequences are filtered to remove redundant sequences prior to alignment, it will make sure those genomes are included in the alignment and in all subsequent analysis.
* cladeSNP was written to be used on a computer cluster with slurm workload manager, but it can be run locally with minimal changes (assuming MAFFT and blastn are installed on your machine).
* The nucleotide positions in the file `clade_associated_SNPs.nucl.txt` are 118 nucleotides smaller than their true position in the SARS-CoV-2 reference genome because the 5' and 3' ends of the genome are trimmed to positions 118 and 29740, respectively, because these regions often assemble poorly.
