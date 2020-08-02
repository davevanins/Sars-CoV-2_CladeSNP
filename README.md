
# SARS-CoV-2 cladeSNP

# Purpose
This pipeline is designed to take the raw SARS-CoV-2 genomes from GISAID, align them to a reference alignment, then identify putative recombinant genomes for manual inspection

# Dependencies

* MAFFT (initially designed an operated using v7.464)
* NCBI Blastn

# Usage
At the beginning of each script is a "user defined variables" section which needs to be updated to reflect the directory of the project, name of the input file downloaded from GISAID, and name of the reference alignment to which the genomes should be aligned.
To replicate the results from our manuscript, we recommend downloading the NCBI SARS-CoV-2 reference genome (accession NC_045512) to align all GISAID sequences against, but you can use any fasta formatted alignment. 

Scripts should be run in the following order:

`01_genome_rename.py`
* Genomes uploaded to GISAID often have different standards for how to handle ambiguous nucleotides, as DNA vs RNA sequences, etc... so this script reformats all sequences in the input file and 

`02_trim_genomes.py`
* This script takes the blast output file and uses that information to consistently cut each genome and generate a set of unique sequences to align to the reference alignment. It finishes by submitting a MAFFT alignment job.

`03_alignment_parse_cladeSNPs.py`
* This script generates the final output for manually identifying putative recombinant genomes. After searching each genome for the clade-defining SNPs and comparing them to the SNP profiles of each clade, it outputs all sequences that have more than 1 difference in its clade-defining SNP profile relative to the nearest clade.

#Interpreting the output
The output lists each query genome on a separate line, along with the number of clade-defining SNP differences that sequence is to the nearest clade. Following that information is a visual depiction of similarities and differences with the SNP profiles of each clade. When a clade shares a nucleotide at a position, a '-' is entered, but when they differ the nucleotide from the clade SNP profile is written. Recombinants are identified by scanning for sequences that appear to be reassortments of two specific clades.

For examples of what validated recombinant sequences look like in the output file, we recommend include the following sequences in your analysis for reference:
* EPI_ISL_468407
* EPI_ISL_452334
* EPI_ISL_475584
* EPI_ISL_454983
* EPI_ISL_464547

# Generating clade-defining SNPs for clades you chose
`alignment_parse_groupSNPs.py` 
* The script takes a reference alignment and a text file that contains a key for what clade each genome in the reference alignment is part of in the format: "genome \t cladeID \n". The output is a file containing the SNPs that differentiate the clades, and a separate file indicating the prevalence of those SNPs in each clade

# Notes
* If there are any sequences which you wish to exclude from your analysis, make a file called "genomes_to_exclude.txt" and put it in the working directory with the GISAID 'ESP_ISL_####' accession numbers, each listed on a new line separately.
* Likewise, if there are genomes you want to be sure are included, make a file called "genomes_to_include.txt" and put it in the working directory. 
* cladeSNP was written to be used on a computer cluster with slurm workload manager, but it can be run locally with minimal edits (assuming MAFFT and blastn are installed on your machine).
