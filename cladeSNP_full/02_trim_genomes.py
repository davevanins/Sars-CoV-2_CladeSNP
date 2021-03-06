import os
###################################   USER DEFINED VARIABLES   ###################################

project_dir = "./gisaid_2020/"
blast_dir = project_dir + "blast/"
output_prefix = 'gisaid_hcov19'

genome_oneline_filename = output_prefix+".fa"
reference_alignment_filename = project_dir+"hcov19.reference_alignment.trim.subset.fa"

flanking_seq_filename = "cov2_flank_seq.fasta"
max_N = 0.1
min_length = 29610
max_length = 29660

slurm_prefix = "#!/bin/bash\n#SBATCH -N 1\n#SBATCH -n 9\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=6000\n#SBATCH --time=48:00:00\ncd "+project_dir+"\n"

############################################ FUNCTIONS ############################################

############################################   MAIN   #############################################
accessions_to_include = []
try:
	include_infile = open(project_dir+"genomes_to_include.txt","r")
	for line in include_infile:
		line = line.strip()
		accessions_to_include.append(line)
	include_infile.close()
except:
	pass

flank_position_dict = {}
blast_infile = open(blast_dir+output_prefix+".flank.blastn.txt","r")
for line in blast_infile:
	#"0-qseqid 1-sseqid 2-pident 3-evalue 4-qlen 5-slen 6-length 7-qstart 8-qend 9-sstart 10-send"
	line = line.strip().split("\t")
	accession = line[0]
	side = line[1]
	pid = float(line[2])
	align_len = int(line[6])
	qstart = int(line[7])
	qend = int(line[8])
	if side == "front":
		cut_site =  qstart
	elif side == "back":
		cut_site = qend
	try:
		prev_site = flank_position_dict[accession][side]
		flank_position_dict[accession][side] = min(prev_site,cut_site)
	except:
		try:
			flank_position_dict[accession][side] = cut_site
		except:
			flank_position_dict[accession] = {}
			flank_position_dict[accession][side] = cut_site
blast_infile.close()


seq_dict = {}
name_dict = {}
genome_infile = open(project_dir+genome_oneline_filename,"r")
accessions_to_include_seq_dict = {} 
quality_outfile = open("genome_info.txt","w")
for line in genome_infile:
	line = line.strip()
	if line[0] == ">":
		accession = line[1:len(line)]
	else:
		seq_in = line
		try:
			front_cut_site = max(flank_position_dict[accession]["front"]-1,0)
			back_cut_site = flank_position_dict[accession]["back"]
			seq_out = seq_in[front_cut_site:back_cut_site]

			len_before = len(seq_out)
			len_after = len(seq_out.replace("N",""))
			Ncount = len_before-len_after
			N_prop = float(Ncount)/float(len_before)
			
			if len(seq_out) >= min_length and len(seq_out) <= max_length:
				quality_outfile.write(accession+"\t"+str(len_before)+"\t"+str(Ncount)+"\n")
				try:
					ref_iso = seq_dict[seq_out]
					name_dict[ref_iso].append(accession)
					if accession in accessions_to_include:
						accessions_to_include_seq_dict[seq_out] = accession
				except:
					seq_dict[seq_out] = accession
					name_dict[accession] = []
					name_dict[accession].append(accession)
		except:
			pass
genome_infile.close()
quality_outfile.close()

##Write unique sequences
unique_sequences_fasta = project_dir+output_prefix+".filtered.unique.fa"
unique_sequences_names = project_dir+output_prefix+".filtered.names"
alignment_outname = project_dir+output_prefix+".filtered.unique.align.fa"

unique_genome_outfile = open(unique_sequences_fasta,"w")
for seq in seq_dict:
	unique_genome_outfile.write(">"+seq_dict[seq]+"\n"+seq+"\n")
for seq in accessions_to_include_seq_dict:
	unique_genome_outfile.write(">"+accessions_to_include_seq_dict[seq]+"\n"+seq+"\n")
unique_genome_outfile.close()
del seq_dict
del accessions_to_include_seq_dict

unique_names_outfile = open(unique_sequences_names,"w")
for ref_accession in name_dict:
	names = name_dict[ref_accession]
	out_line = ref_accession+"\t"+ref_accession
	for i in range(0,len(names)):
		iso = names[i]
		if iso != ref_accession:
			out_line += ","+iso
	out_line += "\n"
	unique_names_outfile.write(out_line)
unique_names_outfile.close()
del name_dict

command = "mafft --thread -9 --nomemsave --keeplength --add "+unique_sequences_fasta+" --reorder "+reference_alignment_filename+" > "+alignment_outname
# os.system(command)
shell_out = open(project_dir+"mafft_align.sh","w")
shell_out.write(slurm_prefix+command+"\n")
shell_out.write("python 03_alignment_parse_cladeSNPs.py\n")

shell_out.close()
os.system("sbatch "+project_dir+"mafft_align.sh")
