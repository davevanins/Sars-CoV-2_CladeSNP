import os
###################################   USER DEFINED VARIABLES   ###################################

project_dir = "./gisaid_2020/"
blast_dir = project_dir + "blast/"
output_prefix = 'gisaid_hcov19'

input_name = 'gisaid_2020_all' #The file downloaded from GISAID.org, assumed to be a .fasta file
reference_alignment_filename = project_dir+"hcov19.reference_alignment.trim.subset.fa"
flanking_seq_filename = "cov2_flank_seq.fasta" #contains the sequence that will determine where each raw genome will be cut to remove the beginning and end of each genome, which has very inconsistent quality

outfile_name = output_prefix+".fa"

max_N = 0.1
min_length = 29500
max_length = 31000

slurm_prefix = "#!/bin/bash\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=2000\n#SBATCH --time=12:00:00\ncd "+project_dir+"\n"

############################################ FUNCTIONS ############################################

def nucleotide_character_edit(seqin,ambiguous_nucleotide_list):
	seqin = seqin.replace("a","A")
	seqin = seqin.replace("c","C")
	seqin = seqin.replace("g","G")
	seqin = seqin.replace("t","T")
	seqin = seqin.replace("u","T")
	seqin = seqin.replace("U","T")
	seqin = seqin.replace("n","N")
	for code in ambiguous_nucleotide_list:
		seqin = seqin.replace(code,"N")
	return seqin

############################################   MAIN   #############################################
accessions_to_exclude = []
try:
	exclude_infile = open(project_dir+"genomes_to_exclude.txt","r")
	for line in exclude_infile:
		line = line.strip()
		accessions_to_exclude.append(line)
	exclude_infile.close()
except:
	pass
accessions_to_include = []
accessions_to_include_found = {}
try:
	include_infile = open(project_dir+"genomes_to_include.txt","r")
	for line in include_infile:
		line = line.strip()
		accessions_to_include.append(line)
		accessions_to_include_found[line] = False
	include_infile.close()
	print(accessions_to_include_found)
except:
	pass

accession_name_dict = {}
try:
	firstline = True
	metadatatsv = open(project_dir+"metadata.tsv","r")
	for line in metadatatsv:
		if firstline == True:
			firstline = False
		else:
			line = line.strip().split("\t")
			seqname = line[0]
			accession = line[2]
			accession_name_dict[seqname] = accession
	metadatatsv.close()
except:
	pass


reference_alignment = open(reference_alignment_filename,"r")
for line in reference_alignment:
	line = line.strip()
	if line[0] == ">":
		accession = line[1:len(line)]
		accessions_to_exclude.append(accession)
reference_alignment.close()
accessions_to_exclude = sorted(list(set(accessions_to_exclude)))

ambiguous_nucleotides = ['B','b','D','d','H','h','K','k','M','m','R','r','S','s','V','v','W','w','Y','y']
genomes_infile = open(project_dir+input_name+".fasta","r")
genomes_included = open(project_dir+input_name+".included.txt","w")
accession_list = []
skip = False
seq_dict = {}
for line in genomes_infile:
	line = line.strip()
	if len(line) >0:
		if line[0] == ">":
			try:
				accession = line.split("|")[1]
			except:
				accession = accession_name_dict[line[1:len(line)]]
			if accession in accessions_to_exclude or accession in accession_list:
				skip = True
			else:
				skip = False
				accession_list.append(accession)
				genomes_included.write(accession+"\n")
		elif skip == False:
			line = nucleotide_character_edit(line,ambiguous_nucleotides)
			try:
				seq_dict[accession] += line
			except:
				seq_dict[accession] = line
genomes_infile.close()
genomes_included.close()
print("Found "+str(len(accession_list))+" non-redundant genomes")
del accession_name_dict

# quality_outfile = open("Ncount.txt","w")
removed_counter = 0
removed_genomes = []
genomes_outfile = open(project_dir+outfile_name,"w")
for accession in seq_dict:
	seq = seq_dict[accession]

	seq = seq.replace("N"," ")
	seq = seq.strip()
	seq = seq.replace(" ","N")
	seq = seq.replace("-","")

	len_before = len(seq)
	len_after = len(seq.replace("N",""))
	Ncount = len_before-len_after
	N_prop = float(Ncount)/float(len_before)
	if (len_before >= min_length and len_before <= max_length and N_prop <= max_N) or accession in accessions_to_include:
		genomes_outfile.write(">"+accession+"\n"+seq+"\n")
		if accession in accessions_to_include:
			accessions_to_include_found[accession] = True
	else:
		removed_counter += 1
		removed_genomes.append(accession)
print("Removed "+str(removed_counter)+" genomes for unusual length and excessive Ns")
print(removed_genomes)
genomes_outfile.close()

if not os.path.exists(blast_dir):
	os.makedirs(blast_dir)

command = "blastn -query "+project_dir+outfile_name+" -subject "+project_dir+flanking_seq_filename+' -evalue 1E-9 -outfmt "6 qseqid sseqid pident evalue qlen slen length qstart qend sstart send" -out '+blast_dir+output_prefix+".flank.blastn.txt"
# os.system(command)

shell_out = open(project_dir+"blastn.sh","w")
shell_out.write(slurm_prefix+command+"\n")
shell_out.write("python 02_trim_genomes.py\n")
shell_out.close()
os.system("sbatch "+project_dir+"blastn.sh")
