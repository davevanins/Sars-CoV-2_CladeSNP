import os
###################################   USER DEFINED VARIABLES   ###################################
args = sys.argv
input_name = args[1]

project_dir = "./"
output_prefix = input_name.split(".f")[0]
blast_template_dir = project_dir+"blast_template_files/"

flanking_seq_filename = "cov2_flank_seq.fasta"
rename_outfile_name = project_dir+output_prefix+".rename.fa"
trim_outfile_name = project_dir+output_prefix+".trim.fa"
outfile_name = output_prefix+".fa"

blast_dir = project_dir + "blast/"
if not os.path.exists(blast_dir):
	os.makedirs(blast_dir)

clade_diff_snp_loc_filename = project_dir+"clade_associated_SNPs.nucl.txt"
exclude_locs_for_dist = []#if you wish to include positions 28881-3 as one site, change this list to: [28764,28765]

min_output_length = 29610
max_output_length = 29660

max_N_prop = 0.01
accessions_to_include = [] #add sequence names to this list to override N and length cutoffs. Optional, only for development or hypothesis testing

############################################ FUNCTIONS ############################################

def nucleotide_character_edit(seqin,ambiguous_nucleotide_list):
	seqin = seqin.replace("a","A")
	seqin = seqin.replace("c","C")
	seqin = seqin.replace("g","G")
	seqin = seqin.replace("t","T")
	seqin = seqin.replace("u","T")
	seqin = seqin.replace("U","T")
	seqin = seqin.replace("n","N")
	seqin = seqin.replace(".","")
	for code in ambiguous_nucleotide_list:
		seqin = seqin.replace(code,"N")
	return seqin

def check_parent_clades(clade_profiles,SNP_locations,input_SNPs):
	clades = []
	potential_parents = []
	for clade in clade_profiles:
		clades.append(clade)
	clade_pairs = []
	for i in range(0,len(clades)):
		clade1 = clades[i]
		for j in range(0,len(clades)):
			clade2 = clades[j]
			if clade1 != clade2:
				pair_f = clade1+","+clade2
				pair_r = clade2+","+clade1
				if pair_f not in clade_pairs and pair_r not in clade_pairs:
					clade_pairs.append(pair_f)
	for pair in clade_pairs:
		clade1 = pair.split(",")[0]
		clade2 = pair.split(",")[1]
		loc1 = []
		loc2 = []
		all_loc = []
		conflict = False
		for k in range(0,len(SNP_locations)):
			loc = SNP_locations[k]
			base1 = clade_profiles[clade1][loc]
			base2 = clade_profiles[clade2][loc]
			q_base = input_SNPs[loc]
			if base1 != q_base and base2 != q_base:
				conflict = True
			if base1 != base2:
				if base1 == q_base:
					loc1.append(loc)
					all_loc.append(loc)
				elif base2 == q_base:
					loc2.append(loc)
					all_loc.append(loc)
		if conflict == False and len(loc1)>0 and len(loc2)>0:
			#count number of breaks to explain recombination	
			break_count = 0
			if min(loc1) < min(loc2):
				current_parent = clade1
			else:
				current_parent = clade2

			all_loc = sorted(all_loc)
			for i in range(0,len(all_loc)):
				loc = all_loc[i]
				if loc in loc1:
					if current_parent != clade1:
						current_parent = clade1
						break_count += 1
				elif loc in loc2:
					if current_parent != clade2:
						current_parent = clade2
						break_count += 1
			tup = (break_count,len(all_loc), pair, all_loc) #len(all_loc) is included for sorting purposes
			potential_parents.append(tup)
	potential_parents = sorted(potential_parents) #sorted so predicted parents with (first) lowest number of predicted breakpoints required to explain sequence with (second) fewest total distinguishing cdSNPs are ordered first
	return potential_parents

############################################   MAIN   #############################################
##Screen genomes to replace ambiguous nucleotides with 'N'
ambiguous_nucleotides = ['B','b','D','d','H','h','K','k','M','m','R','r','S','s','V','v','W','w','Y','y']
genomes_infile = open(project_dir+input_name,"r")
skip = False
seq_dict = {}
used_list = []
for line in genomes_infile:
	line = line.strip()
	if len(line) >0:
		if line[0] == ">":
			accession = line.split("|")[1]
			if accession in used_list:
				skip = True
			else:
				skip = False
				used_list.append(accession)
		elif skip == False:
			line = nucleotide_character_edit(line,ambiguous_nucleotides)
			try:
				seq_dict[accession] += line
			except:
				seq_dict[accession] = line
genomes_infile.close()

##Remove any streches of N's from head or tail of genome (many genomes on GISAID have dozens of N's before and after the assembled content to make a consistent length)
temp_seq_dict = {}
genomes_outfile = open(rename_outfile_name,"w")
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
	genomes_outfile.write(">"+accession+"\n"+seq+"\n")
	temp_seq_dict[accession] = seq
genomes_outfile.close()
del seq_dict

##make blastn command to find location to trim the sequence
command = "blastn -query "+rename_outfile_name+" -subject "+blast_template_dir+flanking_seq_filename+' -evalue 1E-9 -outfmt "6 qseqid sseqid pident evalue qlen slen length qstart qend sstart send" -out '+blast_dir+output_prefix+".trim.blastn.txt"
os.system(command)

##parse the blast output and trim the genomes to the 118 and 29740 sites
flank_position_dict = {}
#"0-qseqid 1-sseqid 2-pident 3-evalue 4-qlen 5-slen 6-length 7-qstart 8-qend 9-sstart 10-send"
blast_infile = open(blast_dir+output_prefix+".trim.blastn.txt","r")
for line in blast_infile:
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

accession_list = []
quality_outfile = open(project_dir+output_prefix+".genome_info.txt","w")
genomes_outfile = open(trim_outfile_name,"w")
trim_seq_dict = {}
for accession in temp_seq_dict:
	seq_in = temp_seq_dict[accession]
	try:
		front_cut_site = max(flank_position_dict[accession]["front"]-1,0)
		back_cut_site = flank_position_dict[accession]["back"]
		seq_out = seq_in[front_cut_site:back_cut_site]

		len_before = len(seq_out)
		len_after = len(seq_out.replace("N",""))
		Ncount = len_before-len_after
		N_prop = float(Ncount)/float(len_before)

		if len(seq_out) >= min_output_length and len(seq_out) <= max_output_length:
			quality_outfile.write(accession+"\t"+str(len_before)+"\t"+str(Ncount)+"\n")
			trim_seq_dict[seq_out] = accession
			accession_list.append(accession)
			genomes_outfile.write(">"+accession+"\n"+seq_out+"\n")
	except:
		pass
quality_outfile.close()
genomes_outfile.close()
print("Including "+str(len(accession_list))+" non-redundant genomes")

##store 14 clade cdSNP profiles
clade_diff_snp_loc_file = open(clade_diff_snp_loc_filename,"r")
clade_diff_snps = {}
clade_list = []
diff_snp_locs = []
first_line = True
for line in clade_diff_snp_loc_file:
	line = line.strip().split("\t")
	if first_line == True:
		first_line = False
		clade_list = line
	else:
		loc = int(line[0])
		# clade_diff_snps[loc] = {}
		diff_snp_locs.append(loc)
		for j in range(1,len(line)):
			clade = clade_list[j-1]
			nt = line[j]
			try:
				clade_diff_snps[clade][loc] = nt
			except:
				clade_diff_snps[clade] = {}
				clade_diff_snps[clade][loc] = nt
clade_diff_snp_loc_file.close()
print(diff_snp_locs)

##perform blastn to find the locations of the cdSNPs in each genome
for loc in diff_snp_locs:
	command = "blastn -query "+blast_template_dir+"SNP_"+str(loc)+".fasta -subject "+trim_outfile_name+' -evalue 1E-5 -outfmt "6 qseqid sseqid pident evalue qlen slen length qstart qend sstart send sseq" -out '+blast_dir+output_prefix+".SNP_"+str(loc)+".blastn.txt"
	os.system(command)


##find cdSNPs for all query genomes
query_snp_dict = {}
query_list = []
eval_dict = {}
for loc in diff_snp_locs:
	blast_infile = open(blast_dir+output_prefix+".SNP_"+str(loc)+".blastn.txt","r")
	#"0-qseqid 1-sseqid 2-pident 3-evalue 4-qlen 5-slen 6-length 7-qstart 8-qend 9-sstart 10-send 11-sseq"
	for line in blast_infile:
		line = line.strip().split("\t")
		SNP_ID = line[0]
		loc = int(line[0].split("_")[1])
		accession = line[1]
		pid = float(line[2])
		evalue = float(line[3])
		qlen =int(line[4])
		align_len = int(line[6])
		qstart = int(line[7])
		qend = int(line[8])
		sseq = line[11]

		cdSNP_loc = (qlen/2)-qstart
		if cdSNP_loc >0:
			cdSNP = sseq[cdSNP_loc]

			query_list.append(accession)
			try:
				try:
					prev_evalue = eval_dict[accession][loc]
					if evalue < prev_evalue:
						query_snp_dict[accession][loc] = cdSNP
						eval_dict[accession][loc] = evalue
				except:
					query_snp_dict[accession][loc] = cdSNP
					eval_dict[accession][loc] = evalue
			except:
				query_snp_dict[accession] = {}
				query_snp_dict[accession][loc] = cdSNP
				eval_dict[accession] = {}
				eval_dict[accession][loc] = evalue
blast_infile.close()
query_list = sorted(list(set(query_list)))
print(str(len(query_list))+" query sequences to screen for recombinants")


parent_clade_outfile = open(project_dir+output_prefix+".putative_recombinant_info.txt","w") #output all possible parent clade pairs
parent_clade_outfile.write("accession\tparent_clade1\tparent_clade2\tclade1_dist\tclade2_dist\tnum_supporting_SNPs\trecombination_break_points\n")
minDist_outfile = open(project_dir+output_prefix+".recombinants.min_diffSNP_distance.txt","w") #output for genomes that pass recombinant screen
minDist_all_outfile = open(project_dir+output_prefix+".all_genomes.min_cdSNP_distance.txt","w") #minimum distance to nearest clade for all input genomes
zeroDist_outfile = open(project_dir+output_prefix+".zero_cdSNP_distance.clade_assignment.txt","w") #output only those genomes that match the 14 clades' cdSNP profiles perfectly

for clade in clade_list:
	minDist_outfile.write(clade+"\t")
	for x in range(0,len(diff_snp_locs)):
		loc = diff_snp_locs[x]
		s_nt = clade_diff_snps[clade][loc]#[clade]
		minDist_outfile.write(s_nt)
	minDist_outfile.write("\n")
minDist_outfile.write("\n")
minDist_outfile.write("\t\t\t\t")
for k in range(0,len(clade_list)):
		clade = clade_list[k]
		minDist_outfile.write("\t"+clade)
minDist_outfile.write("\n")

dist_dict = {}
flagged_accessions = []
for i in range(0,len(query_list)):
	accession = query_list[i]
	dist_dict[accession] = {}
	dist_list = []
	N_found = False
	for k in range(0,len(clade_list)):
		clade = clade_list[k]
		dist = 0
		for j in range(0,len(diff_snp_locs)):
			loc = diff_snp_locs[j]
			if loc not in exclude_locs_for_dist:
				s_nt = clade_diff_snps[clade][loc]
				try:
					q_nt = query_snp_dict[accession][loc]
				except:
					q_nt = "N"
				if q_nt != s_nt and q_nt != "-" and q_nt != "N":
					dist += 1
				elif q_nt == "N":
					N_found = True
		dist_dict[accession][clade] = dist
		tup = (dist,clade)
		dist_list.append(tup)

	SNP_profile = ''
	for j in range(0,len(diff_snp_locs)):
		loc = diff_snp_locs[j]
		try:
			q_nt = query_snp_dict[accession][loc]
		except:
			q_nt = "N"
		SNP_profile += q_nt

	dist_list = sorted(dist_list)
	min_dist_val = dist_list[0][0]
	min_dist_clade = dist_list[0][1]

	minDist_all_outfile.write(accession+"\t"+str(min_dist_val)+"\t"+str(min_dist_clade)+"\t"+SNP_profile+"\n")

	if (min_dist_val >= 2 and N_found == False) or accession in accessions_to_include:

		parent_clades = sorted(check_parent_clades(clade_diff_snps,diff_snp_locs,query_snp_dict[accession]))
		min_break_count = 99
		if parent_clades != []:
			for p in range(0,len(parent_clades)):
				tup = parent_clades[p] #(break_count,len(all_loc), pair, all_loc)
				clade1 = tup[2].split(",")[0]
				clade2 = tup[2].split(",")[1]
				clade_dist1 = dist_dict[accession][clade1]
				clade_dist2 = dist_dict[accession][clade2]
				break_count = tup[0]
				num_supporting_SNPs = tup[1]
				parent_clade_outfile.write(accession+"\t"+clade1+"\t"+clade2+"\t"+str(clade_dist1)+"\t"+str(clade_dist2)+"\t"+str(num_supporting_SNPs)+"\t"+str(break_count)+"\n")
				if break_count <= min_break_count:
					min_break_count = break_count
		elif parent_clades == []:
			min_break_count = 'nan'
		
		if min_break_count != 'nan':
			flagged_accessions.append(accession)
			minDist_outfile_string = ''
			minDist_outfile_string += str(min_dist_val)+"\t"+str(min_dist_clade)+"\t"+str(min_break_count)+"\t"
			for x in range(0,len(diff_snp_locs)):
				loc = diff_snp_locs[x]
				q_nt = query_snp_dict[accession][loc]
				minDist_outfile_string += q_nt
			for k in range(0,len(clade_list)):
				clade = clade_list[k]
				gaped_string = ''
				for x in range(0,len(diff_snp_locs)):
					loc = diff_snp_locs[x]
					q_nt = query_snp_dict[accession][loc]
					s_nt = clade_diff_snps[clade][loc]
					if q_nt != s_nt:
						gaped_string += s_nt
					else:
						gaped_string += '-'
				minDist_outfile_string += '\t'+gaped_string
			minDist_outfile_string += '\n'
			try:
				name_list = name_dict[accession]
				for b in range(0,len(name_list)):
					accession_name = name_list[b]
					minDist_outfile.write(accession_name+"\t"+minDist_outfile_string)
			except:
				minDist_outfile.write(accession+"\t"+minDist_outfile_string)

	elif min_dist_val == 0 and N_found == False:
		try:
			name_list = name_dict[accession]
			for b in range(0,len(name_list)):
				accession_name = name_list[b]
				zeroDist_outfile.write(accession_name+"\t"+str(min_dist_clade)+"\n")
		except:
			zeroDist_outfile.write(accession+"\t"+str(min_dist_clade)+"\n")
minDist_outfile.close()
minDist_all_outfile.close()
zeroDist_outfile.close()

outfile = open(project_dir+output_prefix+".seqs_to_characterize.fasta","w")
for accession in flagged_accessions:
	outfile.write(">"+accession+"\n"+trim_seq_dict[accession]+"\n")
outfile.close()
