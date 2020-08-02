import os
###################################   USER DEFINED VARIABLES   ###################################

project_dir = "./gisaid_2020/"
blast_dir = project_dir + "blast/"
output_prefix = 'gisaid_hcov19'

trimmmed_alignment_filename = project_dir+output_prefix+".filtered.unique.align.fa"
clade_diff_snp_loc_filename = project_dir+"clade_associated_SNPs.nucl.txt"


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

strain_to_group_dict = {}
group_to_strain_dict = {}
group_list = []
group_infile = open(project_dir+"ref_align_groups.txt","r")
for line in group_infile:
	line = line.strip().split("\t")
	strain_to_group_dict[line[0]] = line[1]
	try:
		group_to_strain_dict[line[1]].append(line[0])
	except:
		group_to_strain_dict[line[1]] = []
		group_to_strain_dict[line[1]].append(line[0])
		group_list.append(line[1])
group_infile.close()
group_list = sorted(list(set(group_list)))

msa_infile = open(trimmmed_alignment_filename,"r")
msa_dict = {}
accession_list = []
for line in msa_infile:
	line = line.strip()
	if len(line) >0:
		if line[0] == ">":
			accession = line[1:len(line)]
			accession_list.append(accession)
		else:
			line = line.replace("a","A")
			line = line.replace("c","C")
			line = line.replace("g","G")
			line = line.replace("t","T")
			line = line.replace("n","N")
			try:
				msa_dict[accession] += line
			except:
				msa_dict[accession] = line
msa_in_length = len(msa_dict[accession])
print("Done storing MSA")
print(str(msa_in_length)+" columns")

name_dict = {}
name_infile = open(project_dir+output_prefix+".filtered.names","r")
for line in name_infile:
	line = line.strip().split("\t")
	name_list = line[1].split(",")
	for i in range(0,len(name_list)):
		try:
			name_dict[line[0]].append(name_list[i])
		except:
			name_dict[line[0]] = []
			name_dict[line[0]].append(name_list[i])
name_infile.close()

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
		clade_diff_snps[loc] = {}
		diff_snp_locs.append(loc)
		for j in range(1,len(line)):
			clade = clade_list[j-1]
			nt = line[j]
			clade_diff_snps[loc][clade] = nt
clade_diff_snp_loc_file.close()
print(diff_snp_locs)

query_snp_dict = {}
query_list = []
for i in range(0,len(diff_snp_locs)):
	loc = diff_snp_locs[i]
	query_snp_dict[loc] = {}
	for accession in msa_dict:
		nt = msa_dict[accession][loc]
		query_snp_dict[loc][accession] = nt
		query_list.append(accession)
query_list = sorted(list(set(query_list)))
print(str(len(query_list))+" query sequences to search")

minDist_outfile = open(project_dir+"min_diffSNP_distance.txt","w")
minDist_all_outfile = open(project_dir+"min_diffSNP_distance.all.txt","w")
zeroDist_outfile = open(project_dir+"zero_diffSNP_distance.clade_assignment.txt","w")
for clade in clade_list:
	minDist_outfile.write(clade+"\t")
	for x in range(0,len(diff_snp_locs)):
		loc = diff_snp_locs[x]
		s_nt = clade_diff_snps[loc][clade]
		minDist_outfile.write(s_nt)
	minDist_outfile.write("\n")
minDist_outfile.write("\t\t\t")

minDist_outfile.write("SeqID\tmin_dist_to_clade\tclade\n")
minDist_all_outfile.write("\n")

flagged_accessions = []
for i in range(0,len(query_list)):
	accession = query_list[i]
	dist_list = []
	N_found = False
	for k in range(0,len(clade_list)):
		clade = clade_list[k]
		dist = 0
		for j in range(0,len(diff_snp_locs)):
			loc = diff_snp_locs[j]
			q_nt = query_snp_dict[loc][accession]
			s_nt = clade_diff_snps[loc][clade]
			if q_nt != s_nt and q_nt != "-" and q_nt != "N":
				dist += 1
			elif q_nt == "N":
				N_found = True
		tup = (dist,clade)
		dist_list.append(tup)

	dist_list = sorted(dist_list)
	min_dist_val = dist_list[0][0]
	min_dist_clade = dist_list[0][1]

	minDist_all_outfile.write(accession+"\t"+str(min_dist_val)+"\t"+str(min_dist_clade)+"\n")

	if (min_dist_val >= 2 and N_found == False) or accession in accessions_to_include:
		flagged_accessions.append(accession)
		minDist_outfile.write(accession+"\t"+str(min_dist_val)+"\t"+str(min_dist_clade)+"\t")
		for x in range(0,len(diff_snp_locs)):
			loc = diff_snp_locs[x]
			q_nt = query_snp_dict[loc][accession]
			minDist_outfile.write(q_nt)
		for k in range(0,len(clade_list)):
			clade = clade_list[k]
			gaped_string = ''
			for x in range(0,len(diff_snp_locs)):
				loc = diff_snp_locs[x]
				q_nt = query_snp_dict[loc][accession]
				# minDist_outfile.write(q_nt)
				s_nt = clade_diff_snps[loc][clade]
				if q_nt != s_nt:
					gaped_string += s_nt
				else:
					gaped_string += '-'
			minDist_outfile.write('\t'+gaped_string)
		minDist_outfile.write("\n")
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

ref_strains = []
outfile = open("seqs_to_characterize.fasta","w")
for accession in flagged_accessions:
	outfile.write(">"+accession+"\n"+msa_dict[accession]+"\n")
outfile.close()
