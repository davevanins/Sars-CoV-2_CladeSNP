import os
import random
###################################   USER DEFINED VARIABLES   ###################################

project_dir = "./gisaid_2020/"
blast_dir = project_dir + "blast/"
output_prefix = 'gisaid_hcov19'

output_prefix = 'gisaid_hcov19'
trimmmed_alignment_filename = project_dir+output_prefix+".filtered.unique.align.fa"
clade_diff_snp_loc_filename = project_dir+"clade_associated_SNPs.nucl.txt"

exclude_locs_for_dist = []#if you wish to include positions 28881-3 as one site, change this list to: [28764,28765]
accessions_to_include = [] #this list acts as an override to include low quality sequences in final outputs. Optional.

############################################ FUNCTIONS ############################################

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
			tup = (break_count,len(all_loc), pair, all_loc)
			potential_parents.append(tup)
	return potential_parents

############################################   MAIN   #############################################
try:
	include_infile = open(project_dir+"genomes_to_include.txt","r")
	for line in include_infile:
		line = line.strip()
		accessions_to_include.append(line)
	include_infile.close()
except:
	pass

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

query_snp_dict = {}
query_list = []
for i in range(0,len(diff_snp_locs)):
	loc = diff_snp_locs[i]
	for accession in msa_dict:
		nt = msa_dict[accession][loc]
		try:
			query_snp_dict[accession][loc] = nt
		except:
			query_snp_dict[accession] = {}
			query_snp_dict[accession][loc] = nt
		query_list.append(accession)
query_list = sorted(list(set(query_list)))
print(str(len(query_list))+" query sequences to search")



parent_clade_outfile = open(project_dir+"putative_recombinant_info.txt","w")
parent_clade_outfile.write("accession\tparent_clade1\tparent_clade2\ttclade1_dist\tclade2_dist\tnum_supporting_SNPs\trecombination_break_points\n")
minDist_outfile = open(project_dir+"min_diffSNP_distance.txt","w")
cladeDist_outfile = open(project_dir+"clade_diffSNP_distance.txt","w")
minDist_all_outfile = open(project_dir+"min_diffSNP_distance.all.txt","w")
zeroDist_outfile = open(project_dir+"zero_diffSNP_distance.clade_assignment.txt","w")

for clade in clade_list:
	minDist_outfile.write(clade+"\t")
	for x in range(0,len(diff_snp_locs)):
		loc = diff_snp_locs[x]
		s_nt = clade_diff_snps[clade][loc]
		minDist_outfile.write(s_nt)
	minDist_outfile.write("\n")
minDist_outfile.write("\n")
minDist_outfile.write("\t\t\t\t")
for k in range(0,len(clade_list)):
		clade = clade_list[k]
		minDist_outfile.write("\t"+clade)
		cladeDist_outfile.write("\t"+clade)
minDist_outfile.write("\n")
cladeDist_outfile.write("\n")

dist_dict = {}
flagged_accessions = []
for i in range(0,len(query_list)):
	accession = query_list[i]
	dist_dict[accession] = {}
	cladeDist_outfile.write(accession)

	dist_list = []
	N_found = False
	for k in range(0,len(clade_list)):
		clade = clade_list[k]
		dist = 0
		for j in range(0,len(diff_snp_locs)):
			loc = diff_snp_locs[j]
			if loc not in exclude_locs_for_dist:
				q_nt = query_snp_dict[accession][loc]
				s_nt = clade_diff_snps[clade][loc]
				if q_nt != s_nt and q_nt != "-" and q_nt != "N":
					dist += 1
				elif q_nt == "N":
					N_found = True
		dist_dict[accession][clade] = dist
		tup = (dist,clade)
		dist_list.append(tup)
		cladeDist_outfile.write("\t"+str(dist))

	SNP_profile = ''
	for j in range(0,len(diff_snp_locs)):
		loc = diff_snp_locs[j]
		q_nt = query_snp_dict[accession][loc]
		SNP_profile += q_nt
	
	cladeDist_outfile.write("\n")

	dist_list = sorted(dist_list)
	min_dist_val = dist_list[0][0]
	min_dist_clade = dist_list[0][1]

	minDist_all_outfile.write(accession+"\t"+str(min_dist_val)+"\t"+str(min_dist_clade)+"\t"+SNP_profile+"\n")

	if (min_dist_val >= 2 and N_found == False) or accession in accessions_to_include:
		flagged_accessions.append(accession)

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
cladeDist_outfile.close()
minDist_all_outfile.close()
zeroDist_outfile.close()

ref_strains = []
outfile = open("seqs_to_characterize.fasta","w")
for accession in flagged_accessions:
	outfile.write(">"+accession+"\n"+msa_dict[accession]+"\n")
outfile.close()
