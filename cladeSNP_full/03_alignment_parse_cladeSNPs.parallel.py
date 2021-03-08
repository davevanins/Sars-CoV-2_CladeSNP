import os
import sys
import random
args = sys.argv
num_alignments = int(args[1])
###################################   USER DEFINED VARIABLES   ###################################

project_dir = "/nobackup1b/users/davevan/SARS_COV2/cladeSNP/gisaid_02_16_2021/"
blast_dir = project_dir + "blast/"
output_prefix = 'gisaid_hcov19'

trimmmed_alignment_filename = project_dir+output_prefix+".filtered.unique.align."
names_file = project_dir+output_prefix+".filtered.names"

clade_diff_snp_loc_filename = project_dir+"clade_associated_SNPs.nucl.txt"

exclude_locs_for_dist = []#if you wish to include positions 28881-3 as one site, change this list to: [28764,28765]

min_cdSNP_threshold = 2 #minimum cdSNP distance from nearest clade to be considered

pplacer_reference_name = "hcov_rep.refpkg" #to use pplacer a maximum likelihood reference tree with alignment needs to be generated using "taxit create" 
comparisons_per_job = 600 #this value is for the number of sequences to place on the reference tree using pplacer

slurm_prefix = "#!/bin/bash\n#SBATCH -N 1\n#SBATCH -n 1\n#SBATCH -p sched_mit_chisholm,sched_mit_hill,newnodes\n#SBATCH --mem=40000\n#SBATCH --time=12:00:00\nsource activate python3env\ncd "+project_dir+"\n"

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

			break_site_list = []
			all_loc = sorted(all_loc)
			for i in range(0,len(all_loc)):
				loc = all_loc[i]
				if loc in loc1:
					if current_parent != clade1:
						current_parent = clade1
						break_count += 1
						break_site = int(((all_loc[i]-all_loc[i-1])/2.0)+all_loc[i-1])
						break_site_list.append(break_site)
						# print(str(all_loc[i-1])+"\t"+str(all_loc[i])+"\t"+str(break_site))
				if loc in loc2:
					if current_parent != clade2:
						current_parent = clade2
						break_count += 1
						break_site = int(((all_loc[i]-all_loc[i-1])/2.0)+all_loc[i-1])
						break_site_list.append(break_site)
			break_site_list = sorted(break_site_list)
			tup = (break_count,len(all_loc), pair, all_loc,break_site_list)
			potential_parents.append(tup)
	return potential_parents

def split_seq_by_parent_clade(break_loc_list,full_genome_seq):
	max_loc_stop = 29623
	break_loc_list = list(set(break_loc_list))
	break_loc_list.append(max_loc_stop)
	break_loc_list = sorted(break_loc_list)
	p1_seq = ''
	p2_seq = ''
	loc_start = 0
	clade = 1
	for j in range(0,len(break_loc_list)):
		loc_stop = break_loc_list[j]
		seq = full_genome_seq[loc_start:loc_stop]
		gap = "-"*len(seq)
		if clade == 1:
			p1_seq += seq
			p2_seq += gap
			clade = 2
		elif clade == 2:
			p2_seq += seq
			p1_seq += gap
			clade = 1
		loc_start = loc_stop
	return p1_seq,p2_seq

############################################   MAIN   #############################################
accessions_to_include = [] #this list acts as an override to include low quality sequences in final outputs. Optional.
try:
	include_infile = open(project_dir+"genomes_to_include.txt","r")
	for line in include_infile:
		line = line.strip()
		accessions_to_include.append(line)
	include_infile.close()
except:
	pass

msa_dict = {}
accession_list = []
for alignment_num in range(0,num_alignments):
	msa_infile = open(trimmmed_alignment_filename+str(alignment_num)+".fa","r")
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
name_infile = open(names_file,"r")
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
minDist_outfile.write("accession\tdist_to_nearest_clade\tcdSNP_profile")
for k in range(0,len(clade_list)):
		clade = clade_list[k]
		minDist_outfile.write("\t"+clade)
		cladeDist_outfile.write("\t"+clade)
minDist_outfile.write("\n")
cladeDist_outfile.write("\n")

split_seq_record = {}
dist_dict = {}
pplace_gap_seq_dict = {}
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

	if (min_dist_val >= min_cdSNP_threshold and N_found == False) or accession in accessions_to_include:

		parent_clades = sorted(check_parent_clades(clade_diff_snps,diff_snp_locs,query_snp_dict[accession]))
		if parent_clades != [] and len(msa_dict[accession]) >= 29623:
			flagged_accessions.append(accession)

			accession_full_seq = msa_dict[accession][0:29623]
			pplace_gap_seq_dict[accession+"_full"] = accession_full_seq
			for p in range(0,len(parent_clades)):
				tup = parent_clades[p] #(break_count,len(all_loc), pair, all_loc,break_site_list)
				break_count = tup[0]
				num_supporting_SNPs = tup[1]
				parent_clade_pair = tup[2]
				clade1 = parent_clade_pair.split(",")[0]
				clade2 = parent_clade_pair.split(",")[1]
				clade_dist1 = dist_dict[accession][clade1]
				clade_dist2 = dist_dict[accession][clade2]
				break_sites = tup[4]
				parent_clade_outfile.write(accession+"_s"+str(p)+"\t"+clade1+"\t"+clade2+"\t"+str(clade_dist1)+"\t"+str(clade_dist2)+"\t"+str(num_supporting_SNPs)+"\t"+str(break_count)+"\n")
		
				split_seq_p1,split_seq_p2 = split_seq_by_parent_clade(break_sites,accession_full_seq)
				pplace_gap_seq_dict[accession+"_s"+str(p)+"_p1"] = split_seq_p1
				pplace_gap_seq_dict[accession+"_s"+str(p)+"_p2"] = split_seq_p2

				try:
					split_seq_record[accession] += 1
				except:
					split_seq_record[accession] = 1


			minDist_outfile_string = ''
			minDist_outfile_string += str(min_dist_val)+"\t"
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
outfile = open(project_dir+"putative_recombinant_seqs.fasta","w")
for accession in flagged_accessions:
	outfile.write(">"+accession+"\n"+msa_dict[accession][0:29623]+"\n")
outfile.close()

a = -1
s = -1
for sequence_name in pplace_gap_seq_dict:
	sequence_string = pplace_gap_seq_dict[sequence_name]
	a += 1
	if a%comparisons_per_job == 0:
		s += 1
		print("s: "+str(s)+"\ta: "+str(a))
		sequence_outfile_filename = "recomb_split."+str(s)+".fasta"
		sequence_outfile = open(project_dir+sequence_outfile_filename,"w")
		sequence_outfile.write(">"+sequence_name+"\n"+sequence_string+"\n")
		
		command = "pplacer -p -c "+pplacer_reference_name+" "+sequence_outfile_filename+"\n"

		sh = open(project_dir+"pplace_run."+str(s)+".sh","w")
		sh.write(slurm_prefix+command+"\n")
		sh.close()
	else:
		sequence_outfile.write(">"+sequence_name+"\n"+sequence_string+"\n")
sequence_outfile.close()

for i in range(0,s+1):
	command = "sbatch pplace_run."+str(i)+".sh"
	os.system(command)
	
