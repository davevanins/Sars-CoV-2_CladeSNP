import os
###################################   USER DEFINED VARIABLES   ###################################

project_dir = "./gisaid_2020/"
output_prefix = 'gisaid_hcov19'

input_alignment_filename = "hcov19.reference_alignment.trimmed.fa"
clade_key_filename = "ref_align_groups.txt"

min_focus_prop = 0.95
max_other_prop = 0.05

############################################ FUNCTIONS ############################################


############################################   MAIN   #############################################
strain_to_group_dict = {}
group_to_strain_dict = {}
group_list = []
group_infile = open(project_dir+clade_key_filename,"r")
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


msa_infile = open(project_dir+input_alignment_filename,"r")
msa_dict = {}
accession_list = []
skip_accesion = True
temp_accession = ''
for line in msa_infile:
	line = line.strip()
	if len(line) >0:
		if line[0] == ">":
			accession = line[1:len(line)]
			try:
				group = strain_to_group_dict[accession]
				accession_list.append(accession)
				skip_accesion = False
				temp_accession = accession
			except:
				skip_accesion = True
		elif skip_accesion == False:
			line = line.replace("a","A")
			line = line.replace("c","C")
			line = line.replace("g","G")
			line = line.replace("t","T")
			line = line.replace("n","N")
			try:
				msa_dict[accession] += line
			except:
				msa_dict[accession] = line
print(len(msa_dict))
msa_in_length = len(msa_dict[temp_accession])
print("Done storing MSA")
print(str(msa_in_length)+" columns")

nt_outfile = open(project_dir+"clade_associated_SNPs.low_prop.nucl.txt","w")
prop_outfile = open(project_dir+"clade_associated_SNPs.low_prop.prop.txt","w")
for z in range(0,len(group_list)):
	group = group_list[z]
	nt_outfile.write("\t"+group)
	prop_outfile.write("\t"+group)
nt_outfile.write("\n")
prop_outfile.write("\n")
print(group_list)

column_stats = open(project_dir+output_prefix+".align.trim.column_stats.txt","w")
column_stats.write("\tA\tT\tC\tG\n")

nt_list = ['A','C','G','T']
nt_dict = {}
group_count = {}
indicator_sites = []
for i in range(0,msa_in_length):
	valid_character_count = 0
	nt_dict[i] = {}
	group_count[i] = {}
	for accession in msa_dict:
		try:
			group = strain_to_group_dict[accession]
			nt = msa_dict[accession][i]
			if nt != "N" and nt != "-":
				try:
					group_count[i][group] += 1
				except:
					group_count[i][group] = 1
				try:
					nt_dict[i][group][nt] += 1
				except:
					# nt_dict[i][group] = {'A':0,'C':0,'G':0,'T':0,'N':0,'-':0}
					nt_dict[i][group] = {'A':0,'C':0,'G':0,'T':0}
					nt_dict[i][group][nt] += 1
			if nt != "-":
				valid_character_count += 1
		except:
			pass
	prop_valid = float(valid_character_count)/float(len(accession_list))
	if prop_valid >= 0.95:
		include_column = False
		for nt in nt_list:
			focus_include = False
			other_include = False
			for group in nt_dict[i]:
				nt_prop = float(nt_dict[i][group][nt])/float(group_count[i][group])
				# prop_dict[group] = nt_prop
				if nt_prop >= min_focus_prop:
					focus_include = True
				elif nt_prop <= max_other_prop:
					other_include = True
			if focus_include == True and other_include == True:
				include_column = True
		if include_column == True:
			major_nt_dict = {}
			for group in group_list:
				max_count = 0
				max_nt = ''
				for nt in nt_list:
					if nt_dict[i][group][nt] > max_count:
						max_nt = nt
						max_count = nt_dict[i][group][nt]
				major_nt_dict[group] = max_nt
			indicator_sites.append(i)
			nt_outfile.write(str(i))
			prop_outfile.write(str(i))
			for z in range(0,len(group_list)):
				group = group_list[z]
				major_nt = major_nt_dict[group]
				nt_prop = float(nt_dict[i][group][major_nt])/float(group_count[i][group])
				nt_outfile.write("\t"+str(major_nt))
				prop_outfile.write("\t"+str(nt_prop))
			nt_outfile.write("\n")
			prop_outfile.write("\n")
	nt_dict_sum = {'A':0,'C':0,'G':0,'T':0}
	for group in nt_dict[i]:
		for nt in nt_list:
			nt_count = nt_dict[i][group][nt]
			try:
				nt_dict_sum[nt] += nt_count
			except:
				nt_dict_sum[nt] = nt_count
	column_stats.write(str(i)+"\t"+str(nt_dict_sum["A"])+"\t"+str(nt_dict_sum["T"])+"\t"+str(nt_dict_sum["C"])+"\t"+str(nt_dict_sum["G"])+"\n")
column_stats.close()
nt_outfile.close()
prop_outfile.close()
