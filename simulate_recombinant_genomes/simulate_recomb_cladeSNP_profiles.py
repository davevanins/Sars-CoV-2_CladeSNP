import os
import random
import numpy as np
import scipy.stats
random.seed(2020)
###################################   USER DEFINED VARIABLES   ###################################
project_dir = "./"
clade_diff_snp_loc_filename = project_dir+"clade_associated_SNPs.nucl.txt"

number_of_break_points = 5 #max 100, any number greater than 100 will assign cdSNP identity randomly as if there are an infinite number of breakpoints
min_dist_from_nearest_clade = 2

hist_start = 0
hist_stop = 30000
bin_size = 10 #bases

number_of_replicates = 30
number_of_genomes_to_simulate = 257

exclude_locs_for_dist = []

clade_proportion = {'19A-1':0.0344,'19B-2':0.00957,'19A-2':0.00323,'20A-3':0.225,'19A-4':0.0403,'20B-2':0.363,'20B-1':0.0421,'19B-3':0.0185,'20A-1':0.0336,'19B-4':0.00970,'20C-1':0.184,'20A-2':0.0237,'19A-3':0.00883,'19B-1':0.00469}

############################################ FUNCTIONS ############################################
def mean_confidence_interval(data, confidence=0.95):
	a = 1.0 * np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
	return m, m-h, m+h

def get_random_clade_pairs(clade_prop,num_to_sim):
	clade_pairs = []
	for clade1 in clade_prop:
		for clade2 in clade_prop:
			if clade1 != clade2:
				prop1 = clade_prop[clade1]
				prop2 = clade_prop[clade2]
				prob = prop1*prop2
				add_pairs = int(prob*num_to_sim)
				for i in range(0,add_pairs):
					val = random.randint(0,10000000000)
					tup = (val,clade1,clade2)
					clade_pairs.append(tup)
	clade_pairs = sorted(clade_pairs)
	return clade_pairs

def get_random_locations(num_breaks):
	rand_positions = []
	for i in range(0,num_breaks):
		novel = False
		while novel == False:
			val = random.randint(0,29741)
			if val not in rand_positions:
				rand_positions.append(val)
				novel = True
	rand_positions.append(0)
	rand_positions.append(29741)
	rand_positions = sorted(rand_positions)
	return rand_positions

def count_clade_dist(clade_list,cladeSNP_locations,clade_diff_snps,query_snp_dict):
	dist_list = []
	for k in range(0,len(clade_list)):
		clade = clade_list[k]
		dist = 0
		for j in range(0,len(cladeSNP_locations)):
			loc = cladeSNP_locations[j]
			q_nt = query_snp_dict[loc]
			s_nt = clade_diff_snps[clade][loc]
			if q_nt != s_nt:
				dist += 1
		tup = (dist,clade)
		dist_list.append(tup)
		dist_list = sorted(dist_list)
	return dist_list

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
		all_supporting_locs = []
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
					all_supporting_locs.append(loc)
				elif base2 == q_base:
					loc2.append(loc)
					all_supporting_locs.append(loc)

		if conflict == False and len(loc1)>0 and len(loc2)>0:#(len(loc1)+len(loc1)) > 0:
			
			#count number of breaks to explain recombination	
			break_count = 0
			if min(loc1) < min(loc2):
				current_parent = clade1
			else:
				current_parent = clade2

			all_supporting_locs = sorted(all_supporting_locs)
			for i in range(0,len(all_supporting_locs)):
				loc = all_supporting_locs[i]
				if loc in loc1:
					if current_parent != clade1:
						current_parent = clade1
						break_count += 1
				elif loc in loc2:
					if current_parent != clade2:
						current_parent = clade2
						break_count += 1
			number_of_suporting_cladeSNPs = len(all_supporting_locs)
			tup = (break_count,number_of_suporting_cladeSNPs,pair,all_supporting_locs)
			potential_parents.append(tup)
	potential_parents = sorted(potential_parents)
	return potential_parents


############################################   MAIN   #############################################
if number_of_break_points > 100:
	breakpoint_naming_value = "infinite"
else:
	breakpoint_naming_value = str(number_of_break_points)

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
		loc = int(line[0])+118
		if loc not in exclude_locs_for_dist:
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


### Calculate real overlapping ranges
real_cdSNP_profile_dict = {}
infile = open("min_diffSNP_distance.txt","r")
first_line = True
for line in infile:
	if first_line == True:
		first_line = False
	else:
		line = line.strip().split("\t")
		accession = line[0]
		cdSNP_profile = line[4]
		for i in range(0,len(cdSNP_profile)):
			loc = diff_snp_locs[i]
			cdSNP = cdSNP_profile[i]
			try:
				real_cdSNP_profile_dict[accession][loc] = cdSNP
			except:
				real_cdSNP_profile_dict[accession] = {}
				real_cdSNP_profile_dict[accession][loc] = cdSNP
infile.close()

real_bin_count = {}
for accession in real_cdSNP_profile_dict:
	clade_dist_list = sorted(count_clade_dist(clade_list,diff_snp_locs,clade_diff_snps,real_cdSNP_profile_dict[accession]))
	parent_clades = check_parent_clades(clade_diff_snps,diff_snp_locs,real_cdSNP_profile_dict[accession]) #(break_count,number_of_suporting_cladeSNPs,pair,all_supporting_locs)

	if parent_clades != []:
		min_clade_dist = clade_dist_list[0][0]
		min_observed_breakpoints = parent_clades[0][0]
		if min_clade_dist >= min_dist_from_nearest_clade:
			SNP_profile = ''
			for p in range(0,len(diff_snp_locs)):
				loc = diff_snp_locs[p]
				SNP_profile += real_cdSNP_profile_dict[accession][loc]

			parents = parent_clades[0][2].split(",")
			clade1_obs = parents[0]
			clade2_obs = parents[1]
			diff_locs = sorted(parent_clades[0][3])
			prev_parent = ""
			breaks = []
			for i in range(0,len(diff_locs)):
				loc = diff_locs[i]
				base1 = clade_diff_snps[clade1_obs][loc]
				base2 = clade_diff_snps[clade2_obs][loc]
				qbase = real_cdSNP_profile_dict[accession][loc]
				if prev_parent == "":
					if qbase == base1:
						prev_parent = "one"
					elif qbase == base2:
						prev_parent = "two"
					else:
						print("ERROR")
				else:
					if qbase == base1:
						if prev_parent == "one":
							pass
						elif prev_parent == "two":
							tup = (diff_locs[i-1],diff_locs[i])
							breaks.append(tup)
							prev_parent = "two"
					elif qbase == base1:
						if prev_parent == "two":
							pass
						elif prev_parent == "one":
							tup = (diff_locs[i-1],diff_locs[i])
							breaks.append(tup)
							prev_parent = "one"
			start = 0
			for stop in range(hist_start+bin_size,hist_stop+1,bin_size):
				bin = str(start)+"\t"+str(stop)
				try:
					temp = real_bin_count[bin]
				except:
					real_bin_count[bin] = 0
				for point in breaks:
					if start > point[0] >= stop:
						real_bin_count[bin] += 1
					elif start > point[1] >= stop:
						real_bin_count[bin] += 1
					elif start > point[0] and stop < point[1]:
						real_bin_count[bin] += 1
				start = stop

real_bin_count_list = []
for bin in real_bin_count:
	real_bin_count_list.append(real_bin_count[bin])
median_real_count = np.median(real_bin_count_list)
rel_real_bin_count = {}
for bin in real_bin_count:
	rel_real_bin_count[bin] = float(real_bin_count[bin])



### Simulate recombinant genomes
simulated_genome_summary_outfile = open("simulated_genome_summary.breakpoints-"+str(breakpoint_naming_value)+"_"+str(number_of_replicates)+".txt","w")
simulated_genome_summary_outfile.write("clade1\tclade2\tcdSNP_profile\tmin_dist_val\tmin_observed_breakpoints\n")

bin_count = {}
rel_bin_count = {}
passing_recombinant_count = 0
rand_clade_pairs = get_random_clade_pairs(clade_proportion,number_of_genomes_to_simulate*100)

for rep in range(0,number_of_replicates):
	bin_count[rep] = {}
	rel_bin_count[rep] = {}
	passing_recombinant_count = 0
	while passing_recombinant_count < number_of_genomes_to_simulate:
		for num in range(0,len(rand_clade_pairs)):
			if passing_recombinant_count < number_of_genomes_to_simulate:
				clade1 = rand_clade_pairs[num][1]
				clade2 = rand_clade_pairs[num][2]

				pair = clade1+","+clade2
				pair_r = clade2+","+clade1

				if number_of_break_points <= 100:
					break_points = get_random_locations(number_of_break_points)
				elif number_of_break_points > 100:
					break_points = []

				query_snp_dict = {}
				c1_locs = []
				c2_locs = []
				if random.randint(0,1) == 1: #randomize which clade the first recombination block orginates from
					current_clade = clade1
				else:
					current_clade = clade2

				if number_of_break_points <= 100:
					for k in range(0,len(break_points)-1): #go through recombination blocks and assign SNPs to each parental clade
						start = break_points[k]
						stop = break_points[k+1]
						for p in range(0,len(diff_snp_locs)):
							loc = diff_snp_locs[p]
							if start <= loc < stop:
								if current_clade == clade1:
									c1_locs.append(loc)
									query_snp_dict[loc] = clade_diff_snps[clade1][loc]
								else:
									c2_locs.append(loc)
									query_snp_dict[loc] = clade_diff_snps[clade2][loc]

						if current_clade == clade1:
							current_clade = clade2
						elif current_clade == clade2:
							current_clade = clade1
				elif number_of_break_points > 100:
					for p in range(0,len(diff_snp_locs)):
						if random.randint(0,1) == 1: #randomize which clade the first recombination block orginates from
							current_clade = clade1
						else:
							current_clade = clade2
						loc = diff_snp_locs[p]
						if current_clade == clade1:
							c1_locs.append(loc)
							query_snp_dict[loc] = clade_diff_snps[clade1][loc]
						else:
							c2_locs.append(loc)
							query_snp_dict[loc] = clade_diff_snps[clade2][loc]

				clade_dist_list = sorted(count_clade_dist(clade_list,diff_snp_locs,clade_diff_snps,query_snp_dict))
				parent_clades = check_parent_clades(clade_diff_snps,diff_snp_locs,query_snp_dict) #(break_count,number_of_suporting_cladeSNPs,pair,all_supporting_locs)

				if parent_clades != []:
					min_clade_dist = clade_dist_list[0][0]
					min_observed_breakpoints = parent_clades[0][0]
					if min_clade_dist >= min_dist_from_nearest_clade:
						passing_recombinant_count += 1
						SNP_profile = ''
						for p in range(0,len(diff_snp_locs)):
							loc = diff_snp_locs[p]
							SNP_profile += query_snp_dict[loc]

						parents = parent_clades[0][2].split(",")
						clade1_obs = parents[0]
						clade2_obs = parents[1]
						diff_locs = sorted(parent_clades[0][3])
						prev_parent = ""
						breaks = []
						for i in range(0,len(diff_locs)):
							loc = diff_locs[i]
							base1 = clade_diff_snps[clade1_obs][loc]
							base2 = clade_diff_snps[clade2_obs][loc]
							qbase = query_snp_dict[loc]
							if prev_parent == "":
								if qbase == base1:
									prev_parent = "one"
								elif qbase == base2:
									prev_parent = "two"
								else:
									print("ERROR")
							else:
								if qbase == base1:
									if prev_parent == "one":
										pass
									elif prev_parent == "two":
										tup = (diff_locs[i-1],diff_locs[i])
										breaks.append(tup)
										prev_parent = "two"
								elif qbase == base1:
									if prev_parent == "two":
										pass
									elif prev_parent == "one":
										tup = (diff_locs[i-1],diff_locs[i])
										breaks.append(tup)
										prev_parent = "one"

						start = 0
						for stop in range(hist_start+bin_size,hist_stop+1,bin_size):
							bin = str(start)+"\t"+str(stop)
							try:
								temp = bin_count[rep][bin]
							except:
								bin_count[rep][bin] = 0
							for point in breaks:
								if start > point[0] >= stop:
									bin_count[rep][bin] += 1
								elif start > point[1] >= stop:
									bin_count[rep][bin] += 1
								elif start > point[0] and stop < point[1]:
									bin_count[rep][bin] += 1
							start = stop

						simulated_genome_summary_outfile.write(clade1_obs+"\t"+clade2_obs+"\t"+SNP_profile+"\t"+str(min_clade_dist)+"\t"+str(min_observed_breakpoints)+"\n")

		bin_count_list = []
		for bin in bin_count[rep]:
			bin_count_list.append(bin_count[rep][bin])
		sim_median_count = np.median(bin_count_list)
		for bin in real_bin_count:
			try:
				rel_bin_count[rep][bin] = rel_real_bin_count[bin]/(float(bin_count[rep][bin]))
			except:
				rel_bin_count[rep][bin] = 0.0
		if passing_recombinant_count < number_of_genomes_to_simulate:
			print("Ran out of clade pairs! Cannot complete analysis")
simulated_genome_summary_outfile.close()

range_outfile = open("ranges.break_points-"+str(breakpoint_naming_value)+"_"+str(number_of_replicates)+".txt","w")
start = 0
for stop in range(hist_start+bin_size,hist_stop+1,bin_size):
	bin = str(start)+"\t"+str(stop)
	rel_bin_count_list = []
	for rep in range(0,number_of_replicates):
		rel_bin_count_list.append(rel_bin_count[rep][bin])
	mean_val, lower_conf, upper_conf = mean_confidence_interval(rel_bin_count_list)
	range_outfile.write(str(bin)+"\t"+str(mean_val)+"\t"+str(lower_conf)+"\t"+str(upper_conf)+"\n")
	start = stop
range_outfile.close()
