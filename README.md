
# SARS-CoV-2 CladeSNP tools
This repository contains all of the code to replicate the analsysis of VanInsberghe, Neish, Lowen, and Koelle (2020) "Identification of SARS-CoV-2 recombinant genomes".


# Contents

`cladeSNP_full`
* This pipeline is the main analysis tool used in our manuscript to identify recombinant genomes.

`cladeSNP_blast`
* This is a lightweight version of the full pipeline used in our analysis that is suitable for identifying recombinant genomes it is considerably faster and uses a fraction of the resources as the cladeSNP_full pipeline. 

`simulate_recombinant_genomes`
* This section is used to simulate recombinant genomes to test if there are hotspots for recombination breakpoints in the SARS-CoV-2 genome.

`calculate_ceiling`
* This section is used to calculate the theoretical maximum number of circulating recombinant viruses in a geographic region.
