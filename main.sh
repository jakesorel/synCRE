#!/bin/bash


# #################
# #Set up environment
# ##################
# # cd setup
# # ./setup.sh
# # cd .. 

# ##################
# #Make lookup table
# ###################
# # cd reference
# # ./get_gene_info.sh
# # cd .. 
# eval "$(conda shell.bash hook)"
# source activate synCRE
# # python lookup_table/make_alias_dict.py
# # python lookup_table/make_lookup.py
# # echo "Made lookup table"

# ##################
# #Cluster RNA expression
# ###################
# python expression_clustering/RNA_seq_cluster.py
# echo "Clustered RNAseq expression"


##################
#Intersect archetype list with bed_file
###################
python motif_finder/intersect_archetypes.py
python motif_finder/beds_by_expression.py


###################
#Plot results
##################
eval "$(conda shell.bash hook)"
source activate pygenometracks
python genome_plotting/make_ini.py
./results/genome_plots/run_all.sh