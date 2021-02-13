#!/bin/bash


#################
#Set up environment
##################
# cd setup
# ./setup.sh
# cd .. 

##################
#Make lookup table
###################
# cd reference
# ./get_gene_info.sh
# cd .. 
eval "$(conda shell.bash hook)"
source activate synCRE
python lookup_table/make_alias_dict.py
python lookup_table/make_lookup.py
echo "Made lookup table"

##################
#Cluster RNA expression
###################
python expression_clustering/RNA_seq_cluster.py
echo "Clustered RNAseq expression"




