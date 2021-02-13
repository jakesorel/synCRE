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
# python lookup_table/make_alias_dict.py
# python lookup_table/make_lookup.py
# echo "Made lookup table"

##################
#Cluster RNA expression
###################
python expression_clustering/RNA_seq_cluster.py
echo "Clustered RNAseq expression"


##################
#Intersect archetype list with bed_file
###################

ml BEDTools/2.26.0-foss-2016b
bedtools intersect -wb -a enhancer_beds/o2e33_oosterveen.bed -b archetypes/mm10.archetype_motifs.v1.0.chr16.bed > results/o2e33_oosterveen_archetypes.bed




