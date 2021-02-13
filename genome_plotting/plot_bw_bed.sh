#!/bin/bash

# #
# for i in $(seq 1 286)
# do
#    filename=results/o2e33_oosterveen_archetype_$i.bed
#    echo ${filename}
#    python reformat_functions/clean_bed.py "${filename}"
# done

eval "$(conda shell.bash hook)"
source activate pygenometracks

for i in $(seq 0 7)
do
    python genome_plotting/make_ini.py "$i"
    pyGenomeTracks --tracks genome_plotting/tracks.ini --region chr16:91191337-91193823 -o genome_plotting/bigwig_with_genes_${i}.pdf
    python reformat_functions/compile_archetype_bed.py "$i"
    python genome_plotting/make_ini_cluster.py "$i"
    pyGenomeTracks --tracks genome_plotting/tracks.ini --region chr16:91191337-91193823 -o genome_plotting/bigwig_with_genes_merge_${i}.pdf
done
source activate synCRE

# # bed_names=""
# # cat clustered_genes/archetypes/archetypes_for_cluster_0.txt | while read i
# # do
# #     filename=results/o2e33_oosterveen_archetype_${i}_clean.bed
# #     echo ${filename}
# #     bed_names=${filename}
# # done
# # echo "-------------"
# # echo ${bed_names}
# #
# ##sortBed -i results/o2e33_oosterveen_archetypes.bed > results/o2e33_oosterveen_archetypes_sorted.bed
# ##python reformat_functions/clean_bed.py "results/o2e33_oosterveen_archetypes_sorted.bed"
# #make_tracks_file --trackFiles Nishi/Gli3FLAG_IP_R1.mLb.clN.bigWig ${bed_names} -o tracks.ini
# pyGenomeTracks --tracks tracks.ini --region chr16:91191337-91193823 -o bigwig_with_genes.png

