import sys
import numpy as np
bigwigs = """
[Gli3FLAG_IP_R1.mLb.clN]
file=Nishi/Gli3FLAG_IP_R1.mLb.clN.bigWig
title=Nishi Gli3
color = #666666
min_value = 0
#max_value = auto
height = 1.5
number of bins = 500
nans to zeros = True
show data range = yes
file_type = bigwig

[Olig2_IP_R1.mLb.clN.bigWig]
file=Nishi/Olig2_IP_R1.mLb.clN.bigWig
title=Nishi Olig2
color = #666666
min_value = 0
#max_value = auto
height = 1.5
number of bins = 500
nans to zeros = True
show data range = yes
file_type = bigwig

[Nkx22_IP_R1.mLb.clN.bigWig]
file=Nishi/Nkx22_IP_R1.mLb.clN.bigWig
title=Nishi Nkx22
color = #666666
min_value = 0
height = 1.5
number of bins = 500
nans to zeros = True
show data range = yes
file_type = bigwig

[ATAC_D5_500_M_R1.mLb.clN.bigWig]
file=ATAC/D5_500_M_R1.mLb.clN.bigWig
title=ATAC D5 500 M
color = #666666
min_value = 0
height = 1.5
number of bins = 500
nans to zeros = True
show data range = yes
file_type = bigwig
\n
"""

archetype_template = """
[archetype %i]
file=results/o2e33_oosterveen_archetype_%i_clean.bed
title=A%i
color = darkblue
# height = 0.75
# line_width = 0.5
gene_rows = 2
labels = off
file_type = bed
fontsize = 10
style = UCSC
"""

def generate_archetype_script(i):
    out = archetype_template%(i,i,i)
    return out

foot = """[x-axis]
[spacer]"""

RNA_cluster = int(sys.argv[1])
archetypes = np.loadtxt("clustered_genes/archetypes/archetypes_for_cluster_%d.txt"%RNA_cluster,dtype=np.int64)

f = open('genome_plotting/tracks.ini', 'w')
f.write(bigwigs)
for i in archetypes:
    f.write(generate_archetype_script(i))
f.write(foot)
f.close()  # you can omit in most cases as the destructor will call it
