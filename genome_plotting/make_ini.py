import sys
import numpy as np
import pandas as pd

bigwig_template = """
[%s]
file=%s
title=%s
color = %s
min_value = 0
#max_value = auto
height = %.3f
number of bins = 500
nans to zeros = True
show data range = yes
file_type = bigwig
"""



archetype_template = """
[%s]
file=%s
title=%s
color = %s
height = %.3f
# line_width = 0.5
gene_rows = %d
labels = %s
file_type = bed
fontsize = 10
style = UCSC
"""



foot = """[x-axis]
[spacer]"""

def make_bigwig(name,dir,color="#666",height=1.5):
    return bigwig_template%(name,dir,name,color,height)

def make_bed(name,dir,color="darkblue",height=0.75,gene_rows=2,labels="off"):
    out = archetype_template%(name,dir,name,color,height,gene_rows,labels)
    return out

bigwigs = pd.read_csv("reference/bigwig_files.txt",sep="\t",header=None)
bigwigs.columns = ["name","dir"]

RNA_cluster = int(sys.argv[1])
archetypes = np.loadtxt("clustered_genes/archetypes/archetypes_for_cluster_%d.txt"%RNA_cluster,dtype=np.int64)

f = open('genome_plotting/tracks.ini', 'w')
f.write(bigwigs)
for i in archetypes:
    f.write(generate_archetype_script(i))
f.write(foot)
f.close()  # you can omit in most cases as the destructor will call it
