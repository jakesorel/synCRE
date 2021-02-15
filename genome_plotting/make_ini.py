import numpy as np
import pandas as pd
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *
from pybedtools import BedTool


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



archetype_template_rows = """
[%s]
file=%s
title=%s
color = %s
height = %.3f
# line_width = 0.5
# gene_rows = %d
labels = %s
file_type = bed
fontsize = 10
style = UCSC
"""

archetype_template_height = """
[%s]
file=%s
title=%s
color = %s
height = %.3f
# line_width = 0.5
# gene_rows = %d
labels = %s
file_type = bed
fontsize = 10
style = UCSC
"""



foot = """[x-axis]
[spacer]"""

runline_template = """
pyGenomeTracks --tracks results/genome_plots/%s/config_files/%s/%s.ini --region %s:%d-%d -o results/genome_plots/%s/plots/%s/%s.pdf
"""

def make_bigwig(name,dir,color="#666",height=1.5):
    return bigwig_template%(name,dir,name,color,height)

def make_bed(name,dir,color="darkblue",height=1.5,gene_rows=None,labels="off"):
    if gene_rows is None:
        out = archetype_template_height%(name,dir,name,color,height,labels)
    else:
        out = archetype_template_rows%(name,dir,name,color,gene_rows,labels)
    return out

def make_runline(eCREname, cat,filename,chr,start,end):
    out = runline_template%(eCREname,cat,filename,chr,start,end,eCREname,cat,filename)
    return out
make_directory("results/genome_plots")
bigwigs = pd.read_csv("reference/bigwig_files.txt",sep="\t",header=None)
bigwigs.columns = ["name","dir"]

eCRE_names = os.listdir("reference/eCRE_locs")
eCRE_names = [name.split(".bed")[0] for name in eCRE_names]




for name in eCRE_names:
    make_directory("results/genome_plots/%s"%name)
    make_directory("results/genome_plots/%s/run_files"%name)
    g = open("results/genome_plots/%s/run_files/run.sh"%name,'w')
    make_directory("results/genome_plots/%s/config_files"%name)
    make_directory("results/genome_plots/%s/config_files/all"%name)
    make_directory("results/genome_plots/%s/config_files/by_cluster"%name)
    make_directory("results/genome_plots/%s/config_files/by_cluster_merge"%name)

    make_directory("results/genome_plots/%s/plots"%name)
    make_directory("results/genome_plots/%s/plots/all"%name)
    make_directory("results/genome_plots/%s/plots/by_cluster"%name)
    make_directory("results/genome_plots/%s/plots/by_cluster_merge"%name)

    # bed_dir = "results/archetype_beds/%s"%name
    e_chrom,e_start,e_end = BedTool("reference/eCRE_locs/%s.bed"%name).to_dataframe().values.ravel()


    # ##Config file for all
    # archetypes = np.loadtxt("results/archetype_beds/%s/by_cluster/archetypes_for_cluster_%d.txt"%RNA_cluster,dtype=np.int64)

    cat = "all"
    filename = "all"
    f = open('results/genome_plots/%s/config_files/%s/%s.ini'%(name,cat,filename), 'w')
    for bwname,bwdir in bigwigs.values:
        print(bwdir)
        f.write(make_bigwig(bwname,bwdir))
    f.write(make_bed("All archetypes","results/archetype_beds/%s/all_beds/%s_archetypes_clean.bed"%(name,name)))
    f.write(foot)
    f.close()  # you can omit in most cases as the destructor will call it
    g.write(make_runline(eCREname = name,cat=cat,filename=filename,chr=e_chrom,start=e_start,end=e_end))
    #
    # cat = "merge"
    # archetype_files = np.loadtxt("results/archetype_beds/%s/by_cluster"%name)
    # for archetype_file in archetype_files:
    #     cluster_no = int((archetype_file.split(".bed")[0]).split("cluster_")[1])
    #     cluster_name = "RNA_cluster_%d"%cluster_no
    #     f = open('results/genome_plots/%s/config_files/%s/%s.ini' % (name, cat, cluster_name), 'w')
    #     make_bed("results/archetype_beds/%s/by_cluster/%s"%(name,archetype_file), dir)
    # eCRE_names = os.listdir("reference/eCRE_locs")

    g.close()
