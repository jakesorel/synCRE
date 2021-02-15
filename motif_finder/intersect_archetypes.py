from pybedtools import BedTool
import os, sys
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *


make_directory("results/archetype_beds")
eCRE_names = os.listdir("reference/eCRE_locs")
motifs = pd.read_excel("reference/motif_annotations.xlsx", 0,engine='openpyxl')

for name in eCRE_names:
    eCRE = BedTool("reference/eCRE_locs/%s"%name)
    chromosome = eCRE.to_dataframe()["chrom"].values[0]
    archetypes = BedTool("reference/archetypes/mm10.archetype_motifs.v1.0.%s.bed"%chromosome)
    eCREname = name.split(".bed")[0]
    make_directory("results/archetype_beds/%s"%eCREname)
    make_directory("results/archetype_beds/%s/all_beds"%eCREname)
    make_directory("results/archetype_beds/%s/by_archetype"%eCREname)
    archetype_intersect = eCRE.intersect(archetypes,wb=True).saveas("results/archetype_beds/%s/all_beds/%s_archetypes.bed"%(eCREname,eCREname))
    print("Ran intersect on %s"%name)
    df = archetype_intersect.to_dataframe()
    out_df = pd.DataFrame(df[df.columns[3:9]])
    out_df.to_csv("results/archetype_beds/%s/all_beds/%s_archetypes_clean.bed"%(eCREname,eCREname),sep="\t",index=False,header=False)
    print("Saved clean version")
    bed_names = out_df[out_df.columns[3]].values
    for i, motif in enumerate(motifs["Name"].values):
        cid = motifs.iloc[i]["Cluster_ID"]
        bed_ids = np.where(bed_names == motif)[0]
        out_df.iloc[bed_ids].to_csv("results/archetype_beds/%s/by_archetype/%s_archetype_%d.bed"%(eCREname,eCREname,i), sep="\t", header=False, index=False)
    print("Saved sorted by archetypes")

