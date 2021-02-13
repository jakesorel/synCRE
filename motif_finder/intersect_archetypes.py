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
for name in eCRE_names:
    eCRE = BedTool("reference/eCRE_locs/%s"%name)
    chromosome = eCRE.to_dataframe()["chrom"].values[0]
    archetypes = BedTool("reference/archetypes/mm10.archetype_motifs.v1.0.%s.bed"%chromosome)
    eCREname = name.split(".bed")[0]
    make_directory("results/archetype_beds/%s"%eCREname)
    make_directory("results/archetype_beds/%s/all_beds"%eCREname)
    archetype_intersect = eCRE.intersect(archetypes,wb=True).saveas("results/archetype_beds/%s/all_beds/%s_archetypes.bed"%(eCREname,eCREname))
    print("Ran intersect on %s"%name)
    out_df = pd.DataFrame(np.array([archetype_intersect[i] for i in range(3,9)]).T)
    out_df.to_csv("results/archetype_beds/%s/all_beds/%s_archetypes_clean.bed"%(eCREname,eCREname),sep="\t",index=False,header=False)
    print("Saved clean version")
