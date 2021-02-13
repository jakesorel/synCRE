from pybedtools import BedTool
from reformat_functions.functions import *
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)

make_directory("results/archetype_beds")
eCRE_names = os.listdir("reference/eCRE_locs")
for name in eCRE_names:
    eCRE = BedTool("reference/eCRE_locs/%s"%name)
    chromosome = eCRE.to_dataframe()["chrom"].values[0]
    archetypes = BedTool("reference/archetypes/mm10.archetype_motifs.v1.0.%s.bed"%chromosome)
    eCREname = name.split(".bed")[0]
    make_directory("results/archetype_beds/%s"%eCREname)
    make_directory("results/archetype_beds/%s/all_beds"%eCREname)
    eCRE.intersect(archetypes,wb=True).saveas("results/archetype_beds/%s/all_beds/%s_archetypes.bed"%(eCREname,eCREname))
    print("Ran intersect on %s"%name)
