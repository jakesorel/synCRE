import pandas as pd
import numpy as np
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *

n_RNA_cluster = len(os.listdir("results/expression/archetypes"))
eCRE_names = os.listdir("reference/eCRE_locs")
for eCREname in eCRE_names:
    eCREname = eCREname.split(".bed")[0]
    make_directory("results/archetype_beds/%s/by_cluster" % (eCREname))
    for RNA_cluster in range(n_RNA_cluster):
        archetypes = np.loadtxt("results/expression/archetypes/archetypes_for_cluster_%d.txt"%RNA_cluster,dtype=np.int64)
        print(archetypes)
        first = False

        for i in archetypes:
            print("results/archetype_beds/%s/by_archetype/%s_archetype_%d.bed"%(eCREname,eCREname,i))
            try:
                if first is False:

                    df = pd.read_csv("results/archetype_beds/%s/by_archetype/%s_archetype_%d.bed"%(eCREname,eCREname,i), sep="\t", header=None)
                    first = True
                else:

                    df = pd.concat([df,pd.read_csv("results/archetype_beds/%s/by_archetype/%s_archetype_%d.bed"%(eCREname,eCREname,i),sep="\t",header=None)])
            except:
                a = 1
        df.to_csv("results/archetype_beds/%s/by_cluster/%s_cluster_%d.bed"%(eCREname,eCREname,i),sep="\t",index=False,header=False)
