import pandas as pd
import numpy as np
import sys
from shutil import copyfile
RNA_cluster = int(sys.argv[1])
archetypes = np.loadtxt("clustered_genes/archetypes/archetypes_for_cluster_%d.txt"%RNA_cluster,dtype=np.int64)
first = False
for i in archetypes:
    try:
        if first is False:
            df = pd.read_csv("results/o2e33_oosterveen_archetype_%i_clean.bed" % i, sep="\t", header=None)
            first = True
        else:
            df = pd.concat([df,pd.read_csv("results/o2e33_oosterveen_archetype_%i_clean.bed"%i,sep="\t",header=None)])
    except:
        a = 1
df.to_csv("results/o2e33_oosterveen_cluster_%d.bed"%RNA_cluster,sep="\t",index=False,header=False)
