import pandas as pd
import numpy as np
import sys
from shutil import copyfile
file_name = sys.argv[1]
new_name = file_name.split(".bed")[0] + "_clean.bed"
try:
    df = pd.read_csv(file_name,sep="\t",header=None)
    out_df = pd.DataFrame(np.array([df[i] for i in range(3,9)]).T)
    out_df.to_csv(new_name,sep="\t",index=False,header=False)
except:
    copyfile(file_name,new_name)