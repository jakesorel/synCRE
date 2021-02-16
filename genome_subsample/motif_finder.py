import os, sys
import pandas as pd
import numpy as np
import shutil
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *

make_directory("results/motifs")

runline = """
eval "$(conda shell.bash hook)" \n
source activate moods \n
cd pkgs/moods/scripts \n
pwd \n
python2 moods-dna.py  \
--sep ";" -s ../../../results/fasta/by_eCRE/%s.fa --p-value %.6f \
--lo-bg 2.977e-01 2.023e-01 2.023e-01 2.977e-01 \
-m ../../../reference/motifs/pmf/* -o ../../../results/motifs/%s_p=%.6f.csv
"""

def make_runline(eCRE,pval):
    return runline%(eCRE,pval,eCRE,pval)
eCRE_names = [name.split(".bed")[0] for name in os.listdir("reference/eCRE_locs")]
eCRE_names = ["Olig2","Pax6","Nkx2-2"]
p_vals = [0.001,0.0005,0.0001,0.00005]
for p_val in p_vals:
    for name in eCRE_names:
        os.system(make_runline(name,p_val))
