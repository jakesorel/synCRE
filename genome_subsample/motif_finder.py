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
"python2 moods-dna.py  \
--sep ";" -s ../../../results/fasta/by_eCRE/%s.fa --p-value %.6f \
--lo-bg 2.977e-01 2.023e-01 2.023e-01 2.977e-01 \
-m ../../../reference/motifs/pmf/* -o ../../../results/motifs/%s_p=%.6f.csv"
"""

def make_runline(eCRE,pval):
    return runline%(eCRE,pval,eCRE,pval)
eCRE_names = [name.split(".bed")[0] for name in os.listdir("reference/eCRE_locs")]
os.system("echo pwd")
for name in eCRE_names:
    os.system(make_runline(name,0.0001))
