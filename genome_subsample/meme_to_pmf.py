"""
This script parses the set of ".meme" files in the meme_files directory.

It takes each pmf matrix and saves it to a separate file in the "pmf" directory. Filenames are assigned by the motif-name.
"""

from io import StringIO
import pandas as pd
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *

make_directory("reference/motifs/pmf")

meme_files = os.listdir("reference/motifs/meme_files")
for memefile in meme_files:
    pmf_on = False
    for i, line in enumerate(open("reference/motifs/meme_files/%s"%memefile).readlines()):
        if "MOTIF" in line:
            motif_name = line.split("MOTIF ")[1].split("\n")[0]
            motifs = []
        if (("URL" in line)or(line == "\n"))and(pmf_on==True):
            pmf_on = False
            df = pd.DataFrame(motifs).transpose()
            df.to_csv("reference/motifs/pmf/%s.pmf"%motif_name,sep="\t",header=False,index=False)
        if pmf_on is True:
            motif = pd.read_csv(StringIO(line),delim_whitespace=True,header=None).values.ravel()
            motifs.append(motif)
        if "letter-probability matrix" in line:
            pmf_on = True

# jaspar_files_all = os.listdir("jaspar")
# jaspar_files = []
# for file in jaspar_files_all:
#     if ".jaspar" in file:
#         jaspar_files.append(file)
# for filenm in jaspar_files:
#     expt_name = filenm.split(".jaspar")[0] + ".pmf"
#     f = open("jaspar/%s"%expt_name,"w")
#     for i, line in enumerate(open("jaspar/%s"%filenm).readlines()):
#         if i!=0:
#             f.write(line[6:-2])
#             f.write("\n")
#     f.close()
#     os.remove("jaspar/%s"%filenm)
#
