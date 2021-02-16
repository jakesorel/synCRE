import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *
# lookup = open_dict("lookup_table/lookup_table")
motif_table = pd.read_excel("reference/motif_annotations.xlsx", 1, engine='openpyxl')
lookup_names = dict(zip(motif_table["Motif"],motif_table["Cluster_ID"]))
motifcsvs = os.listdir("results/motifs/raw")
make_directory("results/motifs/bed")
for motifcsv in motifcsvs:
    if ".csv" in motifcsv:
        df = pd.read_csv("results/motifs/raw/%s"%motifcsv,sep=";",header=None)
        locs = df[0]
        names = df[1].values
        names = [name.split(".pmf")[0].split(" ")[0] for name in names]
        motif_ids = []
        for name in names:
            try:
                motif_ids.append(lookup_names[name])
            except KeyError:
                motif_ids.append(-1)
        chrom,start,end = [],[],[]
        for loc in locs:
            chrm,pos = loc.split(":")
            st,en = pos.split("-")
            st,en = int(st),int(en)
            chrom.append(chrm)
            start.append(st)
            end.append(en)
        start,end = np.array(start),np.array(end)
        begin = df[2].values
        start += begin
        end += begin
        bed_df = pd.DataFrame([chrom,start,end,motif_ids]).transpose()
        bed_df = bed_df.loc[bed_df[3] != -1]
        bed_df.to_csv("results/motifs/bed/%s.bed"%(motifcsv.split(".csv")[0]),sep="\t",header=None,index=None)

