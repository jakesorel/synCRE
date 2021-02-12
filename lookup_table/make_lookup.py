#!/usr/bin/python
"""
This script makes a lookup file to get a cluster id from a gene-name.

Uses the aliases.json file which converts all names to the standard names

Also uses manual_additions.csv file, where the remaining unknown gene names are manually curated

Saves the lookup table to a csv and json file.
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *

def split_doubles(motif_names,cluster_ids):
    motif_names_old,cluster_ids_old = motif_names.copy(),cluster_ids.copy()
    motif_names,cluster_ids = [],[]
    for i, name in enumerate(motif_names_old):
        if len(name.split("+")) >1:
            names = name.split("+")
            for name in names:
                motif_names.append(name)
                cluster_ids.append(cluster_ids_old[i])
        else:
            motif_names.append(name)
            cluster_ids.append(cluster_ids_old[i])
    return motif_names,cluster_ids

def fix_Znf(motif_names):
    for i, motif in enumerate(motif_names):
        if ("Zn" in motif)&("Znf" not in motif):
            motif_names[i] = "Znf"+motif.split("Zn")[1]
    return motif_names


def fix_manual_additions(motif_names,manual_addition_file="lookup_table/manual_additions.csv"):
    manual_df = pd.read_csv(manual_addition_file,index_col=0)
    manual_dict = dict(zip(manual_df.Name, manual_df.New_name))
    orig_names = manual_df.Name.values
    for i, motif in enumerate(motif_names):
        if motif in orig_names:
            motif_names[i] = manual_dict[motif]
            print(motif,manual_dict[motif])
    return motif_names


def de_duplicate(motif_names,cluster_ids):
    motif_names_old,cluster_ids_old = np.array(motif_names),np.array(cluster_ids)
    motif_names,cluster_ids = [],[]
    for cluster_id in np.unique(cluster_ids_old):
        id_mask = cluster_id == cluster_ids_old
        unique_names = np.unique(motif_names_old[id_mask])
        for name in unique_names:
            motif_names.append(name)
            cluster_ids.append(cluster_id)
    return motif_names,cluster_ids

if __name__ == "__main__":

    ##Open reference datasets
    aliases = open_dict("lookup_table/aliases")
    true_names = aliases.values()
    TF_names = list(pd.read_csv("reference/TF_list.csv",header=None)[0].values)
    TF_names = standardize_names(TF_names,aliases,true_names)

    ##Compare with the RNA-seq list
    RNA_names = list(pd.read_csv(open("reference/RNA_seq.txt").read(),index_col=0).index)
    RNA_names = standardize_names(RNA_names,aliases,true_names)
    print("these TFs are not found in the RNAseq dataset",set(TF_names).difference(set(RNA_names)))


    ##Open archetype dataset
    archetype_clusters = pd.read_excel("reference/motif_annotations.xlsx", 0,engine='openpyxl')
    motifs = pd.read_excel("reference/motif_annotations.xlsx", 1,engine='openpyxl')

    cluster_ids = motifs["Cluster_ID"].values
    motif_names = motifs["Motif"]

    motif_names,cluster_ids = split_doubles(motif_names,cluster_ids)
    for i, name in enumerate(motif_names):
        motif_names[i] = name.split("_")[0]
        motif_names[i] = motif_names[i].split(".")[0].capitalize()


    motif_names = fix_Znf(motif_names)

    motif_names = standardize_names(motif_names,aliases,true_names)

    motif_names = fix_manual_additions(motif_names)

    motif_names,cluster_ids = de_duplicate(motif_names,cluster_ids)

    make_directory("lookup_table/not_referenced")

    ###Save lookup table
    # save_csv({"Cluster_ID":cluster_ids,"Motif_name":motif_names},"lookup_table/lookup_table.csv")
    save_dict(dict(zip(motif_names,[int(cid) for cid in cluster_ids])),"lookup_table/lookup_table")

    ###Save not_referenced
    save_csv(sorted(list(set(motif_names).difference(set(RNA_names)))),"lookup_table/not_referenced/motifs_not_in_RNA.csv")
    save_csv(sorted(list(set(motif_names).difference(set(TF_names)))),"lookup_table/not_referenced/motifs_not_in_TF_list.csv")
    save_csv(sorted(list(set(TF_names).difference(set(motif_names)))),"lookup_table/not_referenced/TFs_not_in_motif_list.csv")


    # save_csv({"Original_name":list(pd.read_csv("expression/RNA_normCounts_filter1.csv",index_col=0).index),"New_name":RNA_names},"lookup_table/rename RNA-seq.csv")
