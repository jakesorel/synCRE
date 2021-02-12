#!/usr/bin/python
"""
This script generates a dictionary of gene name aliases and saves it to a json file found within the directory lookup_table
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *

all_gene_list = pd.read_csv("reference/gene_info",sep="\t")
official_name = capitalize_list(list(all_gene_list["Symbol"]))
all_names = []
for i, synonym in enumerate(all_gene_list["Synonyms"].values):
    if synonym == "-":
        all_names.append([official_name[i]])
    else:
        syn_vals = capitalize_list(synonym.split("|"))
        all_names.append([official_name[i]]+syn_vals)

aliases = {}
for i, names in enumerate(all_names):
    for name in names:
        aliases[name] = official_name[i]

save_dict(aliases,"lookup_table/aliases")
