import os
import pandas as pd
import json
import gzip, shutil

def make_directory(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def save_csv(data,filename,**kwargs):
    gene_df = pd.DataFrame(data)
    gene_df.to_csv(filename,**kwargs)

def capitalize_list(list):
    return [val.capitalize() for val in list]

def save_dict(dictionary,filename):
    with open('%s.json'%filename, 'w') as fp:
        json.dump(dictionary, fp)

def open_dict(filename):
    with open('%s.json'%filename, 'r') as fp:
        dictionary = json.load(fp)
    return dictionary

def standardize_name(name,aliases,true_names=None):
    if true_names is None:
        true_names = aliases.values()
    if name not in true_names:
        try:
            return aliases[name]
        except KeyError:
            return "*" + name
    else:
        return name

def standardize_names(names,aliases,true_names=None):
    return [standardize_name(name, aliases,true_names) for name in names]


def unzip(file,outfile):
    with gzip.open(file, 'r') as f_in, open(outfile, 'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)