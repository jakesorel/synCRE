import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
os.chdir("..")
from reformat_functions.functions import *



def get_gene_id(df,gene):
    return np.nonzero(df.index.values == gene)[0][0]

def z_score(expr):
    return (expr - np.nanmean(expr))/np.nanstd(expr)


def plot_expression_profile(zmap,title,ax,cmap=plt.cm.viridis,vmin=-2,vmax=2):
    ax.set_title(title)
    ax.imshow(np.flip(zmap,axis=0),vmin=vmin,vmax=vmax,cmap=cmap)
    ax.set_yticks(np.arange(5))
    ax.set_yticklabels(np.flip(np.array(["NMP","p3","pMN","p2","p1"])))
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(["D3","D4","D5","D6"],rotation=90)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmax=vmax, vmin=vmin))
    sm._A = []
    cl = plt.colorbar(sm, ax=ax, pad=0.05, fraction=0.085, aspect=10, orientation="vertical")
    cl.set_label("z-score")

def plot_expression_profile_gene(df,gene,expr_mat,ax,cmap=plt.cm.viridis,vmin=-2,vmax=2):
    gene_id = get_gene_id(df,gene)
    gene_expr = expr_mat[:,:,:,:,gene_id]
    gene_expr_z = z_score(gene_expr)
    zmap = np.nanmean(gene_expr_z,axis=(2,3))
    plot_expression_profile(zmap, gene, ax, cmap=cmap, vmin=vmin, vmax=vmax)


def kmeans(expr_mat_z):
    empty_mask = ~np.isnan(expr_mat_z.reshape(-1, expr_mat_z.shape[-1])).all(axis=1)
    flat_expr_mat_z = expr_mat_z.reshape(-1, expr_mat_z.shape[-1])[empty_mask].T
    kmeans = KMeans(n_clusters=8, random_state=0).fit(flat_expr_mat_z)
    gene_states = kmeans.labels_
    return kmeans,gene_states

if __name__ == "__main__":

    ##Set up directories
    make_directory("results")
    make_directory("results/expression")
    make_directory("results/expression/plots")
    make_directory("results/expression/clusters")
    make_directory("results/expression/all_genes")
    make_directory("results/expression/archetypes")

    ##Load data
    aliases = open_dict("lookup_table/aliases")
    true_names = aliases.values()
    lookup = open_dict("lookup_table/lookup_table")
    gene_names = list(lookup.keys())
    # gene_names = list(pd.read_csv("reference/TF_list.csv",header=None).values.ravel())
    # df = pd.read_csv("reference/RNA_seq/RNA_normCounts_filter1.csv",index_col=0)
    df = pd.read_csv(open("reference/RNA_seq.txt").read(),index_col=0)

    #Standardise names
    df.index = standardize_names(df.index,aliases,true_names)

    #Filter by TFs represented in lookup table
    col_names = list(df.columns)
    genes = list(df.index.values)
    TF_mask = np.zeros(len(genes),dtype=np.bool)
    TF_mask = [gene in gene_names for gene in genes]
    df = df.iloc[np.nonzero(TF_mask)]
    TF_genes = list(df.index.values)

    #Generate matrix from RNA-seq dataframe
    conditions = col_names
    condition_split = np.array([condition.split(".") for condition in conditions])
    condition_split[:,0] = np.array([condition.split("D")[1] for condition in condition_split[:,0]])
    condition_split[:,3] = np.array([condition.split("R")[1] for condition in condition_split[:,3]])
    condition_split[condition_split[:,2]=="NMP",2],condition_split[condition_split[:,2]=="3",2],condition_split[condition_split[:,2]=="M",2],condition_split[condition_split[:,2]=="2",2],condition_split[condition_split[:,2]=="1",2] = 0,10,20,30,40
    condition_split = condition_split.astype(np.float)
    condition_split = np.column_stack((condition_split[:,2],condition_split[:,0],condition_split[:,1],condition_split[:,3]))

    pos_mat,day_mat,conc_mat,rep_mat = np.meshgrid(np.array((0,10,20,30,40)),np.array((3,4,5,6)),np.array((0,10,100,500)),np.array((1,2,3)),indexing="ij")
    expr_mat = np.ones((pos_mat.shape[0],pos_mat.shape[1],pos_mat.shape[2],pos_mat.shape[3],df.shape[0]),dtype=np.float64)*np.nan
    pos_id,day_id,conc_id,rep_id = np.zeros(condition_split.shape[0],dtype=np.int64),np.zeros(condition_split.shape[0],dtype=np.int64),np.zeros(condition_split.shape[0],dtype=np.int64),np.zeros(condition_split.shape[0],dtype=np.int64)
    for i in range(condition_split.shape[0]):
        a,b,c,d = np.nonzero((condition_split[i,0] == pos_mat)*(condition_split[i,1]==day_mat)*(conc_mat==condition_split[i,2])*(rep_mat==condition_split[i,3]))
        pos_id[i] = a[0]
        day_id[i] = b[0]
        conc_id[i] = c[0]
        rep_id[i] = d[0]
        expr_mat[a[0],b[0],c[0],d[0]] = df[conditions[i]]


    ##Plot some candidate genes
    fig, ax = plt.subplots(4,2,sharey=True,sharex=True)
    plot_genes = ["Nkx2-2","Olig2","Pax6","Irx3","Gli3","Sox2","Cdx1","Cdx2"]
    ax = ax.ravel()
    for i, gene in enumerate(plot_genes):
        im = plot_expression_profile_gene(df,gene,expr_mat,ax[i])
    fig.subplots_adjust(hspace=0.5,wspace=0)
    fig.savefig("results/expression/plots/candidate_genes.pdf")

    ##Collapse expression by SAG concentration and repeat, leaving FACS gate and time
    mean_expr_mat = np.nanmean(expr_mat,axis=(2,3))

    ##Remove genes whose max expression is below a threshold (10th percentile of non-zero expression levels)
    percentile_thresh = 10
    thresh = np.nanpercentile(mean_expr_mat[~np.isnan(mean_expr_mat)*(mean_expr_mat!=0)],percentile_thresh)
    print("thresh = ",thresh)
    expressed_mask = (mean_expr_mat>percentile_thresh).any(axis=(0,1))
    expr_mat_reduced = mean_expr_mat[:,:,expressed_mask]

    ##Calculate the fold change (min/max) and remove those that are above thresh (0.25)
    fold_change = np.nanmin(expr_mat_reduced,axis=(0,1))/np.nanmax(expr_mat_reduced,axis=(0,1))
    fc_thresh = 0.25
    changed_mask = fold_change < fc_thresh
    expr_mat_reduced = expr_mat_reduced[:,:,changed_mask]

    ##Calculate z-score across conditions
    expr_mat_z = np.dstack([z_score(expr_mat_reduced[:,:,i]) for i in range(expr_mat_reduced.shape[-1])])

    ##Perform k-means clustering
    kmeans,gene_states = kmeans(expr_mat_z)

    ##Plot k-means clustering
    fig, ax = plt.subplots(4,2,sharex=True,sharey=True)
    ax = ax.ravel()
    for j,i in enumerate(np.unique(gene_states)):
        av_expr = np.nanmean(expr_mat_z[:,:,gene_states == i],axis=-1)
        plot_expression_profile(av_expr,"Cluster %d"%i,ax[j])
    fig.subplots_adjust(hspace=0.5,wspace=0)
    fig.savefig("results/expression/plots/expression_clusters.pdf")

    ##Save lists of genes: all genes, and by cluster
    save_csv(list(np.array(TF_genes)[expressed_mask][changed_mask]),
             "results/expression/all_genes/filtered_genes.txt",header=False,index=False)
    for i in np.unique(gene_states):
        save_csv(list(np.array(TF_genes)[expressed_mask][changed_mask][gene_states == i]),"results/expression/clusters/cluster %d.txt"%i,header=False,index=False)

    ##Save unique archetypes for each expression cluster
    for i in np.unique(gene_states):
        gene_list = list(np.array(TF_genes)[expressed_mask][changed_mask][gene_states == i])
        cid_list = []
        for gene in gene_list:
            try:
                cid_list.append(lookup[gene])
            except:
                a = 1
        save_csv(np.array(np.unique(cid_list)),"results/expression/archetypes/archetypes_for_cluster_%d.txt"%i,header=False,index=False)
