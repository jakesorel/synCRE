import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb

motifs = pd.read_csv("results/archetypes_intersected.csv",index_col=0)
df = pd.read_csv("expression/RNA_normCounts_filter1.csv",index_col=0)
lookup = pd.read_csv("expression/motif_names_manual.csv")
cluster_ids = np.arange(8,dtype=np.int64)
for cluster_id in cluster_ids:
    cluster = list(pd.read_csv("clustered_genes/cluster %d.csv"%cluster_id,index_col=0,header=0)["0"].values)

    ###Find archetype ids
    ar_ids = []
    cid_list = []
    for gene in cluster:
        cid = lookup.loc[lookup["Motif"] == gene]["Cluster_Ids"].values
        for cidd in cid:
            cid_list.append(cidd)
        if cid.size == 0:
            cid = -1
        ar_ids.append(cid)

    cid_list = sorted(list(set(cid_list)))
    # locs = []
    # for cid in cid_list:
    #     lcs = motifs.loc[motifs["Cluster_ID"]==cid]["Location"]
    #     for lc in lcs:
    #         locs.append(lc)
    #
    # start,end = np.zeros(len(locs),dtype=np.int64),np.zeros(len(locs),np.int64)
    # for i, loc in enumerate(locs):
    #     s,e = loc.split(",")
    #     start[i],end[i] = int(s),int(e)
    #
    # nstart,nend = start - start.min(),end-start.min()
    #
    # genome = np.zeros((nend.max()+1))
    # for i in range(start.size):
    #     genome[nstart[i]:nend[i]+1]  = 1
    # #
    # # plt.imshow(genome.reshape(-1,1),aspect="auto")
    # # plt.show()

    expr = []
    ecids = []
    for i, cid in enumerate(cid_list):
        genes = list(lookup.loc[lookup["Cluster_Ids"]==cid]["Motif"].values)
        exprs = np.zeros((len(genes),65))
        for j, gene in enumerate(genes):
            ids = np.where(np.array(df.index == gene))[0]
            if ids.size>1:
                print("no")
            for id in ids:
                exprs[j] = df.iloc[id].values
        ciD = np.ones_like(exprs,dtype=np.int64)*cid
        expr.append()
    fig, ax = plt.subplots(figsize=(10,3))
    expr_df = pd.DataFrame(expr.T)
    expr_df.columns = cid_list
    expr_df = expr_df.melt()
    expr_df.columns = ["Archetype ID","RNA-seq expression"]
    sb.boxplot(data=expr_df,x = "Archetype ID",y="RNA-seq expression",ax = ax,order = np.array(cid_list)[np.argsort(np.median(expr,axis=1))])
    ax.set_xticks(np.arange(len(cid_list)))
    ax.set_xticklabels(labels = np.array(cid_list)[np.argsort(np.median(expr,axis=1))],rotation=90,fontsize=6)
    ax.set_title("Cluster %d"%cluster_id)
    fig.tight_layout()
    fig.savefig("clustered_genes/expression/cluster %d.pdf"%cluster_id)

