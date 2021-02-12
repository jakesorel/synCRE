import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import SimpSOM as sps
import seaborn as sb
from sklearn.cluster import KMeans
from functions import *

# gene_names = pd.read_csv("expression/lookup_motif_gene_archetype_20-12-08.csv")
gene_names = pd.read_csv("expression/TF_list.csv",header=None)
gene_names.columns = ["mouse_genename"]
# gene_names = gene_names.iloc[np.array([~np.isnan(name) if type(name) is float else True for name in gene_names["mouse_genename"].values])]
df = pd.read_csv("expression/RNA_normCounts_filter1.csv",index_col=0)

col_names = list(df.columns)

genes = list(df.index.values)
TF_mask = np.zeros(len(genes),dtype=np.bool)
TF_mask = [gene in gene_names["mouse_genename"].values for gene in genes]
def check_reference_list(gene_list, TF_reference_list):
    print(np.array([gene in gene_list for gene in TF_reference_list]).mean()*100,"%")

check_reference_list(genes,gene_names["mouse_genename"].values)

df = df.iloc[np.nonzero(TF_mask)]
TF_genes = list(df.index.values)


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


def get_gene_id(df,gene):
    return np.nonzero(df.index.values == gene)[0][0]

def z_score(expr):
    return (expr - np.nanmean(expr))/np.nanstd(expr)

def lfc_score(expr):
    max_expr = np.nanmax(expr)
    return np.log(expr/max_expr + 1)

def plot_expression_profile(df,gene,expr_mat,ax):
    # fig, ax = plt.subplots()
    gene_id = get_gene_id(df,gene)
    gene_expr = expr_mat[:,:,:,:,gene_id]
    gene_expr_z = z_score(gene_expr)
    ax.set_title(gene)
    ax.imshow(np.flip(np.nanmean(gene_expr_z,axis=(2,3)),axis=0))
    ax.set_yticks(np.arange(5))
    ax.set_yticklabels(np.flip(np.array(["NMP","p3","pMN","p2","p1"])))
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(["D3","D4","D5","D6"])
    # fig.show()

fig, ax = plt.subplots(3,2,sharey=True,sharex=True)
plot_genes = ["Nkx2-2","Olig2","Pax6","Irx3","Sox4","Sox2"]
ax = ax.ravel()
for i, gene in enumerate(plot_genes):
    plot_expression_profile(df,gene,expr_mat,ax[i])
fig.show()


mean_expr_mat = np.nanmean(expr_mat,axis=(2,3))

percentile_thresh = 10
thresh = np.nanpercentile(mean_expr_mat[~np.isnan(mean_expr_mat)*(mean_expr_mat!=0)],percentile_thresh)
print("thresh = ",thresh)
expressed_mask = (mean_expr_mat>percentile_thresh).any(axis=(0,1))

expr_mat_reduced = mean_expr_mat[:,:,expressed_mask]

fold_change = np.nanmin(expr_mat_reduced,axis=(0,1))/np.nanmax(expr_mat_reduced,axis=(0,1))
#
# plt.hist(fold_change)
# plt.show()
changed_mask = fold_change < 0.25

expr_mat_reduced = expr_mat_reduced[:,:,changed_mask]
expr_mat_z = np.dstack([z_score(expr_mat_reduced[:,:,i]) for i in range(expr_mat_reduced.shape[-1])])


empty_mask = ~np.isnan(expr_mat_z.reshape(-1,expr_mat_z.shape[-1])).all(axis=1)
flat_expr_mat_z = expr_mat_z.reshape(-1,expr_mat_z.shape[-1])[empty_mask].T


from sklearn.cluster import KMeans

kmeans = KMeans(n_clusters=12,random_state=0).fit(flat_expr_mat_z)

Sum_of_squared_distances = []
K = range(1,20)
for k in K:
    km = KMeans(n_clusters=k,random_state= 0)
    km = km.fit(flat_expr_mat_z)
    Sum_of_squared_distances.append(km.inertia_)
plt.plot(Sum_of_squared_distances)
plt.show()

kmeans = KMeans(n_clusters=8,random_state=0).fit(flat_expr_mat_z)

# net = sps.somNet(4,4,flat_expr_mat_z,PBC=True,n_jobs=-1)
# net.train(0.01,10000)
# net.save('filename_weights')
#
# prj=np.array(net.project(flat_expr_mat_z))
#
# prj = prj[:,0] + 1j*prj[:,1]
# unique_states = np.unique(prj)
# gene_states = np.zeros((prj.shape[0]),dtype=np.int64)
# for i, prj_i in enumerate(prj):
#     gene_states[i] = np.nonzero(prj_i == unique_states)[0]

gene_states = kmeans.labels_
fig, ax = plt.subplots(2,4,sharex=True,sharey=True)
ax = ax.ravel()
for j,i in enumerate(np.unique(gene_states)):
    av_expr = np.nanmean(expr_mat_z[:,:,gene_states == i],axis=-1)
    ax[j].set_title(i)
    ax[j].imshow(np.flip(av_expr,axis=0),vmin=-2,vmax=2,cmap=plt.cm.plasma)
    ax[j].set_yticks(np.arange(5))
    ax[j].set_yticklabels(np.flip(np.array(["NMP","p3","pMN","p2","p1"])))
    ax[j].set_xticks(np.arange(4))
    ax[j].set_xticklabels(["D3","D4","D5","D6"])
fig.savefig("clustered_genes/heatmap.pdf",dpi=300)
#
# fig,ax = plt.subplots(2,4,sharex=True,sharey=True)
# ax = ax.ravel()
# for i, axx in enumerate(ax):
#     axx.set_title(i)
#     av_expr = expr_mat_z[:,:,gene_states == i]
#     # z_scores = np.row_stack([z_score(np.nanmean(expr_mat_reduced[:,:,i],axis=1)) for i in range(expr_mat_reduced.shape[-1])])
#     # axx.violinplot(z_scores)
#     axx.violinplot(np.nanmean(av_expr,axis=1).reshape(-1,av_expr.shape[-1]).T)
#     axx.set_xticks(np.arange(1,6))
#     axx.plot((0.5,5.5),(0,0),linestyle="--",color="grey")
#     axx.set_xticklabels(np.array(["NMP","p3","pMN","p2","p1"]),rotation=90)
#     axx.set(xlim=(0.5,5.5))
# fig.savefig("clustered_genes/violinplot.pdf",dpi=300)
#
# gene_states[np.where(np.array(TF_genes)[expressed_mask][changed_mask] == "Cdx2")[0][0]]

for i in np.unique(gene_states):
    gene_df = pd.DataFrame(list(np.array(TF_genes)[expressed_mask][changed_mask][gene_states == i]))
    gene_df.to_csv("clustered_genes/cluster %d.csv"%i)

make_directory("clustered_genes/archetypes")

lookup = open_dict("lookup_table/lookup_table")
for i in np.unique(gene_states):
    gene_list = list(np.array(TF_genes)[expressed_mask][changed_mask][gene_states == i])
    cid_list = []
    for gene in gene_list:
        try:
            cid_list.append(lookup[gene])
        except:
            a = 1
    save_csv(np.array(np.unique(cid_list)),"clustered_genes/archetypes/archetypes_for_cluster_%d.txt"%i,header=False,index=False)
