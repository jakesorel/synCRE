from synCRE import *
# ##############################
# # Make the lookup table
# ##############################
# # lkp = Lookup(RNA_seq_file="reference/RNA_seq.txt")
#
# ##############################
# # Cluster RNA expression
# ##############################
# expr = Expression(RNA_seq_file="reference/RNA_seq_local.txt")
# # expr.run_all()
#
# ##############################
# # Find motifs within each eCRE
# ##############################
mtf = Motif_Finder()
# # mtf.make_pmf()
# # mtf.sample_eCRE_sequence()
# # mtf.find_motifs(4)
# # mtf.make_truth_matrix()
# #
# required_dicts = {"Olig2": {"Olig2": 1, "Nkx2-2": 1, "Gli3": 1,"Sox2":1},
#                   "Nkx2-2": {"Olig2": 1, "Gli3": 1,"Sox2":0},
#                   "Pax6": {"Olig2": 1, "Nkx2-2": 1,"Sox2":1}}
# mtf.filter_motifs_by_hit(required_dicts)
# mtf.plot_motif_distributions()
mtf.motif_to_bed()
# # # mtf.motifs_to_bedgraph()
mtf.motifs_by_archetype(collapse=True)
mtf.motifs_by_cluster(make_bedgraph=False)
mtf.collapse_all_bed()
mtf.collapse_bed_relevant_clusters(clusters=[0,3,4,5])
# print("motifs found")
#
##############################
# Plot data
##############################
eCREs = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
eCREs = os.listdir("results/motifs/bed")

max_val_dict={"Olig2":25,"Nkx2-2":20,"Pax6":20}

for eCRE in eCREs:
    if ".bed" in eCRE:
        eCRE = eCRE.split(".bed")[0]
        plot = GenomePlot(eCRE,plot_constructs=False,plot_bw=False,plot_genes=False,plot_phylo=False)
        # plot.ini_all_motifs()
        # plot.ini_atac(max_value=max_val_dict[eCRE])
        # plot.ini_by_cluster()
        plot.ini_relevant_clusters(split_clusters=False)
        plot.ini_by_cluster_merge(together=True)
        # plot.ini_by_candidate()
        plot.make_plots(parallel=True,suppress=False,width=10)
        print("""
        ################################################
        Plots for %s complete
        ################################################
        """%eCRE)
#
#
# """
# To do:
#
# -- sort bed files - DONE
# -- include genes in plots -- DONE
# -- add colours to the plots (via the input file)
# -- include conservation in plots -- DONE
# -- check lookup class works
# -- clean directory
# -- documentation
# -- don't regenerate lookup if exists already, with over-ride
#
# """
#
# self = expr
#
# self.motif_annotations = "reference/motif_annotations.xlsx"
#
#
# def make_motif_name_dict(self):
#     """
#
#     :return:
#     """
#     self.motifs = pd.read_excel(self.motif_annotations, 1, engine='openpyxl')
#     self.motif_name_dict = dict(zip(self.motifs["Motif"],self.motifs["Cluster_ID"]))
# make_motif_name_dict(self)
#
# fp_dfs = []
# for file in os.listdir("reference/TOBIAS"):
#     df = pd.read_excel("reference/TOBIAS/%s"%file,0, engine='openpyxl')
#     df = df[list([df.columns[2]]) + list(df.columns[6:56:2])]
#     cluster_ids = np.ones(df.shape[0],dtype=np.int64)*-1
#     for i, idd in enumerate(df["motif_id"]):
#         try:
#             cluster_ids[i] = self.motif_name_dict[idd]
#         except:
#             print(idd, "not in reference")
#     df["Cluster_ID"] = cluster_ids
#     fp_dfs.append(df)
#
# fp_df = pd.concat(fp_dfs)
#
# # fp_df = pd.read_csv("reference/TOBIAS/bindetect_results.txt",sep="\t")
#
# fp_df = fp_df.loc[fp_df["Cluster_ID"]!=-1] #filter unmatched
#
# arc_fp_df = []
# for i in range(1,287):
#     arc_fp_df.append(fp_df.loc[fp_df["Cluster_ID"]==i].max(axis=0)) #take the mean value for a given archetype
#
# fp_df = pd.concat(arc_fp_df,axis=1).transpose()
# fp_df.index = np.arange(1,287)
#
# fp_columns = fp_df.columns
#
# pos_range, day_range, conc_range = np.array(("NMP", "3", "M", "2", "1")), np.array((3, 4, 5, 6)),np.array((0, 10, 100, 500))
# fp_mat = np.ones((pos_range.shape[0], day_range.shape[0], conc_range.shape[0], fp_df.shape[0]),
#                    dtype=np.float64) * np.nan
# for i, pos in enumerate(pos_range):
#     for j, day in enumerate(day_range):
#         for k, conc in enumerate(conc_range):
#             name = "D%d_%d_%s_mean_score"%(day,conc,pos)
#             if name in fp_columns:
#                 fp_mat[i,j,k] = fp_df[name]
#
# fp_mat_mean = np.nanmean(fp_mat,axis=2)
# fp_cv = np.array([np.nanstd(fp_mat_mean[:,:,i])/np.nanmean(fp_mat_mean[:,:,i]) for i in range(fp_mat_mean.shape[2])])
# thresh = 0.005
# # fp_mat_mean = fp_mat_mean.T[fp_cv>thresh].T
# fp_z = np.dstack((self.z_score(fp_mat_mean[:,:,i]) for i in range(fp_mat_mean.shape[2])))
#
# empty_mask = ~np.isnan(
#     fp_z.reshape(-1, fp_z.shape[-1])).all(axis=1)
# flat_fp_z = fp_z.reshape(-1, fp_z.shape[-1])[
#     empty_mask].T
#
# self.n_clusters = 11
# self.fp_kmeans = KMeans(n_clusters=self.n_clusters, random_state=0).fit(flat_fp_z)
# self.fp_gene_states = self.fp_kmeans.labels_
#
# fig, ax = plt.subplots(int(np.ceil(self.n_clusters/2)), 2, sharex=True, sharey=True)
# ax = ax.ravel()
# cmap=plt.cm.viridis
# vmin=-2
# vmax=2
# for j, i in enumerate(np.unique(self.fp_gene_states)):
#     av_expr = np.nanmean(fp_z[:, :, self.fp_gene_states == i], axis=-1)
#     self.plot_expression_profile(av_expr, "Cluster %d" % i, ax[j], cmap=cmap, vmin=vmin, vmax=vmax)
# fig.subplots_adjust(hspace=0.5, wspace=0)
# fig.show()
#
#
# fp_mat_mean = np.nanmean(fp_mat,axis=(1,2))
# fp_z = np.column_stack((self.z_score(fp_mat_mean[:,i]) for i in range(fp_mat_mean.shape[1])))
#
# import seaborn as sb
# sb.clustermap(fp_z,row_cluster=False,z_score=False)
# plt.show()
#
# make_directory("results/footprinting")
# make_directory("results/footprinting/z_by_archetype")
# for i in range(286):
#     fig, ax = plt.subplots()
#     cmap=plt.cm.viridis
#     vmin=-2
#     vmax=2
#     self.plot_expression_profile(fp_z[:,:,i], "Archetype %d" % i, ax, cmap=cmap, vmin=vmin, vmax=vmax)
#     fig.subplots_adjust(hspace=0.5, wspace=0)
#     fig.savefig("results/footprinting/z_by_archetype/%d.pdf"%i)
#
