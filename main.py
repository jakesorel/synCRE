from synCRE import *

##############################
# Make the lookup table
##############################
# lkp = Lookup(RNA_seq_file="reference/RNA_seq.txt")

##############################
# Cluster RNA expression
##############################
# expr = Expression(RNA_seq_file="reference/RNA_seq_local.txt")
# expr.run_all()

##############################
# Find motifs within each eCRE
##############################
mtf = Motif_Finder()
mtf.make_pmf()
mtf.sample_eCRE_sequence()
# mtf.find_motifs([0.0001])
# mtf.motif_to_bed()
# # mtf.motifs_to_bedgraph()
# mtf.motifs_by_archetype(collapse=True)
# # mtf.motifs_by_cluster(make_bedgraph=False)
# mtf.collapse_all_bed()
# mtf.collapse_bed_relevant_clusters(clusters=[1,3,5,6])
# #
# # bed = BedTool("results/motifs/by_archetype/Nkx2-2_p=0.001000/archetype_1.bed")
# # bedmerge = bed.merge().saveas("results/motifs/by_archetype/Nkx2-2_p=0.001000/archetype_1_merge.bed")
# #
# # bed.to_dataframe().shape
# # bedmerge.to_dataframe().shape
# # bed = BedTool("results/motifs/bed/Pax6_p=0.000100.bed")
# # start = bed.to_dataframe()["start"].values
# # print((start[1:] - start[:-1]).min())
# print("motifs found")
#
# ##############################
# # Plot data
# ##############################
# # eCREs_all = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
# # eCREs = []
# # for eCRE in eCREs_all:
# #     if "p=0.000100" in eCRE:
# #         eCREs.append(eCRE)
# eCREs = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
#
# # for eCRE in eCREs:
# eCRE = "Nkx2-2_p=0.000100"
# # eCRE = eCREs[0]
# plot = GenomePlot(eCRE,plot_constructs=True,plot_bw=False,plot_genes=False,plot_phylo=False)
# # plot.ini_all_motifs()
# # plot.ini_by_cluster()
# plot.ini_relevant_clusters()
# # plot.ini_by_cluster_merge()
# # plot.ini_by_candidate()
# plot.make_plots(parallel=True,suppress=True)
# print("""
# ################################################
# Plots for %s complete
# ################################################
# """%eCRE)
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
