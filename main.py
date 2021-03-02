from synCRE import *
# ##############################
# # Make the lookup table
# ##############################
# # lkp = Lookup(RNA_seq_file="reference/RNA_seq.txt")
#
# ##############################
# # Cluster RNA expression
# ##############################
# # expr = Expression(RNA_seq_file="reference/RNA_seq_local.txt")
# # expr.run_all()
#
# ##############################
# # Find motifs within each eCRE
# ##############################
mtf = Motif_Finder()
# # mtf.make_pmf()
# # mtf.sample_eCRE_sequence()
# # mtf.find_motifs(0.001)
#
# required_dicts = {"Olig2": {"Olig2": 1, "Nkx2-2": 1, "Gli3": 1},
#                   "Nkx2-2": {"Olig2": 1, "Nkx2-2": 1, "Gli3": 1},
#                   "Pax6": {"Olig2": 1, "Nkx2-2": 1}}
# mtf.filter_motifs_by_hit(required_dicts)
# mtf.plot_motif_distributions()
# mtf.motif_to_bed()
# # # mtf.motifs_to_bedgraph()
# mtf.motifs_by_archetype(collapse=True)
# mtf.motifs_by_cluster(make_bedgraph=False)
# mtf.collapse_all_bed()
mtf.collapse_bed_relevant_clusters(clusters=[0,3,4,5])
# print("motifs found")
#
##############################
# Plot data
##############################
eCREs = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
eCREs = os.listdir("results/motifs/bed")

for eCRE in eCREs:
    if ".bed" in eCRE:
        eCRE = eCRE.split(".bed")[0]
        plot = GenomePlot(eCRE,plot_constructs=True,plot_bw=True,plot_genes=True,plot_phylo=True)
        # plot.ini_all_motifs()
        # plot.ini_by_cluster()
        plot.ini_relevant_clusters()
        # plot.ini_by_cluster_merge()
        # plot.ini_by_candidate()
        plot.make_plots(parallel=True,suppress=True)
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
