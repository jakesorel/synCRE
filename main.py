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
# mtf = Motif_Finder()
# mtf.make_pmf()
# mtf.sample_eCRE_sequence()
# mtf.find_motifs([0.001])
# mtf.motif_to_bed()
# mtf.motifs_to_bedgraph()
# mtf.motifs_by_archetype()
# mtf.motifs_by_cluster(make_bedgraph=False)

# bed = BedTool("results/motifs/bed/Pax6_p=0.000100.bed")
# start = bed.to_dataframe()["start"].values
# print((start[1:] - start[:-1]).min())
print("motifs found")

##############################
# Plot data
##############################
eCREs_all = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
eCREs = []
for eCRE in eCREs_all:
    if "p=0.000100" in eCRE:
        eCREs.append(eCRE)
for eCRE in eCREs:
    # eCRE = "Nkx2-2_p=0.001000"
    # eCRE = eCREs[0]
    plot = GenomePlot(eCRE)
    plot.ini_all_motifs()
    # plot.ini_by_cluster()
    # plot.ini_by_cluster_merge()
    # plot.ini_by_candidate()
    plot.make_plots(parallel=True,suppress=False)
    print("""
################################################
Plots for %s complete
################################################
    """%eCRE)


"""
To do:

-- sort bed files - DONE
-- include genes in plots -- DONE
-- add colours to the plots (via the input file)
-- include conservation in plots -- DONE
-- check lookup class works
-- clean directory
-- documentation
-- don't regenerate lookup if exists already, with over-ride

"""
