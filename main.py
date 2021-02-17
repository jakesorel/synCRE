from synCRE import *

lkp = Lookup(RNA_seq_file="reference/RNA_seq.txt")

# expr = Expression(RNA_seq_file="reference/RNA_seq_local.txt")
# expr.run_all()

# mtf = Motif_Finder()
# mtf.make_pmf()
# mtf.sample_eCRE_sequence()
# mtf.find_motifs([0.001])
# mtf.motif_to_bed()
# mtf.motifs_by_archetype()
# mtf.motifs_by_cluster()

print("motifs found")

eCREs_all = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
eCREs = []
for eCRE in eCREs_all:
    if "p=0.000100" in eCRE:
        eCREs.append(eCRE)
for eCRE in eCREs:
    # eCRE = "Nkx2-2_p=0.001000"
    # eCRE = eCREs[0]
    plot = GenomePlot(eCRE)
    plot.ini_all_motifs(phylo=True)
    # plot.ini_by_cluster()
    plot.ini_by_cluster_merge()
    plot.ini_by_candidate()
    plot.make_plots(parallel=False,suppress=True)
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

"""
