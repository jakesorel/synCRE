from synCRE import *

# expr = Expression(RNA_seq_file="reference/RNA_seq_local.txt")
# expr.run_all()

mtf = Motif_Finder()
# mtf.make_pmf()
# mtf.sample_eCRE_sequence()
# mtf.find_motifs([0.001])
mtf.motif_to_bed()
mtf.motifs_by_archetype()
mtf.motifs_by_cluster()

print("motifs found")

eCREs = [name.split(".bed")[0] for name in os.listdir("results/motifs/bed")]
for eCRE in eCREs:
    plot = GenomePlot(eCRE)
    plot.ini_all_motifs()
    plot.ini_by_cluster()
    plot.ini_by_cluster_merge()
    plot.ini_by_candidate()
    plot.make_plots(parallel=True,suppress=True)
    print("""
################################################
Plots for %s complete
################################################
    """%eCRE)


