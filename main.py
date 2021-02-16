from synCRE import *

# expr = Expression(RNA_seq_file="reference/RNA_seq_local.txt")
# expr.run_all()
#
mtf = Motif_Finder()
# mtf.make_pmf()
# mtf.sample_eCRE_sequence()
# mtf.find_motifs([0.001])
mtf.motif_to_bed()
mtf.motifs_by_archetype()
mtf.motifs_by_cluster()

plot = GenomePlot("Olig2_p=0.000050")
plot.ini_all_motifs()
plot.ini_by_cluster()

