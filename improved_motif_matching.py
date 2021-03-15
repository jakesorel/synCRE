from pybedtools import BedTool
import os
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
def make_directory(dir):
    """

    :param dir:
    :return:
    """
    if not os.path.exists(dir):
        os.mkdir(dir)

make_directory("results/improved_motif")

ATAC_bed = BedTool("/camp/lab/briscoej/working/Joaquina/hpc_camp/ATAC_FullFiltered/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed")
Olig2_bed = BedTool("/camp/lab/briscoej/working/Joaquina/hpc_camp/public_data/Nishi2015_PMID26293298/results/bwa/mergedLibrary/macs/broadPeak/Olig2_IP_R2_peaks.gappedPeak")

genome_dir = '/camp/lab/luscomben/home/shared/ref/genomes/mouse/gencode_GRCm38/GRCm38.releaseM25.primary_assembly.genome.fa'

make_directory("results/fasta/ATAC")

ATAC_bed.sequence(genome_dir, fo="results/fasta/ATAC/mm10_ATAC.fa")

runline = """
eval "$(conda shell.bash hook)" \n
source activate moods \n
cd pkgs/moods/scripts \n
pwd \n
python2 moods-dna.py  \
--sep ";" -s ../../../results/fasta/ATAC/mm10_ATAC.fa --threshold 4 \
--lo-bg 2.977e-01 2.023e-01 2.023e-01 2.977e-01 \
-m ../../../reference/motifs/pmf/OLIG2_MOUSE.H11MO.0.A.pmf -o ../../../results/improved_motif/raw.csv
"""
os.system(runline)

motif_table = pd.read_excel("reference/motif_annotations.xlsx", 1, engine='openpyxl')

lookup_names = dict(zip(motif_table["Motif"], motif_table["Cluster_ID"]))

df = pd.read_csv("results/improved_motif/raw.csv", sep=";", header=None)
seq_len = [len(string) for string in df[5]]
locs = df[0]
names = df[1].values
names = [name.split(".pmf")[0].split(" ")[0] for name in names]
motif_ids = []
for name in names:
    try:
        motif_ids.append(lookup_names[name])
    except KeyError:
        motif_ids.append(-1)
chrom, start, end = [], [], []
for loc in locs:
    chrm, pos = loc.split(":")
    st, en = pos.split("-")
    st, en = int(st), int(en)
    chrom.append(chrm)
    start.append(st)
    end.append(en)
start, end = np.array(start), np.array(end)
start += df[2].values
end = start + seq_len
f = open("results/improved_motif/raw.bed","w")
for c,s,e,m,sc in zip(chrom, start, end, motif_ids,df[4]):
    f.write("%s\t%s\t%s\t%s\t%s\n"%(c,s,e,m,sc))
f.close()

f = open("reference/ChiP/Olig2_IP_R1_peaks.gappedPeak")
g = open("reference/ChiP/Olig_narrow.bed","w")
for line in f.readlines():
    line_dict = dict(zip(["Chrom", "ChromStart","ChromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts","signalValue","pValue","qValue"],line.split("\n")[0].split("\t")))
    chrom = line_dict["Chrom"]
    fstart = int(line_dict["ChromStart"])
    starts = line_dict["blockStarts"].split(",")
    sizes = line_dict["blockSizes"].split(",")
    for start,size in zip(starts,sizes):
        g.write("%s\t%d\t%d\n"%(chrom,fstart+int(start),fstart+int(start)+int(size)))
f.close()
g.close()

Olig2_bed = BedTool("reference/ChiP/Olig2_IP_R1_peaks.broadPeak")
Olig2_bed = BedTool("reference/ChiP/Olig_narrow.bed")

# Olig2_df = Olig2_bed.to_dataframe()
# Olig2_df.columns = "Chrom", "ChromStart","ChromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts","signalValue","pValue","qValue"
motif_bed = BedTool("results/improved_motif/raw.bed")
motif_bed.intersect(Olig2_bed,c=True,wa=True).saveas("results/improved_motif/chip_intersect.bed")

int_df = pd.read_csv("results/improved_motif/chip_intersect.bed",header=None,sep="\t")
int_df.columns = "chrom","start","end","motif","score","chip"
df_sample = int_df
chip_mask = df_sample["chip"]==1

fig, ax = plt.subplots()
sb.kdeplot(data = df_sample.loc[chip_mask],x="score",label="Olig chip+",ax=ax)
sb.kdeplot(data = df_sample.loc[~chip_mask],x="score",label="Olig chip-",ax=ax)
ax.legend()
# ax.set(yscale="log")
fig.savefig("results/improved_motif/density_plot.pdf")

score_g_6 = df_sample["score"]>6
from scipy.stats import gaussian_kde
kde_p = gaussian_kde(df_sample.loc[chip_mask*score_g_6]["score"].values)
kde_n = gaussian_kde(df_sample.loc[(~chip_mask)*score_g_6]["score"].values)

fig, ax = plt.subplots()
s_span = np.linspace(6,14)
ax.plot(s_span,kde_p(s_span),label="Olig2 chip+")
ax.plot(s_span,kde_n(s_span),label="Olig2 chip-")
ax.set(yscale="log",xlabel="score",ylabel="density")
fig.savefig("results/improved_motif/density_plot_log.pdf")
