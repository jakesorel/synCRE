from pybedtools import BedTool
import pandas as pd
consensus_ATAC = pd.read_csv("reference/CRE_info/consensus_peaks.mLb.clN.annotatePeaks.txt",sep="\t",skiprows=[0],header=None)
consensus_ATAC = consensus_ATAC[consensus_ATAC.columns[1:4]].to_csv("reference/eCRE_locs/consensus_peaks.bed",index=False,header=False,sep="\t")
consensus = BedTool("reference/CRE_info/consensus_peaks.bed")
oosterveen = BedTool("reference/CRE_info/oosterveen_mm10.bed")
oosterveendf = BedTool("reference/CRE_info/oosterveen_mm10.bed").to_dataframe()

locs, names = [],[]
for i in range(oosterveendf.shape[0]):
    loc = consensus.intersect(oosterveen.at([i]),u=True).to_dataframe()
    name = oosterveendf["name"].values[i]
    locs.append(loc)
    names.append(name)
# inters = consensus.intersect(oosterveen,wa=True,u=True).to_dataframe()


for i in range(len(locs)):
    # pd.DataFrame([eCRE[inters.columns[0:3]].values]).to_csv("../../reference/eCRE_locs/%s.bed"%eCRE["name"],header=None,index=None,sep="\t")
    locs[i].to_csv("reference/eCRE_locs/%s.bed"%names[i],header=None,index=None,sep="\t")


