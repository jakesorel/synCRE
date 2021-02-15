from pybedtools import BedTool
import pandas as pd
consensus_ATAC = pd.read_csv("../../reference/CRE_info/consensus_peaks.mLb.clN.annotatePeaks.txt",sep="\t",skiprows=[0],header=None)
consensus_ATAC = consensus_ATAC[consensus_ATAC.columns[1:4]].to_csv("../../reference/eCRE_locs/consensus_peaks.bed",index=False,header=False,sep="\t")
consensus = BedTool("../../reference/CRE_info/consensus_peaks.bed")
oosterveen = BedTool("../../reference/CRE_info/oosterveen_mm10.bed")
inters = oosterveen.intersect(consensus,wb=True).to_dataframe()

for i in range(inters.shape[0]):
    eCRE = inters.iloc[i]
    pd.DataFrame([eCRE[inters.columns[0:3]].values]).to_csv("../../reference/eCRE_locs/%s.bed"%eCRE["name"],header=None,index=None,sep="\t")
    # pd.DataFrame([eCRE[inters.columns[0:3]].values]).to_csv("reference/eCRE_locs/%s.bed"%eCRE["name"],header=None,index=None,sep="\t")
