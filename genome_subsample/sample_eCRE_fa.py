# from pybedtools import BedTool,example_filename
# import os, sys
# import pandas as pd
# import numpy as np
# currentdir = os.path.dirname(os.path.realpath(__file__))
# parentdir = os.path.dirname(currentdir)
# sys.path.append(parentdir)
# from reformat_functions.functions import *
#
# #
# # chr16fa = "/camp/lab/luscomben/home/shared/ref/genomes/mouse/gencode_GRCm38/chr_sep_files/chr16.fa.gz"
# # chr16fa = "chr16.fa.gz"
# # unzip(chr16fa,"out.fa")
# # fa_file = example_filename('%s/out.fa'%currentdir)
# # olig2fa = BedTool("reference/eCRE_locs/Olig2.bed").sequence("out.fa",fo="sample.fa")
# # print(open(olig2fa.seqfn).read())
# #

from pybedtools import BedTool,example_filename
import os, sys
import pandas as pd
import numpy as np
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from reformat_functions.functions import *

make_directory("results/fasta")
make_directory("results/fasta/scrap")
make_directory("results/fasta/by_eCRE")

bed_files = os.listdir("reference/eCRE_locs")
for bed_file in bed_files:
    eCRE_name = bed_file.split(".bed")[0]
    eCRE = BedTool("reference/eCRE_locs/%s"%bed_file)
    chrom = eCRE.to_dataframe()["chrom"].values[0]
    chrfa = "/camp/lab/luscomben/home/shared/ref/genomes/mouse/gencode_GRCm38/chr_sep_files/%s.fa.gz"%chrom
    # chr16fa = "chr16.fa.gz"
    unzip(chrfa,"results/fasta/scrap/out.fa")
    try:
        eCRE.sequence("results/fasta/scrap/out.fa",fo="results/fasta/by_eCRE/%s.fa"%eCRE_name)
    except RuntimeWarning:
        pass
    os.remove("results/fasta/scrap/out.fa")
    os.remove("results/fasta/scrap/out.fai")
    print(eCRE_name," complete")
