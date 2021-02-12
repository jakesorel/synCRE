#!/bin/bash
wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz
gzip -d Mus_musculus.gene_info.gz
mv Mus_musculus.gene_info gene_info