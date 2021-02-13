#!/usr/bin/bash
ml purge
ml Anaconda2/2018.12
set -euo pipefail
conda env create -f env.yml
conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7