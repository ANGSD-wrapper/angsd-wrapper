#!/bin/bash

set -e
set -u
set -o pipefail

#   A simple script to hold variables for the Site Frequency Spectrum

#   Name the project
PROJECT=

#   Define a list of samples
SAMPLE_LIST=

#   Define a list of inbreeding coefficients
SAMPLE_INBREEDING=

#   Where do we put the outfiles?
OUTDIR=

#   Define the region being looked at
#       Optional, but ANGSD is expensive to run without specifying regions to look at
REGIONS=

#   Where is ANGSD?
ANGSD_DIR=

#   Genotypes Parameters
#       Listed below are the defaults, please modify for your samples
DO_MAJORMINOR=1
UNIQUE_ONLY=0
MIN_MAPQ=30
MIN_BASEQUAL=20
GT_LIKELIHOOD=2
DO_GENO=7
DO_POST=1
POST_CUTOFF=0.95
DO_MAF=2
SNP_PVAL=1e-6
MIN_IND=1
N_CORES=32