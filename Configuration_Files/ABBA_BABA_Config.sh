#!/bin/bash

set -e
set -u
set -o pipefail

#   A simple script to hold variables for ABBA BABA

#   Name of the project
PROJECT=

#   Define a list of samples
SAMPLE_LIST=

OUTGROUP=nope

#   Ancestral sequence
ANC_SEQ=

#   WHere do we put the outfiles?
OUTDIR=

#   Define the region being looked at
#       Optional, but ANGSD is expensive to run without specifying regions to look at
REGIONS=

#   Where is ANGSD?
ANGSD_DIR=

#   ABBA BABA Parameters
#       Listed below are the defaults, please modify for your samples
DO_CONSENSUS=0
DO_FASTA=2
DO_COUNTS=1
DO_ABBABABA=1
UNIQUE_ONLY=0
MIN_BASEQUAL=20
BAQ=1
MIN_IND=1
MIN_MAPQ=30
N_CORES=32
CHECK_BAM_HEADERS=0
BLOCKSIZE=1000