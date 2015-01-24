#!/bin/bash
#
###############################################################################
#
# "make_indexes.sh" SuperMikeII script 
# created by Jean P. Elbers
# jean.elbers@gmail.com
# last edited 22 Jan 2015
#
###############################################################################
#
#PBS -q single
#PBS -A hpc_startup_jelber2
#PBS -l nodes=1:ppn=4
#PBS -l walltime=04:00:00
#PBS -o /work/jelber2/reference
#PBS -j oe
#PBS -N makes_indexes

# Let's mark the time things get started with a date-time stamp.

date

# Set some handy environment variables.

export WORK_DIR=/work/jelber2/reference

# Make sure the WORK_DIR exists:

mkdir -p $WORK_DIR

# Creates BWA index for P.turtle genome with prefix 
# GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic

cd $WORK_DIR
~/bin/bwa-0.7.12/bwa index -a bwtsw \
-p GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna

# Makes GATK dictionary then fasta index

java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar \
CreateSequenceDictionary \
R=GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
O=GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.dict

~/bin/samtools-1.1/samtools faidx \
GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna

# Build a genome (.stidx) file:

python ~/bin/stampy-1.0.23/stampy.py --noparseNCBI \
-G GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna

# Build a hash (.sthash) file:

python ~/bin/stampy-1.0.23/stampy.py \
-g GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
-H GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic
