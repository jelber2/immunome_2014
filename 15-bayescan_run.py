#! /usr/bin/env python

# PBS cluster job submission in Python
# Run Bayescan to look for SNPs under selection
# By Jean P. Elbers
# jelber2@lsu.edu
# Last modified 18 Dec 2014
###############################################################################
Usage = """

13-bayescan_run.py - version 1.1
STEPS:
1.Runs BayeScan
        ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
        /work/jelber2/immunome/bayescan-beagle/bayes_input.txt \
        -snp \
        -d low_freq_snps.txt \ #optional
        -od . \
        -o Sample \
        -threads 16

File Info
InDir = /work/jelber2/immunome/bayescan-beagle
Input Files = bayescan_input.txt
OutDir = InDir
Output Files = bayescan_all_loci or 
               bayescan_no_loci_with_low_freq_minor_alleles


Usage (execute following code in InDir):

~/scripts/immunome/bayescan_run.py bayescan_no_loci_with_low_freq_minor_alleles

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else: # if you do enter files, then run the program
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome/bayescan-beagle"
    os.chdir(InDir)
    for InFileName in FileList: # do the following steps for each file in the inputstream
        Sample = InFileName
        # Customize your options here
        Queue = "workq"
        Allocation = "hpc_gopo01"
        Processors = "nodes=1:ppn=16"
        WallTime = "14:00:00"
        LogOut = "/work/jelber2/immunome/bayescan-beagle"
        LogMerge = "oe"
        JobName = "%s" % (Sample)
        Command ="""
        ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
        /work/jelber2/immunome/bayescan-beagle/bayes_input.txt \
        -snp \
        -d low_freq_snps.txt \
        -od . \
        -o %s \
        -threads 16""" % (Sample)

        JobString = """
        #!/bin/bash
        #PBS -q %s
        #PBS -A %s
        #PBS -l %s
        #PBS -l walltime=%s
        #PBS -o %s
        #PBS -j %s
        #PBS -N %s

        cd %s
        %s\n""" % (Queue, Allocation, Processors, WallTime, LogOut, LogMerge, JobName, InDir, Command)

        #Create pipe to qsub
        proc = subprocess.Popen(['qsub'], shell=True,
          stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        (child_stdout, child_stdin) = (proc.stdout, proc.stdin)

        #Print JobString
        jobname = proc.communicate(JobString)[0]
        print JobString
        print jobname
