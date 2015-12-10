#! /usr/bin/env python

# PBS cluster job submission in Python
# Use structure on msat loci
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 25 Feb 2015
###############################################################################
Usage = """

16-structure-msats.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/combined/popgen-msats
1.Run structure
        ~/bin/structure/structure_kernel_src/structure \
        -m Sample \
        -e ~/bin/structure/structure_kernel_src/extraparams

File Info:
InDir = /work/jelber2/immunome_2014/combined/popgen-msats
Input Files = mainparams.test.0*
OutDir = InDir
Output Files = structure-results-001-k1_f
               structure-results-001-k2_f
               etc.


Usage (execute following code in InDir):

~/scripts/immunome_2014/16-structure.py mainparams.test.0*

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/combined/popgen-msats"
    OutDir = InDir
    # Customize your options here
    for InFileName in FileList:
        FilePrefix = "mainparams.test."
        Samplepre = InFileName.replace(FilePrefix,'') # creates Sample string
        Sample = Samplepre.replace(".","-") # creates Sample string
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=1"
        WallTime = "12:00:00"
        LogOut = InDir
        LogMerge = "oe"
        JobName = "run-structure-%s" % (Sample)
        Command ="""
        ~/bin/structure/structure_kernel_src/structure \
        -m %s \
        -e ~/bin/structure/structure_kernel_src/extraparams""" % (InFileName)

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
        JobName = proc.communicate(JobString)[0]
        print JobString
        print JobName
