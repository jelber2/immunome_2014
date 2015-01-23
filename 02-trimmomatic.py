#! /usr/bin/env python

# PBS cluster job submission in Python
# Trimmomatic for quality and adapter trimming
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

02-trimmomatic.py - version 1.0

Command:
1.Uses Trimmomatic for quality and adapter trimming
    java -Xmx8g-jar ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar \
    PE \
    -threads 4 \
    -phred33 \
    Sample-R1.fastq.gz \
    Sample-R2.fastq.gz \
    OutDir/Sample-R1-paired.trim.fastq.gz \
    OutDir/Sample-R1-single.trim.fastq.gz \
    OutDir/Sample-R2-paired.trim.fastq.gz \
    OutDir/Sample-R2-single.trim.fastq.gz \
    ILLUMINACLIP:/home/jelber2/bin/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:5 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:40

2.Uses bbmerge (part of BBMAP suite) to merge overlapping paired end reads
  and also outputs reads that do not overlap
    cd OutDir
    ~/bin/bbmap-34.33/bbmerge.sh \
    -Xmx8g \
    t=4 \
    in1=Sample-R1-paired.trim.fastq.gz \
    in2=Sample-R2-paired.trim.fastq.gz \
    out=Sample-merged.trim.fastq.gz \
    outu1=Sample-R1-unmerged.trim.fastq.gz \
    outu2=Sample-R2-unmerged.trim.fastq.gz

3.Combines singleton reads into one file
    cat Sample-R1-single.trim.fastq.gz Sample-R2-single.trim.fastq.gz > Sample-singles.trim.fastq.gz

4.Combines singleton and merged files into one file
    cat Sample-singles.trim.fastq.gz Sample-merged.trim.fastq.gz > Sample-singlesANDmerged.trim.fastq.gz


Directory info:
InDir = /work/jelber2/immunome_2014/run1/fastq/
Input Files = *-R1.fastq.gz
              *-R2.fastq.gz
OutDir = /work/jelber2/immunome_2014/run1/trimmed-data
Important Output Files = Sample-R1-unmerged.trim.fastq.gz
                         Sample-R2-unmerged.trim.fastq.gz
                         Sample-singlesANDmerged.trim.fastq.gz


Usage (execute following code in InDir):

~/scripts/immunome_2014/02-trimmomatic.py *.fastq.gz

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/run1/fastq/"
    OutDir = "/work/jelber2/immunome_2014/run1/trimmed-data"
    NewFolderName = "trimmed-data"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(NewFolderName):
        os.mkdir(NewFolderName) # if NewFolderName does not exist, then make it
    os.chdir(InDir)
    R1 = re.compile(r"\w+\-R1\.fastq\.gz") # search string for obtaining R1 files
    FileList = [f for f in FileList if R1.match(f)] # keeps only R1 files in FileList
    for InFileName in FileList: # so trimmomatic only does command once for each R1 and R2 pair
        FileSuffix = "-R1.fastq.gz" # string to remove from InFileName
        Sample = InFileName.replace(FileSuffix,'') # creates Sample string
        # Customize your options here
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=4"
        WallTime = "04:00:00"
        LogOut = OutDir
        LogMerge = "oe"
        JobName = "Trimmomatic-%s" % (Sample)
        Command = """
        java -jar ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar \
        PE \
        -threads 4 \
        -phred33 \
        %s-R1.fastq.gz \
        %s-R2.fastq.gz \
        %s/%s-R1-paired.trim.fastq.gz \
        %s/%s-R1-single.trim.fastq.gz \
        %s/%s-R2-paired.trim.fastq.gz \
        %s/%s-R2-single.trim.fastq.gz \
        ILLUMINACLIP:/home/jelber2/bin/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:5 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:40

        cd %s

        ~/bin/bbmap-34.33/bbmerge.sh \
        -Xmx8g \
        t=4 \
        in1=%s-R1-paired.trim.fastq.gz \
        in2=%s-R2-paired.trim.fastq.gz \
        out=%s-merged.trim.fastq.gz \
        outu1=%s-R1-unmerged.trim.fastq.gz \
        outu2=%s-R2-unmerged.trim.fastq.gz

        cat %s-R1-single.trim.fastq.gz %s-R2-single.trim.fastq.gz > %s-singles.trim.fastq.gz

        cat %s-singles.trim.fastq.gz %s-merged.trim.fastq.gz > %s-singlesANDmerged.trim.fastq.gz""" % \
        (Sample, Sample, OutDir, Sample, OutDir, Sample, OutDir, Sample, OutDir, Sample,
        OutDir,
        Sample, Sample, Sample, Sample, Sample, 
        Sample, Sample, Sample,
        Sample, Sample, Sample)

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
