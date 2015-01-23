#! /usr/bin/env python

# PBS cluster job submission in Python
# Clean, sort, add Read groups, and Mark duplicates, realign around indels
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

05a-clean_sort_addRG.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/run1/stampy-alignment
1.Uses samtools merge to combine stampy bam files
        ~/bin/samtools-1.1/samtools merge -@ 1\
        Sample.stampy.bam Sample-paired.stampy.bam Sample-singlesANDmerged.stampy.bam

2.Uses samtools flagstat to get alignment metrics on stampy aligned bam file
        ~/bin/samtools-1.1/samtools flagstat \
        Sample.stampy.bam > ../stampy-alignment/Sample.stampy.bam.flagstat

3.Clean the initial stampy BAM file:
        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar CleanSam \
        I=Sample.stampy.bam \
        O=../clean-sort-addRG/Sample-CL.bam

4. Add read groups and sort:
        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar AddOrReplaceReadGroups \
        I=../clean-sort-addRG/Sample-CL.bam \
        O=../clean-sort-addRG/Sample-CL-RG.bam \
        SORT_ORDER=coordinate \
        RGPL=illumina \
        RGPU=barcode \
        RGLB=Lib1 \
        RGID=Sample \
        RGSM=Sample \
        VALIDATION_STRINGENCY=LENIENT


Directory info:

(1)/work/jelber2/immunome_2014/run1/stampy-alignment
(2)                                /clean-sort-addRG

InDir = /work/jelber2/immunome_2014/run1/stampy-alignment
Input Files = *-singlesANDmerged.stampy.bam
              *-paired.stampy.bam

Usage (execute following code in InDir):

~/scripts/immunome_2014/05a-clean_sort_addRG.py *.stampy.bam

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/run1/stampy-alignment"
    OutDir1 = "clean-sort-addRG"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, make it
    os.chdir(InDir)
    Paired = re.compile(r"\w+\-paired\.stampy\.bam") # search string for obtaining paired sam files
    FileList = [f for f in FileList if Paired.match(f)] # keeps only paired files in FileList
    for InFileName in FileList: # so samtools grabs only the file names (i.e., Samples)
        FileSuffix = "-paired.stampy.bam" # string to remove from InFileName
        Sample = InFileName.replace(FileSuffix,'') # creates Sample string
        # Customize your job options here
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=4"
        WallTime = "00:30:00"
        LogOut = "/work/jelber2/immunome_2014/run1/clean-sort-addRG"
        LogMerge = "oe"
        JobName = "clean-sort-addRG-%s" % (Sample)
        Command ="""
        ~/bin/samtools-1.1/samtools merge -@ 1\
        %s.stampy.bam %s-paired.stampy.bam %s-singlesANDmerged.stampy.bam

        ~/bin/samtools-1.1/samtools flagstat \
        %s.stampy.bam > %s.stampy.bam.flagstat

        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar CleanSam \
        I=%s.stampy.bam \
        O=../clean-sort-addRG/%s-CL.bam

        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar AddOrReplaceReadGroups \
        I=../clean-sort-addRG/%s-CL.bam \
        O=../clean-sort-addRG/%s-CL-RG.bam \
        SORT_ORDER=coordinate \
        RGPL=illumina \
        RGPU=barcode \
        RGLB=Lib1 \
        RGID=%s_9Sep2014 \
        RGSM=%s \
        VALIDATION_STRINGENCY=LENIENT""" % \
        (Sample, Sample, Sample,
        Sample, Sample,
        Sample, Sample,
        Sample, Sample, Sample, Sample)

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

