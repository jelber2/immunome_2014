#! /usr/bin/env python

# PBS cluster job submission in Python
# Clean, sort, add Read groups, and Mark duplicates, realign around indels
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

05b-clean_sort_addRG_markdup_realign.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/run1/clean-sort-addRG/
1.Uses samtools merge to combine stampy bam files
        ~/bin/samtools-1.1/samtools merge \
        /work/jelber2/immunome_2014/combined/merged-bams/Sample.bam \
        /work/jelber2/immunome_2014/run1/clean-sort-addRG/Sample-CL-RG.bam \
        /work/jelber2/immunome_2014/run2/clean-sort-addRG/Sample-CL-RG.bam

2.Uses samtools flagstat to get alignment metrics on stampy aligned bam file
        cd /work/jelber2/immunome_2014/combined/merged-bams/
        ~/bin/samtools-1.1/samtools flagstat \
        ../merged-bams/Sample.bam > ../merged-bams/Sample.bam.flagstat

3.Clean the initial stampy BAM file:
        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar CleanSam \
        I=../merged-bams/Sample.bam \
        O=../clean-sort-addRG-markdup/Sample-CL2.bam

4.Mark PCR duplicates and optical duplicates:
        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar MarkDuplicates \
        I=../clean-sort-addRG-markdup/Sample-CL2.bam \
        O=../clean-sort-addRG-markdup/Sample-CL2-MD.bam \
        METRICS_FILE=../clean-sort-addRG-markdup/Sample-CL2-MD.metrics \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
        CREATE_INDEX=true \
        ASSUME_SORTED=false \
        REMOVE_DUPLICATES=false

5.Find INDEL regions within individual BAM files
        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I ../clean-sort-addRG-markdup/Sample-CL2-MD.bam \
        --minReadsAtLocus 4 \
        -o ../realign-around-indels/Sample.merged.intervals

6.Realign the BAM based on indel intervals:
        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I ../clean-sort-addRG-markdup/Sample-CL2-MD.bam \
        -targetIntervals ../realign-around-indels/Sample.merged.intervals \
        -LOD 3.0 \
        -o ../realign-around-indels/Sample-realigned.bam

Directory info:

(1)/work/jelber2/immunome_2014/run1/clean-sort-addRG
(2)/work/jelber2/immunome_2014/combined/merged-bams
(3)                                    /clean-sort-addRG-markdup
(4)                                    /realign-around-indels

InDir = /work/jelber2/immunome_2014/combined/
Input Files = Sample-CL-RG.bam


Usage (execute following code in InDir):

~/scripts/immunome_2014/05b-clean_sort_addRG_markdup_realign *-CL-RG.bam

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/combined/"
    OutDir1 = "merged-bams"
    OutDir2 = "clean-sort-addRG-markdup"
    OutDir3 = "realign-around-indels"
    os.chdir(InDir)
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, make it
    if not os.path.exists(OutDir2):
        os.mkdir(OutDir2) # if OutDir2 does not exist, make it
    if not os.path.exists(OutDir3):
        os.mkdir(OutDir3) # if OutDir3 does not exist, make it
    os.chdir(InDir)
    for InFileName in FileList: # so samtools grabs only the file names (i.e., Samples)
        FileSuffix = "-CL-RG.bam" # string to remove from InFileName
        Sample = InFileName.replace(FileSuffix,'') # creates Sample string
        # Customize your job options here
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=4"
        WallTime = "04:00:00"
        LogOut = "/work/jelber2/immunome_2014/combined/clean-sort-addRG-markdup"
        LogMerge = "oe"
        JobName = "clean-sort-addRG-markdup-realign-%s" % (Sample)
        Command ="""
        ~/bin/samtools-1.1/samtools merge \
        /work/jelber2/immunome_2014/combined/merged-bams/%s.bam \
        /work/jelber2/immunome_2014/run1/clean-sort-addRG/%s-CL-RG.bam \
        /work/jelber2/immunome_2014/run2/clean-sort-addRG/%s-CL-RG.bam

        cd /work/jelber2/immunome_2014/combined/merged-bams/
        ~/bin/samtools-1.1/samtools flagstat \
        ../merged-bams/%s.bam > ../merged-bams/%s.bam.flagstat

        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar CleanSam \
        I=../merged-bams/%s.bam \
        O=../clean-sort-addRG-markdup/%s-CL2.bam

        java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar MarkDuplicates \
        I=../clean-sort-addRG-markdup/%s-CL2.bam \
        O=../clean-sort-addRG-markdup/%s-CL2-MD.bam \
        METRICS_FILE=../clean-sort-addRG-markdup/%s-CL2-MD.metrics \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
        CREATE_INDEX=true \
        ASSUME_SORTED=false \
        REMOVE_DUPLICATES=false

        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I ../clean-sort-addRG-markdup/%s-CL2-MD.bam \
        --minReadsAtLocus 4 \
        -o ../realign-around-indels/%s.merged.intervals

        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I ../clean-sort-addRG-markdup/%s-CL2-MD.bam \
        -targetIntervals ../realign-around-indels/%s.merged.intervals \
        -LOD 3.0 \
        -o ../realign-around-indels/%s-realigned.bam""" % \
        (Sample, Sample, Sample,
        Sample, Sample,
        Sample, Sample,
        Sample, Sample, Sample,
        RefDir, Sample, Sample,
        RefDir, Sample, Sample, Sample)

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

