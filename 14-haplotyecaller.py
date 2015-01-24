#! /usr/bin/env python

# PBS cluster job submission in Python
# Uses GATK-3.3.0 Haplotype Caller for variant calling performed per-sample
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

14-hc_per_sample.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/run1/call-SNPs-recal03
1.Analyze patterns of covariation in the sequence dataset
        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I InDir/Sample-recal03.bam \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -L RefDir/immunome_baits_C_picta-3.0.3.list \
        -stand_call_conf 30 \
        -stand_emit_conf 10 \
        -o Sample-raw-snps-indels.vcf

File Info
InDir = /work/jelber2/immunome_2014/run1/call-SNPs-recal03
Input Files =
       Sample-recal03.bam
OutDir = /work/jelber2/immunome_2014/run1/hc
Output Files = Sample-raw-snps-indels.vcf


Usage (execute following code in InDir):

find . -name '*-recal03.bam' -not -name 'ALL-samples-*' -exec ~/scripts/immunome_2014/14-hc_per_sample.py {} \;

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/run1/call-SNPs-recal03"
    OutDir1 = "hc"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, then make it
    os.chdir(InDir)
    for InFileName in FileList: # do the following steps for each file in the inputstream
        FileSuffix = "-recal03.bam"
        FilePrefix = "./"
        Samplepre = InFileName.replace(FileSuffix,'') # creates Samplepre string
        Sample = Samplepre.replace(FilePrefix,'') # creates Sample string
        # Customize your options here
        Queue = "single"
        Allocation = "hpc_gopo01"
        Processors = "nodes=1:ppn=4"
        WallTime = "01:00:00"
        LogOut = "/work/jelber2/immunome_2014/run1/hc"
        LogMerge = "oe"
        JobName = "hc_per_sample-%s" % (Sample)
        Command ="""
        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T HaplotypeCaller \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I %s/%s-recal03.bam \
        --emitRefConfidence GVCF \
        --variant_index_type LINEAR \
        --variant_index_parameter 128000 \
        -L %s/immunome_baits_C_picta-3.0.3.list \
        -stand_call_conf 30 \
        -stand_emit_conf 10 \
        -o ../hc/%s-raw-snps-indels.vcf""" % \
        (RefDir, InDir, Sample, RefDir, Sample)

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
