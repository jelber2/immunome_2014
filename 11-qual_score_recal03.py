#! /usr/bin/env python

# PBS cluster job submission in Python
# Uses GATK-3.3.0 BaseRecalibrator to recalibrate quality scores
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

11-qual_score_recal03.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/combined/call-SNPs-recal02
1.Analyze patterns of covariation in the sequence dataset
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I Sample-recal02.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -knownSites ALL-samples-recal02-Q30-SNPs.vcf \
    -o ../call-SNPs-recal03/Sample-recal-data.table

2.Do a second pass to analyze covariation remaining after recalibration
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I Sample-recal02.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -knownSites ../call-SNPs-recal02/ALL-samples-recal02-Q30-SNPs.vcf \
    -BQSR ../call-SNPs-recal03/Sample-recal-data.table \
    -o ../call-SNPs-recal03/Sample-post-recal-data.table

3.Generate before/after plots
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T AnalyzeCovariates \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -before ../call-SNPs-recal03/Sample-recal-data.table \
    -after ../call-SNPs-recal03/Sample-post-recal-data.table \
    -plots ../call-SNPs-recal03/Sample-recalibration_plots.pdf

4.Apply the recalibration to your sequence data
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I Sample-recal02.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -BQSR ../call-SNPs-recal03/Sample-recal-data.table \
    -o ../call-SNPs-recal03/Sample-recal03.bam


File Info
InDir = /work/jelber2/immunome_2014/combined/call-SNPs-recal02
Input Files =
       Sample-recal02.bam
       ../call-SNPs-recal02/ALL-samples-recal02-Q30-SNPs.vcf
OutDir = /work/jelber2/immunome_2014/combined/call-SNPs-recal03
Output Files = Sample-recal03.bam


Usage (execute following code in InDir):

find . -name '*-recal02.bam' -not -name 'ALL-samples-*' -exec ~/scripts/immunome_2014/11-qual_score_recal03.py {} \;

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/combined/call-SNPs-recal02"
    OutDir1 = "call-SNPs-recal03"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, then make it
    os.chdir(InDir)
    for InFileName in FileList: # do the following steps for each file in the inputstream
        FileSuffix = "-recal02.bam"
        FilePrefix = "./"
        Samplepre = InFileName.replace(FileSuffix,'') # creates Samplepre string
        Sample = Samplepre.replace(FilePrefix,'') # creates Sample string
        # Customize your options here
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=4"
        WallTime = "01:00:00"
        LogOut = "/work/jelber2/immunome_2014/combined/call-SNPs-recal03"
        LogMerge = "oe"
        JobName = "qual_score_recal03-%s" % (Sample)
        Command ="""
        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I %s-recal02.bam \
        -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
        -knownSites ALL-samples-recal02-Q30-SNPs.vcf \
        -o ../call-SNPs-recal03/%s-recal-data.table

        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I %s-recal02.bam \
        -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
        -knownSites ../call-SNPs-recal02/ALL-samples-recal02-Q30-SNPs.vcf \
        -BQSR ../call-SNPs-recal03/%s-recal-data.table \
        -o ../call-SNPs-recal03/%s-post-recal-data.table

        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T AnalyzeCovariates \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -before ../call-SNPs-recal03/%s-recal-data.table \
        -after ../call-SNPs-recal03/%s-post-recal-data.table \
        -plots ../call-SNPs-recal03/%s-recalibration_plots.pdf

        java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
        -I %s-recal02.bam \
        -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
        -BQSR ../call-SNPs-recal03/%s-recal-data.table \
        -o ../call-SNPs-recal03/%s-recal03.bam""" % \
        (RefDir, Sample, Sample,
        RefDir, Sample, Sample, Sample,
        RefDir, Sample, Sample, Sample,
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
