#! /usr/bin/env python

# PBS cluster job submission in Python
# Use picard to merge realigned-around-indel BAM files then makes index
# Then call SNPs with GATK-3.3.0
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

06-mergeBAM_callSNPs_initial.py - version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/combined/realign-around-indels
1.Merge Bam files
    java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar MergeSamFiles \
    SO=coordinate \
    AS=true \
    CREATE_INDEX=true \
    I=Sample1-realigned.bam \
    I=Sample2-realigned.bam \
    O=../call-SNPs-initial/ALL-samples-initial.bam

2.Call SNPs
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ../call-SNPs-initial/ALL-samples-initial.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ../call-SNPs-initial/ALL-samples-initial-Q30-rawSNPS.vcf

3.Call Indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ../call-SNPs-initial/ALL-samples-initial.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ../call-SNPs-initial/ALL-samples-initial-Q30-indels.vcf

4.Filter SNP calls around indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -V ../call-SNPs-initial/ALL-samples-initial-Q30-rawSNPS.vcf \
    --mask ../call-SNPs-initial/ALL-samples-initial-Q30-indels.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad Validation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 2.0" \
    --filterName "Low Variant Confidence" \
    --genotypeFilterExpression "DP < 10.0" \
    --genotypeFilterName "Low Read Depth Over Sample" \
    --genotypeFilterExpression "GQ < 20.0" \
    --genotypeFilterName "Low GenotypeQuality" \
    -o ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf


File Info:
InDir = /work/jelber2/immunome_2014/combined/realign-around-indels
Input Files = *-realigned.bam
OutDir = /work/jelber2/immunome_2014/combined/call-SNPs-initial
Output Files = ALL-samples-initial.bam
               ALL-samples-initial-Q30-SNPs.vcf


Usage (execute following code in InDir):

~/scripts/immunome_2014/06-mergeBAM_callSNPs_initial.py *-realigned.bam

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    IFileList =[]
    for File in FileList:
        IFile = "I="+File+" \\"
        IFileList.append(IFile)
    IFileListString = '\n'.join(IFileList)
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/combined/realign-around-indels"
    OutDir1 = "call-SNPs-initial"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, make it
    os.chdir(InDir)
    # Customize your options here
    Queue = "single"
    Allocation = "hpc_gopo02"
    Processors = "nodes=1:ppn=4"
    WallTime = "04:00:00"
    LogOut = "/work/jelber2/immunome_2014/combined/call-SNPs-initial"
    LogMerge = "oe"
    JobName = "mergeBAM_callSNPs_initial"
    Command ="""
    java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar MergeSamFiles \
    SO=coordinate \
    AS=true \
    CREATE_INDEX=true \
    %s
    O=../call-SNPs-initial/ALL-samples-initial.bam

    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ../call-SNPs-initial/ALL-samples-initial.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ../call-SNPs-initial/ALL-samples-initial-Q30-rawSNPS.vcf

    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T UnifiedGenotyper \
    -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ../call-SNPs-initial/ALL-samples-initial.bam \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -maxAltAlleles 32 \
    -gt_mode DISCOVERY \
    -glm INDEL \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    -o ../call-SNPs-initial/ALL-samples-initial-Q30-indels.vcf

    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -L /work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -V ../call-SNPs-initial/ALL-samples-initial-Q30-rawSNPS.vcf \
    --mask ../call-SNPs-initial/ALL-samples-initial-Q30-indels.vcf \
    --maskExtension 5 \
    --maskName InDel \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad Validation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --filterExpression "QD < 2.0" \
    --filterName "Low Variant Confidence" \
    --genotypeFilterExpression "DP < 10.0" \
    --genotypeFilterName "Low Read Depth Over Sample" \
    --genotypeFilterExpression "GQ < 20.0" \
    --genotypeFilterName "Low GenotypeQuality" \
    -o ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf""" % \
    (IFileListString,
    RefDir,
    RefDir,
    RefDir)

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
