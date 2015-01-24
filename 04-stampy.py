#! /usr/bin/env python

# PBS cluster job submission in Python
# Stampy alignment to Painted Turtle Genome
# By Jean P. Elbers
# jean.elbers@gmail.com
# Last modified 22 Jan 2015
###############################################################################
Usage = """

04-stampy.py version 1.0

Command:
cd InDir = /work/jelber2/immunome_2014/run1/bwa-alignment
1.Produce only sam header and count number of header lines with wc, also copy
  header file to stampy-alignment directory for bam concatenation later
        ~/bin/samtools-1.1/samtools view -SH \
        -o Sample.header.bwa.sam Sample.bwa.sam
        export HeadLength_i="$(wc -l < Sample.header.bwa.sam)"
        cp Sample.header.bwa.sam ../stampy-alignment/

2.Use awk to remove header based on its length from sam file
        awk -v n=$HeadLength_i '1 <= NR && NR <= n {next} {print}' Sample.bwa.sam > Sample.noheader.bwa.sam

3.Split the sam file without a header into 16 equal parts
        export NumReads_i="$(wc -l < Sample.noheader.bwa.sam)"
        export SplitReads_i="$(echo $[(NumReads_i+15)/16])"
        split -d -l $SplitReads_i Sample.noheader.bwa.sam Sample-

4.Add the header back to each split
        ~/bin/parallel-20150122/src/parallel 'cat Sample.header.bwa.sam {} > {}.sam' ::: Sample-*

5.Convert the split sam files to bam format and remove sam in filename
        ~/bin/parallel-20150122/src/parallel '~/bin/samtools-1.1/samtools view -Sb \
        {} > {}.bwa.bam' ::: Sample-*.sam
        ~/bin/parallel-20150122/src/parallel 'rename .sam "" {}' ::: Sample-*.sam.bwa.bam

6.Map the bam files using stampy
        ~/bin/parallel-20150122/src/parallel 'python ~/bin/stampy-1.0.23/stampy.py \
        -g RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
        -h RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
        -f sam \
        -t 1 \
        --bamkeepgoodreads \
        -M {} > ../stampy-alignment/{}.stampy.sam' ::: Sample-*.bwa.bam

7.Remove all but original sam file associated with the particular sample from
  the bwa-alignment directory
        find . -name 'Sample*' -not -name 'Sample.bwa.sam' -type f -delete

8.Rename bwa.bam from Sample-##.bwa.bam.stampy.sam
        cd ../stampy-alignment/
        ~/bin/parallel-20150122/src/parallel 'rename .bwa.bam "" {}' ::: Sample-*.bwa.bam.stampy.sam

9.Convert sam2bam then remove sam from filename
        ~/bin/parallel-20150122/src/parallel '~/bin/samtools-1.1/samtools view \
        -bS \
        -t RefDir/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fai {} > {}.bam' ::: Sample-*.stampy.sam
        ~/bin/parallel-20150122/src/parallel 'rename .sam "" {}' ::: Sample-*stampy.sam.bam

10.Concatenate bam files
        ~/bin/samtools-1.1/samtools cat \
        -h Sample.header.sam \
        -o Sample.stampy.bam \
        Sample-00.stampy.bam \
        Sample-01.stampy.bam \
        Sample-02.stampy.bam \
        Sample-03.stampy.bam \
        Sample-04.stampy.bam \
        Sample-05.stampy.bam \
        Sample-06.stampy.bam \
        Sample-07.stampy.bam \
        Sample-08.stampy.bam \
        Sample-09.stampy.bam \
        Sample-10.stampy.bam \
        Sample-11.stampy.bam \
        Sample-12.stampy.bam \
        Sample-13.stampy.bam \
        Sample-14.stampy.bam \
        Sample-15.stampy.bam

11.Uses samtools flagstat to get alignment metrics on concatenated aligned sam file
        ~/bin/samtools-1.1/samtools flagstat \
        Sample.stampy.bam > Sample.stampy.bam.flagstat

12.Remove sam files and bam splits associated with the particular sample
   in stampy-alignment directory
        rm Sample*.sam
        rm Sample-*.stampy.bam


Directory info:
InDir = /work/jelber2/immunome_2014/run1/bwa-alignment
Input Files = Sample-singlesANDmerged.bwa.sam
              Sample-paired.bwa.sam
OutDir = /work/jelber2/immunome_2014/run1/stampy-alignment
Output Files = Sample-singlesANDmerged.stampy.bam
               Sample-singlesANDmerged.stampy.bam.flagstat
               Sample-paired.stampy.bam
               Sample-paired.stampy.bam.flagstat


Usage (execute following code in InDir):

~/scripts/immunome_2014/04-stampy.py *.bwa.sam

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome_2014/run1/bwa-alignment"
    OutDir = "/work/jelber2/immunome_2014/run1/stampy-alignment"
    NewFolderName = "stampy-alignment"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(NewFolderName):
        os.mkdir(NewFolderName) # if NewFolderName does not exist, then make it
    os.chdir(InDir)
    i=0
    for InFileName in FileList: # so stampy processes all bam files
        i+=1
        FileSuffix = ".bwa.sam" # string to remove suffix from InFileName
        Sample = InFileName.replace(FileSuffix,'') # creates Sample (note: Sample = Tort1-paired and/or Tort1-singlesANDmerged not Tort1!!)
        # Customize your job options here
        Queue = "workq"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=16"
        WallTime = "08:00:00"
        LogOut = OutDir
        LogMerge = "oe"
        JobName = "stampy-%s" % (Sample)
        Command ="""
        ~/bin/samtools-1.1/samtools view -SH \
        -o %s.header.bwa.sam %s.bwa.sam
        export HeadLength_%d="$(wc -l < %s.header.bwa.sam)"
        cp %s.header.bwa.sam ../stampy-alignment/

        awk -v n=$HeadLength_%d '1 <= NR && NR <= n {next} {print}' %s.bwa.sam > %s.noheader.bwa.sam

        export NumReads_%d="$(wc -l < %s.noheader.bwa.sam)"
        export SplitReads_%d="$(echo $[(NumReads_%d+15)/16])"
        split -d -l $SplitReads_%d %s.noheader.bwa.sam %s-

        ~/bin/parallel-20150122/src/parallel 'cat %s.header.bwa.sam {} > {}.sam' ::: %s-*

        ~/bin/parallel-20150122/src/parallel '~/bin/samtools-1.1/samtools view -Sb \
        {} > {}.bwa.bam' ::: %s-*.sam
        ~/bin/parallel-20150122/src/parallel 'rename .sam "" {}' ::: %s-*.sam.bwa.bam

        ~/bin/parallel-20150122/src/parallel 'python ~/bin/stampy-1.0.23/stampy.py \
        -g %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
        -h %s/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic \
        -f sam \
        -t 1 \
        --bamkeepgoodreads \
        -M {} > ../stampy-alignment/{}.stampy.sam' ::: %s-*.bwa.bam

        find . -name '%s*' -not -name '%s.bwa.sam' -type f -delete

        cd ../stampy-alignment/
        ~/bin/parallel-20150122/src/parallel 'rename .bwa.bam "" {}' ::: %s-*.bwa.bam.stampy.sam

        ~/bin/parallel-20150122/src/parallel '~/bin/samtools-1.1/samtools view \
        -bS {} > {}.bam' ::: %s-*.stampy.sam
        ~/bin/parallel-20150122/src/parallel 'rename .sam "" {}' ::: %s-*stampy.sam.bam

        ~/bin/samtools-1.1/samtools cat \
        -h %s.header.bwa.sam \
        -o %s.stampy.bam \
        %s-00.stampy.bam \
        %s-01.stampy.bam \
        %s-02.stampy.bam \
        %s-03.stampy.bam \
        %s-04.stampy.bam \
        %s-05.stampy.bam \
        %s-06.stampy.bam \
        %s-07.stampy.bam \
        %s-08.stampy.bam \
        %s-09.stampy.bam \
        %s-10.stampy.bam \
        %s-11.stampy.bam \
        %s-12.stampy.bam \
        %s-13.stampy.bam \
        %s-14.stampy.bam \
        %s-15.stampy.bam

        ~/bin/samtools-1.1/samtools flagstat \
        %s.stampy.bam > %s.stampy.bam.flagstat

        rm %s*.sam
       	rm %s-*.stampy.bam""" %\
        (Sample, Sample, i, Sample, Sample,
        i, Sample, Sample,
        i, Sample, i, i, i, Sample, Sample,
        Sample, Sample,
        Sample, Sample,
        RefDir, RefDir, Sample,
        Sample, Sample,
        Sample,
        Sample, Sample,
        Sample, Sample, Sample, Sample, Sample, Sample, Sample, Sample, Sample,
        Sample, Sample, Sample, Sample, Sample, Sample, Sample, Sample, Sample,
        Sample, Sample,
        Sample, Sample)

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
