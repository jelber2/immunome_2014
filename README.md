#immunome_2014
========
###README.md for analytical methods for following manuscript:
    Elbers, J.P., R.W. Clostio, and S.S. Taylor (2016) Population genetic
    inferences using immune gene SNPs mirror patterns inferred by microsatellites.
    Intended for Molecular Ecology.
####Python scripts for running PBS job submissions on LSU's SuperMikeII cluster.
####Also bash code for analyzing resulting data on a CentOS 6.5 machine.
#####NOTE:
    Microsatellite data used the following population abbreviations:
    FGP, SD, GG, FL
    For simplicity, we renamed population abbreviations to follow the state they
    occurred in, so populations are renamed (in order):
    LA, AL, GA, FL respectively.
========
#STEPS FOR QUALITY CONTROL, MAPPING, & SNP CALLING
##Download fastq.gz.zip files for the two MiSeq runs from BaseSpace
###Copy data to supermikeII
    #run1 = miseq data from 9Sep2014
    rsync --archive --stats --progress /work/jelber2/immunome_2014/run1/analysis_13944931_fastq.zip \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/run1/ -n
    #run2 = miseq data from 15Sep2014
    rsync --archive --stats --progress /work/jelber2/immunome_2014/run2/analysis_14120117_fastq.zip \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/run2/ -n
    # create folder for combined data
    cd /work/jelber2/immunome_2014/
    mkdir combined
###Unzip data (on SuperMikeII)
    #run1 data
    cd /work/jelber2/immunome_2014/run1/
    unzip /work/jelber2/immunome_2014/run1/analysis_13944931_fastq.zip
    cd Data/Intensities/BaseCalls
    mkdir /work/jelber2/immunome_2014/run1/fastq
    mv * /work/jelber2/immunome_2014/run1/fastq
    cd /work/jelber2/immunome_2014/run1/
    rm -r Data/
    #run2 data
    cd /work/jelber2/immunome_2014/run2/
    unzip /work/jelber2/immunome_2014/run2/analysis_14120117_fastq.zip
    cd Data/Intensities/BaseCalls
    mkdir /work/jelber2/immunome_2014/run2/fastq
    mv * /work/jelber2/immunome_2014/run2/fastq
    cd /work/jelber2/immunome_2014/run2/
    rm -r Data/
###Rename the data
    #rename run1 files
    cd /work/jelber2/immunome_2014/run1/fastq/
    #remove Undetermined reads
    rm Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz
    #make file list
    ls *.fastq.gz > convert_base_filenames.txt
    #use regular expressions to create mv command to rename file
    perl -pe "s/(\w+)_(S\d+)_(L001)_(R\d)_(001).(fastq.gz)\n/mv \1_\2_\3_\4_\5.\6 \1-\4.\6\n/" \
    convert_base_filenames.txt > convert_base_filenames2.txt
    #open the file with less and select everything with the mouse then press ctrl+shift+c
    less convert_base_filenames2.txt
    #press q to quit, then paste command into shell with ctrl+shift+v
    #rename run2 files
    cd /work/jelber2/immunome_2014/run2/fastq/
    #remove Undetermined reads
    rm Undetermined_S0_L001_R1_001.fastq.gz Undetermined_S0_L001_R2_001.fastq.gz
    #make file list
    ls *.fastq.gz > convert_base_filenames.txt
    #use regular expressions to create mv command to rename file
    perl -pe "s/(\w+)_(S\d+)_(L001)_(R\d)_(001).(fastq.gz)\n/mv \1_\2_\3_\4_\5.\6 \1-\4.\6\n/" \
    convert_base_filenames.txt > convert_base_filenames2.txt
    #open the file with less and select everything with the mouse then press ctrl+shift+c
    less convert_base_filenames2.txt
    #press q to quit, then paste command into shell with ctrl+shift+v
##Install program and get reference genome
###trimmomatic-0.32
    cd /home/jelber2/bin/
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip
    unzip Trimmomatic-0.32.zip
    mv Trimmomatic-0.32.zip Trimmomatic-0.32
    #PATH=~/home/jelber2/bin/Trimmomatic-0.32/trimmomatic-0.32.jar
###bbmerge-5.4 (part of bbmap-34.33)
    cd /home/jelber2/bin/
    mkdir bbmap-34.33
    cd bbmap-34.33/
    wget http://downloads.sourceforge.net/project/bbmap/BBMap_34.33.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2F%3Fsource%3Ddlp&ts=1421955805&use_mirror=iweb
    mv BBMap_34.33.tar.gz\?r\=http\:%2F%2Fsourceforge.net%2Fprojects%2Fbbmap%2F\?source\=dlp BBMap_34.33.tar.gz
    tar xzf BBMap_34.33.tar.gz
    cd bbmap/
    mv * ..
    cd ..
    rm -r bbmap
    #PATH=~/bin/bbmap-34.33/bbmerge.sh
###bwa-0.7.12
    cd /home/jelber2/bin/
    wget https://github.com/lh3/bwa/archive/0.7.12.tar.gz
    mv 0.7.12 bwa-0.7.12.tar.gz
    tar xzf bwa-0.7.12.tar.gz
    mv bwa-0.7.12.tar.gz bwa-0.7.12
    cd bwa-0.7.12/
    make
    #PATH=~/bin/bwa-0.7.12/bwa
###stampy-1.0.23
    cd /home/jelber2/bin/
    wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
    tar xzf Stampy-latest.tgz
    cd stampy-1.0.23
    make
    #PATH=~/bin/stampy-1.0.23/stampy.py
###java jre1.7.0
    #had to download using firefox on my Centos machine
    #saved in /home/jelber2/bin/
    rsync --stats --archive --progress /home/jelber2/bin/jre-7-linux-x64.tar.gz jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    tar xzf jre-7-linux-x64.tar.gz
    mv jre-7-linux-x64.tar.gz jre1.7.0
    #add  PATH += /home/jelber2/bin/jre1.7.0/bin  to .soft file
    nano ~/.soft
    #then resoft
    #resoft
    #PATH=~/bin/jre1.7.0/bin
###picard-1.128
    #on my Centos machine
    cd /home/jelber2/bin/
    wget https://github.com/broadinstitute/picard/releases/download/1.128/picard-tools-1.128.zip
    rsync --stats --archive --progress /home/jelber2/bin/picard-tools-1.128.zip jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    unzip picard-tools-1.128.zip
    mv picard-tools-1.128.zip picard-tools-1.128
    #PATH=~/bin/picard-tools-1.128/picard.jar
###GATK-3.3.0
    #had to download using firefox on my Centos machine
    #saved in /home/jelber2/bin/GATK-3.3.0
    rsync --stats --archive --progress /home/jelber2/bin/GATK-3.3.0/ jelber2@mike.hpc.lsu.edu:/home/jelber2/bin/GATK-3.3.0/ -n
    #switched to SuperMikeII
    cd /home/jelber2/bin/
    cd GATK-3.3.0
    tar xjf GenomeAnalysisTK-3.3-0.tar.bz2
    #PATH=~/bin/GATK-3.3.0/GenomeAnalysisTK.jar
###samtools-1.1
    cd /home/jelber2/bin/
    wget http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F1.1%2F&ts=1421967581&use_mirror=softlayer-dal
    tar xjf samtools-1.1.tar.bz2 
    mv samtools-1.1.tar.bz2 samtools-1.1
    cd samtools-1.1
    make
    nano ~/.soft #add the following line to .soft file using nano
    PATH += /home/jelber2/bin/samtools-1.1/
###parallel-20150122
    cd /home/jelber2/bin/
    wget ftp://ftp.gnu.org/gnu/parallel/parallel-20150122.tar.bz2
    tar xjf parallel-20150122.tar.bz2
    mv parallel-20150122.tar.bz2 parallel-20150122
    #PATH=~/bin/parallel-20150122/src/parallel
###Got painted turtle reference genome
    cd /work/jelber2/reference/
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
###Ran 01-make_indexes.sh on SuperMikeII to make indexes
    qsub ~/scripts/immunome_2014/01-make_indexes.sh
##Quality Control
###Ran 02-trimmomatic.py on fastq.gz in:
    # /work/jelber2/immunome_2014/run1/fastq/
    cd /work/jelber2/immunome_2014/run1/fastq/
    ~/scripts/immunome_2014/02-trimmomatic.py *.fastq.gz
    # /work/jelber2/immunome_2014/run2/fastq/
    cd /work/jelber2/immunome_2014/run2/fastq/
    #changes run1 to run2
    perl -pe "s/run1/run2/g" ~/scripts/immunome_2014/02-trimmomatic.py > \
    ~/scripts/immunome_2014/02-trimmomatic-run2.py
    python ~/scripts/immunome_2014/02-trimmomatic-run2.py *.fastq.gz
##Mapping
###Ran 03-bwa.py on trim.fastq.gz on trim.fastq.gz in:
    # /work/jelber2/immunome_2014/run1/trimmed-data/
    cd /work/jelber2/immunome_2014/run1/trimmed-data/
    ~/scripts/immunome_2014/03-bwa.py *.trim.fastq.gz
    # /work/jelber2/immunome_2014/run2/trimmed-data/
    cd /work/jelber2/immunome_2014/run2/trimmed-data/
    #changes run1 to run2
    perl -pe "s/run1/run2/g" ~/scripts/immunome_2014/03-bwa.py > \
    ~/scripts/immunome_2014/03-bwa-run2.py
    python ~/scripts/immunome_2014/03-bwa-run2.py *.trim.fastq.gz
###Ran 04-stampy.py on bwa.sam in:
    # /work/jelber2/immunome_2014/run1/bwa-alignment/
    cd /work/jelber2/immunome_2014/run1/bwa-alignment/
    ~/scripts/immunome_2014/04-stampy.py *.bwa.sam
    # /work/jelber2/immunome_2014/run2/bwa-alignment/
    cd /work/jelber2/immunome_2014/run2/bwa-alignment/
    #changes run1 to run2
    perl -pe "s/run1/run2/g" ~/scripts/immunome_2014/04-stampy.py > \
    ~/scripts/immunome_2014/04-stampy-run2.py
    python ~/scripts/immunome_2014/04-stampy-run2.py *.bwa.sam
##SNP Calling
###Ran 05a-clean_sort_addRG.py stampy.bam in:
    # /work/jelber2/immunome_2014/run1/stampy-alignment/
    cd /work/jelber2/immunome_2014/run1/stampy-alignment/
    ~/scripts/immunome_2014/05a-clean_sort_addRG.py *.stampy.bam
    # /work/jelber2/immunome_2014/run2/stampy-alignment/
    cd /work/jelber2/immunome_2014/run2/stampy-alignment/
    #changes run1 to run2 and RGID=%s_9Sep2014 to RGID=%s_15Sep2014
    perl -pe "s/run1/run2/g" ~/scripts/immunome_2014/05a-clean_sort_addRG.py |
    perl -pe "s/RGID=%s_9Sep2014/RGID=%s_15Sep2014/g" \
    > ~/scripts/immunome_2014/05a-clean_sort_addRG-run2.py
    python ~/scripts/immunome_2014/05a-clean_sort_addRG-run2.py *.stampy.bam
###Ran 05b-clean_sort_addRG_markdup_realign.py
    cd /work/jelber2/immunome_2014/run1/clean-sort-addRG/
    ~/scripts/immunome_2014/05b-clean_sort_addRG_markdup_realign.py *-CL-RG.bam
###Need to create interval list to call SNPs in the immunome_baits target region
    cd /work/jelber2/reference/
    java -Xmx8g -jar ~/bin/picard-tools-1.128/picard.jar BedToIntervalList \
    I=immunome_baits_C_picta-3.0.3.bed \
    SEQUENCE_DICTIONARY=GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.dict \
    O=immunome_baits_C_picta-3.0.3.interval.list
###Ran 06-mergeBAM_callSNPs_initial.py
    cd /work/jelber2/immunome_2014/combined/realign-around-indels/
    ~/scripts/immunome_2014/06-mergeBAM_callSNPs_initial.py *-realigned.bam
###Need to get gsalib for R in order to produce recalibration plots
    #see link below for further details
    #http://gatkforums.broadinstitute.org/discussion/1244/what-is-a-gatkreport
    #on SuperMikeII
    R
    install.packages("gsalib")
    #YOU WILL GET THE FOLLOWING
    #Warning in install.packages("gsalib") :
    #'lib = "/home/packages/R/2.15.1/gcc-4.4.6/lib64/R/library"' is not writable
    #Would you like to use a personal library instead?  (y/n) type y
    #Would you like to create a personal library
    #~/R/x86_64-unknown-linux-gnu-library/2.15
    #to install packages into?  (y/n) type y
    #quit R
    q()
    ##Hmm, the above trick didn't work!
    ##Will have to generate the plots on my Centos machine
###Ran 07-qual_score_recal01.py
    cd /work/jelber2/immunome_2014/combined/realign-around-indels/
    #excludes file ALL-samples-realigned.bam
    find . -name '*-realigned.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_2014/07-qual_score_recal01.py {} \;
###Ran 08-mergeBAM_callSNPs_recal01.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal01/
    ~/scripts/immunome_2014/08-mergeBAM_callSNPs_recal01.py *-recal01.bam
###Ran 09-qual_score_recal02.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal01/
    #excludes file ALL-samples-recal01.bam
    find . -name '*-recal01.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_2014/09-qual_score_recal02.py {} \;
###Ran 10-mergeBAM_callSNPs_recal02.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal02/
    ~/scripts/immunome_2014/10-mergeBAM_callSNPs_recal02.py *-recal02.bam
###Ran 11-qual_score_recal03.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal02/
    #excludes file ALL-samples-recal02.bam
    find . -name '*-recal02.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_2014/11-qual_score_recal03.py {} \;
###Ran 12-mergeBAM_callSNPs_recal03.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal03/
    ~/scripts/immunome_2014/12-mergeBAM_callSNPs_recal03.py *-recal03.bam
###Calculate seq_metrics.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal03/
    java -Xmx4g -jar ~/bin/picard-tools-1.128/picard.jar CalculateHsMetrics \
    BAIT_INTERVALS=/work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    BAIT_SET_NAME=Immunome \
    TARGET_INTERVALS=/work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    METRIC_ACCUMULATION_LEVEL=SAMPLE \
    R=/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    I=ALL-samples-recal03.bam \
    O=ALL-samples-recal03-baitsonly.hsmetrics.txt
###plot coverage for each sample
####run bedtools coverage on all bam files, then keep only lines with 'all' on them
    #FROM http://gettinggeneticsdone.blogspot.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
    #ALSO FROM https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir plot_coverage
    cd plot_coverage/
    #command below takes 5-10 minutes
    while read i
    do
    ~/bin/bedtools-2.22.1/bin/bedtools coverage -abam ../call-SNPs-recal03/$i-recal03.bam \
    -b /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.bed \
    -hist | grep ^all > $i.baitcoverage.all.txt
    done < ../call-SNPs-recal03/samplelist
###Need to use featureCounts to summarize number of genes, reads per gene, etc
####Get Subread
    #featureCounts is part of the Subread package http://bioinf.wehi.edu.au/featureCounts/
    cd ~/bin/
    wget http://downloads.sourceforge.net/project/subread/subread-1.4.6/subread-1.4.6-Linux-x86_64.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.4.6%2F&ts=1423446986&use_mirror=iweb
    mv subread-1.4.6-Linux-x86_64.tar.gz?r=http:%2F%2Fsourceforge.net%2Fprojects%2Fsubread%2Ffiles%2Fsubread-1.4.6%2F subread-1.4.6-Linux-x86_64.tar.gz
    tar xzf subread-1.4.6-Linux-x86_64.tar.gz
    mv subread-1.4.6-Linux-x86_64.tar.gz subread-1.4.6-Linux-x86_64
####Get genometools-1.5.4 to annotate introns in gff file
    cd ~/bin/
    wget http://genometools.org/pub/genometools-1.5.4.tar.gz
    tar xzf genometools-1.5.4.tar.gz
    mv genometools-1.5.4.tar.gz genometools-1.5.4
    cd genometools-1.5.4
    #on MacOSX
    make
    #on CentOS
    #install ruby first
    #become superuser
    su
    yum install ruby.x86_64
    #stop being a super user
    exit
    #make the executable using 64bit mode but without cairo
    make 64bit=yes cairo=no
    #test the install - will take a long time (>30min?)
    make 64bit=yes cairo=no test
#####Use genometools to get introns
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-seqs/
    ~/bin/genometools-1.5.4/bin/gt gff3 -addintrons yes \
    -o /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff
####Need to prefilter the GFF file for immune genes
    cd /media/immunome_2014/work/jelber2/reference/
    #intersect the gff file
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    -b immunome_baits_C_picta-3.0.3.bed \
    > GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gff
####Convert GFF file of gene annotations to GTF
    cd /media/immunome_2014/work/jelber2/reference/
    perl -pe "s/\S+=GeneID:(\d+).+/gene_id \"\1\";/g" \
    GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gff \
    > GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf
#####run subread
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/featureCounts/
######all samples at once
    #at the gene level
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o ALL.gene -F GTF -T 2 --ignoreDup \
    ../call-SNPs-recal03/ALL-samples-recal03.bam
    #at exon level
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o ALL.exon -F GTF -T 2 -f --ignoreDup \
    ../call-SNPs-recal03/ALL-samples-recal03.bam
######each sample separately
    #at the gene level
    while read i
    do
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o $i.gene -F GTF -T 2 --ignoreDup ../call-SNPs-recal03/$i-recal03.bam
    done < ../call-SNPs-recal03/samplelist
    #at the exon level
    while read i
    do
    ~/bin/subread-1.4.6-Linux-x86_64/bin/featureCounts \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic_immunome_baits.gtf \
    -o $i.exon -F GTF -T 2 -f --ignoreDup ../call-SNPs-recal03/$i-recal03.bam
    done < ../call-SNPs-recal03/samplelist
####write an R function to do the following
#####count number of possible immune genes
    wc -l AL102.gene
    #total genes = 632 (after subtracting 2 header lines)
#####count number of possible immune gene exons
    wc -l AL102.exon
    #total exons = 37275 (after subtracting 2 header lines)
#####how many different immune genes were captured
    grep -Pv "0$" ALL.gene | wc -l
    #558 (after subtracting 2 header lines)
#####how many different immune genes exons were captured
    grep -Pv "0$" ALL.exon |sort | cut -f 1| uniq -c | perl -pe "s/ +/\t/g" | perl -ane '$sum += $F[0]; END {print $sum}'
    #4430 (after subtracting 2 header lines)
#####count number of genes per sample
    while read i
    do
    test=$(grep -Pv "0$" $i.gene | wc -l)
    echo -e $i'\t'$test > $i.gene.count
    done < ../call-SNPs-recal03/samplelist
    cat *.gene.count > gene.counts.per.sample
#####count number of exons per sample
    while read i
    do
    grep -Pv "0$" $i.exon |sort | cut -f 1| uniq -c | \
    perl -pe "s/ +/\t/g" | \
    perl -ane '$sum += $F[0]; END {print $sum; print "\n"}' > $i.exon.count
    done < ../call-SNPs-recal03/samplelist
    cat *.exon.count > exon.counts.per.sample
###Ran 14-haplotypecaller.py
    cd /work/jelber2/immunome_2014/combined/call-SNPs-recal03/
    #excludes file ALL-samples-recal03.bam
    find . -name '*-recal03.bam' -not -name 'ALL-samples-*' \
    -exec ~/scripts/immunome_2014/14-haplotypecaller.py {} \;
###Ran GenotypeGVCFs to perform joint genotyping
    #on Cenots machine
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/hc/
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -stand_call_conf 30 \
    -stand_emit_conf 10 \
    --max_alternate_alleles 32 \
    --variant AL102-raw-snps-indels.vcf \
    --variant AL103-raw-snps-indels.vcf \
    --variant AL106-raw-snps-indels.vcf \
    --variant AL108-raw-snps-indels.vcf \
    --variant FL846-raw-snps-indels.vcf \
    --variant FL855-raw-snps-indels.vcf \
    --variant FL857-raw-snps-indels.vcf \
    --variant FL880-raw-snps-indels.vcf \
    --variant GG1044-raw-snps-indels.vcf \
    --variant GG1435-raw-snps-indels.vcf \
    --variant GG1835-raw-snps-indels.vcf \
    --variant GG462-raw-snps-indels.vcf \
    --variant LA62-raw-snps-indels.vcf \
    --variant LA66-raw-snps-indels.vcf \
    --variant LA77-raw-snps-indels.vcf \
    --variant LA78-raw-snps-indels.vcf \
    -o ALL-samples-raw-snps-indels.vcf
###Added expressions to filter variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -L /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    -V ALL-samples-raw-snps-indels.vcf \
    --clusterWindowSize 10 \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
    --filterName "Bad_Validation" \
    --filterExpression "QUAL < 30.0" \
    --filterName "LowQual" \
    --genotypeFilterExpression "DP < 10.0" \
    --genotypeFilterName "Low_Read_Depth_Over_Sample" \
    --genotypeFilterExpression "GQ < 20.0" \
    --genotypeFilterName "Low_GenotypeQuality" \
    -o ALL-samples-Q30-snps-indels.vcf
####Got only Indel variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-indels.vcf \
    -selectType INDEL
####Got only SNP variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-snps.vcf \
    -selectType SNP
###Got "Truthing" SNPs for Variant Quality Score Recalibration
####First got concordant SNPs between HaplotypeCaller and UnifiedGenotyper
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    --variant ALL-samples-Q30-snps.vcf \
    --concordance ../call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
    -o concordant-snps-HCvsUG.vcf
#####Next got only SNPs passing filters
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    --variant concordant-snps-HCvsUG.vcf \
    -o concordant-snps-HCvsUG-PASS.vcf \
    --excludeFiltered
####Second got concordant Indels between HaplotypeCaller and UnifiedGenotyper
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T SelectVariants \
    --variant ../hc/ALL-samples-Q30-indels.vcf \
    --concordance ../call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
    -o ../hc/concordant-indels-HCvsUG.vcf
#####Next got only indels passing filters
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    --variant ../hc/concordant-indels-HCvsUG.vcf \
    -o ../hc/concordant-indels-HCvsUG-PASS.vcf \
    --excludeFiltered
#####Finally replaced 'PASS' with 'INDEL'
    perl -pe "s/PASS/INDEL/g" \
    ../hc/concordant-indels-HCvsUG-PASS.vcf \
    > ../hc/concordant-indels-HCvsUG-PASS-renamed-INDEL.vcf
####Ran Variant Recalibrator
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir vqsr
    cd vqsr
#####Recalibrated snps
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 ../hc/concordant-snps-HCvsUG-PASS.vcf \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
    -recalFile VQSR-snps.recal \
    -mode SNP \
    -tranchesFile VQSR-snps.tranches \
    -rscriptFile VQSR-snps.plots.R \
    --maxGaussians 4
#####Recalibrated indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 ../hc/concordant-indels-HCvsUG-PASS.vcf \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
    -recalFile VQSR-indels.recal \
    -mode INDEL \
    -tranchesFile VQSR-indels.tranches \
    -rscriptFile VQSR-indels.plots.R \
    --maxGaussians 2
    #COULD NOT GET TO WORK with --maxGaussians 4 or 3, presumably because the
    #indel dataset (~200 indels) is too small.
####Applied the recalibration on snps
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    --ts_filter_level 99.5 \
    -tranchesFile VQSR-snps.tranches \
    -recalFile VQSR-snps.recal \
    -o ALL-samples-Q30-snps-recal.vcf
####Applied the recalibration on indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    --ts_filter_level 99.0 \
    -tranchesFile VQSR-indels.tranches \
    -recalFile VQSR-indels.recal \
    -o ALL-samples-Q30-indels-recal.vcf
###Needed to use beagle to improve SNPs (using Linkage Disequilibrium) called by Haplotype Caller
####Downloaded beagle
    cd ~/bin
    wget http://faculty.washington.edu/browning/beagle/beagle.r1398.jar
####Ran beagle on snps and indels separately
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir beagle
    cd beagle
    #snps
    java -Xmx8000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_2014/combined/vqsr/ALL-samples-Q30-snps-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_2014/combined/beagle/ALL-samples-Q30-snps-recal-beagle
    #indels
    java -Xmx8000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_2014/combined/vqsr/ALL-samples-Q30-indels-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_2014/combined/beagle/ALL-samples-Q30-indels-recal-beagle
========
#STEPS FOR POPGEN
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir popgen
    mkdir popgen-msat
##For SNPs
    cd popgen
###Get only polymorphic loci
####Downloaded vcftools
    cd /home/jelber2/bin/
    wget http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1411515317&use_mirror=superb-dca2
    tar -xzf vcftools_0.1.12b.tar.gz 
    mv vcftools_0.1.12b.tar.gz vcftools_0.1.12b
    cd vcftools_0.1.12b/
    nano ~/.soft #add the following two lines to using nano.soft file
    PATH+=/home/jelber2/bin/tabix-0.2.6
    PERL5LIB = /home/jelber2/bin/vcftools_0.1.12b/perl
    resoft #to refresh soft file
    cd /home/jelber2/bin/vcftools_0.1.12b/
    make #compile vcftools
    # Path to vcftools executable
    /home/jelber2/bin/vcftools_0.1.12b/bin/vcftools
####Remove loci with AF=1
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/
    #removes loci with AF=1
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle.vcf.gz | grep -v "AF=1" \
    > ALL-samples-Q30-snps-recal-beagle2.vcf
    #still leaves some unwanted non-polymorphic SNPs?
    #try calculating allele frequencies then
    ~/bin/vcftools_0.1.12b/bin/vcftools \
    --vcf ALL-samples-Q30-snps-recal-beagle2.vcf \
    --freq --out ALL-samples-Q30-snps-recal-beagle
    #get only non-polymorphic loci
    grep ":1" ALL-samples-Q30-snps-recal-beagle.frq | cut -f 1-2 > nonpolymorphicloci
    grep -v "\s2\s" ALL-samples-Q30-snps-recal-beagle.frq | cut -f 1-2 > multiallelicloci
    cat nonpolymorphicloci multiallelicloci > multiallelic_or_nonpolymorphicloci
    #calculate number of di-,tri-,tetra-allelic loci
    cut -f 3 ALL-samples-Q30-snps-recal-beagle.frq | sort | uniq -c
    #    di = 20947 (includes non-polymorphic loci = 3046)
    #   tri = 758
    # tetra = 7
    #
    #filter out multiallelic_or_nonpolymorphicloci
    #might take a few minutes
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/remove\tlocus\t\n/" ALL-samples-Q30-snps-recal-beagle2.vcf
    done < multiallelic_or_nonpolymorphicloci
    grep -v "remove" ALL-samples-Q30-snps-recal-beagle2.vcf \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
###Check loci for linkage disequilibrium and Hardy-Weinberg Equilibrium
####Had to make populations.txt file and population-specific files for vcf filtering
    # e.g., format (\t=tab)
    #Sample\tpopulation1
    #populations.txt, so you have only torts from each population in a file, one per line
    #now get population-specific files
    grep "AL" populations.txt > Alabama
    grep "GG" populations.txt > Georgia
    grep "FL" populations.txt > Florida
    grep "LA" populations.txt > Louisiana
####Hardy-Weinberg Equilibrium test
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.AL --keep Alabama
#####R function to count number of sites out of HWE (i.e., p_HWE < 0.05 after fdr correction)
    R
    ALhwe <-read.table(file ="hwe.AL.hwe", header = TRUE)
    ALhwe.fdr <- p.adjust(p = ALhwe$P_HWE, method = "fdr", n = length(ALhwe$P_HWE))
    summary(ALhwe.fdr)
    LAhwe <-read.table(file ="hwe.LA.hwe", header = TRUE)
    LAhwe.fdr <- p.adjust(p = LAhwe$P_HWE, method = "fdr", n = length(LAhwe$P_HWE))
    summary(LAhwe.fdr)
    FLhwe <-read.table(file ="hwe.FL.hwe", header = TRUE)
    FLhwe.fdr <- p.adjust(p = FLhwe$P_HWE, method = "fdr", n = length(FLhwe$P_HWE))
    summary(FLhwe.fdr)
    GGhwe <-read.table(file ="hwe.GG.hwe", header = TRUE)
    GGhwe.fdr <- p.adjust(p = GGhwe$P_HWE, method = "fdr", n = length(GGhwe$P_HWE))
    summary(GGhwe.fdr)
####outputs linkage disequilibrium pvalues
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --geno-chisq --out geno.AL --keep Alabama
#####R function to count number of sites out of linkage equi (i.e., PVAL < 0.05 after fdr correction)
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/popgen")
    ALld = na.omit(read.table(file ="geno.AL.geno.chisq", header = TRUE))
    ALld.fdr=p.adjust(p = ALld$PVAL, method = "fdr", n = length(ALld$PVAL))
    summary(ALld.fdr)
    #summary(ALld$PVAL)
    GGld = na.omit(read.table(file ="geno.GG.geno.chisq", header = TRUE))
    GGld.fdr=p.adjust(p = GGld$PVAL, method = "fdr", n = length(GGld$PVAL))
    summary(GGld.fdr)
    #summary(GGld$PVAL)
    FLld = na.omit(read.table(file ="geno.FL.geno.chisq", header = TRUE))
    FLld.fdr=p.adjust(p = FLld$PVAL, method = "fdr", n = length(FLld$PVAL))
    summary(FLld.fdr)
    #summary(FLld$PVAL)
    LAld = na.omit(read.table(file ="geno.LA.geno.chisq", header = TRUE))
    LAld.fdr=p.adjust(p = LAld$PVAL, method = "fdr", n = length(LAld$PVAL))
    summary(LAld.fdr)
    #summary(LAld$PVAL)
###Calculate nucleotide diversity (pi per site)
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.AL --keep Alabama
###Calculate Tajima's D
    #need to figure out "correct" window size
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.120.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 1200 --out tajimaD.1200.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 2400 --out tajimaD.2400.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.AL --keep Alabama
###Convert VCF file to structure, fstat, and smartpca input
####VCF to Structure
#####use PGDSpider
    #get the program
    cd ~/bin/
    wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.4.zip
    unzip PGDSpider_2.0.7.4.zip
    mv PGDSpider_2.0.7.4.zip PGDSpider_2.0.7.4
    #use following command to generate spider.conf.xml and spid file in ../popgen directory
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/structure-input.txt \
    -outputformat STRUCTURE
    #edit spider.conf.xml
    nano ~/bin/PGDSpider_2.0.7.4/spider.conf.xml #to add path to samtools
    #change <entry key="PathBcftools"></entry>
    #to <entry key="PathBcftools">/home/jelber2/bin/samtools-0.1.19/bcftools/bcftools</entry>
    #change <entry key="PathSamtools"></entry>
    #to <entry key="PathSamtools">/home/jelber2/bin/samtools-0.1.19/samtools</entry>
    #save and exit
    #edit the spid file
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/template_VCF_STRUCTURE.spid
    #contents after editing, minus the leading spaces
        # spid-file generated: Wed Jan 14 18:47:06 CST 2015
        # VCF Parser questions
        PARSER_FORMAT=VCF
        # Do you want to include a file with population definitions?
        VCF_PARSER_POP_QUESTION=true
        # Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
        VCF_PARSER_REGION_QUESTION=
        # What is the ploidy of the data?
        VCF_PARSER_PLOIDY_QUESTION=DIPLOID
        # Only output following individuals (ind1, ind2, ind4, ...):
        VCF_PARSER_IND_QUESTION=
        # Output genotypes as missing if the read depth of a position for the sample is below:
        VCF_PARSER_READ_QUESTION=
        # Take most likely genotype if "PL" or "GL" is given in the genotype field?
        VCF_PARSER_PL_QUESTION=
        # Do you want to exclude loci with only missing data?
        VCF_PARSER_EXC_MISSING_LOCI_QUESTION=
        # Select population definition file:
        VCF_PARSER_POP_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/populations.txt
        # Only output SNPs with a phred-scaled quality of at least:
        VCF_PARSER_QUAL_QUESTION=
        # Do you want to include non-polymorphic SNPs?
        VCF_PARSER_MONOMORPHIC_QUESTION=false
        # Output genotypes as missing if the phred-scale genotype quality is below:
        VCF_PARSER_GTQUAL_QUESTION=
        #
        # STRUCTURE Writer questions
        WRITER_FORMAT=STRUCTURE
        # Save more specific fastSTRUCTURE format?
        STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
        # Specify the locus/locus combination you want to write to the STRUCTURE file:
        STRUCTURE_WRITER_LOCUS_COMBINATION_QUESTION=
        # Specify which data type should be included in the STRUCTURE file  (STRUCTURE can only analyze one data type per file):
        STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
        # Do you want to include inter-marker distances?
        STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
    #saved as vcf2structure.spid
#####Do VCF to STRUCTURE file conversion
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/structure-input.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2structure.spid > structure-input.log
####VCF to GENEPOP
#####create the spid file
    #replace the STRUCTURE section of vcf2structure.spid with the following for GENEPOP
    #minus the leading spaces
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2structure.spid
        # GENEPOP Writer questions
        WRITER_FORMAT=GENEPOP
        # Specify the locus/locus combination you want to write to the GENEPOP file:
        GENEPOP_WRITER_LOCUS_COMBINATION_QUESTION=
        # Specify which data type should be included in the GENEPOP file  (GENEPOP can only analyze one data type per file):
        GENEPOP_WRITER_DATA_TYPE_QUESTION=SNP
    #saved as vcf2genepop.spid
#####Do VCF to GENEPOP file conversion
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/genepop-input.txt \
    -outputformat GENEPOP \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2genepop.spid > genepop-input.log
####VCF to FSTAT
#####create the spid file
    #replace the STRUCTURE section of vcf2structure.spid with the following for FSTAT
    #minus the leading spaces
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2structure.spid
        # FSTAT Writer questions
        WRITER_FORMAT=FSTAT
        # Specify which data type should be included in the FSTAT file  (FSTAT can only analyze one data type per file):
        FSTAT_WRITER_DATA_TYPE_QUESTION=SNP
        # Save label file
        FSTAT_WRITER_LABEL_FILE_QUESTION=
        # Do you want to save an additional file with labels (population names)?
        FSTAT_WRITER_INCLUDE_LABEL_QUESTION=false
        # Specify the locus/locus combination you want to write to the FSTAT file:
        FSTAT_WRITER_LOCUS_COMBINATION_QUESTION=
    #saved as vcf2fstat.spid
#####Do VCF to FSTAT file conversion
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/fstat-input.txt \
    -outputformat FSTAT \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2fstat.spid > fstat-input.log
##Look at population structure with structure
###Installed structure from source, to run on SuperMikeII
####Get source
    cd ~/bin/
    mkdir structure
    cd structure
    wget http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/structure_kernel_source.tar.gz
    tar xzf structure_kernel_source.tar.gz
    cd structure_kernel_src/
####compile
    make
####used the following setting for mainparams file created using Windows front-end version
        #define OUTFILE /work/jelber2/immunome_2014/combined/popgen/structure-results
        #define INFILE /work/jelber2/immunome_2014/combined/popgen/structure-input.txt
        #define NUMINDS 16
        #define NUMLOCI 17901
        #define LABEL 1
        #define POPDATA 1
        #define POPFLAG 0
        #define LOCDATA 0
        #define PHENOTYPE 0
        #define MARKERNAMES 1
        #define MAPDISTANCES 0
        #define ONEROWPERIND 0
        #define PHASEINFO 0
        #define PHASED 0
        #define RECESSIVEALLELES 0
        #define EXTRACOLS 0
        #define MISSING 
        #define PLOIDY 2
        #define MAXPOPS 1
        #define BURNIN 50000
        #define NUMREPS 100000
        #define NOADMIX 0
        #define LINKAGE 0
        #define USEPOPINFO 0
        #define LOCPRIOR 0
        #define INFERALPHA 1
        #define ALPHA 1.0
        #define POPALPHAS 0
        #define UNIFPRIORALPHA 1
        #define ALPHAMAX 10.0
        #define ALPHAPROPSD 0.025
        #define FREQSCORR 1
        #define ONEFST 0
        #define FPRIORMEAN 0.01
        #define FPRIORSD 0.05
        #define INFERLAMBDA 0
        #define LAMBDA 1.0
        #define COMPUTEPROB 1
        #define PFROMPOPFLAGONLY 0
        #define ANCESTDIST 0 
        #define STARTATPOPINFO 0
        #define METROFREQ 10
        #define UPDATEFREQ 1
####Ran structure using the following command for k1,k2,k3,k4,k5
    cd /work/jelber2/immunome_2014/combined/popgen/
    ~/bin/structure/structure_kernel_src/structure \
    -m mainparams.test.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    #etc.
    #note we used default settings for extraparams
    #(i.e., the correlated allele frequency and the admixture ancestry models)
    #implemented on SuperMike II using /home/jelber2/scripts/immunome_2014/16-structure.py
####Used STRUCTURE HARVESTER web v0.6.94 to select best K values
####Used CLUMPAK web to visualize population assignments
    http://clumpak.tau.ac.il/
###Calculated basic genetic stats and performed pca with r package hierfstat
    #analyses detailed in hierfstat.r
##For msats
###Create input files
####Original
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats
    #original file = 11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-alldata-jpe-2015-02-25.arp
####Same msats populations as for SNPs
    #kept only LA,SD,GG,FL populations, named file as:
    #11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp
####Same msats populations and individuals as for SNPs
    #used only same individuals (minus GG population, which did not contain GG1044, GG1435, GG1835)
    #so used random selection of 3 others from within GG, named file as:
    #11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-same-jpe-2015-02-25.arp
    #used following code to randomly choose 3 torts from GG
    grep -Po "GG\d+" \
    11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp | \
    sort -R | head -n 3 > GG.3random
####Same populations for msats and SNPs but 4 random individuals per pop
    grep -Po "LA\d+" \
    11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp | \
    sort -R | head -n 4 > LA.4random
    grep -Po "SD\d+" \
    11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp | \
    sort -R | head -n 4 > SD.4random
    grep -Po "GG\d+" \
    11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp | \
    sort -R | head -n 4 > GG.4random
    grep -Po "FL\d+" \
    11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp | \
    sort -R | head -n 4 > FL.4random
###Get only polymorphic loci
    #they are all polymorphic
###Check loci for linkage disequilibrium and Hardy-Weinberg Equilibrium
    #see Rachel's paper that the loci are in linkage equilibrium and HWE
###Convert msat loci from ARLEQUIN to FSTAT,STRUCTURE FORMATS
####ARLEQUIN to FSTAT
    #first created template_SPID file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp \
    -inputformat ARLEQUIN \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/fstat-msat-input.txt \
    -outputformat FSTAT
    #edited template_ARLEQUIN_FSTAT.spid
    #saved the followin as arlequin2fstat.spid (minus the leading spaces)
        # spid-file generated: Wed Feb 25 15:09:01 CST 2015
        # Arlequin Parser questions
        PARSER_FORMAT=ARLEQUIN
        #
        # FSTAT Writer questions
        WRITER_FORMAT=FSTAT
        # Specify which data type should be included in the FSTAT file  (FSTAT can only analyze one data type per file):
        FSTAT_WRITER_DATA_TYPE_QUESTION=MICROSAT
        # Save label file
        FSTAT_WRITER_LABEL_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/fstat-msats-populations.txt
        # Do you want to save an additional file with labels (population names)?
        FSTAT_WRITER_INCLUDE_LABEL_QUESTION=true
        # Specify the locus/locus combination you want to write to the FSTAT file:
        FSTAT_WRITER_LOCUS_COMBINATION_QUESTION=
    #convert the file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/fstat-msats-input.txt \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/arlequin2fstat.spid > arlequin2fstat.log
####ARLEQUIN to STRUCTURE
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp \
    -inputformat ARLEQUIN \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/structure-msats-input.txt \
    -outputformat STRUCTURE
    #replace the FSTAT section of arlequin2fstat.spid with the following for STRUCTURE
    #minus the leading spaces
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/arlequin2fstat.spid
        # STRUCTURE Writer questions
        WRITER_FORMAT=STRUCTURE
        # Save more specific fastSTRUCTURE format?
        STRUCTURE_WRITER_FAST_FORMAT_QUESTION=false
        # Specify the locus/locus combination you want to write to the STRUCTURE file:
        STRUCTURE_WRITER_LOCUS_COMBINATION_QUESTION=
        # Specify which data type should be included in the STRUCTURE file  (STRUCTURE can only analyze one data type per file):
        STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
        # Do you want to include inter-marker distances?
        STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
    #saved files as arlequin2structure.spid
    #convert the file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.4/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/11.03.09-RWC-arlequin-inpt-for-animal-conserv-paper-LA_SD_GG_FLdata-jpe-2015-02-25.arp \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/structure-msats-input.txt \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/arlequin2structure.spid > arlequin2structure.log
    #add loci names
    sed -i 1i"\\\t\t1 2 3 4 5 6 7 8 9 10" structure-msats-input.txt
###Look at population structure with structure
    #copy mainparams files
    cp ../popgen/mainparams.test.0* .
    #modify the directory and number of loci with perl
    perl -pi -e "s/popgen/popgen-msats/g" mainparams.test.0*
    perl -pi -e "s/structure-input/structure-msats-input/g" mainparams.test.0*
    perl -pi -e "s/structure-results-/structure-results-msats-/g" mainparams.test.0*
    perl -pi -e "s/17901/10/g" mainparams.test.0*
    perl -pi -e "s/\/work/\/media\/immunome_2014\/work/g" mainparams.test.0*
    perl -pi -e "s/16/101/g" mainparams.test.0*
####Ran structure using the following command for k1,k2,k3,k4
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen-msats/
    ls mainparams.test.0* > params
    while read i
    do
    ~/bin/structure/structure_kernel_src/structure \
    -m $i \
    -e ~/bin/structure/structure_kernel_src/extraparams
    done < params
####Used STRUCTURE HARVESTER web v0.6.94 to select best K values
####Used CLUMPAK web to visualize population assignments
    http://clumpak.tau.ac.il/
###Calculated basic genetic stats with r package hierfstat and also pca
    #analyses detailed in hierfstat-msats.r
========
#STEPS FOR GETTING NAMES OF GENES WITH POLYMORPHIC SNPs
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir snp_genes
    cd snp_genes
    # use bedtools to intersect gene annotations with polymorphic snps
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff \
    -b ../popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf | \
    # use awk to get only "genes"
    awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' - > snp_genes.gff
    # use perl to get only gene names, then sort them, and get only unique names
    perl -pe "s/\S+=GeneID:(\d+).+/\1/g" snp_genes.txt |sort | uniq > snp_genes2.txt
    #use R to get the descriptions
    #the following is implemented in:
    # C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/immunome_2014/snps_gene_names.R
    R
    library("genomes")
    setwd("C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/immunome_2014/")
    setwd("/Users/jelbers/Documents/Documents/LSU/Dissertation/Manuscripts/immunome_2014/")
    input <- read.table("snps-genes.names.count.txt")
    Sys.setenv(email="jelber2@lsu.edu")
    #
    output <- genomes::efetch(id = input$V2,
                              "gene",
                              "gb",
                              "xml")
    output2 <- output[grepl("<Gene-ref_desc>.+</Gene-ref_desc>",
                            output,
                            perl=TRUE)]
    output3 <- sub("\\s+<Gene-ref_desc>(.+)</Gene-ref_desc>",
                   "\\1",
                   output2,
                   perl=TRUE)
    allgenes <- sort(unique(output3))
    MHC <- allgenes[agrepl("histo",allgenes)]
    TLR <- allgenes[grepl("toll",allgenes)]
    q()
========
#STEPS FOR VARIANT ANNOTATION
##Download Tools First
###Downloaded snpEff
    #ideally want to know if variants will affect protein structure and possibly immune gene function
    cd /work/jelber2/reference
    wget http://iweb.dl.sourceforge.net/project/snpeff/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
###Added Chrysemys_picta_bellii-3.0.3 to snpEff.config using nano
    cd /media/immunome_2014/work/jelber2/reference/snpEff
    nano snpEff.config # added the following four lines after the Capsella_rubella_v1.0 entry (remove 4 spaces on left if cut and pasting)
    # Chrysemys_picta_bellii-3.0.3
    Chrysemys_picta_bellii-3.0.3.genome : western painted turtle
    	Chrysemys_picta_bellii-3.0.3.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/
    	Chrysemys_picta_bellii-3.0.3.M.codonTable : Standard
###Created data directory for Chrysemys_picta_bellii-3.0.3 genome
    cd /media/immunome_2014/work/jelber2/reference/snpEff
    mkdir data
    cd data
    mkdir Chrysemys_picta_bellii-3.0.3
    cd Chrysemys_picta_bellii-3.0.3
    # downloaded FASTA file
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    # snpEff requires genome.fa file to be called "sequences.fa"
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz sequences.fa.gz
    # have to unzip sequences.fa.gz
    gunzip sequences.fa.gz
    # downloaded gff3 file (i.e., gene annotation file)
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    # snpEff requires gene annotation file be called "genes.gff"
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz genes.gff.gz
    # unzipped genes.gff.gz
    gunzip genes.gff.gz
    # download protein sequences
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa.gz
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa.gz protein.fa.gz
    gunzip protein.fa.gz
###Built snpEff database for Chrysemys_picta_bellii-3.0.3
    cd /media/immunome_2014/work/jelber2/reference/snpEff
    # used snpEff_build.py script to implement command below, which took < 30 minutes
    java -jar -Xmx8g /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar build -gff3 -v Chrysemys_picta_bellii-3.0.3 2>&1 | tee Chrysemys_picta_bellii-3.0.3.build
###Downloaded bcftools
    cd ~/bin/
    git clone --branch=develop git://github.com/samtools/htslib.git
    git clone --branch=develop git://github.com/samtools/bcftools.git
    cd bcftools; make
#Need to look for protein altering variants shared by samples in the same population
###Split vcf file from GATK for snpEff
    #snpEff needs ALL-samples*.vcf file split by sample (i.e., into Sample1.vcf, Sample2.vcf)
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/call-SNPs-recal03/
    ls *-recal03.bam | grep -Po '^\w+'| sort -u | grep -v 'ALL' > samplelist
    mkdir ../split-vcfs
    cd ../split-vcfs/
    cp ../popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf .
    #compress snps and index with tabix
    ~/bin/samtools-1.1/htslib-1.1/bgzip ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz
    #split files
    #code to split each vcf file
    while read i
    do
    ~/bin/bcftools/bcftools view -s $i ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf.gz > $i-snps.vcf
    done < ../call-SNPs-recal03/samplelist
###Ran snpEff on each split vcf file
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/split-vcfs/
    # command below to run snpEff on all samples in samplelist
    # not implemented on SuperMikeII b/c process was < 15 min
    while read i
    do
    java -Xmx8g -jar /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    $i-snps.vcf > $i-snps-snpeff.vcf
    mv snpEff_genes.txt $i-snps-snpeff-genes.txt
    mv snpEff_summary.html $i-snps-snpeff-summary.html
    done < ../call-SNPs-recal03/samplelist
###Ran VariantAnnotator on each snpeff file
    while read i
    do
    rm $i-snps.vcf.idx
    rm $i-snps-snpeff.vcf.idx
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -A SnpEff \
    --variant $i-snps.vcf \
    --snpEffFile $i-snps-snpeff.vcf \
    -L $i-snps.vcf \
    -o $i-snps-annotated.vcf
    done < ../call-SNPs-recal03/samplelist
###Merge split, annotated vcfs
    #compress then index split snp files
    while read i
    do
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f $i-snps-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-snps-annotated.vcf.gz
    done < ../call-SNPs-recal03/samplelist
    #merge snp files and index them
    ~/bin/bcftools/bcftools merge \
    -o ALL-samples-snps-annotated.vcf \
    -O v -m none \
     ../split-vcfs/*-snps-annotated.vcf.gz
========
