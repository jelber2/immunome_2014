#immunome_2014
========
###Python scripts for running PBS job submissions on LSU's SuperMikeII cluster.
###The jobs in this repository are for analyzing Illumina NGS reads from a target enrichment experiment.
###Probes target painted turtle immune response genes and thus capture the "immunome" of other chelonians.

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
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/popgen")
    ALhwe <-read.table(file ="hwe.AL.hwe", header = TRUE)
    ALhwde <- nrow(ALhwe[ALhwe$P_HWE < 0.05, ])
    ALhwde
    FLhwe <-read.table(file ="hwe.FL.hwe", header = TRUE)
    FLhwde <- nrow(FLhwe[FLhwe$P_HWE < 0.05, ])
    FLhwde
    GGhwe <-read.table(file ="hwe.GG.hwe", header = TRUE)
    GGhwde <- nrow(GGhwe[GGhwe$P_HWE < 0.05, ])
    GGhwde
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
####Ran structure using the following command for k1,k2,k3,k4
    cd /work/jelber2/immunome_2014/combined/popgen/
    ~/bin/structure/structure_kernel_src/structure \
    -m mainparams.test.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    #etc.
    #note we used default settings for extraparams
    #(i.e., the correlated allele frequency and the admixture ancestry models)
    #implemented on SuperMike II using /home/jelber2/scripts/immunome_2014/16-structure.py
####Used STRUCTURE HARVESTER weeb v0.6.94 to select best K values
####Used CLUMPP v1.1.2b to average data from multiple runs
    #get CLUMPP
    wget http://web.stanford.edu/group/rosenberglab/software/CLUMPP_Linux64.1.1.2.tar.gz
    tar xzf CLUMPP_Linux64.1.1.2.tar.gz 
    mv CLUMPP_Linux64.1.1.2.tar.gz CLUMPP_Linux64.1.1.2/.
    #used CLUMPP
    #
    #
####Used DISTRUCT v1.1 to visualize population assignments.
    #get DISTRUCT
    cd ~/bin/
    wget http://web.stanford.edu/group/rosenberglab/software/distruct1.1.tar.gz
    tar xzf distruct1.1.tar.gz 
    mv distruct1.1.tar.gz distruct1.1/.
    #used DISTRUCT
    #
    #
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
#STEPS FOR LOOKING FOR SNPs UNDER SELECTION
##Download tools
###Download BayeScan
    cd ~/bin/
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
    unzip BayeScan2.1.zip
    mv BayeScan2.1.zip BayeScan2.1
    cd BayeScan2.1/
    cd binaries/
    chmod u+x BayeScan2.1_linux64bits # makes the file executable
    # Path to BayeScan
    /home/jelber2/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
###Download Simple Fool's Guide (SFG) to RNA-seq scripts to convert vcf file to BayeScan input format
    cd ~/scripts/immunome_2014/
    mkdir fromSFG
    cd fromSFG
    wget http://sfg.stanford.edu/Scripts.zip
    unzip Scripts.zip 
    mv Scripts\ for\ SFG/ Scripts_for_SFG
####Run make_bayescan_input.py from SFG
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir bayescan
    cd bayescan
#####1.Add Genotype Qualities to ALL-samples-Q30-snps-recal-beagle.vcf.gz
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle.vcf.gz | perl -pe "s/(GT:DS:GP)/\1:GQ/" \
    > ALL-samples-Q30-snps-recal-beagle-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" \
    ALL-samples-Q30-snps-recal-beagle-fixed.vcf \
    > ALL-samples-Q30-snps-recal-beagle-fixed2.vcf
#####2.Get only polymorphic SNP loci
    #remove loci with AF=1
    grep -v "AF=1" ALL-samples-Q30-snps-recal-beagle-fixed2.vcf > ALL-samples-Q30-snps-recal-beagle-fixed3.vcf
    #filter out multiallelic_or_nonpolymorphicloci
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/remove\tlocus\t\n/" ALL-samples-Q30-snps-recal-beagle-fixed3.vcf
    done < ../popgen/multiallelic_or_nonpolymorphicloci
    grep -v "remove" ALL-samples-Q30-snps-recal-beagle-fixed3.vcf \
    > ALL-samples-Q30-snps-recal-beagle-fixed2-polymorphic.vcf
#####3.Ran make_bayescan_input.py
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_snps.txt, used_snp_genos.txt
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/
    #snps
    python ~/scripts/immunome_2014/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py \
     ../bayescan/ALL-samples-Q30-snps-recal-beagle-fixed2-polymorphic.vcf \
    ../popgen/populations.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.snps
    mv low_freq_snps.txt low_freq_snps.txt.snps
    mv population-info.txt population-info.txt.snps
    mv snpkey.txt snpkey.txt.snps
    mv used_snp_genos.txt used_snp_genos.txt.snps
    #copy files to SuperMikeII
    rsync --stats --progress --archive \
    /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/combined/bayescan/ -n
###Ran 15-bayescan_run.py on SuperMikeII
    #snps
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome_2014/combined/bayescan/bayes_input.txt.snps
    -snp \
    -d low_freq_snps.txt.snps \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.snps \
    -threads 16
###View bayescan results
    #initiate R in the terminal
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/")
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R.r")
    #snps
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_snps_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.snps_fst.txt", FDR=0.05)
    #save the candidate loci to a text file
    write(noMAF_snps_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_snps.txt", ncolumns= 1,append= FALSE)
    q()
###View bayescan results in IGV
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.snps snpkey.txt.snps2
    #code to create IGV batch file for noMAF loci
    while read i
    do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.snps2
    done < noMAF_loci_FDR_0.05_outlier_snps.txt
    grep 'goto' snpkey.txt.snps2 > noMAF_loci_FDR_0.05_outlier_snps_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome_2014/combined/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #open noMAF_loci_FDR_0.05_outlier_snps_igv.txt
###Filter annotated VCF file by outlier snps
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz > ALL-samples-snps-annotated2.vcf
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.05_outlier_snps_igv.txt > noMAF_loci_FDR_0.05_outlier_snps_vcf.txt
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_SNP\t\2\n/" ALL-samples-snps-annotated2.vcf
    done < noMAF_loci_FDR_0.05_outlier_snps_vcf.txt
    grep 'OUTLIER_SNP\|^#' ALL-samples-snps-annotated2.vcf | grep -v "contig" > ALL-samples-outlier-snps.vcf
    # get only gene names, note that some SNPs are intergenic
    grep -v "#" ALL-samples-outlier-snps.vcf | \
    perl -pe "s/.+SNPEFF_GENE_NAME=(\w+);.+\n/\1\n/" | \
    perl -pe "s/NW_.+\n/intergenic\n/g" | \
    sort | uniq -c | \
    perl -pe "s/( )+/\t/g" > ALL-samples-outlier-snps-gene-names.txt
    # how many SNPs are under selection?
    perl -ane '$sum += $F[0]; END {print $sum; print "\n"}' ALL-samples-outlier-snps-gene-names.txt
    # 66
    # how many genes have SNPs under selection
    grep -v "intergenic" ALL-samples-outlier-snps-gene-names.txt | wc -l
    # 31
    # how many SNPs are intergenic
    grep "intergenic" ALL-samples-outlier-snps-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 8
    # how many SNPs are contained in genes
    grep -v "intergenic" ALL-samples-outlier-snps-gene-names.txt | perl -ane '$sum += $F[0]; END {print $sum; print "\n"}'
    # 58
========








OLDER SCRIPTS NOT USED IN ANALYSIS
##Get FASTA read depth 20, min length 60
####Get bedtools2.22.1
    cd ~/bin/
    wget https://github.com/arq5x/bedtools2/releases/download/v2.22.1/bedtools-2.22.1.tar.gz
    tar xzf bedtools-2.22.1.tar.gz
    mv bedtools2 bedtools-2.22.1
    mv bedtools-2.22.1.tar.gz bedtools-2.22.1
    cd bedtools-2.22.1
    make
###Get intervals at least read depth 20
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir fasta-intervals
    mkdir fasta-seqs
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-intervals/
    #Use GATK Callableloci
    while read i
    do
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T CallableLoci \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -I ../call-SNPs-recal03/$i-recal03.bam \
    -L /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.interval.list \
    --summary $i.callableloci.summary \
    --minBaseQuality 20 \
    --minMappingQuality 10 \
    --minDepth 20 \
    --minDepthForLowMAPQ 10 \
    --format BED \
    --out $i.callableloci
    #Use grep to get only "CALLABLE" intervals
    #Use bedtools merge to get contiguous regions
    grep "CALLABLE" $i.callableloci | ~/bin/bedtools-2.22.1/bin/bedtools merge > $i.callableloci.cont.bed
    #Use awk to get interval lengths
    awk -v OFS='\t' '{a=$3-$2;print $1,$2,$3,a;}' $i.callableloci.cont.bed > $i.callableloci.bylength.bed
    done < ../call-SNPs-recal03/samplelist
###Get only regions at least 60bp long using R
    R
    #set the working directory
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-intervals/")
    #get the desired files
    print(files <- list.files(pattern=".bylength.bed$"))
    #convert files to samples (i.e., AL102 instead of AL102.callableloci.bylength.bed)
    print(samples <- gsub("prefixToTrash-0|\\.callableloci\\.bylength\\.bed",
          "", files, perl=TRUE), sep="")
    #for each sample:
    #1.Create the string to name the input file
    #2.Read in the file
    #3.Get only intervals of length 60 or greater
    #4.Create the string to name the output file
    #5.Write the output file
    for (i in samples){
      filein <- paste(i, ".callableloci.bylength.bed", sep="")
      i.callable <- read.table(filein)
      i.callablesixty = i.callable[i.callable$V4 > 59, ]
      fileout <- paste(i, ".callableloci.depth20.len60.bed", sep="")
      write.table(i.callablesixty,
                  file= fileout,
                  append = FALSE, quote = FALSE, sep = "\t",
                  eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                  col.names = FALSE)
    }
    q()
###Use bedtools then grep to get regions shared among all samples
    ~/bin/bedtools-2.22.1/bin/bedtools multiinter -header -i *.callableloci.depth20.len60.bed > ALLsamples
    #use grep to get intervals shared by all samples (n total = 16)
    grep -P "\t16\t" ALLsamples | cut -f 1-3 > ALLsamples.callableloci.bed
    #get only intervals that are the same length for all samples
    ~/bin/bedtools-2.22.1/bin/bedtools intersect -a ALLsamples.callableloci.bed -b *.callableloci.depth20.len60.bed -f 1.0 -r -u > ALLsamples.callableloci.samelength.bed
    #use awk to calculate lengths
    awk -v OFS='\t' '{a=$3-$2;print $1,$2,$3,a;}' ALLsamples.callableloci.samelength.bed > ALLsamples.callableloci.bylength.bed
###Use R to get intervals >= 60 bp
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-intervals/")
    callable <- read.table("ALLsamples.callableloci.bylength.bed")
    callablesixty = callable[callable$V4 > 59, ]
    write.table(callablesixty,
                file= "ALLsamples.callableloci.depth20.len60.bed",
                append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE)
    q()
###use awk to convert from 0-based to 1-based positions then use perl to convert format for samtools region
    awk -v OFS='\t' '{a=$2+1;print $1,a,$3,$4;}' ALLsamples.callableloci.depth20.len60.bed | 
    perl -pe "s/(\w+_\w+\.\d)\t(\d+)\t(\d+)\t\d+/\1:\2-\3/" > ../fasta-seqs/loci2filterbams.txt
###use GATK's FastaAlternateReferenceMaker to output consensus.fa with degenerate bases
    #first combine Q30 snps and indels after beagle genotype imputation
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/beagle/
    #concatentate snps and indels with GATK
    java -Xmx4g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    --variant ALL-samples-Q30-snps-recal-beagle.vcf \
    --variant ALL-samples-Q30-indels-recal-beagle.vcf \
    -o ALL-samples-Q30-snps-indels-recal-beagle.vcf \
    --assumeIdenticalSamples \
    -genotypeMergeOptions UNSORTED
    #updated files with rsync
    rsync --stats --progress --archive /media/immunome_2014/work/jelber2/immunome_2014/combined/beagle/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/combined/beagle/ -n
    rsync --stats --progress --archive /media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-seqs/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/combined/fasta-seqs/ -n
    #ran the following on SuperMikeII using create_fasta.sh
    #using 16 cores for at most 32 hours - actually took 6.5 hrs
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-seqs/
    ~/bin/parallel-20150122/src/parallel \
    'while read i
    do
    java -Xmx2g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -T FastaAlternateReferenceMaker \
    -o $i.{}.fa \
    -L {} \
    --variant ../beagle/ALL-samples-Q30-snps-indels-recal-beagle.vcf \
    -IUPAC $i \
    --lineWidth 10000
    perl -pi -e "s/>.+\\n/>$i.{}\\n/" $i.{}.fa
    done < ../call-SNPs-recal03/samplelist
    cat *.{}.fa > {}.fa
    rm *.{}.fa' :::: loci2filterbams.txt
    #get files from SuperMikeII
    rsync --stats --progress --archive jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/combined/fasta-seqs/ \
    /media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-seqs/ -n
###get PHASE
    cd ~/bin/
    #note had to manually download and place in ~/bin/ because had to enter a password for download
    wget http://stephenslab.uchicago.edu/phase/phasecode/phase.2.1.1.linux.tar.gz
    tar xzf phase.2.1.1.linux.tar.gz 
    mv phase.2.1.1.linux.tar.gz phase.2.1.1.linux
###get SeqPHASE to generate PHASE input file and process PHASE output
    cd ~/bin/
    mkdir seqphase
    cd seqphase/
    wget http://seqphase.mpg.de/seqphase/seqphase2014.zip
    unzip seqphase2014.zip
###get muscle
    cd ~/bin/
    mkdir muscle-3.8.31
    cd muscle-3.8.31/
    wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
    tar xzf muscle3.8.31_i86linux64.tar.gz
###run muscle,SeqPHASE,PHASE,then SeqPHASE using gnu parallel
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/fasta-seqs/
    ~/bin/parallel-20150122/src/parallel \
    '~/bin/muscle-3.8.31/muscle3.8.31_i86linux64 -in {}.fa -out {}.fa.aln
    ~/bin/seqphase/seqphase1.pl -1 {}.fa.aln -p {}
    ~/bin/phase.2.1.1.linux/PHASE {}.inp {}.out
    ~/bin/seqphase/seqphase2.pl -c {}.const -i {}.out_pairs -o {}.fa.phased' :::: loci2filterbams.txt
###Results
    #How many regions are there among the 16 samples that are at least 60bp and
    #have a read depth of 20 reads per individual?
    ls *.fa | wc -l
    #1917
    #How many of these regions are polymorphic, and thus phaseable?
    ls *.fa.phased | wc -l
    #912
    #How many genes and exons do these 912 regions represent
    ls *.fa.phased | perl -pe "s/(\w+_\d+\.\d):(\d+)-(\d+)\.fa\.phased\n/\1\t\2\t\3\n/g" |
    awk -v OFS='\t' '{a=$2-1;print $1,a,$3;}' - | sort -k 1,1 -k2,2n > ALLsamples.callableloci.depth20.len60.poly.bed
    #use bedtools to intersect the gff and bed file
    ~/bin/bedtools-2.22.1/bin/bedtools intersect \
    -a /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.introns \
    -b ALLsamples.callableloci.depth20.len60.poly.bed -wao > ALLsamples.callableloci.depth20.len60.poly.bed
    grep -v "RefSeq" test | grep -Pv ".\t-1\t-1\t0" | \
    grep -Pv "Gnomon\tmRNA" | grep -Pv "Gnomon\tCDS" | \
    grep -Pv "Gnomon\ttranscript" > test2
##Calculate dn/ds (aka - ka/ks) ratios
    #use seqinr package
    #http://www.inside-r.org/packages/cran/seqinr/docs/kaks
    read.alignment(file = file, format=FASTA)
    






##Use Analysis of Next-generation Sequencing Data (ANGSD)
###Get angsd
    cd ~/bin/
    wget http://popgen.dk/software/download/angsd/angsd0.614.tar.gz
    tar xfz angsd0.614.tar.gz
    cd angsd0.614
    make
###Get ngsF
    cd ~/bin/
    git clone https://github.com/fgvieira/ngsF.git
    make
###Get ngsPopGen
    cd ~/bin/
    git clone https://github.com/mfumagalli/ngsPopGen.git
    cd ngsPopGen
    make
###Get NgsAdmix
    cd ~/bin/
    mkdir ngsadmix
    cd ngsadmix/
    wget popgen.dk/software/NGSadmix/ngsadmix32.cpp
    g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix
###First need to estimate genotype likelihoods (GL) using ANGSD
    #make a directory for the files
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir angsd
    cd angsd
    #makes a file called pop.bamfile.list of the recalibrated bam files
    ls /media/immunome_2014/work/jelber2/immunome_2014/combined/call-SNPs-recal03/AL*.bam | grep -v "ALL" > AL.bamfile.list
    ls /media/immunome_2014/work/jelber2/immunome_2014/combined/call-SNPs-recal03/FL*.bam > FL.bamfile.list
    ls /media/immunome_2014/work/jelber2/immunome_2014/combined/call-SNPs-recal03/GG*.bam > GG.bamfile.list
    ls /media/immunome_2014/work/jelber2/immunome_2014/combined/call-SNPs-recal03/LA*.bam > LA.bamfile.list
    cat AL.bamfile.list FL.bamfile.list GG.bamfile.list LA.bamfile.list > ALL.bamfile.list
    #estimate GL for each population separately
    #options
    #-GL 2 = use GATK GL model
    #-doGlf 3 = output format for ngsF (binary 3 times likelihood)
    #-doMajorMinor 1 = infer major and minor allele from GL
    #-doMaf 2 = frequency (fixed major unknown minor, which is estimated)
    #-SNP_pval 2e-6 = likelihood ratio test for SNP, this pvalue corresponds to 29.233671 likelihood units
    #-minInd 4 = discard the sites where we don't have data from -minInd individuals
    # note
    # use nano to create the file bamfilelist, containing one population (e.g., AL, FL, GG,LA) per line
    #
###Estimate genotype likelihoods
    while read i
    do
    ~/bin/angsd0.614/angsd \
    -GL 2 \
    -out $i.genolike -nThreads 2 \
    -doGlf 3 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 2e-6 \
    -minInd 4 \
    -bam $i.bamfile.list
    done < bamfilelist
    #determine the number of n_sites with R
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/angsd")
    ALdata <-read.delim("AL.genolike.mafs.gz", header=T)
    FLdata <-read.delim("FL.genolike.mafs.gz", header=T)
    GGdata <-read.delim("GG.genolike.mafs.gz", header=T)
    LAdata <-read.delim("LA.genolike.mafs.gz", header=T)
    nrow(ALdata)
    #21882
    nrow(FLdata)
    #17809
    nrow(GGdata)
    #18510
    nrow(LAdata)
    #15850
    #unzip *genolike.glf.gz because ngsF can't read the compressed form
    gunzip *.genolike.glf.gz
###Use ngsF
    ~/bin/ngsF/ngsF -n_ind 4 -n_threads 2 -init_values r -min_epsilon 1e-9 \
    -n_sites 21882 -glf AL.genolike.glf -out AL.outputF
    ~/bin/ngsF/ngsF -n_ind 4 -n_threads 2 -init_values r -min_epsilon 1e-9 \
    -n_sites 17809 -glf FL.genolike.glf -out FL.outputF
    ~/bin/ngsF/ngsF -n_ind 4 -n_threads 2 -init_values r -min_epsilon 1e-9 \
    -n_sites 18510 -glf GG.genolike.glf -out GG.outputF
    ~/bin/ngsF/ngsF -n_ind 4 -n_threads 2 -init_values r -min_epsilon 1e-9 \
    -n_sites 15850 -glf LA.genolike.glf -out LA.outputF
###Calculate popgen statistics
####1.Estimate the folded site allele frequency (Saf) likelihood
    #options
    #-doSaf 2 = use inbreeding coefficients from ngsF (-indF outputF) so
    # we don't have to assume HWE
    #-anc = estimate ancestral allele using reference genome
    #-fold 1 = fold the Saf
    # first fix FASTA index if you need to
    #~/bin/samtools-0.1.19/samtools faidx /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna
    #now estimate folded Saf
    while read i
    do
    ~/bin/angsd0.614/angsd -bam $i.bamfile.list -out $i.outFold -nThreads 2 \
    -doSaf 2 -indF $i.outputF -doMAF 2 -doMajorMinor 1 \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -GL 2 \
    -fold 1
    done < bamfilelist
####2.Obtain the maximum likelihood estimate of the folded SFS
    # 4 is number of individuals
    #-P 2 is the number of threads/cores
    while read i
    do
    ~/bin/angsd0.614/misc/realSFS $i.outFold.saf 4 -P 2 > $i.outFold.sfs
    done < bamfilelist
####3.Calculate the thetas
    #-doThetas 1 = calculate thetas
    #-pest outFold.sfs = our estimate of the site frequency spectrum
    while read i
    do
    ~/bin/angsd0.614/angsd -bam $i.bamfile.list -out $i.outFold \
    -doThetas 1 \
    -doSaf 2 -indF $i.outputF -doMAF 2 -doMajorMinor 1 \
    -pest $i.outFold.sfs \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -GL 2 -fold 1
    done < bamfilelist
####4.Estimate Tajimas D
    #create a binary version of thete.thetas.gz
    while read i
    do
    ~/bin/angsd0.614/misc/thetaStat make_bed $i.outFold.thetas.gz
    done < bamfilelist
    #calculate Tajimas D with fixed sliding window
    while read i
    do
    ~/bin/angsd0.614/misc/thetaStat do_stat $i.outFold.thetas.gz -nChr 4 -win 50000 -step 10000
    mv $i.outFold.thetas.gz.pestPG $i.fixwin.outFold.thetas.gz.pestPG
    done < bamfilelist
    #calculate Tajimas D with "smart" sliding window
    while read i
    do
    ~/bin/angsd0.614/misc/thetaStat do_stat $i.outFold.thetas.gz -nChr 4
    mv $i.outFold.thetas.gz.pestPG $i.varwin.outFold.thetas.gz.pestPG
    done < bamfilelist
    #note only values for tw (Watterson theta) Tajima (Tajima's D) are meaningful
    #use awk to output only useful columns and then grep to remove regions without SNPs
    #fixed window
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' AL.fixwin.outFold.thetas.gz.pestPG | grep -v "0$" > AL.fixwin.thetas
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' FL.fixwin.outFold.thetas.gz.pestPG | grep -v "0$" > FL.fixwin.thetas
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' GG.fixwin.outFold.thetas.gz.pestPG | grep -v "0$" > GG.fixwin.thetas
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' LA.fixwin.outFold.thetas.gz.pestPG | grep -v "0$" > LA.fixwin.thetas
    #variable window
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' AL.varwin.outFold.thetas.gz.pestPG > AL.varwin.thetas
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' FL.varwin.outFold.thetas.gz.pestPG > FL.varwin.thetas
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' GG.varwin.outFold.thetas.gz.pestPG > GG.varwin.thetas
    awk -v OFS='\t' '{print $1,$2,$3,$4,$9,$14}' LA.varwin.outFold.thetas.gz.pestPG > LA.varwin.thetas
###Use NgsAdmix
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir ngsadmix
    cd ngsadmix/
    ~/bin/angsd0.614/angsd -GL 2 -out genolike -nThreads 2 -doGlf 2 \
    -doMajorMinor 1  -doMaf 2 -SNP_pval 2e-6 -bam ../angsd/ALL.bamfile.list \
    -minMapQ 30 -minQ 20
    #remove all loci with 'nan'
    gunzip -c genolike.beagle.gz | grep -v 'nan' | gzip > genolike.beagle2.gz 
    #run the program
    ~/bin/ngsadmix/NGSadmix -likes genolike.beagle2.gz -K 4 -P 2 -o myoutfiles \
    -minMaf 0.05 -minInd 14 -tol 2e-9 -tolLike50 2e-9
    #plot the data in R
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/ngsadmix")
    admix<-t(as.matrix(read.table("myoutfiles.qopt")))
    barplot(admix,col=1:4,space=0,border=NA,xlab="Individuals",ylab="admixture")
###Use ngsPopGen
####FST
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir ngsPopGen
    cd ngsPopGen/
    #following steps from:
    #https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
    #Steps:
    #1.Filter sites using angsd
    #ALL samples
    # -minMaf 0.05 = keep only sites with minor allele freq > 0.05, which 
    #  equal to the frequency of singletons with 10 diploids
    ~/bin/angsd0.614/angsd -b ../angsd/ALL.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 16 \
    -out test.pops.angsd -P 2 -setMinDepth 20 -setMaxDepth 100 -GL 2 -doSaf 1 -doMaf 2 \
    -minMaf 0.05 -doMajorMinor 1
    #AL population
    ~/bin/angsd0.614/angsd -b ../angsd/AL.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 4 \
    -out test.AL.angsd -P 2 -setMinDepth 10 -setMaxDepth 50 -GL 2 -doSaf 1
    #FL population
    ~/bin/angsd0.614/angsd -b ../angsd/FL.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 4 \
    -out test.FL.angsd -P 2 -setMinDepth 10 -setMaxDepth 50 -GL 2 -doSaf 1
    #GG population
    ~/bin/angsd0.614/angsd -b ../angsd/GG.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 4 \
    -out test.GG.angsd -P 2 -setMinDepth 10 -setMaxDepth 50 -GL 2 -doSaf 1
    #LA population
    ~/bin/angsd0.614/angsd -b ../angsd/LA.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -remove_bads -unique_only -minMapQ 30 -minQ 20 -only_proper_pairs 1 -trim 0 -minInd 4 \
    -out test.LA.angsd -P 2 -setMinDepth 10 -setMaxDepth 50 -GL 2 -doSaf 1
    #2.Determine identical sites among all pops, which is required to do
    #  pairwise FST comparisons (ALvsFL, ALvsGG, ALvsLA, FLvsGG, FLvsLA, GGvsLA)
    less test.AL.angsd.saf.pos.gz test.FL.angsd.saf.pos.gz \
    test.GG.angsd.saf.pos.gz test.LA.angsd.saf.pos.gz | \
    sort -S 50% | uniq -d | sort -k1,1 -S 50% > intersectALL.txt
    #index the intersection
    ~/bin/angsd0.614/angsd sites index intersectALL.txt
    #create the rf file
    cut -f1 intersectALL.txt |sort|uniq >intersectALLchrs.txt
    #3.Compute sample allele frequencies at intersection sites
    #index bam files in format *.bam.bai, so I can use the -rf option below
    #first make bamlist
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/call-SNPs-recal03/
    ls *-recal03.bam > bamlist
    #now index bam files (will take a while)
    while read i
    do
    ~/bin/samtools-1.1/samtools index $i
    done < bamlist
    cd ../ngsPopGen/
    #now compute allele frequencies
    ~/bin/angsd0.614/angsd -b ../angsd/AL.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -out test.ALvsFL -P 2 -GL 2 -doSaf 1 -doMaf 2 -doMajorMinor 1 -rf intersectALLchrs.txt -sites intersectALL.txt
    ~/bin/angsd0.614/angsd -b ../angsd/FL.bamfile.list \
    -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna \
    -out test.FLvsAL -P 2 -GL 2 -doSaf 1 -doMaf 2 -doMajorMinor 1 -rf intersectALLchrs.txt -sites intersectALL.txt
    #4.Estimate the 2D-SFS to be used as prior
    N_SITES_ALL=`wc -l intersectALL.txt | cut -f 1 -d " "`
    ~/bin/ngsPopGen/ngs2dSFS -postfiles test.ALvsFL.saf test.FLvsAL.saf -outfile test.ALvsFL.2dsfs -nind 4 4 -nsites $N_SITES_ALL
    #5.calculate FST
    ~/bin/ngsPopGen/ngsFST -postfiles test.ALvsFL.saf test.FLvsAL.saf -priorfile test.ALvsFL.2dsfs -nind 4 4 -nsites $N_SITES_ALL -outfile test.ALvsFL.fst
    #6.Create positions file
    while read i
    do
    grep "$i" intersectALL.txt | cut -f 2 > tmp$i
    done < intersectALLchrs.txt
    cat tmp* > ALLpositions.txt
    rm tmp*
    #R script doesn't seem to work?
    Rscript ~/bin/ngsPopGen/scripts/plotFST.R -i test.ALvsFL.fst -o test.ALvsFL.fst.out -p ALLpositions.txt -w 1 -s 1
####calculate nucleotide diversity, etc.
    #1.Estimate the marginal SFS to be used as priors.
    #8 is the diploid number of chr (ind=4x2=8)
    ~/bin/angsd0.614/misc/realSFS test.ALvsFL.saf 8 -P 2 > test.ALvsFL.sfs
    ~/bin/angsd0.614/misc/realSFS test.FLvsAL.saf 8 -P 2 > test.FLvsAL.sfs
    #From these priors, we calculate the sample allele frequency posterior probabilities for each population.
    ~/bin/angsd0.614/angsd -bam ../angsd/AL.bamfile.list -out test.ALvsFL -doSaf 1 -pest test.ALvsFL.sfs -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna -GL 2 -P 2 -rf intersectALLchrs.txt -sites intersectALL.txt
    ~/bin/angsd0.614/angsd -bam ../angsd/AL.bamfile.list -out test.FLvsAL -doSaf 1 -pest test.FLvsAL.sfs -anc /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna -GL 2 -P 2 -rf intersectALLchrs.txt -sites intersectALL.txt
    N_SITES_ALL=`wc -l intersectALL.txt | cut -f 1 -d " "`
    #calcualte stats
    ~/bin/ngsPopGen/ngsStat -npop 2 -postfiles test.ALvsFL.saf test.FLvsAL.saf -nsites $N_SITES_ALvsFL -nind 4 4 -outfile test.ALvsFL.stat
    #Rscript does not seem to work?
    Rscript ~/bin/ngsPopGen/scripts/plotSS.R -i test.ALvsFL.stat -p ALLpositions.txt -o test.ALvsFL.stat.pdf -n pop1-pop2 -w 1 -s 1
####perform PCA
    #We first calculate genotype posterior probabilities, assuming HWE, using ANGSD.
    ~/bin/angsd0.614/angsd -b ../angsd/ALL.bamfile.list -nInd 16 -doGeno 32 -doPost 1 -out test.ALL -P 2 -rf intersectALLchrs.txt -sites intersectALL.txt -GL 2 -doMajorMinor 1 -doMaf 2
    zcat test.ALL.geno.gz > test.ALL.geno
    #We can now estimate the covariance matrix.
    N_SITES_ALL=`wc -l intersectALL.txt | cut -f 1 -d " "`
    ~/bin/ngsPopGen/ngsCovar -probfile test.ALL.geno -outfile test.ALL.covar -nind 16 -call 0 -nsites $N_SITES_ALL
    #For plotting purposes, we create a dummy PLINK cluster file.
    Rscript -e 'write.table(cbind(seq(1,16),rep(1,16),c(rep("AL",4),rep("FL",4),rep("GG",4),rep("LA",4))), row.names=F, sep=" ", col.names=c("FID","IID","CLUSTER"), file="test.pops.clst", quote=F)'
    #This script will calculate principal components and plot them.
    Rscript ~/bin/ngsPopGen/scripts/plotPCA.R -i test.ALvsFL.covar -c 1-2 -a test.pops.clst -o test.pca.ALL.pdf
    evince test.pca.ALL.pdf
========
#STEPS FOR VARIANT PREDICTION
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
    -o ALL-samples-snps-annotated.vcf.gz \
    -O z -m none \
     ../split-vcfs/*-snps-annotated.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-snps-annotated.vcf.gz
###Get only high quality non-synonymous alleles
####There were none for the indels, so had to remove NON_SYNONYMOUS_CODING
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir venny-data
    cd venny-data
    #script to automate unzipping bgzipped files
    #and get only  with non-synonymous snps
    #NOTE COMMENT CODE FOR GENOTYPE AND UNCOMMENT
    #CODE FOR HAPLOTYPE IF YOU WANT HAPLOTYPES
    #PRECEEDING STEPS ASSUME YOU WANT GENOTYPES
    #snps
    while read i
    do
    zcat ../split-vcfs/$i-snps-annotated.vcf.gz | \
    grep -v '#' | grep 'PASS' | grep 'NON_SYNONYMOUS_CODING' | \
    cut -f 1-2,10 | \
    #code for same haplotypes (i.e., 0/1 =! 1/0)
    #perl -pe "s/(\w+\.\d)\t(\d+)\t(\d\|\d).+\n/\1\t\2\t\3\n/" > $i-snps-haplotype.txt
    #code for same genotypes (i.e., 0/1 = 1/0)
    perl -pe "s/(\w+\.\d)\t(\d+)\t(\d)\|(\d).+\n/\1\t\2\t\3\t\4\n/" | \
    awk -v OFS='\t' '{a=$3+$4;print $1,$2,a;}' - > $i-snps-genotype.txt
    done < ../call-SNPs-recal03/samplelist



####For getting snp alleles shared amongst populations
#####http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Venn-Diagrams



####1.Opened each set of text files (4 files per set) separately with gedit
    #e.g.,
    gedit AL*-snps-genotype.txt
    #opens AL102-snps-genotype.txt AL106-snps-genotype.txt AL103-snps-genotype.txt AL108-snps-genotype.txt
####2.Copy and pasted the contents of each file into Venny
    #Oliveros, J.C. (2007) VENNY. An interactive tool for comparing lists with Venn Diagrams.
    #http://bioinfogp.cnb.csic.es/tools/venny/index.html
####3.Saved the intersection of all 4 samples as pop-shared-snps-genotype.txt (e.g., AL-shared-snps-genotype.txt)
####4.Saved venn diagram image as pop-shared-snps-genotype.png
    #where 'pop' is either LA, AL, GG, or FL
####5.Saved unique parts (i.e., complement) of each Sample as Sample-unique-snps-genotype.txt
####6.Repeated steps 1-4 using the four populations
    gedit LA-shared-snps-genotype.txt AL-shared-snps-genotype.txt GG-shared-snps-genotype.txt FL-shared-snps-genotype.txt
####7.Saved venn diagram as immunome_venny_shared_all_nonsynonymous_snps.png
####8.Saved snp alleles shared by individuals in each pop as:
    LA-only-shared-snps-genotype.txt
    AL-only-shared-snps-genotype.txt
    GG-only-shared-snps-genotype.txt
    FL-only-shared-snps-genotype.txt
####9.Get only snps unique to a population
    cat LA*-unique-snps-genotype.txt > LA-unique-snps-genotype.txt
    cat AL*-unique-snps-genotype.txt > AL-unique-snps-genotype.txt
    cat GG*-unique-snps-genotype.txt > GG-unique-snps-genotype.txt
    cat FL*-unique-snps-genotype.txt > FL-unique-snps-genotype.txt
####10.Copy and pasted the contents of each file into Venny
    gedit LA-unique-snps-genotype.txt AL-unique-snps-genotype.txt GG-unique-snps-genotype.txt FL-unique-snps-genotype.txt
####11.Saved snp alleles unique to each pop as:
    LA-only-unique-snps-genotype.txt
    AL-only-unique-snps-genotype.txt
    GG-only-unique-snps-genotype.txt
    FL-only-unique-snps-genotype.txt
###Used IGV to visualize shared snps
####Get igv
    cd ~/bin
    wget http://www.broadinstitute.org/igv/projects/downloads/IGV_2.3.40.zip
    unzip IGV_2.3.40.zip 
    mv IGV_2.3.40.zip IGV_2.3.40
    cd IGV_2.3.40/
####Make igv .genome file from C_picta genome
    cd /work/jelber2/reference
    #get FASTA(FNA) file for C_picta-3.0.3
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    #get gff gene annotation for C_picta-3.0.3
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    #unzip the two files
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    #modify igv.sh and change 2000m to 4000m, so you have at least 4GB ram
    nano ~/bin/IGV_2.3.40/igv.sh
    #run igv
    ~/bin/IGV_2.3.40/igv.sh
    #Click on the 'Genomes' tab on the main menu then select, 'Create .genome File...'
    #Under 'Unique identifier' type 'C_picta-3.0.3'
    #Under 'Descriptive name' type 'C_picta-3.0.3'
    #Under 'FASTA file' select 'GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna'
    #Under 'Gene file' select 'GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff'
    #Click 'Ok'
    #Save as GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.genome in the 
    # /media/immunome_2014/work/jelber2/reference/ directory
####Make directory for igv batch files to view snps
    cd /work/jelber2/immunome
    mkdir igv-snps
    cd igv-snps
####Make target files to filter vcfs
    cut -f 1-2 ../venny-data/LA-only-unique-snps-genotype.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > LA-only_nonysn_unique_snps.txt
    cut -f 1-2 ../venny-data/LA-only-shared-snps-genotype.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > LA-only_nonysn_snps.txt
    cut -f 1-2 ../venny-data/AL-only-shared-snps-genotype.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > AL-only_nonysn_snps.txt
    cut -f 1-2 ../venny-data/GG-only-shared-snps-genotype.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > GG-only_nonysn_snps.txt
    cut -f 1-2 ../venny-data/FL-only-shared-snps-genotype.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > FL-only_nonysn_snps.txt
    #split files into 50 sites
    split -l 100 -d LA-only_nonysn_unique_snps.txt LA-only_nonysn_unique_snps.txt-
    split -l 50 -d LA-only_nonysn_snps.txt LA-only_nonysn_snps.txt-
    split -l 50 -d AL-only_nonysn_snps.txt AL-only_nonysn_snps.txt-
    split -l 50 -d GG-only_nonysn_snps.txt GG-only_nonysn_snps.txt-
    split -l 50 -d FL-only_nonysn_snps.txt FL-only_nonysn_snps.txt-
####Run IGV
    #loading the ALLsamples-annotated.vcf file
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome_2014/combined/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #Once loaded, Click on 'Tools' tab in the main menu
    #Select the batch file LA-only_nonysn_snps_target.txt-00 in the 
    #/work/jelber2/immunome_2014/combined/igv-snps/ directory
    #let IGV do it's thing (it will take a while)
    #repeat for each number (i.e., 01,02,03) and population (AL,GG,FL)
========
#STEPS FOR LOOKING FOR SNPs UNDER SELECTION
##Download tools
###Download BayeScan
    cd ~/bin/
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
    unzip BayeScan2.1.zip
    mv BayeScan2.1.zip BayeScan2.1
    cd BayeScan2.1/
    cd binaries/
    chmod u+x BayeScan2.1_linux64bits # makes the file executable
    # Path to BayeScan
    /home/jelber2/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
###Download Simple Fool's Guide (SFG) to RNA-seq scripts to convert vcf file to BayeScan input format
    cd ~/scripts/immunome_2014/
    mkdir fromSFG
    cd fromSFG
    wget http://sfg.stanford.edu/Scripts.zip
    unzip Scripts.zip 
    mv Scripts\ for\ SFG/ Scripts_for_SFG
####Run make_bayescan_input.py from SFG
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/
    mkdir bayescan
    cd bayescan
#####1.Add Genotype Qualities to ALL-samples-Q30-snps-recal-beagle.vcf.gz
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle.vcf.gz | perl -pe "s/(GT:DS:GP)/\1:GQ/" \
    > ALL-samples-Q30-snps-recal-beagle-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" \
    ALL-samples-Q30-snps-recal-beagle-fixed.vcf \
    > ALL-samples-Q30-snps-recal-beagle-fixed2.vcf
#####2.Get only polymorphic SNP loci
    #remove loci with AF=1
    grep -v "AF=1" ALL-samples-Q30-snps-recal-beagle-fixed2.vcf > ALL-samples-Q30-snps-recal-beagle-fixed3.vcf
    #filter out multiallelic_or_nonpolymorphicloci
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/remove\tlocus\t\n/" ALL-samples-Q30-snps-recal-beagle-fixed3.vcf
    done < ../popgen/multiallelic_or_nonpolymorphicloci
    grep -v "remove" ALL-samples-Q30-snps-recal-beagle-fixed3.vcf \
    > ALL-samples-Q30-snps-recal-beagle-fixed2-polymorphic.vcf
#####3.Ran make_bayescan_input.py
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_snps.txt, used_snp_genos.txt
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/
    #snps
    python ~/scripts/immunome_2014/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py \
     ../bayescan/ALL-samples-Q30-snps-recal-beagle-fixed2-polymorphic.vcf \
    ../popgen/populations.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.snps
    mv low_freq_snps.txt low_freq_snps.txt.snps
    mv population-info.txt population-info.txt.snps
    mv snpkey.txt snpkey.txt.snps
    mv used_snp_genos.txt used_snp_genos.txt.snps
    #copy files to SuperMikeII
    rsync --stats --progress --archive \
    /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ \
    jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome_2014/combined/bayescan/ -n
###Ran 15-bayescan_run.py on SuperMikeII
    #snps
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome_2014/combined/bayescan/bayes_input.txt.snps
    -snp \
    -d low_freq_snps.txt.snps \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.snps \
    -threads 16
###View bayescan results
    #initiate R in the terminal
    R
    setwd("/media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/")
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R.r")
    #snps
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_snps_results <- plot_bayescan("bayescan_no_loci_with_low_freq_minor_alleles.snps_fst.txt", FDR=0.05)
    #save the candidate loci to a text file
    write(noMAF_snps_results$outliers, file= "noMAF_loci_FDR_0.05_outlier_snps.txt", ncolumns= 1,append= FALSE)
    q()
###View bayescan results in IGV
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.snps snpkey.txt.snps2
    #code to create IGV batch file for noMAF loci
    while read i
    do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.snps2
    done < noMAF_loci_FDR_0.05_outlier_snps.txt
    grep 'goto' snpkey.txt.snps2 > noMAF_loci_FDR_0.05_outlier_snps_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome_2014/combined/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #open noMAF_loci_FDR_0.05_outlier_snps_igv.txt
###Filter annotated VCF file by outlier snps
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/
    zcat ../split-vcfs/ALL-samples-snps-annotated.vcf.gz > ALL-samples-snps-annotated2.vcf
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.05_outlier_snps_igv.txt > noMAF_loci_FDR_0.05_outlier_snps_vcf.txt
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_SNP\t\2\n/" ALL-samples-snps-annotated2.vcf
    done < noMAF_loci_FDR_0.05_outlier_snps_vcf.txt
    grep 'OUTLIER_SNP\|^#' ALL-samples-snps-annotated2.vcf > ALL-samples-outlier-snps.vcf
========



OLDER POPGEN
###Get plink-1.90b3c
    cd ~/bin/
    mkdir plink-1.90b3c
    cd plink-1.90b3c
    wget https://www.cog-genomics.org/static/bin/plink150202/plink_linux_x86_64.zip
    unzip plink_linux_x86_64.zip
####Use plink
    #what does 50 5 and 0.5 options mean?
    #a) consider a window of 50 SNPs (50)
    #b) calculate LD between each pair of SNPs in the window
    #c) remove one of a pair of SNPs if the LD is greater than 0.5 (0.5)
    #d) shift the window 5 SNPs forward and repeat the procedure. (5)
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/
    #make plink formatted population file using awk to put column 1 in place of col 2
    cp ../bayescan/populations.txt .
    awk -v OFS='\t' '{print $2,$1}' populations.txt > plink_populations_ALL
    #use grep to get only samples for each population
    grep "AL" plink_populations_ALL > plink_populations_AL
    grep "GG" plink_populations_ALL > plink_populations_GG
    grep "FL" plink_populations_ALL > plink_populations_FL
    grep "LA" plink_populations_ALL > plink_populations_LA
    #now time assess pairwise linkage disequilibrium for all loci within each pop separately
    ~/bin/plink-1.90b3c/plink --ped plink-input.ped --map plink-input.map \
    --indep-pairwise 50 5 0.5 --allow-extra-chr --keep plink_populations_AL
    mv plink.prune.in plink.AL.prune.in
    mv plink.prune.out plink.AL.prune.out
    mv plink.log plink.AL.log
    ~/bin/plink-1.90b3c/plink --ped plink-input.ped --map plink-input.map \
    --indep-pairwise 50 5 0.5 --allow-extra-chr --keep plink_populations_GG
    mv plink.prune.in plink.GG.prune.in
    mv plink.prune.out plink.GG.prune.out
    mv plink.log plink.GG.log
    ~/bin/plink-1.90b3c/plink --ped plink-input.ped --map plink-input.map \
    --indep-pairwise 50 5 0.5 --allow-extra-chr --keep plink_populations_FL
    mv plink.prune.in plink.FL.prune.in
    mv plink.prune.out plink.FL.prune.out
    mv plink.log plink.FL.log
    ~/bin/plink-1.90b3c/plink --ped plink-input.ped --map plink-input.map \
    --indep-pairwise 50 5 0.5 --allow-extra-chr --keep plink_populations_LA
    mv plink.prune.in plink.LA.prune.in
    mv plink.prune.out plink.LA.prune.out
    mv plink.log plink.LA.log
    #combine the lists of SNPs to exlcuded because not in linkage equilibrium
    sort plink.??.prune.out | uniq -d > plink.prune.out
    #now remove the loci in linkage disequilibrium
    ~/bin/plink-1.90b3c/plink --ped plink-input.ped --map plink-input.map \
    --exclude plink.prune.out --recode --allow-extra-chr
    mv plink.map plink-output.map
    mv plink.ped plink-output.ped

##Use ARLEQUIN
###run PGDSpider to create ARLEQUIN output file
    #first create template spid file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/arlequin-input.txt \
    -outputformat ARLEQUIN
    #edit the vcf2structure.spid file
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2structure.spid
    #the VCF portion will be the same as the vcf2structure.spid file's
    #need to remove portion about STRUCTURE and replace with ARLEQUIN options
    #pasted code below without leading spaces
        # Arlequin Writer questions
        WRITER_FORMAT=ARLEQUIN
        # Specify the locus/locus combination you want to write to the Arlequin file:
        ARLEQUIN_WRITER_LOCUS_COMBINATION_QUESTION=
        # Specify the DNA locus you want to write to the Arlequin file or write "CONCAT" for concatenation:
        ARLEQUIN_WRITER_ONCATENATE_QUESTION=
        # Specify which data type should be included in the Arlequin file  (Arlequin can only analyze one data type per file):
        ARLEQUIN_WRITER_DATA_TYPE_QUESTION=SNP
    #saved file as vcf2arlequin.spid
###run PGDSpider with spid file to create arlequin-input.txt
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/arlequin-input.arp \
    -outputformat ARLEQUIN \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2arlequin.spid
###run Arlequin
    #easier to use graphical user interface
##Use Genepop
###Get Genepop
    cd ~/bin/
    mkdir genepop-4.3
    cd genepop-4.3/
    wget http://kimura.univ-montp2.fr/%7Erousset/sources.tar.gz
    tar xzf sources.tar.gz
    g++ -DNO_MODULES -o Genepop GenepopS.cpp -O3
###run PGDSpider to create GENEPOP input file
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/
    #first create template spid file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/genepop-input.txt \
    -outputformat GENEPOP
    #edit the vcf2structure.spid file
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2structure.spid
    #the VCF portion will be the same as the vcf2structure.spid file's
    #need to remove portion about STRUCTURE and replace with GENEPOP options
    #pasted code below without leading spaces
        # GENEPOP Writer questions
        WRITER_FORMAT=GENEPOP
        # Specify the locus/locus combination you want to write to the GENEPOP file:
        GENEPOP_WRITER_LOCUS_COMBINATION_QUESTION=
        # Specify which data type should be included in the GENEPOP file  (GENEPOP can only analyze one data type per file):
        GENEPOP_WRITER_DATA_TYPE_QUESTION=SNP
    #saved file as vcf2genepop.spid
###run PGDSpider with spid file to create genepop-input.txt
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/genepop-input.txt \
    -outputformat GENEPOP \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2genepop.spid
    #edit genepop-input.txt to put a name on first line
    nano genepop-input.txt
###run PGDSpider to create PED plink-input.ped
    cd /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/
    #first create template spid file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/plink-input.ped \
    -outputformat PED
    #edit the vcf2structure.spid file
    nano /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2structure.spid
    #the VCF portion will be the same as the vcf2structure.spid file's
    #need to remove portion about STRUCTURE and replace with PED options
    #pasted code below without leading spaces
        # PED Writer questions
        WRITER_FORMAT=PED
        # Do you want to save an additional MAP file with loci information?
        PED_WRITER_MAP_QUESTION=true
        # Replacement character for allele encoded as 0 (0 encodes for missing data in PED):
        PED_WRITER_ZERO_CHAR_QUESTION=
        # Save MAP file
        PED_WRITER_MAP_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/plink-input.map
        # Specify the locus/locus combination you want to write to the PED file:
        PED_WRITER_LOCUS_COMBINATION_QUESTION=
    #saved file as vcf2ped.spid
###run PGDSpider with spid file to create plink-input.ped
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/bayescan/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/plink-input.ped \
    -outputformat PED \
    -spid /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/vcf2ped.spid

###run Genepop
    ~/bin/genepop-4.3/Genepop
    #select data as
    /media/immunome_2014/work/jelber2/immunome_2014/combined/popgen/genepop-input.txt
    #Run per locus diveristy option 5, then option 1
    #Run per locus diversity option 5, then option 2

