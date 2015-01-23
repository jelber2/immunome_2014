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
##Install program and get refernce genome
###trimmomatic-32
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
    #PATH=~/bin/samtools-1.1/samtools
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
    /work/jelber2/immunome_2014/run1/fastq
    /work/jelber2/immunome_2014/run2/fastq (change directory in script to run2)
##Mapping
###Ran 03-bwa.py on trim.fastq.gz on trim.fastq.gz in:
    /work/jelber2/immunome_2014/run1/trimmed-data
    /work/jelber2/immunome_2014/run2/trimmed-data (change directory in script to run2)
###Ran 04-stampy.py on bwa.sam in:
    /work/jelber2/immunome_2014/run1/bwa-alignment
    /work/jelber2/immunome_2014/run2/bwa-alignment (change directory in script to run2)
##SNP Calling
###Ran 05a-clean_sort_addRG.py stampy.bam in:
    /work/jelber2/immunome_2014/run1/stampy-alignment (make RGID=%s_9Sep2014, and directory to run1)
    /work/jelber2/immunome_2014/run2/stampy-alignment (make RGID=%s_15Sep2014, and directory to run2)
###Ran 05b-clean_sort_addRG_markdup_realign.py
###Ran 06-mergeBAM_callSNPs_initial.py
###Ran 07-qual_score_recal01.py
###Ran 08-mergeBAM_callSNPs_recal01.py
###Ran 09-qual_score_recal02.py
###Ran 10-mergeBAM_callSNPs_recal02.py
###Ran 11-qual_score_recal03.py
###Ran 12-mergeBAM_callSNPs_recal03.py
###Ran 13-seq_metrics.py
###Ran 14-haplotypecaller.py
###Ran GenotypeGVCFs to perform joint genotyping
    cd /work/jelber2/immunome_2014/run1/hc/
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T GenotypeGVCFs \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
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
###Add expressions to filter variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantFiltration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -L /media/immunome_2014/work/jelber2/reference/immunome_baits_C_picta-3.0.3.list \
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
####Get only Indel variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-indels.vcf \
    -selectType INDEL
####Get only SNP variants
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -T SelectVariants \
    -V ALL-samples-Q30-snps-indels.vcf \
    -o ALL-samples-Q30-snps.vcf \
    -selectType SNP
###Get "Truthing" SNPs for Variant Quality Score Recalibration
####First get concordant SNPs between HaplotypeCaller and UnifiedGenotyper
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -T SelectVariants \
    --variant ALL-samples-Q30-snps.vcf \
    --concordance ../call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
    -o concordant-snps-HCvsUG.vcf
#####Now get only SNPs passing filters
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    --variant concordant-snps-HCvsUG.vcf \
    -o concordant-snps-HCvsUG-PASS.vcf \
    --excludeFiltered
####Second get concordant Indels between HaplotypeCaller and UnifiedGenotyper
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -T SelectVariants \
    --variant ../hc/ALL-samples-Q30-indels.vcf \
    --concordance ../call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
    -o ../hc/concordant-indels-HCvsUG.vcf
#####Now get only indels passing filters
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    --variant ../hc/concordant-indels-HCvsUG.vcf \
    -o ../hc/concordant-indels-HCvsUG-PASS.vcf \
    --excludeFiltered
#####Finally replace 'PASS' with 'INDEL'
    perl -pe "s/PASS/INDEL/g" \
    ../hc/concordant-indels-HCvsUG-PASS.vcf \
    > ../hc/concordant-indels-HCvsUG-PASS-renamed-INDEL.vcf
####Run Variant Recalibrator
    cd /media/immunome_2014/work/jelber2/immunome_2014/run1/
    mkdir vqsr
    cd vqsr
#####Recalibrate snps
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 ../hc/concordant-snps-HCvsUG-PASS.vcf \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP -an InbreedingCoeff \
    -recalFile VQSR-snps.recal \
    -mode SNP \
    -tranchesFile VQSR-snps.tranches \
    -rscriptFile VQSR-snps.plots.R \
    --maxGaussians 4
#####Recalibrate indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantRecalibrator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    -resource:concordantSet,VCF,known=true,training=true,truth=true,prior=10.0 ../hc/concordant-indels-HCvsUG-PASS.vcf \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
    -recalFile VQSR-indels.recal \
    -mode INDEL \
    -tranchesFile VQSR-indels.tranches \
    -rscriptFile VQSR-indels.plots.R \
    --maxGaussians 3
    #COULD NOT GET TO WORK with --maxGaussians 4, presumably because the indel 
    #indel dataset (~200 indels) is too small.
####Applying the recalibration on snps
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -input ../hc/ALL-samples-Q30-snps.vcf \
    --ts_filter_level 99.5 \
    -tranchesFile VQSR-snps.tranches \
    -recalFile VQSR-snps.recal \
    -o ALL-samples-Q30-snps-recal.vcf
####Applying the recalibration on indels
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T ApplyRecalibration \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -input ../hc/ALL-samples-Q30-indels.vcf \
    --ts_filter_level 99.0 \
    -tranchesFile VQSR-indels.tranches \
    -recalFile VQSR-indels.recal \
    -o ALL-samples-Q30-indels-recal.vcf
###Need to use beagle to improve SNPs (using Linkage Disequilibrium) called by Unified Genotyper
####Downloaded beagle
    cd ~/bin
    wget http://faculty.washington.edu/browning/beagle/beagle.r1398.jar
####Ran beagle on snps and indels separately
    cd /media/immunome_2014/work/jelber2/immunome_2014/run1/
    mkdir beagle
    cd beagle
    #snps
    java -Xmx8000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_2014/run1/vqsr/ALL-samples-Q30-snps-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_2014/run1/beagle/ALL-samples-Q30-snps-recal-beagle
    #indels
    java -Xmx8000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/media/immunome_2014/work/jelber2/immunome_2014/run1/vqsr/ALL-samples-Q30-indels-recal.vcf\
    nthreads=2 \
    out=/media/immunome_2014/work/jelber2/immunome_2014/run1/beagle/ALL-samples-Q30-indels-recal-beagle
========
#STEPS FOR VARIANT PREDICTION
##Download Tools First
###Downloaded snpEff
    #ideally want to know if variants will affect protein structure and possibly immune gene function
    cd /work/jelber2/reference
    wget http://iweb.dl.sourceforge.net/project/snpeff/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
###Added Chrysemys_picta_bellii-3.0.3 to snpEff.config using nano
    cd /work/jelber2/reference/snpEff
    nano snpEff.config # added the following four lines after the Capsella_rubella_v1.0 entry (remove 4 spaces on left if cut and pasting)
    # Chrysemys_picta_bellii-3.0.3
    Chrysemys_picta_bellii-3.0.3.genome : western painted turtle
    	Chrysemys_picta_bellii-3.0.3.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/
    	Chrysemys_picta_bellii-3.0.3.M.codonTable : Standard
###Created data directory for Chrysemys_picta_bellii-3.0.3 genome
    cd /work/jelber2/reference/snpEff
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
    cd /work/jelber2/reference/snpEff
    # used snpEff_build.py script to implement command below, which took < 30 minutes
    java -jar -Xmx8g /work/jelber2/reference/snpEff/snpEff.jar build -gff3 -v Chrysemys_picta_bellii-3.0.3 2>&1 | tee Chrysemys_picta_bellii-3.0.3.build
###Downloaded vcftools
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
###Downloaded bcftools
    cd /home/jelber2/bin
    git clone --branch=develop git://github.com/samtools/htslib.git
    git clone --branch=develop git://github.com/samtools/bcftools.git
    cd bcftools; make
#Need to look for protein altering variants shared by samples in the same population
###Split vcf file from GATK for snpEff
    #snpEff needs ALL-samples*.vcf file split by sample (i.e., into Sample1.vcf, Sample2.vcf)
    cd /work/jelber2/immunome/call-SNPs-recal03/
    ls *-recal03.bam | grep -Po '^\w+'| sort -u | grep -v 'ALL' > samplelist
    mkdir ../split-vcfs
    cd ../split-vcfs
    cp ../call-SNPs-recal03/samplelist .
    zcat ../beagle/ALL-samples-Q30-snps-recal-beagle.vcf.gz > ALL-samples-Q30-snps-recal-beagle.vcf
    zcat ../beagle/ALL-samples-Q30-indels-recal-beagle.vcf.gz > ALL-samples-Q30-indels-recal-beagle.vcf
    #compress snps and index with tabix
    ~/bin/samtools-1.1/htslib-1.1/bgzip ALL-samples-Q30-snps-recal-beagle.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-Q30-snps-recal-beagle.vcf.gz
    #compress indels and index with tabix
    ~/bin/samtools-1.1/htslib-1.1/bgzip ALL-samples-Q30-indels-recal-beagle.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-Q30-indels-recal-beagle.vcf.gz
    #code to split each snp vcf file
    while read i
    do
    ~/bin/bcftools/bcftools view -s $i ALL-samples-Q30-snps-recal-beagle.vcf.gz -O v -o $i-snps.vcf
    done < samplelist
    #code to split each indel vcf file
    while read i
    do
    ~/bin/bcftools/bcftools view -s $i ALL-samples-Q30-indels-recal-beagle.vcf.gz -O v -o $i-indels.vcf
    done < samplelist
###Ran snpEff on each split vcf file
    #snps
    cd /work/jelber2/immunome/split-vcfs/
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
    done < samplelist
    #indels
    while read i
    do
    java -Xmx8g -jar /media/immunome_2014/work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    $i-indels.vcf > $i-indels-snpeff.vcf
    mv snpEff_genes.txt $i-indels-snpeff-genes.txt
    mv snpEff_summary.html $i-indels-snpeff-summary.html
    done < samplelist
###Ran VariantAnnotator on each snpeff file
    #snps
    while read i
    do
    rm $i-snps.vcf.idx
    rm $i-snps-snpeff.vcf.idx
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -A SnpEff \
    --variant $i-snps.vcf \
    --snpEffFile $i-snps-snpeff.vcf \
    -L $i-snps.vcf \
    -o $i-snps-annotated.vcf
    done < samplelist
    #indels
    while read i
    do
    rm $i-indels.vcf.idx
    rm $i-indels-snpeff.vcf.idx
    java -Xmx8g -jar ~/bin/GATK-3.3.0/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /media/immunome_2014/work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -A SnpEff \
    --variant $i-indels.vcf \
    --snpEffFile $i-indels-snpeff.vcf \
    -L $i-indels.vcf \
    -o $i-indels-annotated.vcf
    done < samplelist
###Merge split, annotated vcfs
    #compress then index split snp files
    while read i
    do
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f $i-snps-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-snps-annotated.vcf.gz
    done < samplelist
    #compress then index split indel files
    while read i
    do
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f $i-indels-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-indels-annotated.vcf.gz
    done < samplelist
    #merge snp files and index them
    ~/bin/bcftools/bcftools merge -f PASS -o ALL-samples-snps-annotated.vcf.gz -O z -m both ../split-vcfs/*-snps-annotated.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-snps-annotated.vcf.gz
    #merge indel files and index them
    ~/bin/bcftools/bcftools merge -f PASS -o ALL-samples-indels-annotated.vcf.gz -O z -m both ../split-vcfs/*-indels-annotated.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-indels-annotated.vcf.gz
###Get only high quality non-synonymous alleles
####There were none for the indels, so had to remove NON_SYNONYMOUS_CODING
    cd /work/jelber2/immunome
    mkdir /venny-data
    cd venny-data
    cp ../split-vcfs/samplelist .
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
    perl -pe "s/(\w+\.\d)\t(\d+)\t(\d)\|(\d).+\n/\1\t\2\t\3\t\4\n/" | awk -v OFS='\t' '{a=$3+$4;print $1,$2,a;}' - > $i-snps-genotype.txt
    done < samplelist
    #indels 
    while read i
    do
    zcat ../split-vcfs/$i-indels-annotated.vcf.gz | \
    grep -v '#' | grep 'PASS' | \
    cut -f 1-2,10 | \
    #code for same haplotypes (i.e., 0/1 =! 1/0)
    #perl -pe "s/(\w+\.\d)\t(\d+)\t(\d\|\d).+\n/\1\t\2\t\3\n/" > $i-indels-haplotype.txt
    #code for same genotypes (i.e., 0/1 = 1/0)
    perl -pe "s/(\w+\.\d)\t(\d+)\t(\d)\|(\d).+\n/\1\t\2\t\3\t\4\n/" | awk -v OFS='\t' '{a=$3+$4;print $1,$2,a;}' - > $i-indels-genotype.txt
    done < samplelist
####For getting snps alleles shared amongst populations
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
####12.Repeated steps 1-11 for indels
###Used IGV to visualize shared snps/indels
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
    # /work/jelber2/reference/ directory
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
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #Once loaded, Click on 'Tools' tab in the main menu
    #Select the batch file LA-only_nonysn_snps_target.txt-00 in the 
    #/work/jelber2/immunome/igv-snps/ directory
    #let IGV do it's thing (it will take a while)
    #repeat for each number (i.e., 01,02,03) and population (AL,GG,FL)
========
#STEPS FOR LOOKING FOR SNPs UNDER SELECTION
##Download tools
###Download BayeScan
    cd /home/jelber2/bin
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
    unzip BayeScan2.1.zip
    mv BayeScan2.1.zip BayeScan2.1
    cd BayeScan2.1/
    cd binaries/
    chmod u+x BayeScan2.1_linux64bits # makes the file executable
    # Path to BayeScan
    /home/jelber2/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
###Download Simple Fool's Guide (SFG) to RNA-seq scripts to convert vcf file to BayeScan input format
    cd ~/scripts/immunome/
    mkdir fromSFG
    cd fromSFG
    wget http://sfg.stanford.edu/Scripts.zip
    unzip Scripts.zip 
    mv Scripts\ for\ SFG/ Scripts_for_SFG
####Run make_bayescan_input.py from SFG
    cd /work/jelber2/immunome
    mkdir bayescan-beagle
    cd bayescan-beagle
#####1.Had to make populations.txt file
    # e.g., format (\t=tab)
    #Sample\tpopulation1
#####2.Add Genotype Qualities to ALL-samples-Q30-snps-recal-beagle.vcf.gz
    cd /work/jelber2/immunome/split-vcfs/
    zcat ALL-samples-Q30-snps-recal-beagle.vcf.gz > ALL-samples-Q30-snps-recal-beagle.vcf
    perl -pe "s/(GT:DS:GP)/\1:GQ/" ALL-samples-Q30-snps-recal-beagle.vcf > ALL-samples-Q30-snps-recal-beagle-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" ALL-samples-Q30-snps-recal-beagle-fixed.vcf > ALL-samples-Q30-snps-recal-beagle-fixed2.vcf
#####3.Add Genotype Qualities to ALL-samples-Q30-indels-recal-beagle.vcf.gz
    zcat ALL-samples-Q30-indels-recal-beagle.vcf.gz > ALL-samples-Q30-indels-recal-beagle.vcf
    perl -pe "s/(GT:DS:GP)/\1:GQ/" ALL-samples-Q30-indels-recal-beagle.vcf > ALL-samples-Q30-indels-recal-beagle-fixed.vcf
    perl -pe "s/(\d\|\d:\d:\d,\d,\d)/\1:30/g" ALL-samples-Q30-indels-recal-beagle-fixed.vcf > ALL-samples-Q30-indels-recal-beagle-fixed2.vcf
#####4.Ran make_bayescan_input.py
    #30 = min genotype quality
    #4 = min number of good quality genotype required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotype were used
    #> = creates a file so you know the values for each population
    #output = bayes_input.tx, snpkey.txt, low_freq_snps.txt, used_snp_genos.txt
    cd /work/jelber2/immunome/bayescan-beagle
    #snps
    python ~/scripts/immunome/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py ../split-vcfs/ALL-samples-Q30-snps-recal-beagle-fixed2.vcf populations.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.snps
    mv low_freq_snps.txt low_freq_snps.txt.snps
    mv population-info.txt population-info.txt.snps
    mv snpkey.txt snpkey.txt.snps
    mv used_snp_genos.txt used_snp_genos.txt.snps
    #indels
    python ~/scripts/immunome/fromSFG/Scripts_for_SFG/make_bayescan_input_using_phased_data.py ../split-vcfs/ALL-samples-Q30-indels-recal-beagle-fixed2.vcf populations.txt 30 4 1 1 > population-info.txt
    mv bayes_input.txt bayes_input.txt.indels
    mv low_freq_snps.txt low_freq_snps.txt.indels
    mv population-info.txt population-info.txt.indels
    mv snpkey.txt snpkey.txt.indels
    mv used_snp_genos.txt used_snp_genos.txt.indels
    #copy files to SuperMikeII
    rsync --stats --progress --archive /work/jelber2/immunome/bayescan-beagle/ jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome/bayescan-beagle/ -n
###Ran 13-bayescan_run.py on SuperMikeII
    #snps
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome/bayescan-beagle/bayes_input.txt.snps
    -snp \
    -d low_freq_snps.txt.snps \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.snps \
    -threads 16
    #indels
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome/bayescan-beagle/bayes_input.txt.indels
    -snp \
    -d low_freq_snps.txt.indels \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles.indels \
    -threads 16
###View bayescan results
    #initiate R in the terminal
    R
    #source the plot_R.r script from Bayescan
    source("/home/jelber2/bin/BayeScan2.1/R functions/plot_R.r")
    #snps
    #plot fst values without minor alleles below minor allele frequency of 1 copy
    noMAF_snps_results <- plot_bayescan("/media/immunome_2014/work/jelber2/immunome/bayescan-beagle/bayescan_no_loci_with_low_freq_minor_alleles_snps_fst.txt", FDR=0.01)
    #save the candidate loci to a text file
    write(noMAF_snps_results$outliers, file= "/media/immunome_2014/work/jelber2/immunome/bayescan-beagle/noMAF_loci_FDR_0.01_outlier_snps.txt", ncolumns= 1,append= FALSE)
    #indels
    noMAF_indels_results <- plot_bayescan("/media/immunome_2014/work/jelber2/immunome/bayescan-beagle/bayescan_no_loci_with_low_freq_minor_alleles_indels_fst.txt", FDR=0.05)
    #save the candidate loci to a text file
    write(noMAF_indels_results$outliers, file= "/media/immunome_2014/work/jelber2/immunome/bayescan-beagle/noMAF_loci_FDR_0.05_outlier_indels.txt", ncolumns= 1,append= FALSE)
###View bayescan results in IGV
    #snps
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.snps snpkey.txt.snps2
    #code to create IGV batch file for noMAF loci
    while read i
    do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.snps2
    done < noMAF_loci_FDR_0.01_outlier_snps.txt
    grep 'goto' snpkey.txt.snps2 > noMAF_loci_FDR_0.01_outlier_snps_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-snps-annotated.vcf.gz
    #open noMAF_loci_FDR_0.01_outlier_snps_igv.txt
    #indels
    #create a copy of snpkey.txt, so it can be modified
    cp snpkey.txt.indels snpkey.txt.indels2
    #code to create IGV batch file for noMAF loci
    while read i
    do
    perl -pi -e "s/^$i\t(.+)_(.+)\n/goto \1:\2\n/" snpkey.txt.indels2
    done < noMAF_loci_FDR_0.05_outlier_indels.txt
    grep 'goto' snpkey.txt.indels2 > noMAF_loci_FDR_0.05_outlier_indels_igv.txt
    #view in IGV
    ~/bin/IGV_2.3.40/igv.sh /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-indels-annotated.vcf.gz
    #open noMAF_loci_FDR_0.05_outlier_indels_igv.txt
###Filter annotated VCF file by outlier snps and indels
    cd /media/immunome_2014/work/jelber2/immunome/bayescan-beagle
    #snps
    zcat /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-snps-annotated.vcf.gz > ALL-samples-snps-annotated.vcf
    cp ALL-samples-snps-annotated.vcf ALL-samples-snps-annotated.vcf2
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.01_outlier_snps_igv.txt > noMAF_loci_FDR_0.01_outlier_snps_vcf.txt
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_SNP\t\2\n/" ALL-samples-snps-annotated.vcf2
    done < noMAF_loci_FDR_0.01_outlier_snps_vcf.txt
    grep 'OUTLIER_SNP\|^#' ALL-samples-snps-annotated.vcf2 > ALL-samples-outlier-snps.vcf
    #indels
    zcat /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-indels-annotated.vcf.gz > ALL-samples-indels-annotated.vcf
    cp ALL-samples-indels-annotated.vcf ALL-samples-indels-annotated.vcf2
    perl -pe "s/goto (\w+\.\d):(\d+)\n/\1\t\2\n/" noMAF_loci_FDR_0.05_outlier_indels_igv.txt > noMAF_loci_FDR_0.05_outlier_indels_vcf.txt
    while read i
    do
    perl -li -e $i
    perl -pi -e "s/(^$i)\t\.\t(.+)\n/\1\tOUTLIER_INDEL\t\2\n/" ALL-samples-indels-annotated.vcf2
    done < noMAF_loci_FDR_0.05_outlier_indels_vcf.txt
    grep 'OUTLIER_INDEL\|^#' ALL-samples-indels-annotated.vcf2 > ALL-samples-outlier-indels.vcf
========
#STEPS FOR POPGEN
    cd /media/immunome_2014/work/jelber2/immunome
    mkdir popgen
    cd popgen
###Get only polymorphic SNP loci
    #vcf file for PGDSpider tool
    grep -v "AF=1" /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-Q30-snps-recal-beagle.vcf \
    > ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf
    #vcf file for vcf2smartpca script
    grep -v "AF=1" /media/immunome_2014/work/jelber2/immunome/split-vcfs/ALL-samples-Q30-snps-recal-beagle-fixed2.vcf \
    > ALL-samples-Q30-snps-recal-beagle-fixed2-polymorphic.vcf
###Ran vcf2smartpca_JPE.py (_JPE= modified script for SFG)
    cd /media/immunome_2014/work/jelber2/immunome/popgen
    cp ../bayescan-beagle/populations.txt .
    python ~/scripts/immunome/fromSFG/Scripts_for_SFG/vcf2smartpca_JPE.py \
    ALL-samples-Q30-snps-recal-beagle-fixed2-polymorphic.vcf smartpca-input 30 populations.txt 
###Get smartpca (part of EigenSoft)
    cd ~/bin
    wget http://cdn1.sph.harvard.edu/wp-content/uploads/sites/181/2014/05/EIG5.0.2.tar.gz
    tar -xzf EIG5.0.2.tar.gz
    cd EIG5.0.2
    cd src
    mk eigenstrat
    mv eigenstrat ../bin
    #Had to do the following to get smartpca to work
    #I have 'liblapack.so.3' not but 'liblapack.so.3gf'
    #Use code below to trick make symbolic link so programs will initiate '3' when looking for '3gf'
    su #to become superuser then enter password
    ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so.3gf
###Use smartpca
    cd /media/immunome_2014/work/jelber2/immunome/popgen
    #on snps
    ~/bin/EIG5.0.2/bin/smartpca -p par.smartpca-input > smartpca-input-log.txt
    perl -pi -e "s/ +/\t/g" smartpca-input_21722loci.evec #converts file from space-delimited to tab-delimited
###Use twstats to look at sig. of eigenvals
    ~/bin/EIG5.0.2/bin/twstats -t ~/bin/EIG5.0.2/POPGEN/twtable -i smartpca-input_21722loci.eval -o smartpca-input_twstats.txt
##Use Structure
###First need to convet vcf file into structure format
    #use PGDSpider
    ##get the program
    cd ~/bin/
    wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.3.zip
    unzip PGDSpider_2.0.7.3.zip
    mv PGDSpider_2.0.7.3.zip PGDSpider_2.0.7.3
    #use following command to generate spider.conf.xml and spid file in ../popgen directory
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome/popgen/structure-input.txt \
    -outputformat STRUCTURE
    #edit spider.conf.xml
    nano spider.conf.xml #to add path to samtools
    #change <entry key="PathBcftools"></entry>
    #to <entry key="PathBcftools">/home/jelber2/bin/samtools-0.1.19/bcftools/bcftools</entry>
    #change <entry key="PathSamtools"></entry>
    #to <entry key="PathSamtools">/home/jelber2/bin/samtools-0.1.19/samtools</entry>
    #save and exit
    #edit the spid file
    nano /media/immunome_2014/work/jelber2/immunome/popgen/template_VCF_STRUCTURE.spid
    #contents after editing
    #saved as vcf2structure.spid
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
        VCF_PARSER_POP_FILE_QUESTION=/media/immunome_2014/work/jelber2/immunome/popgen/populations.txt
        # Only output SNPs with a phred-scaled quality of at least:
        VCF_PARSER_QUAL_QUESTION=
        # Do you want to include non-polymorphic SNPs?
        VCF_PARSER_MONOMORPHIC_QUESTION=true
        # Output genotypes as missing if the phred-scale genotype quality is below:
        VCF_PARSER_GTQUAL_QUESTION=
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
###run PGDSpider with vcf2structure.spid file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar \
    -inputfile /media/immunome_2014/work/jelber2/immunome/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome/popgen/structure-input.txt \
    -outputformat STRUCTURE \
    -spid /media/immunome_2014/work/jelber2/immunome/popgen/vcf2structure.spid
###Installed structure from source, to run on SuperMikeII
####Get source
    cd ~/bin
    mkdir structure
    cd structure
    wget http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/structure_kernel_source.tar.gz
    tar xzf structure_kernel_source.tar.gz
    cd structure_kernel_source
####compile
    make
    #used mainparams file created using Windows front-end version
        #define OUTFILE /work/jelber2/immunome/structure-polymorphic-loci/results
        #define INFILE /work/jelber2/immunome/structure-polymorphic-loci/structure-input.txt
        #define NUMINDS 16
        #define NUMLOCI 21722
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
        #define BURNIN 10000
        #define NUMREPS 20000
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
    #ran the program using the following options
    #implemented in /home/jelber2/scripts/immunome/structure.sh
    ~/bin/structure/structure_kernel_src/structure \
    -m /work/jelber2/immunome/structure-polymorphic-loci/mainparams.test.k1 \
    -e ~/bin/structure/structure_kernel_src/extraparams
    #note for mainparams.test.k1, #define MAXPOPS 1
    #note for mainparams.test.k2, #define MAXPOPS 2
    #etc.
    #note we used default settings for extraparams(i.e., the correlated allele frequency and the admixture ancestry models)
##Use ARLEQUIN
###run PGDSpider to create ARLEQUIN output file
    #first create template spid file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.3/PGDSpider2-cli.jar
    -inputfile /media/immunome_2014/work/jelber2/immunome/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome/popgen/arlequin-input.txt \
    -outputformat ARLEQUIN
    #edit the vcf2structure.spid file
    nano /media/immunome_2014/work/jelber2/immunome/popgen/vcf2structure.spid
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
    -inputfile /media/immunome_2014/work/jelber2/immunome/popgen/ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf \
    -inputformat VCF \
    -outputfile /media/immunome_2014/work/jelber2/immunome/popgen/arlequin-input.arp \
    -outputformat ARLEQUIN \
    -spid /media/immunome_2014/work/jelber2/immunome/popgen/vcf2arlequin.spid
###run Arlequin
    #easier to use graphical user interface
##Use VCFtools to calculate Tajima's D, nucleotide diversity, et.c
    #first make edit populations.txt, so you have only torts from each population in a file, one per line
    #nucleotide diversity (pi per site)
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --site-pi --out pi.AL --keep Alabama
    #inbreeding coefficient "f"
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --het --out het.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --het --out het.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --het --out het.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --het --out het.AL --keep Alabama
    #Hardy-Weinberg Equilibrium test
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --hardy --out hwe.AL --keep Alabama
    #--genor2 calculates the squared correlation coefficient between genotypes encoded as 0, 1 and 2 to represent the number of 
    #non-reference alleles in each individual. This is the same as the LD measure reported by PLINK.
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --genor2 --out genold.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --genor2 --out genold.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --genor2 --out genold.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --genor2 --out genold.AL --keep Alabama
    #singletons
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --singletons --out singletons
    #relatedness
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --relatedness --out relatedness

    #need to figure the correct window size for tajimas d and snp density
    #Tajima's D
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.120.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 1200 --out tajimaD.1200.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 2400 --out tajimaD.2400.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.FL --keep Florida
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.GG --keep Georgia
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.LA --keep Louisiana
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --TajimaD 120 --out tajimaD.AL --keep Alabama
    #SNP density
    ~/bin/vcftools_0.1.12b/bin/vcftools --vcf ALL-samples-Q30-snps-recal-beagle-polymorphic.vcf --SNPdensity 120 --out snpden.FL --keep Florida

##Use Analysis of Next-generation Sequencing Data (ANGSD)
###install the program
    cd ~/bin
    wget http://popgen.dk/software/download/angsd/angsd0.614.tar.gz
    tar xzf angsd0.614.tar.gz
    mv angsd0.614.tar.gz angsd0.614
    cd angsd0.614
    make
###get ngsF to estimate inbreeding coefficients, so don't have to assume HWE
    ~/bin/
    git clone https://github.com/fgvieira/ngsF.git
    cd ngsF/
    make
###First need to estimate genotype likelihoods (GL)
    #make a directory for the files
    cd /media/immunome_2014/work/jelber2/immunome
    mkdir /angsd
    cd angsd
    #makes a file called bam_list of the recalibrated bam files
    ls /media/immunome_2014/work/jelber2/immunome/call-SNPs-recal03/*.bam | \
    grep -v ALL > bamfile.list
    #estimate GL
    #options
    #-GL 2 = use GATK GL model
    #-doGlf 3 = output format for ngsF (binary 3 times likelihood)
    #-doMajorMinor 1 = infer major and minor allele from GL
    #-doMaf 2 = frequency (fixed major unknown minor, which is estimated)
    #-SNP_pval 2e-6 = likelihood ratio test for SNP, this pvalue corresponds to 29.233671 likelihood units
    #-minInd 16 = discard the sites where we don't have data from -minInd individuals
    ~/bin/angsd0.614/angsd \
    -GL 2 \
    -out genolike -nThreads 2 \
    -doGlf 3 \
    -doMajorMinor 1 \
    -doMaf 2 \
    -SNP_pval 2e-6 \
    -minInd 16 \
    -bam bamfile.list
    #determine the number of n_sites with R
    #because wc -l genolike.mafs.gz was not correct?
    #1809 - 1 = 1808
    R
    setwd("/media/immunome_2014/work/jelber2/immunome/angsd/")
    foodata <-read.delim("genolike.mafs.gz", header=T)
    str(foodata)
    #'data.frame':	33131 obs. of  7 variables:
    #n_sites = 33131
    #unzip genolike.glf.gz because ngsF can't read the compressed form
    gunzip genolike.glf.gz
###Use ngsF
    ~/bin/ngsF/ngsF \
    -n_ind 16 \
    -n_threads 2 \
    -init_values r \
    -min_epsilon 1e-9 \
    -n_sites 33131 \
    -glf genolike.glf \
    -out outputF
###Calculate popgen statistics
####1.Estimate the folded site allele frequency (Saf) likelihood
    #options
    #-doSaf = assume HWE
    #-doSaf 2 = use inbreeding coefficients from ngsF (-indF outputF) so
    # we don't have to assume HWE (couldn't get to work), but shouldn't matter
    # because almost all of the inbreeding coefficients were 0 or nearly 0
    #-anc = estimate ancestral allele using reference genome
    #-fold 1 = fold the Saf
    ~/bin/angsd0.614/angsd -bam bamfile.list -out outFold -nThreads 2 \
    -doSaf 1 \
    -indF outputF \
    -anc /work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -GL 2 \
    -fold 1
####2.Obtain the maximum likelihood estimate of the SFS
    #16 is number of individuals
    #-P 2 is the number of threads/cores
    ~/bin/angsd0.614/misc/realSFS outFold.saf 16 -P 2 > outFold.sfs
####3.Calculate the thetas
    #-doThetas 1 = calculate thetas
    #-pest outFold.sfs = our estimate of the site frequency spectrum
    ~/bin/angsd0.614/angsd -bam bamfile.list -out outFold \
    -doThetas 1 \
    -doSaf 1 \
    -pest outFold.sfs \
    -anc /work/jelber2/reference/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3.fna \
    -GL 2 -fold 1
####4.Estimate Tajimas D
    #create a binary version of thete.thetas.gz 
    ~/bin/angsd0.614/misc/thetaStat make_bed outFold.thetas.gz
    #calculate Tajimas D
    ~/bin/angsd0.614/misc/thetaStat do_stat outFold.thetas.gz -nChr 16
    #note only values for tw (Watterson theta) Tajima (Tajima's D) are meaningful

