##R script for calculating population genetic variaton
##statistics for immunome_2014 data
##input data are 17901 polymorphic snp loci

#load the library hierarchial fstat (From developer of Fstat)
library("hierfstat")

#set the working directory to folder containing the fstat data
#note: used PGDSpider to convert VCF of polymorphic snps to genepop and fstat formats

#read in the data
tortdata <-read.fstat.data(fname = "fstat-input.txt")

#calculate allelic richness, with a sample size of 4 (4 x 2 b/c diploid)
tort.ar <-allelic.richness(data = tortdata, min.n = 8)

#calculate basic genetic diversity stats such obs. and exp heterozygosity
tort.div <-basic.stats(data = tortdata)

#calculate pairwise FST values
tort.fst <-pp.fst(tortdata)

#note: population 1 = AL(SD), pop2=FL(FL), pop3=GA(GG), pop4=LA(FGP)
# AL is my abbreviation, while (SD,FL,GG,FGP) are Rachel Clostio's pop abbreviations

#show pairwise FST values
tort.fst$fst.pp

#show lower (ll) and upper (ul) confidence limits for pairwise FST
tortids <- c("AL102","AL103", "AL106", "AL108", "FL846", "FL855", "FL857", "FL880", "GG1044", "GG1435", "GG1835", "GG462", "LA62", "LA66", "LA77", "LA78")

#do a PCA
tortpca<-indpca(dat = tortdata, ind.labels = tortids)