##R script for subsampling SNP data
##statistics for immunome_2014 data
##input data are 17901 polymorphic snp loci

#load the library hierarchial fstat (From developer of Fstat)
library("hierfstat")
#set the working directory to folder containing the fstat data
setwd("C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/immunome_2014/popgen/")
load("C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/immunome_2014/popgen-msats/hierfstat-part-msats.RData") # created using hierfstat-part-msats.R

#read in the data
part.df.snps <-read.fstat.data(fname = "fstat-input.txt")

# get only the population or 1st column of the data
part.pop <- part.df.snps[1]

part.val.list <- c(rep(10,10),
              rep(20,10),
              rep(40,10),
              rep(100,10),
              rep(200,10),
              rep(400,10),
              rep(800,10),
              rep(1600,10),
              rep(3200,10),
              rep(6400,10),
              rep(13200,10),
              rep(17901,1))

# make an empty data.frame
part.df4.snps <- ""

# go through each number of SNPs in val.list
for (i in part.val.list){
  part.df.snps <- part.df.snps[-1] # get rid of 1st column (populations)
  part.df2.snps <- part.df.snps[ ,sample(names(part.df.snps), i)] # randomly sample columns excluding column 1
  part.df2.pop <-cbind(part.pop,part.df2.snps) # add column 1 back as first column
  part.tort.div3 <-basic.stats(data = part.df2.pop) # calculate basic diversity (Ho,He,etc.)
  part.tort.ar3 <-allelic.richness(data = part.df2.pop, min.n = 8) # calculate allelic richness
  part.tort.fst3 <-pp.fst(part.df2.pop) # calculate FST values
  # set up the SNP diversity data 
  part.snp.div3 = cbind("Site"=c("SD","FL","GG","LA"),
                  "AR"=colMeans(part.tort.ar3$Ar),
                  "Ho"=colMeans(part.tort.div3$Ho),
                  "He"=colMeans(part.tort.div3$Hs),
                  "Marker"=rep("snps",4)) #create the data.frame snp.data
  # set up the microsat diversity data
  part.msat.div = cbind("Site"=c("LA","SD","GG","FL"),
                 "AR"=part.msat.AR,
                 "Ho"=part.msat.Ho,
                 "He"=part.msat.Hs,
                 "Marker"=rep("msats",4)) #create the data.frame msat.div
  part.msat.div = as.data.frame(part.msat.div)
  part.div.data3 = rbind(part.msat.div,part.snp.div3)
  part.div.data3$AR <- as.numeric(levels(part.div.data3$AR)[part.div.data3$AR])
  part.div.data3$Ho <- as.numeric(levels(part.div.data3$Ho)[part.div.data3$Ho])
  part.div.data3$He <- as.numeric(levels(part.div.data3$He)[part.div.data3$He])
  part.div.data3.sort <-part.div.data3[with(part.div.data3, order(Marker,Site)), ]
  # set up the FST data
  part.fst.msat <- as.data.frame(cbind("Comparison"=c("LAxSD",
                                                 "LAxGG",
                                                 "SDxGG",
                                                 "LAxFL",
                                                 "SDxFL",
                                                 "GGxFL"), #make column called "Comparison" with 6 values
                                  "FST"=c(na.omit(part.msat.fst)), #make column called "FST", convert FST values into vector, discard na's values, end up with 6 values
                                  "Marker"=c(rep("msats",6)))) #make column called "Marker" with 6 instances of "msats"
  part.fst.msat$FST <- as.numeric(levels(part.fst.msat$FST)[part.fst.msat$FST]) #convert FST values from factor to numeric
  part.fst.snp3 <- as.data.frame(cbind("Comparison"=c("SDxFL",
                                                "SDxGG",
                                                "GGxFL",
                                                "LAxSD",
                                                "LAxFL",
                                                "LAxGG"), #make column called "Comparison" with 6 values
                                 "FST"=c(na.omit(as.vector(part.tort.fst3$fst.pp))), #make column called "FST", convert FST values into vector, discard na's values, end up with 6 values
                                 "Marker"=c(rep("snp",6)))) #make column called "Marker" with 6 instances of "snp"
  part.fst.snp3$FST <- as.numeric(levels(part.fst.snp3$FST))[part.fst.snp3$FST] #convert FST values from factor to numeric
  part.fst.data3 = rbind(part.fst.msat, part.fst.snp3) #combine the data
  part.fst.data3.sort <- part.fst.data3[with(part.fst.data3, order(Marker,Comparison)), ] #sort the 
  # calculate p values for Pearon's r - the correlation coefficient
  part.fst.cor.pval3 <- cor.test(part.fst.data3.sort$FST[1:6], part.fst.data3.sort$FST[7:12], 
                           method='pearson', alternative = "greater",
                           exact = TRUE)
  part.AR.pval3 <- cor.test(part.div.data3.sort$AR[1:4], part.div.data3.sort$AR[5:8],
                      method='pearson', alternative = "greater", 
                      exact = TRUE,conf.level = .95)
  part.Ho.pval3 <- cor.test(part.div.data3.sort$Ho[1:4], part.div.data3.sort$Ho[5:8],
                      method='pearson', alternative = "greater",
                      exact = TRUE,conf.level = .95)
  part.He.pval3 <- cor.test(part.div.data3.sort$He[1:4], part.div.data3.sort$He[5:8],
                      method='pearson', alternative = "greater",
                      exact = TRUE,conf.level = .95)
  # combine the data in one data.frame
  part.df3.snps <- as.data.frame(cbind("AR"  = part.AR.pval3$p.value,
                                  "Ho"  = part.Ho.pval3$p.value,
                                  "He"  = part.He.pval3$p.value,
                                  "Fst" = part.fst.cor.pval3$p.value,
                                  "Size" = i))
  # combine previous data one iteration at a time
  part.df4.snps <- as.data.frame(rbind(part.df4.snps,part.df3.snps))
}

# get rid of the first empty value
part.df4.snps <- part.df4.snps[-1,]
# makes the values numeric or factors
part.df4.snps$AR <- as.numeric(paste(part.df4.snps$AR))
part.df4.snps$Ho <- as.numeric(paste(part.df4.snps$Ho))
part.df4.snps$He <- as.numeric(paste(part.df4.snps$He))
part.df4.snps$Fst <- as.numeric(paste(part.df4.snps$Fst))
part.df4.snps$Size <- as.factor(paste(part.df4.snps$Size))
# add variable number so we can reorder values of x axis later
part.df4.snps$Num <- c(1:110)
