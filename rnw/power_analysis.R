##R script for calculating population genetic variaton
##statistics for immunome_2014 data
##input data are 17901 polymorphic snp loci

#load the library hierarchial fstat (From developer of Fstat)
library("hierfstat")
load("hierfstat_msats_correct_R_workspace.RData") # created using hierfstat-msats-correct.R

#read in the data
df.snps <-read.fstat.data(fname = "fstat-input.txt")

# get only the population or 1st column of the data
pop <- df.snps[1]

val.list <- c(rep(10,10),
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
df4.snps <- ""

# go through each number of SNPs in val.list
for (i in val.list){
  df.snps <- df.snps[-1] # get rid of 1st column (populations)
  df2.snps <- df.snps[ ,sample(names(df.snps), i)] # randomly sample columns excluding column 1
  df2.pop <-cbind(pop,df2.snps) # add column 1 back as first column
  tort.div3 <-basic.stats(data = df2.pop) # calculate basic diversity (Ho,He,etc.)
  tort.ar3 <-allelic.richness(data = df2.pop, min.n = 8) # calculate allelic richness
  tort.fst3 <-pp.fst(df2.pop) # calculate FST values
  # set up the SNP diversity data 
  snp.div3 = cbind("Site"=c("SD","FL","GG","LA"),
                  "AR"=colMeans(tort.ar3$Ar),
                  "Ho"=colMeans(tort.div3$Ho),
                  "He"=colMeans(tort.div3$Hs),
                  "Marker"=rep("snps",4)) #create the data.frame snp.data
  # set up the microsat diversity data
  msat.div = cbind("Site"=c("LA","SD","GG","FL"),
                 "AR"=msat.AR,
                 "Ho"=msat.Ho,
                 "He"=msat.Hs,
                 "Marker"=rep("msats",4)) #create the data.frame msat.div
  msat.div = as.data.frame(msat.div)
  div.data3 = rbind(msat.div,snp.div3)
  div.data3$AR <- as.numeric(levels(div.data3$AR)[div.data3$AR])
  div.data3$Ho <- as.numeric(levels(div.data3$Ho)[div.data3$Ho])
  div.data3$He <- as.numeric(levels(div.data3$He)[div.data3$He])
  div.data3.sort <-div.data3[with(div.data3, order(Marker,Site)), ]
  # set up the FST data
  fst.msat <- as.data.frame(cbind("Comparison"=c("LAxSD",
                                                 "LAxGG",
                                                 "SDxGG",
                                                 "LAxFL",
                                                 "SDxFL",
                                                 "GGxFL"), #make column called "Comparison" with 6 values
                                  "FST"=c(na.omit(msat.fst)), #make column called "FST", convert FST values into vector, discard na's values, end up with 6 values
                                  "Marker"=c(rep("msats",6)))) #make column called "Marker" with 6 instances of "msats"
  fst.msat$FST <- as.numeric(levels(fst.msat$FST)[fst.msat$FST]) #convert FST values from factor to numeric
  fst.snp3 <- as.data.frame(cbind("Comparison"=c("SDxFL",
                                                "SDxGG",
                                                "GGxFL",
                                                "LAxSD",
                                                "LAxFL",
                                                "LAxGG"), #make column called "Comparison" with 6 values
                                 "FST"=c(na.omit(as.vector(tort.fst3$fst.pp))), #make column called "FST", convert FST values into vector, discard na's values, end up with 6 values
                                 "Marker"=c(rep("snp",6)))) #make column called "Marker" with 6 instances of "snp"
  fst.snp3$FST <- as.numeric(levels(fst.snp3$FST))[fst.snp3$FST] #convert FST values from factor to numeric
  fst.data3 = rbind(fst.msat, fst.snp3) #combine the data
  fst.data3.sort <- fst.data3[with(fst.data3, order(Marker,Comparison)), ] #sort the 
  # calculate p values for Pearon's r - the correlation coefficient
  fst.cor.pval3 <- cor.test(fst.data3.sort$FST[1:6], fst.data3.sort$FST[7:12], 
                           method='pearson', alternative = "greater",
                           exact = TRUE)
  AR.pval3 <- cor.test(div.data3.sort$AR[1:4], div.data3.sort$AR[5:8],
                      method='pearson', alternative = "greater", 
                      exact = TRUE,conf.level = .95)
  Ho.pval3 <- cor.test(div.data3.sort$Ho[1:4], div.data3.sort$Ho[5:8],
                      method='pearson', alternative = "greater",
                      exact = TRUE,conf.level = .95)
  He.pval3 <- cor.test(div.data3.sort$He[1:4], div.data3.sort$He[5:8],
                      method='pearson', alternative = "greater",
                      exact = TRUE,conf.level = .95)
  # combine the data in one data.frame
  df3.snps <- as.data.frame(cbind("AR"  = AR.pval3$p.value,
                                  "Ho"  = Ho.pval3$p.value,
                                  "He"  = He.pval3$p.value,
                                  "Fst" = fst.cor.pval3$p.value,
                                  "Size" = i))
  # combine previous data one iteration at a time
  df4.snps <- as.data.frame(rbind(df4.snps,df3.snps))
}

# get rid of the first empty value
df4.snps <- df4.snps[-1,]
# makes the values numeric or factors
df4.snps$AR <- as.numeric(paste(df4.snps$AR))
df4.snps$Ho <- as.numeric(paste(df4.snps$Ho))
df4.snps$He <- as.numeric(paste(df4.snps$He))
df4.snps$Fst <- as.numeric(paste(df4.snps$Fst))
df4.snps$Size <- as.factor(paste(df4.snps$Size))
# add variable number so we can reorder values of x axis later
df4.snps$Num <- c(1:110)
