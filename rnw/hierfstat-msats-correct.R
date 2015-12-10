##R script for calculating population genetic variaton
##statistics for immunome_2014 data
##input data are 10 polymorphic microsat loci

#load the library hierarchial fstat (From developer of Fstat)
library("hierfstat")

# set the working directory to folder containing the fstat data
# note: used PGDSpider to ARlequin to fstat format

# read in the data
msat.data <-read.fstat.data(fname = "fstat-msats-input.txt")


## Allelic richness (AR)
msat.AR <- allelic.richness(data = msat.data,
                            min.n = nrow(msat.data[msat.data$Pop == 4,])*2)

msat.AR.LA.mean <-mean(msat.AR$Ar[,1],na.rm = TRUE)
msat.AR.SD.mean <-mean(msat.AR$Ar[,2],na.rm = TRUE)
msat.AR.GG.mean <-mean(msat.AR$Ar[,3],na.rm = TRUE)
msat.AR.FL.mean <-mean(msat.AR$Ar[,4],na.rm = TRUE)


msat.AR <- c(msat.AR.LA.mean,
                  msat.AR.SD.mean,
                  msat.AR.GG.mean,
                  msat.AR.FL.mean)


## Observed (Ho) and Expected (He/Hs) Heterozygosity
tort.msat.div <- basic.stats(data = msat.data)
msat.Ho <- colMeans(tort.msat.div$Ho, na.rm = TRUE)
msat.Hs <- colMeans(tort.msat.div$Hs, na.rm = TRUE)


## FST
msat.fst <- pp.fst(msat.data) #turns the 2D dataframe into a 1D vector
msat.fst <- as.vector(msat.fst$fst.pp) #turns the 2D dataframe into a 1D vector

## PCA
## enter population labels here
## 1 = LA, 2 = SD, 3 = GG, 4 = FL
msat.ids <- c(rep("LA",nrow(msat.data[msat.data$Pop == 1,])),
                   rep("SD",nrow(msat.data[msat.data$Pop == 2,])),
                   rep("GG",nrow(msat.data[msat.data$Pop == 3,])),
                   rep("FL",nrow(msat.data[msat.data$Pop == 4,])))
# make empty data frames
msat.pca <- indpca(dat = msat.data, ind.labels = msat.ids)
