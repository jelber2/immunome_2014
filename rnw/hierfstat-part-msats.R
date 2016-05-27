##R script for calculating population genetic variaton
##statistics for immunome_2014 data
##input data are 10 polymorphic microsat loci for 16 individuals

#load the library hierarchial fstat (From developer of Fstat)
library("hierfstat")

# set the working directory to folder containing the fstat data
# note: used PGDSpider to ARlequin to fstat format
setwd("C:/Users/jelber2/Dropbox/LSU/Dissertation/Manuscripts/immunome_2014/popgen-msats/")

# read in the data
part.msat.data <-read.fstat.data(fname = "fstat-part-msats-input.txt")


## Allelic richness (AR)
part.msat.AR <- allelic.richness(data = part.msat.data,
                            min.n = nrow(part.msat.data[part.msat.data$Pop == 4,])*2)

part.msat.AR.LA.mean <-mean(part.msat.AR$Ar[,1],na.rm = TRUE)
part.msat.AR.SD.mean <-mean(part.msat.AR$Ar[,2],na.rm = TRUE)
part.msat.AR.GG.mean <-mean(part.msat.AR$Ar[,3],na.rm = TRUE)
part.msat.AR.FL.mean <-mean(part.msat.AR$Ar[,4],na.rm = TRUE)


part.msat.AR <- c(part.msat.AR.LA.mean,
                  part.msat.AR.SD.mean,
                  part.msat.AR.GG.mean,
                  part.msat.AR.FL.mean)


## Observed (Ho) and Expected (He/Hs) Heterozygosity
tort.part.msat.div <- basic.stats(data = part.msat.data)
part.msat.Ho <- colMeans(tort.part.msat.div$Ho, na.rm = TRUE)
part.msat.Hs <- colMeans(tort.part.msat.div$Hs, na.rm = TRUE)


## FST
part.msat.fst <- pp.fst(part.msat.data) #turns the 2D dataframe into a 1D vector
part.msat.fst <- as.vector(part.msat.fst$fst.pp) #turns the 2D dataframe into a 1D vector

## PCA
## enter population labels here
## 1 = LA, 2 = SD, 3 = GG, 4 = FL
part.msat.ids <- c(rep("LA",nrow(part.msat.data[part.msat.data$Pop == 1,])),
                   rep("SD",nrow(part.msat.data[part.msat.data$Pop == 2,])),
                   rep("GG",nrow(part.msat.data[part.msat.data$Pop == 3,])),
                   rep("FL",nrow(part.msat.data[part.msat.data$Pop == 4,])))
# make empty data frames
part.msat.pca <- indpca(dat = part.msat.data, ind.labels = part.msat.ids)
