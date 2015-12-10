library("genomes")
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