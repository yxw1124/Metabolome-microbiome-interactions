##############################################
#####merge metabolome and microbiome data#####
##############################################
rm(list = ls())
library(phyloseq)
library(tidyverse)
load("physeq_ccrc_pruned_clr.Rdata")
sample_lookup_table <- read.table("Data/MetadataCodesCRC.txt", header = TRUE,
                                  colClasses = "character")
Metabolomics <- read.table(file = "Data/CRCMetabolomics_CleanedLogTransformed.txt", sep = "\t", 
                           skip = 0, header = TRUE, comment.char = "", check.names = FALSE, row.names = 1)

otu <- otu_table(physeq_ccrc_pruned_clr)
X_otu <- t(otu)
colnames(X_otu) <- tax_table(physeq_ccrc_pruned_clr)[,6]

c <- Metabolomics$GROUP=="Control"
crc <- Metabolomics$GROUP=="CRC"
Metabolomics_c <- Metabolomics[c,]
Metabolomics_crc <- Metabolomics[crc,]
Metabolomics_ccrc <- rbind(Metabolomics_c,Metabolomics_crc)
Metabolomics_ccrc$samplecode <- row.names(Metabolomics_ccrc)
Metabolomics_ccrc_id <- right_join(Metabolomics_ccrc,sample_lookup_table, by ="samplecode")
Metabolomics_ccrc_id <- na.omit(Metabolomics_ccrc_id)
row.names(Metabolomics_ccrc_id) <- Metabolomics_ccrc_id$sampleid
Metabolomics_ccrc_id <- Metabolomics_ccrc_id %>% select(c(-sampleid,-samplecode))
X_otu <- as.data.frame(X_otu)
omics <- merge(Metabolomics_ccrc_id, X_otu, by = "row.names", all = TRUE)
omics <- omics[-3,]
row.names(omics) <- omics$Row.names
X_omics <- omics[,c(-1,-2,-3)]
X_omics <- as.matrix(X_omics)
X_omics <- apply(X_omics,2,as.numeric)
row.names(X_omics) <- row.names(omics)
Y_omics <- omics$GROUP
Y_omics <- as.factor(Y_omics)
save(X_omics,Y_omics,file = "Omics_genus.RData")

##add class
load("Omics_genus.RData")
X_omics <- as.data.frame(X_omics)
Y_omics <- as.data.frame(Y_omics)
X_omics$Class <- Y_omics
omics <- as.matrix(X_omics)
save(omics,file = "omics_genus_Class.Rdata")





