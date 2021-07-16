########################################
####metabolome multivariate analysis####
########################################

#Delete redundant samples
rm(list = ls())
library(phyloseq)
library(tidyverse)
load("physeq_ccrc_pruned_clr.Rdata")
sample_lookup_table <- read.table("Data/MetadataCodesCRC.txt", header = TRUE,
                                  colClasses = "character")
Metabolomics <- read.table(file = "CRCMetabolomics_CleanedLogTransformed.txt", sep = "\t", 
                           skip = 0, header = TRUE, comment.char = "", check.names = FALSE, row.names = 1)
otu <- otu_table(physeq_ccrc_pruned_clr)
X_otu <- t(otu)
c <- Metabolomics$GROUP=="Control"
crc <- Metabolomics$GROUP=="CRC"
Metabolomics_c <- Metabolomics[c,]
Metabolomics_crc <- Metabolomics[crc,]
Metabolomics_ccrc <- rbind(Metabolomics_c,Metabolomics_crc)
Metabolomics_ccrc$samplecode <- row.names(Metabolomics_ccrc)
Metabolomics_ccrc_id <- right_join(Metabolomics_ccrc,sample_lookup_table, by ="samplecode")
Metabolomics_ccrc_id <- na.omit(Metabolomics_ccrc_id)
row.names(Metabolomics_ccrc_id) <- Metabolomics_ccrc_id$sampleid
X_otu <- as.data.frame(X_otu)
omics <- merge(Metabolomics_ccrc_id, X_otu, by = "row.names", all = TRUE)
omics <- omics[-3,]
row.names(omics) <- omics$Row.names
X_omics <- omics[,c(-1,-2,-3)]

X_Metabolomics <- X_omics[,1:90]
X_Metabolomics <- apply(X_Metabolomics,2,as.numeric)
row.names(X_Metabolomics ) <- row.names(omics)
X_Metabolomics <- as.matrix(X_Metabolomics)
Y_Metabolomics <- omics$GROUP
Y_Metabolomics <- as.factor(Y_Metabolomics)

X_omics <- as.matrix(X_omics)
X_omics <- apply(X_omics,2,as.numeric)
row.names(X_omics) <- row.names(omics)
Y_omics <- omics$GROUP
Y_omics <- as.factor(Y_omics)

save(X_omics,Y_omics, file = "0mics_otu.RData")
save(X_Metabolomics,Y_Metabolomics,file = "Metabolomics.RData")

##PCA & PLS-DA
load("Metabolomics.RData")
library(do)
library(mixOmics)
library(ggplot2)
Y_Metabolomics <- Replace(data=Y_Metabolomics, from="Control",to="C")

MyResult.plsda <- plsda(X_Metabolomics, Y_Metabolomics, ncomp=2)
plotIndiv(MyResult.plsda, ind.names = F, ellipse = TRUE, legend = TRUE,  
          title = 'PLS-DA on metabolome', legend.title = "Class",
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2', star = TRUE)

MyResult.pca <- pca(X_Metabolomics,ncomp=2)

plotIndiv(MyResult.pca, group = Y_Metabolomics, legend = TRUE, 
          ind.names = F, ellipse = TRUE, title = 'PCA on metabolome',
          X.label = 'Component 1', Y.label = 'Component 2', legend.title = "Class")

