##microbiome composition
rm(list = ls())
library(phyloseq)
library(microbiome)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
load("physeq_ccrc_pruned.RData")

# get relative abudance-Genus level
# Averaged by group
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
pseq.genus <- aggregate_rare(physeq_ccrc_pruned, level="Genus", detection = 50, prevalence = 0.25)
pseq.genusrel <- microbiome::transform(pseq.genus, "compositional")
sample_data(pseq.genusrel)$Class <- factor(sample_data(pseq.genusrel)$Class)
p_group1 <- plot_composition(pseq.genusrel, average_by = "Class", otu.sort = "abundance")+
  ggtitle("Relative abundance")+ 
  scale_fill_brewer("Genus") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        plot.title = element_text(size = 12), 
        legend.title=element_text(size = 10),legend.text = element_text(size = 10))+
  scale_fill_manual(values = getPalette(31))+
  labs(fill = "Genus")+
  scale_y_continuous(labels=scales::percent)+
  labs(x = "Class", y = "Relative abundance")
ggsave(p_group1, filename = "~/Desktop/Relative abundance of genus level with group.png", width=6, height=7)

# get relative abudance-Phylum level
# Averaged by group
pseq.phylum <- aggregate_rare(physeq_ccrc_pruned, level="Phylum", detection = 50, 
                              prevalence = 0.25)
pseq.phylumrel <- microbiome::transform(pseq.phylum, "compositional")
sample_data(pseq.phylumrel)$Class <- factor(sample_data(pseq.phylumrel)$Class)
p_group2 <- plot_composition(pseq.phylumrel, average_by = "Class", otu.sort = "abundance")+
  ggtitle("Relative abundance")+ 
  scale_fill_brewer("Phylum") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12),
        plot.title = element_text(size = 12), 
        legend.title=element_text(size = 10),legend.text = element_text(size = 10))+
  scale_fill_brewer("Phylum", palette = "Paired")+
  labs(fill = "Phylum")+
  scale_y_continuous(labels=scales::percent)+
  labs(x = "Class", y = "Relative abundance")
ggsave(p_group2, filename = "Relative abundance of phylum level with group.png", width=6, height=7)

##Alpha-diversity
##Chao1 index
library(ggpubr)
p1 <- plot_richness(physeq_ccrc_pruned, x = "Class", color="Class", measures="Chao1")
plot1 <- p1 + geom_boxplot(data=p1$data, aes(x=Class, y=value, color=NULL), alpha=0.1,  show.legend = F)+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method = "t.test", show.legend = F)+
  geom_point(size=3.5, alpha=0.6)+
  labs(y = "Schao1")
ggsave(plot1, filename = "alpha_plot_chao1.png", width=4, height=6)

##Shannon index
p2 <- plot_richness(physeq_ccrc_pruned, x = "Class", color="Class", measures="Shannon")
plot2 <- p2 + geom_boxplot(data=p2$data, aes(x=Class, y=value, color=NULL), alpha=0.1,  show.legend = F)+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  stat_compare_means(method = "t.test", show.legend = F)+
  geom_point(size=3.5, alpha=0.6)+
  labs(y = "Hshannon")
ggsave(plot2, filename = "alpha_plot_shannon.png", width=4, height=6)


##Beta-diversity/PCoA
##Bray-Curtis distance
rm(list = ls())
library(phyloseq)
library(ggplot2)
load("physeq_ccrc_pruned.Rdata")
#bary
PCoA_1 <- ordinate(physeq_ccrc_pruned, "PCoA", "bray")
p1 <- plot_ordination(physeq_ccrc_pruned, PCoA_1, aes(Axis.1, Axis.2), color="Class")+
  stat_ellipse(level = 0.95, show.legend = F)+
  ggtitle("PCoA on Bray-Curtis distance")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),
        plot.title = element_text(size = 14), 
        legend.title=element_text(size = 14),legend.text = element_text(size = 14))+
  theme_bw()
ggsave(p1, filename = "~/Desktop/PCoA on Bray-Curtis distance.png", width=8, height=6)

##Canberra distance
PCoA_2 <- ordinate(physeq_ccrc_pruned, "PCoA", "canberra")
p2 <- plot_ordination(physeq_ccrc_pruned, PCoA_2, aes(Axis.1, Axis.2), color="Class") +
  stat_ellipse(level = 0.95, show.legend = F)+
  ggtitle("PCoA on Canberra distance")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14),
        plot.title = element_text(size = 14), 
        legend.title=element_text(size = 14),legend.text = element_text(size = 14))+ 
  theme_bw()
ggsave(p2, filename = "~/Desktop/PCoA on Canberra distance.png", width=8, height=6)


