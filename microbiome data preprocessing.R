#Microbiome
#Data preprocessing
rm(list = ls())
library(tidyverse)
library(phyloseq)

mb <- readxl::read_xlsx("Data/SUPPLEMENTARY_TABLE4.xlsx", sheet = 1) %>%
  rename(taxon_ID = ...1)
tax <- readxl::read_xlsx("Data/SUPPLEMENTARY_TABLE4.xlsx", sheet = 2) %>%
  rename(taxon_ID = `Feature ID`)
metadata <- readxl::read_xlsx("Data/SUPPLEMENTARY_TABLE1.xlsx") %>%
  rename(sample_code = Código)
sample_lookup_table <- read.table("Data/MetadataCodesCRC.txt", header = TRUE,
                                  colClasses = "character")

mb <- mb %>% 
  left_join(tax) %>%
  select(-taxon_ID, -Confidence) %>%
  relocate(Kingdom:Species)

tax_mat <- mb %>%
  select(Kingdom:Species) %>%
  rowid_to_column() %>%
  mutate(rowid = paste0("OTU", rowid)) %>%
  column_to_rownames("rowid") %>%
  as.matrix()

otu_mat <- mb %>%
  select(!(Kingdom:Species)) %>%
  rowid_to_column() %>%
  mutate(rowid = paste0("OTU", rowid)) %>%
  column_to_rownames("rowid") %>%
  as.matrix()

# Turn into phyloseq objects
physeq <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE), 
                   tax_table(tax_mat))

# Add sample data
sampledata <- metadata %>%
  select(sample_code, Diagnosis, Class) %>%
  left_join(sample_lookup_table, by = c("sample_code" = "samplecode")) %>%
  select(-sample_code) %>%
  filter(sampleid %in% sample_names(physeq)) %>%
  column_to_rownames("sampleid") %>%
  sample_data()
physeq <- merge_phyloseq(physeq, sampledata)
sample_variables(physeq)

# Collapse to Genus level, keep taxa for which Genus is NA
physeq <- tax_glom(physeq, "Genus", NArm = FALSE)

# Subset samples: consider only C (healthy control) and CRC
physeq_ccrc <- subset_samples(physeq, Class %in% c("C", "CRC"))

# Relative abundances
physeq_ccrc_rel <- transform_sample_counts(physeq_ccrc, function(x) x / sum(x))

# Remove features with relative abundance less than 0.0001 in at least 40% samples
physeq_ccrc_pruned <- physeq_ccrc_rel %>%
  filter_taxa(function(x) sum(x >= 0.0001) >= 0.4*length(x), prune = TRUE)

# Subset features
physeq_ccrc_pruned <- subset_taxa(physeq_ccrc, rownames(tax_table(physeq_ccrc)) %in% rownames(tax_table(physeq_ccrc_pruned)))
save(physeq_ccrc_pruned, file = "physeq_ccrc_pruned.RData")

# Now apply clr transformation 
library(microbiome)
physeq_ccrc_pruned_clr <- microbiome::transform(physeq_ccrc_pruned,'clr')
save(physeq_ccrc_pruned_clr, file = "physeq_ccrc_pruned_clr.RData")



