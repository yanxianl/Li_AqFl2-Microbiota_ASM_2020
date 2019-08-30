# Load packages
library(here) # [CRAN] A replacement for 'file.path()', locating the files relative to the project root
library(tidyverse) # [CRAN] # Easily Install and Load the 'Tidyverse'   
library(qiime2R) # [github::jbisanz/qiime2R] # import qiime2 artifacts to R 
library(phyloseq) # [Bioconductor] Handling and analysis of high-throughput microbiome census data 

# Data import and tidy #########################################################
# metadata
metadata <- read_tsv(here("data", "metadata.tsv"), comment = "#q2") 
metadata <- metadata %>% 
  rename(SampleID = "#SampleID") %>%
  # filter negative controls and mock from feature table
  filter(!SampleType %in% c("Blank-extraction", "Blank-library", "Mock")) %>%
  mutate(SampleType = factor(SampleType, levels = c("REF-DIC", "REF-DIM", "IM-DIC", "IM-DIM")))

# feature table
table <- read_qza(here("data", "qiime2", "table-no-chlo-mito-lowPre-contam-ctrl-with-phyla.qza"))
count_tab <- table$data %>% 
  as.data.frame() 
# taxonomy
taxonomy <- read_qza(here("data", "qiime2", "taxonomy-silva132.qza"))
tax_tab <- taxonomy$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence)

# phylogenetic tree  
tree <- read_qza(here("data", "qiime2", "rooted-tree.qza"))

# alpha-diversity
obs <- read_qza(here("data", "qiime2", "core-metrics-results", "observed_otus_vector.qza"))
shn <- read_qza(here("data", "qiime2", "core-metrics-results", "shannon_vector.qza"))
pd <- read_qza(here("data", "qiime2", "core-metrics-results", "faith_pd_vector.qza"))
evn <- read_qza(here("data", "qiime2", "core-metrics-results", "evenness_vector.qza"))

# beta-diversity: distance metrics
dist_jac <- read_qza(here("data", "qiime2", "core-metrics-results", "jaccard_distance_matrix.qza"))
dist_bc <- read_qza(here("data", "qiime2", "core-metrics-results", "bray_curtis_distance_matrix.qza"))
dist_uwuf <- read_qza(here("data", "qiime2", "core-metrics-results", "unweighted_unifrac_distance_matrix.qza"))
dist_wuf <- read_qza(here("data", "qiime2", "core-metrics-results", "weighted_unifrac_distance_matrix.qza"))
dist_aitchison <- read_qza(here("data", "qiime2", "robust-Aitchison-PCA", "Aitchison-distance.qza"))

# beta-diversity: ordination
ord_jac <- read_qza(here("data", "qiime2", "core-metrics-results", "jaccard_pcoa_results.qza"))
ord_bc <- read_qza(here("data", "qiime2", "core-metrics-results", "bray_curtis_pcoa_results.qza"))
ord_uwuf <- read_qza(here("data", "qiime2", "core-metrics-results", "unweighted_unifrac_pcoa_results.qza"))
ord_wuf <- read_qza(here("data", "qiime2", "core-metrics-results", "weighted_unifrac_pcoa_results.qza"))
ord_aitchison <- read_qza(here("data", "qiime2", "robust-Aitchison-PCA", "rpca-ordination.qza"))

# build a phyloseq object
physeq <- phyloseq(otu_table(as.matrix(count_tab), taxa_are_rows = T),
                   phy_tree(tree$data), 
                   tax_table(as.matrix(tax_tab)), 
                   sample_data(metadata %>% column_to_rownames("SampleID")))

# Export data ##################################################################
save(physeq, file = here("data", "processed", "phyloseq.RData"))

save(obs, shn, pd, evn, file = here("data", "processed", "alpha-diversity.RData"))

save(dist_jac, dist_bc, dist_uwuf, dist_wuf, dist_aitchison, 
     file = here("data", "processed", "beta-diversity_distance.RData"))

save(ord_jac, ord_bc, ord_uwuf, ord_wuf, ord_aitchison, 
     file = here("data", "processed", "beta-diversity_ordination.RData"))

