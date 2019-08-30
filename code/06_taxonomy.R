# Load packages
library(here) # [CRAN] A replacement for 'file.path()', locating the files relative to the project root
library(tidyverse) # [CRAN] Easily install and load the 'Tidyverse'   
library(cowplot)  # [CRAN] Streamlined plot theme and plot annotations for 'ggplot2' 
library(qiime2R) # [github::jbisanz/qiime2R] Import qiime2 artifacts to R 
library(MicrobeR) # [github::jbisanz/MicrobeR] Handy functions for downstream visualization fo microbiome data
library(microbiome) # [Bioconductor] Microbiome analytics
library(phyloseq) # [Bioconductor] Handling and analysis of high-throughput microbiome census data 
library(DT) # [CRAN] An R interface to the DataTables library
library(venn) # [CRAN] Make Venn's diagram

# Load functions
source(here("code", "functions", "make_taxa_barplot.R"))

# Load data 
load(here("data", "processed", "phyloseq.RData"))

# Taxnomic compostion of mock and negative controls ############################
# Taxnomic compostion of mock --------------------------------------------------
# Import and tidy observed mock composition 
mock_obs <- read_qza(here("data", "qiime2", "mock-observed-no-contam-l7-rel.qza"))

mock_obs <- mock_obs$data %>%
  as.data.frame() %>%
  rownames_to_column("tax_obs") %>%
  rename(MOCK1 = "AqFl2-M1", MOCK2 = "AqFl2-M2") %>%
  # prune taxonomy to the lowest level
  mutate(tax_obs = gsub("D_0.*D_6__", "", tax_obs),  
         tax_obs = gsub("D_0.*D_5__|;__", "", tax_obs))

# Import and tidy expected mock composition 
mock_exp <- read_tsv(here("data", "reference", "mock_expected.tsv"))   

mock_exp <- select(mock_exp, -"AqFl2-M2") %>%
  rename(tax_exp = "#OTU ID", MOCK = "AqFl2-M1") %>%
  arrange(tax_exp) %>%
  mutate(tax_exp = gsub("D_0.*D_6__", "", tax_exp)) # prune taxonomy to the lowest level

# Merge tables  
mock_merged <- bind_cols(mock_exp, mock_obs) %>%
  mutate(Gram_staining = c("(G+)", "(G+)", "(G+)", "(G+)", "(G+)", "(G-)", "(G-)", "(G-)"),
         "tax" = paste(tax_exp, Gram_staining, "/", tax_obs)) %>%
  select(tax, MOCK, MOCK1, MOCK2) %>%
  mutate(tax = factor(.$tax, levels = unique(.$tax))) %>% # use the exact order of taxa for plotting
  gather(key = "type", value = "abundance", -tax) %>%
  mutate_if(is.numeric, funs(100*.)) # use percentage as unit for relative abundance

# Make stacked barplot
p_mock <- ggplot(mock_merged, aes(type, abundance)) + 
  geom_bar(aes(fill = tax), stat = "identity") + 
  annotate("segment", x = 1.6, xend = 3.4, y = 102, yend = 102) +
  annotate("text", x = c(1, 2.5), y = c(105, 105), 
           label = c("Expected", "Observed"), size = 5) +
  scale_y_continuous(limits = c(0, 105), breaks = 0:10*10, expand = expand_scale(mult = c(0, 0.05))) +
  labs(x = "Sample ID", y = "Relative abundance (%)", fill = "Taxonomy (expected/observed)") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8, "Paired"),
                  # use italic font for genus/species names
                    labels = c(expression(paste(italic("Bacillus subtilis"), "(G+) / ", italic("Bacillus"))), 
                               expression(paste(italic("Listeria monocytogenes"), "(G+) / ", italic("Listeria monocytogenes"))),
                               expression(paste(italic("Staphylococcus aureus"), "(G+) / ", italic("Staphylococcus aureus"))),
                               expression(paste(italic("Enterococcus faecalis"), "(G+) / ", italic("Enterococcus faecalis"))),
                               expression(paste(italic("Lactobacillus fermentum"), "(G+) / ", italic("Lactobacillus fermentum"))),
                               expression(paste(italic("Escherichia coli"), "(G-) / ", italic("Escherichia-Shigella"))),
                               expression(paste(italic("Salmonella enterica"), "(G-) / ", italic("Salmonella"))),
                               expression(paste(italic("Pseudomonas aeruginosa"), "(G-) / ", italic("Pseudomonas")))
                               )
                 ) +
  theme_cowplot() +
  theme(panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "black"),
        legend.justification = "top", # move legend to top right position (default is on the right)
        legend.text.align = 0) # align legend text to left


# Taxnomic compostion of negative controls -------------------------------------
# Import and tidy table
tab <- read_tsv(here("data", "qiime2", "decontam", "table-no-chlo-mito-lowPre-with-phyla-tax-rel.tsv"), 
                skip = 1)

tab_neg <- tab %>% 
  # sum up features found in the negative controls
  mutate(ctrl = rowSums(.[grep("EB|LB", names(.))], na.rm = TRUE)) %>% 
  # remove features not present in negative controls
  filter(ctrl > 0) %>% 
  select(taxonomy, matches('EB|LB'), ctrl) %>% 
  # the following 15 lines of codes produce the lowest level of taxonomy for each feature
  separate(taxonomy, sep = ";", 
           c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  mutate(Class = ifelse(is.na(Class)|grepl("uncultured|Ambiguous", Class),
                        Phylum, 
                        Class)) %>% 
  mutate(Order = ifelse(is.na(Order)|grepl("uncultured|Ambiguous", Order),
                        Class,
                        Order)) %>%
  mutate(Family = ifelse(is.na(Family)|grepl("uncultured|Ambiguous", Family), 
                         Order, 
                         Family)) %>%
  mutate(Genus = ifelse(is.na(Genus)|grepl("uncultured|Ambiguous", Genus),
                        Family,
                        Genus)) %>%
  mutate(Species = ifelse(is.na(Species)|grepl("uncultured|Ambiguous", Species),
                          Genus, Species)) %>%
  select(-(Kingdom:Genus)) %>%
  # add Greengenes style prefix to taxonomy if species level annotation is not available
  mutate(Species = gsub("D_1", "p", Species), Species = gsub("D_2", "c", Species), 
         Species = gsub("D_3", "o", Species), Species = gsub("D_4", "f", Species), 
         Species = gsub("D_5", "g", Species), Species = gsub("D_6__", "", Species)) %>%
  rename(taxonomy = Species) %>%
  # the following 2 lines of codes collapse features by taxonomy
  group_by(taxonomy) %>%
  summarise_if(is.numeric, funs(sum)) %>%
  # use percentage as unit for relative abundance
  mutate_if(is.numeric, funs(100*.))  %>% 
  # sort features by their relative abundance
  arrange(desc(ctrl)) %>%
  select(-ctrl) %>%
  # use the sorted order of taxa (high to low abundance) for plotting
  mutate(taxonomy = factor(.$taxonomy, levels = unique(.$taxonomy))) 

# Display taxa in negative controls
datatable(tab_neg)
  
# Make stacked barplot
p_neg <- slice(tab_neg, 1:10) %>% # select the top 10 most abundance taxa for plotting
  gather(key = "SampleID", value = "abundance", -taxonomy) %>%
  mutate(SampleID = gsub("AqFl2-EB", "EB", SampleID),
         SampleID = gsub("AqFl2-LB", "LB", SampleID)) %>%
  ggplot(aes(SampleID, abundance)) + 
    geom_bar(aes(fill = taxonomy), stat = "identity") + 
    annotate("segment", x = c(0.6, 4.6), xend = c(4.4, 6.4), 
             y = c(102, 102), yend = c(102,102)) +
    annotate("text", x = c(2.5, 5.5), y = c(105, 105), 
             label = c("Extraction blank", "Library blank"), size = 5) +
    scale_y_continuous(limits = c(0, 105), breaks = 0:10*10, 
                       expand = expand_scale(mult = c(0, 0.05))) +
    labs(x = "Sample ID", y = "Relative abundance (%)", fill = "Taxonomy") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(10, "Paired"),
                      guide = guide_legend(label.theme = element_text(face = "italic"))) +
    theme_cowplot() +
    theme(panel.grid.major.y = element_line(size = 0.5, linetype = 'dashed', colour = "black"),
          legend.justification = "top") # move legend to top right position (default is on the right)

# Assemble plots ---------------------------------------------------------------
plot_grid(p_mock, p_neg, labels = "AUTO", ncol = 2, rel_widths = c(1, 1)) 

ggsave(here("result", "figures", "Figure 1.tiff"), width = 16, height = 8, 
       units = "in", dpi = 300, compression = "lzw")

# Taxnomic compostion of samples ############################################### 
# Extract feature table, taxonomy and metadata from the phyloseq object
count_tab <- as.data.frame(otu_table(physeq)) 
tax_tab <- as.data.frame(tax_table(physeq)) 
metadata <- data.frame(sample_data(physeq), check.names = FALSE)

# Summarize taxonomy at phylum level -------------------------------------------
# Collapse feature table at phylum level
taxa_p <- Summarize.Taxa(count_tab, tax_tab)$Phylum %>%
  # the following 3 lines of codes prune the taxonomy to contain phylum names only
  rownames_to_column("tax") %>%
  mutate(tax = gsub("k__.*p__", "", tax)) %>%
  column_to_rownames("tax") 

# Make taxa barplot
p_phy <- make_taxa_barplot(count_table = taxa_p, metadata = metadata, ntaxa = 10,
                           group = SampleType,  plot_mean = F) +
  labs(fill = "Taxonomy")

# Summarize taxonomy at genus level --------------------------------------------
# Collapse feature table at genus level
taxa_g <- Summarize.Taxa(count_tab, tax_tab)$Genus %>%
  rownames_to_column("tax") %>%
  separate(tax, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %>% 
  mutate(Class = ifelse(Class == "NA"|grepl("uncultured|Ambiguous", Class),
                        Phylum, 
                        Class)) %>% 
  mutate(Order = ifelse(Order == "NA"|grepl("uncultured|Ambiguous", Order),
                        Class,
                        Order)) %>%
  mutate(Family = ifelse(Family == "NA"|grepl("uncultured|Ambiguous", Family), 
                         Order, 
                         Family)) %>%
  mutate(Genus = ifelse(Genus == "NA"|grepl("uncultured|Ambiguous", Genus),
                        Family,
                        Genus)) %>%
  select(-(Kingdom:Family)) %>%
  mutate(Genus = gsub("g__", "", Genus)) %>%
  # the following 4 lines of codes collapse features by taxonomy
  gather("SampleID", "counts", -Genus) %>%
  group_by(SampleID, Genus) %>%
  summarize(counts_colsum = sum(counts)) %>%
  spread("SampleID", "counts_colsum") %>%
  column_to_rownames("Genus") 

# Make taxa barplot showing top 10 genera
p_gen <- make_taxa_barplot(count_table = taxa_g, metadata = metadata, ntaxa = 10, 
                           group = SampleType,  plot_mean = F) +
  labs(fill = "Taxonomy") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Paired"),
                    # use italic font for genus names
                    labels = c("Others", 
                               expression(italic("Lactobacillus")),
                               expression(italic("Actinomyces")),
                               expression(italic("Bacillus")),
                               expression(italic("Photobacterium")),
                               expression(italic("Corynebacterium 1")),
                               "o__Lactobacillales",
                               "f__Spirochaetaceae",
                               expression(italic("Mycoplasma")),
                               expression(italic("Aliivibrio")),
                               expression(italic("Brevinema")))) +
  theme(legend.text.align = 0) 
  
# Assemble taxa barplot --------------------------------------------------------
plot_grid(p_phy, p_gen, labels = "AUTO", ncol = 1, align = 'v', axis = "lr") 

ggsave(here("result", "figures", "Figure 2.tiff"), width = 12, height = 8, 
       units = "in", dpi = 300, compression = "lzw")

# Core microbiome ##############################################################
# Calculate feature prevalence
prev_rDic <- subset_samples(physeq, SampleType == "REF-DIC") %>% prevalence() 
prev_rDim <- subset_samples(physeq, SampleType == "REF-DIM") %>% prevalence() 
prev_iDic <- subset_samples(physeq, SampleType == "IM-DIC") %>% prevalence() 
prev_iDim <- subset_samples(physeq, SampleType == "IM-DIM") %>% prevalence() 

# Gather core features
core_taxa <- cbind.data.frame(prev_rDic, prev_rDim, prev_iDic, prev_iDim) %>%
  rownames_to_column("featureID") %>%
  # get core features based on 80% prevalence threshold
  filter(prev_rDic >= 0.8|prev_rDim >= 0.8|prev_iDic >= 0.8|prev_iDim >= 0.8)

# Core feature taxonomy
rownames_to_column(tax_tab, "featureID") %>%
  inner_join(core_taxa, by = "featureID") %>%
  rename("REF-DIC" = prev_rDic, "REF-DIM" = prev_rDim,
         "IM-DIC" = prev_iDic, "IM-DIM" = prev_iDim) %>%
  mutate(prev_all = rowSums(.[9:12])) %>%
  arrange(desc(prev_all)) %>%
  select(-prev_all) %>%
  datatable()

# Convert feature prevalence to boolean values for plotting
core_taxa_venn <- core_taxa %>% 
  rename("REF-DIC" = prev_rDic, "REF-DIM" = prev_rDim,
         "IM-DIC" = prev_iDic, "IM-DIM" = prev_iDim) %>%
  mutate_if(is.numeric, ~if_else(.x >= 0.8, 1, 0))

# Make Venn diagram
tiff(here("result", "figures", "Figure 3.tiff"), compression = "lzw", 
     units = "in", res = 300, height = 6, width = 6)

venn(core_taxa_venn[2:5], ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 1)

dev.off()

