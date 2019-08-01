library(tidyverse) # [CRAN] # Easily Install and Load the 'Tidyverse'   
library(data.table) # [CRAN] # Extension of `data.frame` 
library(qiime2R) # [github::jbisanz/qiime2R] # qiime2R 
library(MicrobeR) # [not installed on this machine] # not installed on this machine
library(cowplot) # [github::wilkelab/cowplot] # Streamlined Plot Theme and Plot Annotations for 'ggplot2' 
library(phyloseq) # Handling and analysis of high-throughput microbiome census data 
library(microbiome) # Microbiome Analytics 

setwd("C:/Users/yanxianl/OD/AquaFly/AqFl2/Manuscript2/analyses")

# Data import and tidy #########################################################
metadata <- read_tsv("data/metadata.tsv", comment = "#q2") 
metadata <- metadata %>% 
  rename(SampleID = "#SampleID") %>%
  filter(!SampleType %in% c("Blank-extraction", "Blank-library", "Mock")) %>%
  mutate(SampleType = factor(SampleType, levels = c("REF-DIC", "REF-DIM", "IM-DIC", "IM-DIM"))) %>%
  # define color scheme for making the rarefaction plot: brewer.pal(4, "Dark2")
  mutate(color = gsub("REF-DIC", "#1B9E77", SampleType), 
         color = gsub("REF-DIM", "#D95F02", color),
         color = gsub("IM-DIC", "#7570B3", color),
         color = gsub("IM-DIM", "#E7298A", color),
         color = factor(color, levels = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")))

table <- read_qza("data/processed/table-no-chlo-mito-sin-lowPre-contam-with-phyla.qza")
count_tab <- table$data %>% 
  as.data.frame() %>%
  select(-(73:80))

taxonomy <- read_qza("data/processed/taxonomy-silva132.qza")
tax_tab <- taxonomy$data %>% 
  as.data.frame() %>%
  mutate(Taxon = gsub("D_0", "k", Taxon), Taxon = gsub("D_1", "p", Taxon),
         Taxon = gsub("D_2", "c", Taxon), Taxon = gsub("D_3", "o", Taxon),
         Taxon = gsub("D_4", "f", Taxon), Taxon = gsub("D_5", "g", Taxon),
         Taxon = gsub("D_6", "s", Taxon)) %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  select(-Confidence)
  
tree <- read_qza("data/processed/rooted-tree.qza")

physeq <- phyloseq(otu_table(as.matrix(count_tab), taxa_are_rows = T),
                   phy_tree(tree$data), 
                   tax_table(as.matrix(tax_tab)), 
                   sample_data(metadata %>% column_to_rownames("SampleID")))

obs <- read_qza("core-metrics-results-merged/observed_otus_vector.qza")
shn <- read_qza("core-metrics-results-merged/shannon_vector.qza")
pd <- read_qza("core-metrics-results-merged/faith_pd_vector.qza")
evn <- read_qza("core-metrics-results-merged/evenness_vector.qza")

dist_aitchison <- read_qza("deicode/distance.qza")
dist_uwuf <- read_qza("core-metrics-results-merged/unweighted_unifrac_distance_matrix.qza")
dist_wuf <- read_qza("core-metrics-results-merged/weighted_unifrac_distance_matrix.qza")
dist_bc <- read_qza("core-metrics-results-merged/bray_curtis_distance_matrix.qza")

ord_aitchison <- read_qza("deicode/ordination.qza")
ord_uwuf <- read_qza("core-metrics-results-merged/unweighted_unifrac_pcoa_results.qza")
ord_wuf <- read_qza("core-metrics-results-merged/weighted_unifrac_pcoa_results.qza")
ord_bc <- read_qza("core-metrics-results-merged/bray_curtis_pcoa_results.qza")

# Plot taxonomy ################################################################
# Use a modified function to make taxa barplot
source("down_stream_analysis/R_functions/make_barplot.R")

# Summarize taxonomy at phylum level
taxa_p <- Summarize.Taxa(count_tab, tax_tab)$Phylum %>%
  rownames_to_column("Taxon") %>%
  mutate(Taxon = gsub("k__.*p__", "", Taxon)) %>%
  column_to_rownames("Taxon")

# Plot top 10 phyla
make_barplot(table = taxa_p, metadata = metadata, grouping_var = SampleType,
             taxa_to_plot = 10, collapse = F)

ggsave('down_stream_analysis/taxonomy_phylum_top10.tiff', compression = "lzw",
       units = "in", dpi = 300)

# Summarize taxonomy at genus level
taxa_g <- Summarize.Taxa(count_tab, tax_tab)$Genus %>%
  rownames_to_column("Taxon") %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) %>% 
  mutate(Class = ifelse(Class == "NA" | Class == "Ambiguous_taxa" | grepl("uncultured", Class),
                        Phylum, 
                        Class)) %>% 
  mutate(Order = ifelse(Order == "NA" | Order == "Ambiguous_taxa" | grepl("uncultured", Order),
                        Class,
                        Order)) %>%
  mutate(Family = ifelse(Family == "NA" | Family == "Ambiguous_taxa" | grepl("uncultured", Family), 
                         Order, 
                         Family)) %>%
  mutate(Genus = ifelse(Genus == "NA" | Genus == "Ambiguous_taxa" | grepl("uncultured", Genus),
                        Family,
                        Genus)) %>%
  select(-(Kingdom:Family)) %>%
  mutate(Genus = gsub("g__", "", Genus)) %>%
  gather("SampleID", "counts", -Genus) %>%
  group_by(SampleID, Genus) %>%
  summarize(counts_colsum = sum(counts)) %>%
  spread("SampleID", "counts_colsum") %>%
  column_to_rownames("Genus")

# Plot top 10 genera
make_barplot(table = taxa_g, metadata = metadata, grouping_var = SampleType,
             taxa_to_plot = 10, collapse = F)

ggsave('down_stream_analysis/taxonomy_genus_top10.tiff', compression = "lzw",
       units = "in", dpi = 300)

# Core features ##############################################################
# Core features, Venn diagram ------------------------------------------------
# Calculate feature prevalence
prev_pch <- subset_samples(physeq, SampleType == "REF-DIC") %>% prevalence() 
prev_tlpLake <- subset_samples(physeq, SampleType == "REF-DIM") %>% prevalence() 
prev_tlpCage <- subset_samples(physeq, SampleType == "IM-DIC") %>% prevalence() 
prev_tlpPond <- subset_samples(physeq, SampleType == "IM-DIM") %>% prevalence() 

# Gather core features
core_taxa <- cbind.data.frame(prev_pch, prev_tlpLake, prev_tlpCage, prev_tlpPond) %>%
  rownames_to_column("featureID") %>%
  # get core features based on 80% prevalence threshold
  filter(prev_pch >= 0.8 | prev_tlpLake >= 0.8 | prev_tlpCage >= 0.8 | prev_tlpPond >= 0.8)

# Core feature taxonomy
rownames_to_column(tax_tab, "featureID") %>%
  inner_join(core_taxa, by = "featureID") %>%
  rename(prevalence_REF-DIC = prev_pch, prevalence_tilapia_lake = prev_tlpLake,
         prevalence_tilapia_cage = prev_tlpCage, prevalence_tilapia_pond = prev_tlpPond) %>%
  mutate(prev_all = rowSums(.[9:12])) %>%
  arrange(desc(prev_all)) %>%
  select(-prev_all) %>%
  DT::datatable()
    
# Convert feature prevalence to boolean values for plotting
core_taxa_venn <- core_taxa %>% 
  rename(REF-DIC = prev_pch, tilapia_lake = prev_tlpLake, 
         tilapia_cage = prev_tlpCage, tilapia_pond = prev_tlpPond) %>%
  mutate_if(is.numeric, ~if_else(.x >= 0.8, 1, 0))

# Make Venn diagram
library(venn)

tiff('down_stream_analysis/core_features_venn.tiff', compression = "lzw", 
     units = "in", res = 300, height = 6, width = 6)

venn(core_taxa_venn[2:5], ellipse = TRUE, zcolor = "style", cexil = 1, cexsn = 1)
par(oma = c(0,0,1,0))
title(main = "Core features based on 80% prevalence threshold", outer = T)

dev.off()

# Core features, heatmap ----------------------------------------------------- 
# Use a modified function to get best taxonomic annotation for features
source("down_stream_analysis/R_functions/format_to_besthit.R")

physeq_bestHit <- transform(physeq, "compositional") %>%
  format_to_besthit()

tiff('down_stream_analysis/core_features_heatmap.tiff', compression = "lzw", 
     units = "in", res = 300, height = 5, width = 10)

plot_core(physeq_bestHit, 
          plot.type = "heatmap", 
          colours = rev(RColorBrewer::brewer.pal(5, "Spectral")),
          prevalences = seq(0.05, 1, 0.05), 
          min.prevalence = 0.5) + 
          labs(x = "Detection threshold (%)",
               y = "Features", title = NULL) +
          theme_cowplot() +
          theme(legend.key.size = unit(0.7, "cm"),
                legend.text = element_text(size = 10))

dev.off()

# Alpha-diversity ##############################################################
# Rarefaction curves -----------------------------------------------------------
tiff('down_stream_analysis/rarefaction.tiff', compression = "lzw", 
     units = "in", res = 300, height = 6, width = 10)

par(mar = c(5, 4, 1.5, 8), xpd = TRUE) # display legend outside the plotting area
rarecurve(t(count_tab), step = 100, col = as.character(metadata$color), lwd = 2, ylab = "ASVs", label = F)
title(main = "Rarefaction curves")
legend("right", inset = c(-0.2, 0), legend = levels(metadata$SampleType),
       title = "Sample type", col = levels(metadata$color), cex = 0.8, lwd = 3, lty = 1)
par(xpd = FALSE)
abline(v = min(colSums(count_tab))) # highlight the minimum read count

dev.off()

# Alpha-diversity plot ---------------------------------------------------------
# Gather alpha-diversity indices
adiv <- cbind(obs$data, shn$data, pd$data, evn$data) %>%
  as.data.frame() %>%
  rownames_to_column("#SampleID") %>%
  rename(SampleID = "#SampleID") %>%
  left_join(metadata, ., by = "SampleID") %>%
  gather("alpha_diversity", "value", observed_otus:pielou_e) %>%
  mutate(alpha_diversity = factor(alpha_diversity, levels = c("observed_otus", "pielou_e", "shannon", "faith_pd")))

# make plot
adiv_plot <- ggplot(adiv, aes(x = SampleType, y = value)) +
  geom_boxplot(aes(fill = SampleType), outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  facet_wrap(~alpha_diversity, nrow = 1, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.15))) +  
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Alpha-diversity group significance -------------------------------------------
library(lme4)
library(ggResidpanel)
library(influence.ME)
library(pbkrtest)
library(emmeans)
library(ggsignif)

# Split the dataframe by target genes
adiv_spl <- split(adiv, f = adiv$alpha_diversity)

# Based on the exploratory anaylysis, we'll also log transform the response 
# variables to stabilize the mean-variance relationship
adiv_lmer <- lapply(adiv_spl, 
                    function(x) 
                    lmer(log10(value + 1) ~ SampleType + (1|Sampling_site), REML = T, data = x)) 

### Model diagnostics ###
# 1.Residual analysis 
# Create a panel with all plots available for a model fit using lmer 
pdf("down_stream_analysis/alpha_diversity_group_significance_resid_panel.pdf", 
    width = 14, height = 10) 

lapply(
  seq_along(adiv_lmer), 
  function(x) 
  {
    # Extract titles 
    main <- ggdraw() + draw_label(names(adiv_lmer)[x], fontface='bold')
    # Make residual diagnostic plots
    resid_panel <- resid_panel(adiv_lmer[[x]], plots = "all", qqbands = T)
    # Assemble plots
    plot_grid(main, resid_panel, ncol = 1, rel_heights = c(1, 10))
  }
)

dev.off()

# Create a panel of plots of the residuals versus the predictor variables
pdf("down_stream_analysis/alpha_diversity_group_significance_resid_xpanel.pdf", 
    width = 14, height = 10) 

lapply(
  seq_along(adiv_lmer), 
  function(x) 
  {
    # Extract titles 
    main <- ggdraw() + draw_label(names(adiv_lmer)[x], fontface='bold')
    # Make residual diagnostic plots
    resid_xpanel <- resid_xpanel(adiv_lmer[[x]])
    # Assemble plots
    plot_grid(main, resid_xpanel, ncol = 1, rel_heights = c(1, 10))
  }
)

dev.off()

# 2.Influence analysis 
# Based on the residual analysis, influence analysis is not necessary

# Get p values of the fixed effect using parametric bootstrap 
# Refit models using maximal likelihood estimation
adiv_lmerML <- lapply(adiv_spl, 
                      function(x) 
                      lmer(log10(value + 1) ~ SampleType + (1|Sampling_site), REML = F, data = x))  

# Construct reduced models
adiv_lmerML_nofix <- lapply(adiv_lmerML, function(x) update(x, . ~ . -SampleType))

# Create clusters for parallel computing to speed up parametric bootstrap 
nc <- parallel::detectCores() 
clus <- parallel::makeCluster(rep("localhost", nc)) 

set.seed(1910)
# Parametric bootstrap comparisons 
pb <- map2(adiv_lmerML, adiv_lmerML_nofix, ~PBmodcomp(.x, .y, cl = clus))

# Stop the clusters
parallel::stopCluster(clus)

# Get the p values of parametric bootstrap comparisons
p_fixef <- map_df(pb, ~.x$test["PBtest", "p.value"]) %>%
  gather(key = alpha_diversity, value = p_value) 
  
# post hoc test using emmeans
adiv_emms <- lapply(adiv_lmer, 
  function(x) emmeans(x, ~SampleType) %>% contrast("pairwise") %>% summary())

# Gather significant contrasts for p value annotation
ann <- bind_rows(adiv_emms, .id = "alpha_diversity") %>% 
  filter(p.value <= 0.05) %>% 
  separate(contrast, sep = " - ", c("start", "end")) %>% # start/end position of p value labeling on x axis
  mutate(p.value = formatC(p.value, format = "f", digits = 3), # format digits of p values
         p.value = paste("p = ", p.value), # add label "p =" to p values
         y = c(1900, 2300, 2100, 0.8, 8.5, 9.5, 90, 110, 100)) %>% # position of p value label on y axis 
  mutate(alpha_diversity = factor(alpha_diversity, levels = c("observed_otus", "pielou_e", "shannon", "faith_pd")))

# Add the significant p values. The warning about the missing aesthetics can be ignored.
adiv_plot + geom_signif(data = ann, 
                        aes(xmin = start, xmax = end, annotations = p.value, y_position = y),
                        textsize = 4, 
                        manual = T)

ggsave('down_stream_analysis/alpha_diversity_group_significance.tiff', 
       units = "in", dpi = 300, compression = "lzw", width = 10)

# Beta-diversity ###############################################################
# Beta-diversity visualiztion -------------------------------------------------
metadata$SampleID <- as.numeric(metadata$SampleID)

# Bray-curtis distance based PCoA
ord_bc$data$Vectors %>%
  left_join(metadata, by = "SampleID") %>%
  ggplot(aes(x = PC1, y = PC2, color = SampleType)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  labs(x = paste("PC1: ", round(100*ord_bc$data$ProportionExplained[1], digits = 2), "%"),
       y = paste("PC2: ", round(100*ord_bc$data$ProportionExplained[2], digits = 2), "%"),
       title = "Bray-curtis distance based PCoA", color = "Sample type") +
  scale_fill_brewer(palette = "Dark2") +
  theme_cowplot() +
  panel_border(colour = "black")

ord_bc_ps <- ordinate(physeq, "PCoA", "bray")
plot_ordination(physeq, ord_bc_ps, color = "SampleType") + 
  ggtitle("Bray-curtis distance based PCoA") + 
  geom_point(size = 2) +
  theme_classic() + 
  scale_color_brewer("Location", palette = "Set2") +
  stat_ellipse() 

# Weighted-UniFrac distance based PCoA
ord_wuf$data$Vectors %>%
  left_join(metadata, by = "SampleID") %>%
  ggplot(aes(x = PC1, y = PC2, color = SampleType)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  labs(x = paste("PC1: ", round(100*ord_wuf$data$ProportionExplained[1], digits = 2), "%"),
       y = paste("PC2: ", round(100*ord_wuf$data$ProportionExplained[2], digits = 2), "%"),
       title = "Weighted-UniFrac distance based PCoA", color = "Sample type") +
  scale_fill_brewer(palette = "Dark2") +
  theme_cowplot() +
  panel_border(colour = "black") +
  stat_ellipse()

ord_wUF <- ordinate(physeq, "PCoA", "unifrac", weighted = T)
plot_ordination(physeq, ord_wUF, color = "SampleType") + 
  ggtitle("Weighted UniFrac") + 
  geom_point(size = 2) +
  theme_classic() + 
  scale_color_brewer("Location", palette = "Set2") +
  stat_ellipse()

# Unweighted-UniFrac distance based PCoA
ord_uwuf$data$Vectors %>%
  left_join(metadata, by = "SampleID") %>%
  ggplot(aes(x = PC1, y = PC2, color = SampleType)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  labs(x = paste("PC1: ", round(100*ord_uwuf$data$ProportionExplained[1], digits = 2), "%"),
       y = paste("PC2: ", round(100*ord_uwuf$data$ProportionExplained[2], digits = 2), "%"),
       title = "Unweighted-UniFrac distance based PCoA", color = "Sample type") +
  scale_fill_brewer(palette = "Dark2") +
  theme_cowplot() +
  panel_border(colour = "black")

ord_uwuf_ps <- ordinate(physeq, "PCoA", "unifrac", weighted = F)
plot_ordination(ord_uwuf_ps, ord_uwUF, color = "SampleType") + 
  ggtitle("Unweighted UniFrac") + 
  geom_point(size = 2) +
  theme_classic() + 
  scale_color_brewer("Location", palette = "Set2") +
  stat_ellipse()

# Aitchison distance based roboust PCA
ord_aitchison$data$Vectors %>%
  mutate(SampleID = as.numeric(SampleID)) %>%
  left_join(metadata, by = "SampleID") %>%
  ggplot(aes(x = PC1, y = PC2, color = SampleType)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(size = 4) +
  labs(x = paste("PC1: ", round(100*ord_aitchison$data$ProportionExplained[1], digits = 2), "%"),
       y = paste("PC2: ", round(100*ord_aitchison$data$ProportionExplained[2], digits = 2), "%"),
       title = "Aitchison distance based robust PCA ", color = "Sample type") +
  scale_fill_brewer(palette = "Dark2") +
  theme_cowplot() +
  panel_border(colour = "black")

# Testing homogeneity of multivariate dispersions -------------------------------
# Bray-curtis distance ----------------------------------------------
disp_bc <- betadisper(dist_bc$data, metadata$SampleType, type = "median", bias.adjust = TRUE)
disp_bc

## Plot the distances to centroids on the first two PCoA axes 
disp_bc_labs <- paste("PCoA", 1:length(disp_bc$eig), "(", 
                      round(100 * disp_bc$eig / sum(disp_bc$eig), 2), "%)")

tiff('down_stream_analysis/PCoA_bray_curtis.tiff', compression = "lzw", 
     units = "in", res = 300, height = 6, width = 8)

plot(disp_bc, 
     hull = FALSE, 
     ellipse = TRUE, # use 1 sd data ellipse 
     pch = c(16, 16, 16, 16), # set shape to filled cicles for all samples 
     col = levels(metadata$color), # use pre-defined color scheme
     label = F, # turn off labels for controids
     main = "Bray-curtis distance based PCoA", 
     sub = NULL, # turn off subtitle
     xlab = disp_bc_labs[1], # PCoA1
     ylab = disp_bc_labs[2], # PCoA2
     cex = 1,
     lwd = 2) 
abline(h = 0, v = 0, lty = "dotted")
legend("topright", legend = levels(metadata$SampleType), 
       bty = "n", pch = 16, col = levels(metadata$color))

dev.off()

# Boxplot showing the distances to centroid for each group 
boxplot(disp_bc, xlab = "Sample type", notch = TRUE, col = "lightblue")

# Permutaion test
permutest(disp_bc, pairwise = TRUE, permutations = 999)

# Weighted-UniFrac distance -----------------------------------------
disp_wuf <- betadisper(dist_wuf$data, metadata$SampleType, type = "median", bias.adjust = TRUE)
disp_wuf

## Plot the distances to centroids on the first two PCoA axes 
disp_wuf_labs <- paste("PCoA", 1:length(disp_wuf$eig), "(", 
                      round(100 * disp_wuf$eig / sum(disp_wuf$eig), 2), "%)")

tiff('down_stream_analysis/PCoA_weighted_unifrac.tiff', compression = "lzw", 
     units = "in", res = 300, height = 6, width = 8)

plot(disp_wuf, 
     hull = FALSE, 
     ellipse = TRUE, # use 1 sd data ellipse 
     pch = c(16, 16, 16, 16), # set shape to filled cicles for all samples 
     col = levels(metadata$color), # use pre-defined color scheme
     label = F, # turn off labels for controids
     main = "Weighted-UniFrac distance based PCoA", 
     sub = NULL, # turn off subtitle
     xlab = disp_wuf_labs[1], # PCoA1
     ylab = disp_wuf_labs[2], # PCoA2
     cex = 1,
     lwd = 2) 
abline(h = 0, v = 0, lty = "dotted")
legend("topright", legend = levels(metadata$SampleType), 
       bty = "n", pch = 16, col = levels(metadata$color))

dev.off()

# Boxplot showing the distances to centroid for each group 
boxplot(disp_wuf, xlab = "Sample type", notch = TRUE, col = "lightblue")

# Permutaion test
permutest(disp_wuf, pairwise = TRUE, permutations = 999)

# Unweighted-UniFrac distance -----------------------------------------
disp_uwuf <- betadisper(dist_uwuf$data, metadata$SampleType, type = "median", bias.adjust = TRUE)
disp_uwuf

## Plot the distances to centroids on the first two PCoA axes 
disp_uwuf_labs <- paste("PCoA", 1:length(disp_uwuf$eig), "(", 
                        round(100 * disp_uwuf$eig / sum(disp_uwuf$eig), 2), "%)")

tiff('down_stream_analysis/PCoA_unweighted_unifrac.tiff', compression = "lzw", 
     units = "in", res = 300, height = 6, width = 8)

plot(disp_uwuf, 
     hull = FALSE, 
     ellipse = TRUE, # use 1 sd data ellipse 
     pch = c(16, 16, 16, 16), # set shape to filled cicles for all samples 
     col = levels(metadata$color), # use pre-defined color scheme
     label = F, # turn off labels for controids
     main = "Unweighted-UniFrac distance based PCoA", 
     sub = NULL, # turn off subtitle
     xlab = disp_uwuf_labs[1], # PCoA1
     ylab = disp_uwuf_labs[2], # PCoA2
     cex = 1,
     lwd = 2) 
abline(h = 0, v = 0, lty = "dotted")
legend("topright", legend = levels(metadata$SampleType), 
       bty = "n", pch = 16, col = levels(metadata$color))

dev.off()

# Boxplot showing the distances to centroid for each group 
boxplot(disp_uwuf, xlab = "Sample type", notch = TRUE, col = "lightblue")

# Permutaion test
permutest(disp_uwuf, pairwise = TRUE, permutations = 999)

# Aitchison distance ------------------------------------------------
disp_aitchison <- betadisper(dist_aitchison$data, metadata$SampleType, type = "median", bias.adjust = TRUE)
disp_aitchison

## Plot the distances to centroids on the first two PCoA axes 
disp_aitchison_labs <- paste("PCoA", 1:length(disp_aitchison$eig), "(", 
                             round(100 * disp_aitchison$eig / sum(disp_aitchison$eig), 2), "%)")

tiff('down_stream_analysis/PCoA_Aitchison.tiff', compression = "lzw", 
     units = "in", res = 300, height = 6, width = 8)

plot(disp_aitchison, 
     hull = FALSE, 
     ellipse = TRUE, # use 1 sd data ellipse 
     pch = c(16, 16, 16, 16), # set shape to filled cicles for all samples 
     col = levels(metadata$color), # use pre-defined color scheme
     label = F, # turn off labels for controids
     main = "Aitchison distance based PCoA", 
     sub = NULL, # turn off subtitle
     xlab = disp_aitchison_labs[1], # PCoA1
     ylab = disp_aitchison_labs[2], # PCoA2
     cex = 1,
     lwd = 2) 
abline(h = 0, v = 0, lty = "dotted")
legend("topright", legend = levels(metadata$SampleType), 
       bty = "n", pch = 16, col = levels(metadata$color))

dev.off()

# Boxplot showing the distances to centroid for each group 
boxplot(disp_aitchison, xlab = "Sample type", notch = TRUE, col = "lightblue")

# Permutaion test
permutest(disp_aitchison, pairwise = TRUE, permutations = 999)

# Beta-diversity group significance --------------------------------------------
library(vegan) # [CRAN] # Community Ecology Package 
library(BiodiversityR)

# Nested PERMANOVA based on Bray-curtis distance
nested.npmanova(dist_bc$data ~ SampleType + Sampling_site, data = metadata, permutations = 999)


# Nested PERMANOVA based on Weighted-UniFrac distance


# Nested PERMANOVA based on Unweighted-UniFrac distance 


# Nested PERMANOVA based on Aitchison distance


# Define permutation schemes  
perm <- how(complete = TRUE, # The allowed permutations for the current permutation scheme is 720
            within = Within(type = "none"),
            plots = with(metadata, Plots(strata = Tank, type = "free")))

# Run nested PERMANOVA using unweighted/weighted-unifrac and bray-curtis distance

perm_uwuf <- adonis2(dist_uwuf$data ~ SampleType, data = metadata, permutations = perm)

perm_wuf <- adonis2(dist_wuf$data ~ SampleType, data = metadata, permutations = perm)

perm_bc <- adonis2(dist_bc$data ~ SampleType, data = metadata, permutations = perm)

ad1 <- adonis(dist_uwuf$data ~ SampleType, data = metadata, permutations = perm)
ad2 <- adonis2(dist_uwuf$data ~ SampleType, data = metadata, permutations = perm)


# Differential abundance testing -----------------------------------------------
library(ALDEx2)
library(philr)
library(phylofactor)
library(DESeq2)
source("down_stream_analysis/R_functions/ancom2.R")


## Example
ys_data=read.delim("ys_genera.txt",header=T)

otu_test=ys_data[,c(1,5:665)]
map_test=ys_data[,c(1:4)]


comparison_test=ANCOM.main(OTUdat=otu_test,
                           Vardat=map_test[which(map_test$SEX!="unknown"&map_test$COUNTRY!="GAZ:Venezuela"),],
                           adjusted=TRUE,
                           repeated=F,
                           main.var="COUNTRY",
                           adj.formula="SEX+AGE",
                           repeat.var=NULL,
                           longitudinal=F,
                           multcorr=2,
                           sig=0.05,
                           prev.cut=0.90)

comparison_test$W.taxa


