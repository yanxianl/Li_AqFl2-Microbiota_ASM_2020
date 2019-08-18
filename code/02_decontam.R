library(here) # [CRAN] A replacement for 'file.path()', locating the files relative to the project root.
library(tidyverse) # [CRAN] Easily Install and Load the 'Tidyverse' 
library(cowplot) # [CRAN] # Streamlined Plot Theme and Plot Annotations for 'ggplot2' 
library(ggstatsplot) # [CRAN] # 'ggplot2' Based Plots with Statistical Details  

# Set system locale to "Greek" so that greek letters in the metadata can be recognized
Sys.setlocale(category = "LC_ALL", locale = "Greek")

# Import and tidy data #########################################################
# Metadata ---------------------------------------------------------------------
mtd <- read_tsv(here("data", "metadata.tsv"), comment = "#q2:type") %>%
  # Get rid of `#` in the column name
  rename(SampleID = `#SampleID`) 

str(mtd)

# FeatureTable -----------------------------------------------------------------
tab <- read_tsv(here("data", "qiime2", "exported-files", "table-no-chlo-mito-sin-lowPre-with-phyla-tax-rel.tsv"), 
                skip = 1)
str(tab)

# Tidy the featureTable
tab <- tab %>%
  # sum up features found in the negative and positive controls
  mutate(Ctrl = rowSums(.[grep("EB|LB|M", names(.))], na.rm = TRUE)) %>% 
  # remove features not present in negative and positive controls
  filter(Ctrl > 0) %>% 
  # sort features by their relative abundance
  arrange(desc(Ctrl)) %>% 
  # remove the newly created variable "Ctrl"
  select(-Ctrl) %>% 
  # the following 4 lines of codes transpose the dataframe
  rownames_to_column %>% 
  # prevent change of the row order after spreading the dataframe
  mutate(rowname = factor(rowname, levels = unique(rowname))) %>% 
  gather(var, value, -rowname) %>%
  spread(rowname, value) %>%
  # use first row as the column name
  `colnames<-`(.[1, ]) %>% 
  # remove the first row
  slice(-1) %>% 
  rename(SampleID = `#OTU ID`)

# Extract the taxonomy and convert to the long format
tax <- tab %>%
  filter(SampleID == "taxonomy") %>%
  t() %>%
  data.frame() %>%
  # use first row as the column name
  `colnames<-`(.[1, ]) %>% 
  # remove the first row
  slice(-1) %>% 
  # make 80 (number of samples) replicates of each row 
  uncount(nrow(tab)-1) 

# Merge featureTable and metadata ----------------------------------------------
tab_merged <- inner_join(mtd, tab, by = "SampleID") %>%
  arrange(match(SampleType, c("REF-DIC", "REF-DIM", "IM-DIC", "IM-DIM",  
                              "Blank-library", "Blank-extraction", "Mock"))) %>%
  gather(featureID, abundance, 74:ncol(.)) %>%
  # convert "abundance" to a numeric variable as percentage
  mutate(abundance = as.numeric(abundance),
         abundance = abundance * 100) %>% 
  # add taxonomy after the featureID column 
  add_column(tax = tax[, 1], .after = "featureID") %>% 
  as.data.frame()

str(tab_merged)

# Conversions of variable types 
tab_merged$featureID <- factor(tab_merged$featureID, levels = unique(tab_merged$featureID)) 
tab_merged$SampleID <- factor(tab_merged$SampleID, levels = unique(tab_merged$SampleID)) 
tab_merged$SampleType <- factor(tab_merged$SampleType, 
                            levels = c("REF-DIC", "REF-DIM", "IM-DIC", "IM-DIM",
                                       "Blank-library",  "Blank-extraction",
                                       "Mock")) 

# Barplot ######################################################################
# Barplot can visualize the pravelance and abundance of potential contaminant 
# features in the samples and controls

# Split the dataframe by "featureID"
tab_merged_spl <- split(tab_merged, f = tab_merged$featureID)

# Make bar plots for each feature
pdf(here("data", "qiime2", "decontam", "taxa_bar_plots_unstacked.pdf"), width = 16, height = 10) 

lapply(
  seq_along(tab_merged_spl), # add index to the elements in the list 
  function(x) 
  {
    # Make bar plots including samples and negative controls
    p1 <- tab_merged_spl[[x]] %>%
      filter(SampleType != "Mock") %>%
      ggplot(aes(x = SampleID, y = abundance, fill = SampleType)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) + 
        labs(title = unique(tab_merged_spl[[x]][, "tax"]),
             subtitle = paste0("featureID: ", names(tab_merged_spl)[1]),
             y = "relative abundance (%)") +
        scale_fill_brewer(palette = "Dark2", drop = F) +
        guides(fill = guide_legend(nrow = 1)) +
        theme_bw() +
        theme(plot.title = element_text(size = 12),
              axis.title.x = element_text(size = 12),
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title.y = element_text(size = 12),
              legend.position = "bottom")
    
    # Make bar plots including mock only
    p2 <- tab_merged_spl[[x]] %>%
      filter(SampleType == "Mock") %>%
      ggplot(aes(x = SampleID, y = abundance, fill = SampleType)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) + 
      labs(x = "", y = "") +
      scale_fill_brewer(palette = "Dark2", drop = F) +
      guides(fill = guide_legend(nrow = 1)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
    
    # Assemble plots
    plot_grid(p1, p2 + theme(legend.position="none"), 
              nrow = 1, align = 'h', axis = "bt", # align plots with equal y axia height
              rel_widths = c(15, 1))
  }
)

dev.off() 

# Correlation analysis #########################################################
# We'll use Cq values as a proxy of microbial DNA concentration in the samples.
# A clear positive correlation is indicative of contamination.

# Before we proceed with correlation analysis, we'll Check Cq values and see if
# they make sense. Normally, we expect lower Cq values for digesta and mock, and 
# higher Cq values for mucosa and negative controls.
tab_merged %>% 
  # Get Cq values from one of the features
  filter(featureID == unique(.$featureID)[1]) %>%
  drop_na(qPCRCqValue) %>%
  ggplot(aes(x = SampleType, y = qPCRCqValue)) +
  geom_boxplot(aes(fill = SampleType), outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(y = "Cq value", title = "Quantification of 16S rRNA gene copy numbers by qPCR") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave(here("exploratory", "qPCR_Cq_values.pdf"), width = 10, height = 8, units = "in", dpi = 300)

# Correlatin analysis ----------------------------------------------------------
# Exclude samples with zero count of features
tab_merged <- filter(tab_merged, abundance != 0)

# Split the dataframe by featureID
tab_merged_spl <- split(tab_merged, f = tab_merged$featureID)

# Make correlation plots
pdf(here("data", "qiime2", "decontam", "correlation_analysis.pdf"), width = 15, height = 12) 

lapply(
  seq_along(tab_merged_spl), # add index to elements in the list 
  function(x) 
  {
    # Extract taxonomy and featureID as title and subtitle, respectively
    main <- ggdraw() + draw_label(unique(tab_merged_spl[[x]]$tax), fontface='bold')
    sub <- ggdraw() + draw_label(paste0("featureID: ", unique(tab_merged_spl[[x]]$featureID)), fontface='bold')
    
    # Correlation using all samples
    p_all <- tab_merged_spl[[x]]%>% 
      ggplot(aes(x = qPCRCqValue, y = abundance)) +
      geom_point(aes(color = SampleType)) +
      geom_smooth(method = "lm") +
      scale_y_continuous(limits = c(0, NA), 
                         expand = expand_scale(mult = c(0, 0.1))) + 
      labs(x = "Cq value", y = "relative abundance (%)",
           title  = "Correlation with all samples") +
      scale_fill_brewer(palette = "Dark2") +
      guides(colour = guide_legend(nrow = 1)) +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            legend.position = "bottom")
    
    # Corraltion analysis within each sample type
    p_sub <- tab_merged_spl[[x]] %>%
      filter(!SampleType %in% c("Blank-library", "Blank-extraction", "Mock")) %>%
      # make a list column grouped by 'SampleType'
      group_split(SampleType, keep = T) %>% 
      map(function(x)
        # turn off stats if the number of observation is less than 6
        if(sum(!is.na(x$abundance)) < 6){ 
          ggscatterstats(data = x,
                         x = qPCRCqValue, 
                         y = abundance,
                         type = "spearman", # type of test that needs to be run
                         conf.level = 0.95, # confidence level
                         xlab = "Cq value", # label for x axis
                         ylab = "relative abundance (%)", # label for y axis
                         title = paste0("Correlation within the sample type: ", unique(x$SampleType)),
                         k = 3, # no. of decimal places in the results
                         marginal = F, # whether show ggExtra::ggMarginal() plots or not
                         results.subtitle = F, # turn off stats
                         messages = F) # turn off messages and notes
        } else {
          # turn on stats if the number of observation is no less than 6
          ggscatterstats(data = x,
                         x = qPCRCqValue, 
                         y = abundance,
                         type = "spearman", 
                         conf.level = 0.95, 
                         xlab = "Cq value", 
                         ylab = "relative abundance (%)", 
                         title = paste0("Correlation within the sample type: ", unique(x$SampleType)),
                         k = 3, 
                         marginal = F, 
                         messages = F) 
        }
      )
    
    
    # Assemble plots
    panel_1 <- plot_grid(main, sub, p_all, ncol = 1, rel_heights = c(1, 1, 12)) 
    panel_2 <- plot_grid(plotlist = p_sub, ncol = 2, rel_heights = c(1, 1)) 
    plot_grid(panel_1, panel_2, ncol = 1, rel_heights = c(1, 1.2))
    
  }
)

dev.off() 

# List of contaminant sequences ################################################
# Bsed on the evaluation of the taxa bar plots, correlation plots and our prior 
# knowledge of comman contamiant taxa found in our lab, the following features
# are considered as contaminant sequences in the samples:
contam_sample <- tab_merged %>%
  select(featureID, tax) %>%
  distinct() %>%
  filter(row_number() %in% c(1, 2, 4, 5, 22, 28, 30, 34, 35, 36, 38, 39, 41, 42,
                             43, 45, 46, 48, 49, 50, 51, 53, 54, 55, 57))

# The mock comes with a garunteed impurity level of <0.01% (by DNA abundance). 
# Therefore, as long as we observe any alien taxa present at >0.01% in the mock, 
# we can conclude that they are introduced artificially by the workflow. The 
# following features are considered as contaminant features in the mock:
contam_mock <- tab_merged %>%
  select(featureID, tax) %>%
  distinct() %>%
  # the row numbers of the non-contaminant ASVs are used for reverse selection 
  filter(row_number() %in% c(1, 4, 22, 28, 30, 34, 35, 36, 38, 39, 41, 42, 43, 44, 
                             45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57))

# Export featureIDs of the contaminant sequences
write.table(contam_sample, 
            file = here("data", "qiime2", "decontam", "contaminant-features-sample.tsv"), 
            quote = FALSE, 
            sep = '\t',
            row.names = FALSE)

write.table(contam_mock, 
            file = here("data", "qiime2", "decontam", "contaminant-features-mock.tsv"), 
            quote = FALSE, 
            sep = '\t',
            row.names = FALSE)
