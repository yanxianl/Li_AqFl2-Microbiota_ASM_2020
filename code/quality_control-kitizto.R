library(tidyverse) # tidy, manipulate and plot data
library(cowplot) # publication-ready theme for ggplot2 graphics
library(ggstatsplot) # ggplot2 graphics with statistical test details  

# Import and tidy data #########################################################
# Metadata ---------------------------------------------------------------------
mdt <- read_tsv("metadata.tsv") %>% rename(SampleID = `#SampleID`) 
str(mdt)

# FeatureTable -----------------------------------------------------------------
tb <- read_tsv("table-merged-no-chlo-mito-sin-lowPre-with-phyla-taxonomy-relative.tsv", 
               skip = 1)
str(tb)

# Tidy the featureTable
tb <- tb %>%
  # sum up features found in the negative controls
  mutate(neg = `49-N1` + `50-N2`) %>% 
  # remove features not present in negative controls
  filter(neg > 0) %>% 
  # order potential contaminant features by their relative abundance
  arrange(desc(neg)) %>% 
  # remove the newly created variable "neg"
  select(-neg) %>% 
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
taxa <- tb %>%
  filter(SampleID == "taxonomy") %>%
  t() %>%
  data.frame() %>%
  # use first row as the column name
  `colnames<-`(.[1, ]) %>% 
  # remove the first row
  slice(-1) %>% 
  # make 80 (number of samples) replicates of each row 
  uncount(nrow(tb)-1) 

# Merge featureTable and metadata ----------------------------------------------
tb_ful <- inner_join(mdt, tb, by = "SampleID") %>%
  arrange(DNA_extraction_batch) %>%
  gather(featureID, abundance, 18:ncol(.)) %>%
  # convert "abundance" to a numeric variable as percentage
  mutate(abundance = as.numeric(abundance),
         abundance = abundance * 100) %>% 
  # add taxonomy after the featureID column 
  add_column(taxa = taxa[, 1], .after = "featureID") %>% 
  as.data.frame()

str(tb_ful)

# Housekeeping
tb_ful$featureID <- factor(tb_ful$featureID, levels = unique(tb_ful$featureID)) 
tb_ful$SampleID <- factor(tb_ful$SampleID, levels = unique(tb_ful$SampleID)) 

# Correlation between qPCR and Qubit reading ###################################
pdf("Correlation between Cq values and Qubit readings.pdf", 
    width = 8, height = 6)

ggscatterstats(data = filter(tb_ful, featureID == unique(featureID)[1]),
               x = qPCR_Cq, 
               y = Qubit,
               type = "spearman", # type of test that needs to be run
               conf.level = 0.95, # confidence level
               xlab = "Cq value", # label for x axis
               ylab = "Qubit reading (ng/ul)", # label for y axis
               title = "Correlation between Cq values and Qubit readings",
               k = 3, # no. of decimal places in the results
               marginal.type = "density", # type of marginal distribution 
               xfill = "#0072B2", # color fill for x-axis marginal distribution
               yfill = "#009E73", # color fill for y-axis marginal distribution
               xalpha = 0.6, # transparency for x-axis marginal distribution
               yalpha = 0.6, # transparency for y-axis marginal distribution
               point.width.jitter = 0.2, # amount of horizontal jitter for data points
               point.height.jitter = 0.4, # amount of vertical jitter for data points
               messages = F) # turn off messages and notes

dev.off() 

# Barplot ######################################################################
# Barplot can visualize the pravelance and abundance of potential contaminant 
# features in the samples and controls

# Split the dataframe by featureID
tb_ful_spl <- split(tb_ful, f = tb_ful$featureID)

# Make a bar chart for each taxa 
pdf("Relative abundance of features found in negative controls.pdf", 
    width = 16, height = 10) 

lapply(
  seq_along(tb_ful_spl), # add index to the elements in the list 
  function(x) 
  {
    ggplot(tb_ful_spl[[x]], aes(x = SampleID, y = abundance, fill = Sample_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_y_continuous(limits = c(0, NA), expand = expand_scale(mult = c(0, 0.1))) +
      labs(title = unique(tb_ful_spl[[x]][, "taxa"]),
           subtitle = paste0("featureID: ", names(tb_ful_spl)[x]),
           y = "relative abundance (%)") +
      scale_fill_brewer(palette = "Dark2") +
      theme_bw() +
      theme(plot.title = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(size = 12),
            legend.position = "bottom")
  }
)

dev.off() 

# Correlation between Cq values and relative abundance of features #############
# We'll use Cq values as a proxy for microbial DNA concentration in the samples.
# A clear positive correlation is indicative of contaminant features.

# Exclude negative controls and samples with zero count of features 
tb_ful <- filter(tb_ful, abundance != 0 & Sample_type != "Negative_control")

# Split the dataframe by featureID
tb_ful_spl <- split(tb_ful, f = tb_ful$featureID)

# Make correlation plots using ggscatterstats
pdf("Correlation between 16S rRNA gene copy numbers in samples and feature abundance.pdf", 
    width = 14, height = 10) 

lapply(
  seq_along(tb_ful_spl), # add index to elements in the list 
  function(x) 
  {
    # Extract title and subtitle
    main <- ggdraw() + draw_label(unique(tb_ful_spl[[x]]$taxa), fontface='bold')
    sub <- ggdraw() + draw_label(paste0("featureID: ", unique(tb_ful_spl[[x]]$featureID)), fontface='bold')
    
    # Corraltion analysis within each DNA extraction batch
    plist <- tb_ful_spl[[x]] %>%
      group_split(DNA_extraction_batch, keep = T) %>% 
      map(function(x)
        if(sum(!is.na(x$abundance)) < 6){ 
          ggscatterstats(data = x,
                         x = qPCR_Cq, 
                         y = abundance,
                         type = "pearson", # type of test that needs to be run
                         conf.level = 0.95, # confidence level
                         xlab = "Cq value", # label for x axis
                         ylab = "relative abundance (%)", # label for y axis
                         title = paste0("Correlation within DNA extraction batch: ", unique(x$DNA_extraction_batch)),
                         k = 3, # no. of decimal places in the results
                         marginal = F, # whether show ggExtra::ggMarginal() plots or not
                         results.subtitle = F, # turn off stats
                         messages = F) # turn off messages and notes
        } else {
          ggscatterstats(data = x,
                         x = qPCR_Cq, 
                         y = abundance,
                         type = "pearson", 
                         conf.level = 0.95, 
                         xlab = "Cq value", 
                         ylab = "relative abundance (%)", 
                         title = paste0("Correlation within DNA extraction batch: ", unique(x$DNA_extraction_batch)),
                         k = 3, 
                         marginal = F, 
                         messages = F) 
        }
      )
    
    # Corraltion analysis with all samples
    if(sum(!is.na(tb_ful_spl[[x]]$abundance)) < 6){ 
        p <- ggscatterstats(data = tb_ful_spl[[x]],
                            x = qPCR_Cq, 
                            y = abundance,
                            type = "pearson", 
                            conf.level = 0.95, 
                            xlab = "Cq value",
                            ylab = "relative abundance (%)", 
                            title = "Correlation with all samples",
                            k = 3, 
                            marginal = F, 
                            results.subtitle = F, # turn off stats
                            messages = F) 
        } else {
        p <- ggscatterstats(data = tb_ful_spl[[x]],
                            x = qPCR_Cq, 
                            y = abundance,
                            type = "pearson", 
                            conf.level = 0.95, 
                            xlab = "Cq value", 
                            ylab = "relative abundance (%)", 
                            title = "Correlation with all samples",
                            k = 3, 
                            marginal = F, 
                            messages = F) 
        }
    
    # Assemble plots
    plot_grid(main, sub, p, plotlist = plist, ncol = 1, rel_heights = c(1, 1, 10, 10, 10)) # rel_heights values control title margins    
    
  }
)

dev.off() 

# List of contaminant sequences ################################################
# Bsed on the distribution of potential contaminant taxa in the samples and 
# controls, and correlation figures and our prior knowledge 
# knowledge, the following features are considered as contaminants:
contam <- tb_ful %>%
  select(featureID, taxa) %>%
  distinct() %>%
  filter(row_number() %in% c(3, 5, 7, 9, 10, 12, 13, 14, 15, 19, 22, 
                             24, 25, 27, 28, 30, 31, 36, 39, 41, 48))

# Export the featureID of contaminant sequences
write.table(contam, 
            file = 'contaminant-features.tsv', 
            quote = FALSE, 
            sep = '\t',
            row.names = FALSE)
