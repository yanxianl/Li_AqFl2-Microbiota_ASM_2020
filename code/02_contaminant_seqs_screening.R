library(tidyverse) # tidy, manipulate and plot data
library(cowplot) # publication-ready theme for ggplot2 graphics
library(ggstatsplot) # ggplot2 graphics with statistical test details  

# Set system locale to "Greek" so that greek letters can be recognized
Sys.setlocale(category = "LC_ALL", locale = "Greek")

# Import and tidy data #########################################################
# Metadata ---------------------------------------------------------------------
# Skip the 2nd line as it converts all columns to character after importing
col_names <- names(read_tsv("metadata.tsv", n_max = 0))
mdt <- read_tsv("metadata.tsv", col_names = col_names, skip = 2) %>%
  # Get rid of `#` in the column name
  rename(SampleID = `#SampleID`) 

str(mdt)

# FeatureTable -----------------------------------------------------------------
tb <- read_tsv("table-no-chlo-mito-sin-lowPre-with-phyla-tax-rel.tsv", 
               skip = 1)
str(tb)

# Tidy the featureTable
tb <- tb %>%
  # sum up features found in the negative controls
  mutate(neg = rowSums(.[grep("EB|LB", names(.))], na.rm = TRUE)) %>% 
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
  arrange(match(SampleType, c("Blank-library", "Blank-extraction", "Mock",  
                              "REF-DIC", "REF-DIM", "IM100-DIC", "IM100-DIM"))) %>%
  gather(featureID, abundance, 74:ncol(.)) %>%
  # convert "abundance" to a numeric variable as percentage
  mutate(abundance = as.numeric(abundance),
         abundance = abundance * 100) %>% 
  # add taxonomy after the featureID column 
  add_column(taxa = taxa[, 1], .after = "featureID") %>% 
  as.data.frame()

str(tb_ful)

# Some housekeeping conversions of variable types 
tb_ful$featureID <- factor(tb_ful$featureID, levels = unique(tb_ful$featureID)) 
tb_ful$SampleID <- factor(tb_ful$SampleID, levels = unique(tb_ful$SampleID)) 
tb_ful$SampleType <- factor(tb_ful$SampleType, 
                            levels = c("Blank-library",  "Blank-extraction",
                                       "Mock", "REF-DIC", "REF-DIM", "IM100-DIC", 
                                       "IM100-DIM")) 

# Barplot ######################################################################
# Barplot can visualize the pravelance and abundance of potential contaminant 
# features in the samples and controls

# Split the dataframe by "taxa"
tb_ful_spl <- split(tb_ful, f = tb_ful$featureID)

# Make a bar chart for each taxa and save results as a pdf file
pdf(file = "bar_charts.pdf", width = 16, height = 10) 

lapply(
  seq_along(tb_ful_spl), # add index to the elements in the list 
  function(x) 
  {
    ggplot(tb_ful_spl[[x]], aes(x = SampleID, y = abundance, fill = SampleType)) +
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

# Correlation analysis #########################################################
# We'll use Cq values as a proxy for microbial DNA concentration in the samples.
# A clear positive correlation is indicative of contaminant features.

# Before we proceed with correlation analysis, we'll Check Cq values and see if
# they make sense. Normally, we expect lower Cq values for digesta and mock, and 
# higher Cq values for mucosa and negative controls.
tb_ful %>% 
  filter(featureID == unique(.$featureID)[1]) %>%
  drop_na(qPCRCqValue) %>%
  ggplot(aes(x = SampleType, y = qPCRCqValue)) +
  geom_boxplot(aes(fill = SampleType), outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  labs(y = "Cq value", title = "Quantification of 16S rRNA gene by qPCR") +
  theme_cowplot() +
  scale_fill_brewer(palette = "Dark2")

ggsave('Cq_boxPlot.pdf', width = 10, height = 8, units = "in", dpi = 300)

# Correlatin analysis #########################################################
# Exclude samples with zero count of feature
tb_ful <- filter(tb_ful, abundance != 0)

# Split the dataframe by featureID
tb_ful_spl <- split(tb_ful, f = tb_ful$featureID)[1:3]

# Make correlation plots
pdf(file = "correlation_analysis.pdf", width = 14, height = 10) 

lapply(
  seq_along(tb_ful_spl), # add index to elements in the list 
  function(x) 
  {
    # extract taxonomy and featureID for each plot
    taxa <- ggdraw() + draw_label(unique(tb_ful_spl[[x]][, "taxa"]), 
                                  fontface='bold')
    
    featureID <- ggdraw() + draw_label(paste0("featureID: ", names(tb_ful_spl)[x]), 
                                       fontface='bold')
    
    # correlation using all samples
    p1 <- tb_ful_spl[[x]]%>% 
      ggplot(aes(x = qPCRCqValue, y = abundance)) +
      geom_point(aes(color = SampleType)) +
      geom_smooth(method = "lm") +
      scale_y_continuous(limits = c(0, NA), 
                         expand = expand_scale(mult = c(0, 0.1))) + 
      labs(x = "Cq value", y = "relative abundance (%)",
           title  = "Correlation with all samples") +
      scale_fill_brewer(palette = "Dark2") +
      theme_bw() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            legend.position = "bottom")
    
    # Corraltion within sample type
    p2 <- tb_ful_spl[[x]] %>%
      filter(!SampleType %in% c("Blank-library", "Blank-extraction", "Mock")) %>%
      grouped_ggscatterstats(data = .,
                             x = qPCRCqValue, 
                             y = abundance,
                             type = "spearman", # type of test that needs to be run
                             conf.level = 0.95, # confidence level
                             xlab = "Cq value", # label for x axis
                             ylab = "relative abundance (%)", # label for y axis
                             #bf.message = TRUE, # display bayes factor message
                             k = 3, # no. of decimal places in the results
                             marginal.type = "density", # type of marginal distribution to be displayed
                             xfill = "#0072B2", # color fill for x-axis marginal distribution
                             yfill = "#009E73", # color fill for y-axis marginal distribution
                             xalpha = 0.6, # transparency for x-axis marginal distribution
                             yalpha = 0.6, # transparency for y-axis marginal distribution
                             point.width.jitter = 0.2, # amount of horizontal jitter for data points
                             point.height.jitter = 0.4, # amount of vertical jitter for data points
                             grouping.var = SampleType, # grouping variable
                             title.prefix = "SampleType",
                             messages = FALSE, # turn off messages and notes
                             ncol = 2, # number of columns for graph panels
                             title.text = "Correlation within sample type",
                             title.size = 14)
    
    # assemle plots
    plot_grid(taxa, featureID, p1, p2, ncol = 1, rel_heights = c(1, 1, 18, 24)) 
    
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
  filter(!row_number() %in% c(3,4,5,16,17,19,20,23))

# Export the featureID of contaminant sequences
write.table(contam, 
            file = 'contaminant-features.tsv', 
            quote = FALSE, 
            sep = '\t',
            row.names = FALSE)

# Check if over-dominance of contaminant taxa affect diversity estimation ######
library(reshape2)
library(dplyr)
library(ggplot2)

df <- read.csv("C:/Users/yanxianl/OD/AquaFly/AqFl2/AqFl2-alpha-diversity.csv", 
               header = T, na.strings = c(""))
head(df)

# convert numeric variable to string
df <- within(df, {
  SampleID <- as.character(SampleID)
})

# convert to log data format for ggplot2
dfm <- melt(df, id = 1:5)
head(dfm)
# Reorder rows in the desired order
dfm$Diet <- factor(dfm$Diet, levels = c("REF", "IM100"))
dfm$SampleType <- factor(dfm$SampleType, levels = c("REF-DIC", "REF-DIM", "IM100-DIC", "IM100-DIM"))
dfm$variable <- factor(dfm$variable, levels = c("observed_otus", "faith_pd", "shannon", "pielou_e"))

########### plotting box plot######### Without dots or whisker
p_SampleType <- ggplot(dfm, aes(x = Pseudomonas_veronii, y = value)) + 
  geom_point() +
  geom_smooth(method = 'lm', color = 'black', linetype = 'dashed') +
  facet_grid(variable~SampleType, scales = "free_y") + 
  theme(axis.title.x = element_text(size = 20, face = "bold"),
        axis.title.y = element_text(size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  scale_x_continuous(limits = c(0, 1)) +
  #ggtitle("Correlation analysis")
  
  p_SampleOrigin <- ggplot(dfm, aes(x = Pseudomonas_veronii, y = value)) + 
  geom_point(aes(color = Diet)) +
  geom_smooth(method = 'lm', color = 'black', linetype = 'dashed') +
  facet_grid(variable~SampleOrigin, scales = "free_y") + 
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + 
  ggtitle("Correlation analysis")

p_Diet <- ggplot(dfm, aes(x = Pseudomonas_veronii, y = value)) + 
  geom_point(aes(color = SampleOrigin)) +
  geom_smooth(method = 'lm', color = 'black', linetype = 'dashed') +
  facet_grid(variable~Diet, scales = "free_y") + 
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) +  
  ggtitle("Correlation analysis")
