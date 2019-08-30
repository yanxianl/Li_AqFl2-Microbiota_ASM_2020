# The make_taxa_barplot() function is modified from the Microbiome.Barplot() function 
# from the MicrobeR package developed by Jordan Bisanz.

#' Visualize taxonomy with a stacked bar plot
#'
#' @description A wrapper function to create a stacked taxa barplot. 
#'    The most abundant features (defaults to 10, based on rowMeans) will be plotted unless specified. 
#'    Anything of over 10 features will use default coloring which may be difficult to interpret.
#' @param count_table A feature table (ASVs/OTUs) with counts where sample names are column names and feature IDs are in row names
#' @param metadata A sample metadata where sample names are row names
#' @param group A metadata variable to group the samples
#' @param ntaxa The number of features to be displayed on the barplot
#' @param plot_mean whether to show average feature abundance by the grouping variable or not
#' @return A stacked barplot showing relative abundance of features
#' @usage make_taxa_barplot(count_table = otu_table, metadata = metadata, group = var, ntaxa = 10, plot_mean = F)
#' @export

make_taxa_barplot <- function(count_table, metadata, ntaxa, group, plot_mean){
  
  if(missing(ntaxa) & nrow(count_table)>10){ntaxa = 10}
  else if (missing(ntaxa)) {ntaxa = nrow(count_table)}
  
  if(missing(plot_mean)){plot_mean = F}
  
  if(missing(group) & plot_mean == T){stop('The argument "group" must be specified if plot_mean = TRUE')}
  
  count_table <- Make.Percent(count_table)
  count_table <- count_table[order(rowMeans(count_table), decreasing = T), ]
  if(ntaxa < nrow(count_table)){ 
    Others <- colSums(count_table[(ntaxa + 1):nrow(count_table), ])
    count_table <- rbind(count_table[1:ntaxa, ], Others)
  }
  
  forplot <- TidyConvert.ToTibble(count_table, "Taxa") %>% gather(-Taxa, key = "SampleID", value = "Abundance") 
  forplot$Taxa <- factor(forplot$Taxa,levels = rev(unique(forplot$Taxa)))
  
  if(!missing(metadata) & !missing(group)){
    if(TidyConvert.WhatAmI(metadata) == "data.frame" | TidyConvert.WhatAmI(metadata) == "matrix") {metadata <- TidyConvert.ToTibble(metadata, "SampleID")}
    forplot <- inner_join(forplot, metadata, by = "SampleID")
  }
  
  plot <- ggplot(forplot, aes(x = SampleID, y = Abundance, fill = Taxa)) +
            geom_bar(stat = "identity") +
            labs(x = "SampleID", y = "Relative abundance (%)") +
            scale_y_continuous(breaks = 0:10*10, expand = expand_scale(mult = c(0, 0.02))) + 
            theme_cowplot() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  legend.text = element_text(size = 10),
                  legend.justification = "top") 
  
  if(!missing(group) & plot_mean == T){
    group <- enquo(group)
    plot <- group_by(forplot, !!group, Taxa) %>% 
            summarise(Abundance_mean = mean(Abundance)) %>%
            ggplot(aes(x = !!group, y = Abundance_mean, fill = Taxa)) +
              geom_bar(stat = "identity") +
              labs(x = group, y = "Relative abundance (%)") +
              scale_y_continuous(breaks = 0:10*10, expand = expand_scale(mult = c(0, 0.02))) + 
              theme_cowplot() +
              theme(legend.text = element_text(size = 10),
                    legend.justification = "top")  
  } 
  
  if(ntaxa <= 12){
    plot <- plot + scale_fill_brewer(palette = "Paired")
  } else {
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    colors <- getPalette(ntaxa + 1)
    plot <- plot + scale_fill_manual(values = colors)
  }
  
  if(!missing(group) & plot_mean == F){
    group <- enquo(group)
    plot <- plot + facet_grid(cols = vars(!!group), scales = "free", space = "free")
  }
  return(plot)
}
