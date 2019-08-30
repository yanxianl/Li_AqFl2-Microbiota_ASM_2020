# The make_barplot() function is modified from the Microbiome.Barplot() function 
# from the MicrobeR package developed by Jordan Bisanz.

#' Visualize taxonomy with a stacked bar plot
#'
#' @description A wrapper function to create a stacked taxa barplot. 
#'    The most abundant features (defaults to 10, based on rowMeans) will be plotted unless specified. 
#'    Anything of over 10 features will use default coloring which may be difficult to interpret.
#' @param table A feature table (ASVs/OTUs) with counts where sample names are column names and feature IDs are in row names
#' @param metadata A sample metadata where sample names are row names
#' @param grouping_var A metadata variable to group the samples
#' @param taxa_to_plot The number of features to be displayed on the barplot
#' @param collapse whether to show average feature abundance by the grouping variable or not
#' @return A stacked barplot showing relative abundance of features
#' @usage make_barplot(table = table, metadata = metadata, grouping_var = var, taxa_to_plot = 10, collapse = F)
#' @export

make_barplot <- function(table, metadata, grouping_var, taxa_to_plot, collapse){
  
  if(missing(taxa_to_plot) & nrow(table)>10){taxa_to_plot = 10}
  else if (missing(taxa_to_plot)) {taxa_to_plot = nrow(table)}
  
  if(missing(collapse)){collapse = F}
  
  if(missing(grouping_var) & collapse == T){stop('The argument "grouping_var" must be specified if collapse = TRUE')}
  
  table <- Make.Percent(table)
  table <- table[order(rowMeans(table), decreasing = T), ]
  if(taxa_to_plot < nrow(table)){ 
    Others <- colSums(table[(taxa_to_plot + 1):nrow(table), ])
    table <- rbind(table[1:taxa_to_plot, ], Others)
  }
  
  forplot <- TidyConvert.ToTibble(table, "Taxa") %>% gather(-Taxa, key = "SampleID", value = "Abundance") 
  forplot$Taxa <- factor(forplot$Taxa,levels = rev(unique(forplot$Taxa)))
  
  if(!missing(metadata) & !missing(grouping_var)){
    if(TidyConvert.WhatAmI(metadata) == "data.frame" | TidyConvert.WhatAmI(metadata) == "matrix") {metadata <- TidyConvert.ToTibble(metadata, "SampleID")}
    forplot <- inner_join(forplot, metadata, by = "SampleID")
  }
  
  plot <- ggplot(forplot, aes(x = SampleID, y = Abundance, fill = Taxa)) +
            geom_bar(stat = "identity") +
            labs(x = "SampleID", y = "Relative abundance (%)") +
            scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) + 
            theme_cowplot() +
            guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  legend.text = element_text(size = 10),
                  legend.position = "right") 
  
  if(!missing(grouping_var) & collapse == T){
    grouping_var <- enquo(grouping_var)
    plot <- group_by(forplot, !!grouping_var, Taxa) %>% 
            summarise(Abundance_mean = mean(Abundance)) %>%
            ggplot(aes(x = !!grouping_var, y = Abundance_mean, fill = Taxa)) +
              geom_bar(stat = "identity") +
              labs(x = grouping_var, y = "Relative abundance (%)") +
              scale_y_continuous(expand = expand_scale(mult = c(0, 0.02))) + 
              theme_cowplot() +
              guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
              theme(legend.text = element_text(size = 10),
                    legend.position = "right")  
  } 
  
  if(taxa_to_plot <= 12){
    plot <- plot + scale_fill_brewer(palette = "Paired")
  } else {
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    colors <- getPalette(taxa_to_plot + 1)
    plot <- plot + scale_fill_manual(values = colors)
  }
  
  if(!missing(grouping_var) & collapse == F){
    grouping_var <- enquo(grouping_var)
    plot <- plot + facet_grid(cols = vars(!!grouping_var), scales = "free", space = "free")
  }
  return(plot)
}
