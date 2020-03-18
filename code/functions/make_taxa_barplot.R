# The make_taxa_barplot() function is modified from the Microbiome.Barplot() function 
# from the MicrobeR package developed by Jordan Bisanz.

#' Visualize taxonomy with a stacked bar plot
#'
#' @description A wrapper function to create a stacked taxa barplot. 
#'    The most abundant features (defaults to 10, based on rowMeans) will be plotted unless specified. 
#'    Anything of over 10 features may be difficult to interpret.
#' @param table A feature table (ASVs/OTUs) with counts where sample names are column names and taxa names are row names.
#' @param metadata A sample metadata where sample names are row names.
#' @param group_by A metadata column to group the samples.
#' @param ntaxa The number of taxa to be plotted.
#' @param nrow The number of rows for faceting barplots.
#' @param plot_mean whether to show average taxa abundance within each level of the grouping variable or not.
#' @param cluster_sample whether to cluster samples based on the similarity of taxonomic compositions or not.
#' @param sample_label Which metadata column to be used as the x axis tick label. 
#'    The default is the unique sample identifier "SampleID". 
#' @param italize_taxa_name whether to italize taxa names or not. 
#'    Taxa containing two under scores "__" in the name prefix are shown as plain text by default. 
#' @return A stacked barplot showing relative abundance of features.
#' @usage make_taxa_barplot(table = otu_table, metadata = metadata, group_by = var, ntaxa = 10, 
#'                          plot_mean = F, cluster_sample, sample_label, italize_taxa_name).
#' @export

make_taxa_barplot <- function(
  table, 
  metadata, 
  group_by, 
  ntaxa,
  nrow,
  plot_mean, 
  cluster_sample,
  sample_label,
  italize_taxa_name){
  
  if(missing(ntaxa) & nrow(table)>10){
    ntaxa = 10} else if (missing(ntaxa)) {
    ntaxa = nrow(table)
      }
  
  if(missing(plot_mean)){plot_mean = F}
  
  if(missing(group_by) & plot_mean == T){stop('The argument "group_by" must be specified if plot_mean = TRUE')}
  
  table <- Make.Percent(table)
  table <- table[order(rowMeans(table), decreasing = T), ]
  if(ntaxa < nrow(table)){ 
    Others <- colSums(table[(ntaxa + 1):nrow(table), ])
    table <- rbind(table[1:ntaxa, ], Others)
  }
  
  forplot <- TidyConvert.ToTibble(table, "Taxa") %>% gather(-Taxa, key = "SampleID", value = "Abundance") 
  forplot$Taxa <- factor(forplot$Taxa,levels = rev(unique(forplot$Taxa)))
  
  if(!missing(metadata) & !missing(group_by)){
    if(TidyConvert.WhatAmI(metadata) == "data.frame" | TidyConvert.WhatAmI(metadata) == "matrix") {metadata <- TidyConvert.ToTibble(metadata, "SampleID")}
    forplot <- inner_join(forplot, metadata, by = "SampleID")
  }
  
  # hierarchical clustering of samples
  if(missing(cluster_sample)){cluster_sample = F}
  
  if(cluster_sample == T & plot_mean == F){
    if(missing(group_by)){
      order <- select(forplot, Taxa, SampleID, Abundance) %>%
        spread(key = SampleID, value = Abundance) %>%
        column_to_rownames("Taxa") %>%
        t() %>%
        as.matrix() %>%
        dist() %>%
        hclust() %>%
        as.dendrogram() %>%
        labels()
      
    } else {
      group_by <- enquo(group_by)
      
      hclust <- select(forplot, Taxa, SampleID, Abundance, !!group_by) %>%
        group_nest(!!group_by) %>%
        mutate(hclust = map(data, ~{
          spread(.x, key = SampleID, value = Abundance) %>%
            column_to_rownames("Taxa") %>%
            t() %>%
            as.matrix() %>%
            dist() %>%
            hclust() %>%
            as.dendrogram() %>%
            labels()
        }))
      
      order <- unlist(hclust$hclust)
      
      forplot <- arrange(forplot, match(SampleID, order)) %>%
        mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
    }
  }
  
  # make an initial plot
  plot <- ggplot(forplot, aes(x = SampleID, y = Abundance, fill = Taxa)) +
            geom_bar(stat = "identity") +
            labs(x = "SampleID", y = "Relative abundance (%)") +
            scale_y_continuous(breaks = 0:10*10, expand = expansion(mult = c(0, 0.02))) + 
            theme_cowplot() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  #legend.justification = "top",
                  legend.text.align = 0) 
  
  # use group means for plotting if plot_mean = T
  if(!missing(group_by) & plot_mean == T){
    group_by <- enquo(group_by)
    plot <- group_by(forplot, !!group_by, Taxa) %>% 
            summarise(Abundance_mean = mean(Abundance)) %>%
            ggplot(aes(x = !!group_by, y = Abundance_mean, fill = Taxa)) +
              geom_bar(stat = "identity") +
              labs(x = group_by, y = "Relative abundance (%)") +
              scale_y_continuous(breaks = 0:10*10, expand = expansion(mult = c(0, 0.02))) + 
              theme_cowplot() +
              theme(legend.text = element_text(size = 10), 
                    #legend.justification = "top",
                    legend.text.align = 0) 
  } 
  
  # add facets to the plot if group_by = T
  if(!missing(group_by) & plot_mean == F){
    if(missing(nrow)){nrow = 1}
    plot <- plot + 
      facet_wrap(vars(!!group_by), nrow = nrow, scales = "free_x")
  }
  
  # fill color manually
  if(ntaxa <= 11){
    plot <- plot + 
      scale_fill_manual(values = brewer.pal(n = 12, name = "Paired")) 
  } else {
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    colors <- getPalette(ntaxa + 1)
    plot <- plot + scale_fill_manual(values = colors)
  }
  
  if(missing(italize_taxa_name)){italize_taxa_name = F}
  
  # make R expressions for costomizing legend text
  if(italize_taxa_name == T){
    labs <- forplot %>%
      mutate(Taxa = ifelse(grepl("__|Others", Taxa), 
                           paste0("plain(", Taxa, ")"), 
                           paste0("italic(", Taxa, ")")),
             # tilde (~) gets handled as a "space" in R expressions
             Taxa = gsub("\\s+", "~", Taxa))
  }
   
  # define color scheme; use italic font for legend text if italize_taxa_name = T 
  if(ntaxa <= 11){
    if(italize_taxa_name == T){
      plot <- plot + scale_fill_manual(values = brewer.pal(n = 12, name = "Paired"),
                                       labels = parse(text = rev(unique(labs$Taxa))))
    } else {
      plot <- plot + scale_fill_manual(values = brewer.pal(n = 12, name = "Paired"))
    }
    
  } else {
    getPalette <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))
    colors <- getPalette(ntaxa + 1)
    if(italize_taxa_name == T){
      plot <- plot + scale_fill_manual(values = colors, labels = parse(text = rev(unique(labs$Taxa))))
    } else {
      plot <- plot + scale_fill_manual(values = colors)
    }
    
  }
  
  # customize the x-axis label and x-axis tick labels
  if(!missing(sample_label) & plot_mean == F){
    sample_label <- enquo(sample_label)
    xlab <- as.data.frame(forplot) %>%
      distinct(SampleID, !!sample_label) %>%
      mutate(sample_label = as.character(!!sample_label))
    
    plot <- plot + 
      scale_x_discrete(breaks = xlab$SampleID, labels = xlab$sample_label) +
      labs(x = sample_label)
    }
  
  return(plot)
}
