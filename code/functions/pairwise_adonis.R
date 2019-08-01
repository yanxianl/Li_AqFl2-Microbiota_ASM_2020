library(pairwiseAdonis) # [not installed on this machine] # not installed on this machine
metadata$Tank <- as.factor(metadata$Tank)
metadata <- as.data.frame(metadata)

pw <- pairwise.adonis(dist_uwuf$data, metadata$Sample_type, reduce='Ref', p.adjust.m='fdr')

summary(pw)
pairwise.adonis2(dist_uwuf$data~Sample_type, data = metadata, permutations = perm)

co <- combn(unique(as.character(metadata$Sample_type)),2)
matrix.order <- colnames(as.matrix(dist_aitchison$data))
metadata <- metadata %>% arrange(match(SampleID, matrix.order)) # sampleID order in distance matrix and metadata needs to be identical

for(elem in 1:ncol(co)){
  dist.sub <- as.matrix(dist_aitchison$data)[metadata$Sample_type %in% c(as.character(co[1,elem]), as.character(co[2,elem])),
                                             metadata$Sample_type %in% c(as.character(co[1,elem]), as.character(co[2,elem]))]
  
  metadata.sub <- filter(metadata, Sample_type %in% c(as.character(co[1,elem]), as.character(co[2,elem])))
  
  perm <- how(complete = TRUE, # Allow complete enumeration of permutations (720)
              within   = Within(type = "none"),
              plots    = with(metadata.sub, Plots(strata = Tank, type = "free")))                
  
  ad <- adonis2(dist.sub ~ Sample_type, data = metadata.sub, permutations = perm)
  print(ad)
}

ad <- adonis(dist.sub ~ Sample_type, data = metadata.sub, permutations = perm)
pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
Df <- c(Df,ad$aov.tab[1,1])
SumsOfSqs <- c(SumsOfSqs, ad$aov.tab[1,2])
F.Model <- c(F.Model,ad$aov.tab[1,4]);
R2 <- c(R2,ad$aov.tab[1,5]);
p.value <- c(p.value,ad$aov.tab[1,6])

p.adjusted <- p.adjust(p.value,method=p.adjust.m)

sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'
pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)

if(!is.null(reduce)){
  pairw.res <- subset (pairw.res, grepl(reduce,pairs))
  pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
  
  sig = c(rep('',length(pairw.res$p.adjusted)))
  sig[pairw.res$p.adjusted <= 0.1] <-'.'
  sig[pairw.res$p.adjusted <= 0.05] <-'*'
  sig[pairw.res$p.adjusted <= 0.01] <-'**'
  sig[pairw.res$p.adjusted <= 0.001] <-'***'
  pairw.res <- data.frame(pairw.res[,1:7],sig)
}
class(pairw.res) <- c("pwadonis", "data.frame")
return(pairw.res)
} 