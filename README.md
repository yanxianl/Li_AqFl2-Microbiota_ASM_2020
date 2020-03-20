<!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yanxianl/Li_AqFl2-Microbiota_ASM_2020/master?urlpath=rstudio)
<!-- badges: end -->

## Differential response of digesta- and mucosa-associated intestinal microbiota to dietary black soldier fly (*Hermetia illucens*) larvae meal in seawater phase Atlantic salmon (*Salmo salar*)

Knowledge on the nutritional value of black soldier fly (*Hermetia illucens*) as a feed ingredient for fish is accumulating, yet less is known regarding how it may modulate fish intestinal microbiota. In a 16-week seawater feeding trial, Atlantic salmon (*Salmo salar*) were fed either a commercially-relevant reference diet, or an insect meal diet wherein black soldier fly larvae meal was included to replace all the fish meal and most of the pea protein concentrate in the reference diet. The digesta- and mucosa-associated distal intestinal microbiota were profiled by sequencing V1-2 regions of the 16S rRNA gene. Regardless of diet, we observed substantial differences between digesta- and mucosa-associated intestinal microbiota. Microbial richness and diversity were much higher in the digesta than the mucosa. The insect meal diet altered the distal intestinal microbiota assemblage resulting in higher microbial richness and diversity. The diet effect was dependent on the sample origin, with digesta-associated intestinal microbiota showing more pronounced changes than the mucosa-associated microbiota. Multivariate association analysis identified two mucosa-enriched taxa, *Brevinema andersonii* and unclassified Spirochaetaceae, associated with the expression of genes related to immune responses and barrier function in the distal intestine, respectively. The insect meal diet increased the relative abundance of several taxa that have been reported in fish, including *Actinomyces*, *Bacillus*, *Brevibacterium*, *Corynebacterium 1* and *Enterococcus*. In conclusion, the response of digesta- and mucosa-associated intestinal microbiota to dietary inclusion of insect meal was different, with the latter being more resilient to the dietary change.

### Overview
Here's an overview of the file organization in this project.
```
root
├── code                              # all the scripts used for the analysis
│   ├── functions                     # functions for automating tasks
│   ├── utilities                     # utility scripts for miscellaneous tasks
│   ├── 00_setup.bash                 # download raw and reference data for the analysis
│   ├── 01_dada2.Rmd                  # sequence denoising by the dada2 pipeline
│   ├── 02_qiime2_part1.bash          # taxonomic assignment in qiime2
│   ├── 03_preprocessing.Rmd          # feature table filtering    
│   ├── 04_qiime2_part2.bash          # phylogeny and core-metrics-results
│   ├── 05_qiime2R.Rmd                # export qiime2 artifacts into R
│   ├── 06_taxonomy.Rmd               # taxonomic analysis
│   ├── 07_alpha-diversity.Rmd        # alpha-diversity visualization and statistical analysis
│   ├── 08_beta-diversity.Rmd         # beta-diversity visualization and statistical analysis
│   ├── 09_metadata_association.Rmd   # association testing between microbial clades and sample metadata
│   └── README.md
├── data               # all the data, including raw, reference and intermediate data
│   ├── metadata.tsv   # sample metadata
│   ├── raw            # raw data
│   ├── reference      # reference data
│   ├── qPCR           # qPCR assay reports, plat-calibration and Cq values
│   ├── dada2          # outputs from dada2 including the representative sequences and feature table
│   ├── qiime2         # outputs from qiime2
│   ├── preprocessing  # plots for the identification of contaminants; filtered feature table   
│   ├── qiime2R        # .RData containing outputs from qiime2
│   ├── permanova      # input data and results of the PERMANOVA
│   └── maaslin2       # default outputs from the maaslin2 program
├── image   # pictures/photos relevant to the analysis
├── result  # final results published with the paper
│   ├── figures    
│   ├── tables     
│   └── README.md 
├── LICENSE.md  # license for the copyright of scripts published in this repository
└── README.md
```
### How to regenerate this repository

#### Dependencies and locations
* [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) should be located in the default HOME directory of the linux instance.
* [grabseqs (0.7.0)](https://github.com/louiejtaylor/grabseqs) should be installed via the Miniconda3.
* [QIIME2 (2020.2)](https://docs.qiime2.org/2020.2/) should be installed within a Miniconda3 environment named as `qiime2-2020.2`.
  * QIIME2 library: [DEICODE (0.2.3)](https://library.qiime2.org/plugins/deicode/19/) should be installed within the conda environment of qiime2 (`qiime2-2020.2`).
* [Pandoc (1.12.4.2)](https://pandoc.org/index.html) should be located in the user's PATH.
* R (3.6.3) should be located in the user's PATH.
* R packages (name_version): 
  * ape_5.3
  * biomformat_1.12.0 
  * circlize_0.4.8
  * ComplexHeatmap_2.0.0
  * cowplot_1.0.0
  * dada2_1.12.1
  * DT_0.11
  * EMAtools_0.1.3 
  * emmeans_1.4.4
  * factoextra_1.0.6 
  * ggResidpanel_0.3.0
  * ggsignif_0.6.0 
  * ggstatsplot_0.3.1
  * gridExtra_2.3
  * gt_0.1.0 
  * here_0.1
  * knitr_1.27
  * lmerTest_3.1-1
  * lsr_0.5 
  * Maaslin2_0.99.12 
  * MicrobeR_0.3.1
  * microbiome_1.6.0 
  * PerformanceAnalytics_1.5.3
  * philr_1.10.1
  * phyloseq_1.28.0 
  * picante_1.8
  * plotly_4.9.1
  * qiime2R_0.99.13
  * RColorBrewer_1.1-2
  * rlang_0.4.4 
  * rmarkdown_2.1 
  * scales_1.1.0
  * tidyverse_1.3.0
  * vegan_2.5-6
  * venn_1.9
  
#### Running analysis
All the codes should be run from the project's root directory.

1.Download or clone the github repository to the project's root directory.
```bash
# clone the github repository
git clone https://github.com/yanxianl/Li_AqFl2-Microbiota_ASM_2020.git

# delete the following folders which would otherwise cause problems when running `04_qiime2_part2.bash`
rm -rf data/qiime2/core-metrics-results/ data/qiime2/robust-Aitchison-pca/
```
2.Download raw sequence data and reference database/phylogenetic tree for the analysis.
```bash
bash code/00_setup.bash
```
3.Sequence denoising by dada2.
```bash
Rscript -e "rmarkdown::render('code/01_dada2.Rmd')"
```
4.Taxonomic assignment in qiime2.
```bash
bash code/02_qiime2_part1.bash
```
5.Filter the feature table to remove: 1).chloroplast/mitochondria sequences and those without a phylum-level taxonomic annotation;
2).low-prevalence features that only present in one sample; 3).contaminating features.
```bash
Rscript -e "rmarkdown::render('code/03_preprocessing.Rmd')"
```
6.Phylogeny and core-metrics-results in qiime2.
```bash
bash code/04_qiime2_part2.bash
```
7.Export qiime2 outputs into R.
```bash
Rscript -e "rmarkdown::render('code/05_qiime2R.Rmd')"
```
8.Taxonomic analysis.
```bash
Rscript -e "rmarkdown::render('code/06_taxonomy.Rmd')"
```
9.Alpha-diversity visualization and statistical analysis.
```bash
Rscript -e "rmarkdown::render('code/07_alpha-diversity.Rmd')"
```
9.Beta-diversity visualization and statistical analysis.
```bash
Rscript -e "rmarkdown::render('code/08_beta-diversity.Rmd')"
```
10.Association testing between microbial clades and sample metadata.
```bash
Rscript -e "rmarkdown::render('code/09_metadata_association.Rmd')"
```
### To-do list
* Add a driver script to automate all the analysis, e.g., `make ` or `snakemake`.
* Make lightweight RMarkdown files ([03, 05-09]_\*.Rmd) binder-ready, which can be directly run online by clicking the `launch binder` badge located at the top of this README file.
### Acknowledgements
The initial file and directory structure of this project is based on the [template](https://github.com/SchlossLab/new_project/releases/latest) shared by [Dr. Pat Schloss](http://www.schlosslab.org/) to improve the reproducibility of microbiome data analysis. For trainings and tutorials on reproducible data analysis in microbiome research, check the [*Riffomonas*](http://www.riffomonas.org/) project.
