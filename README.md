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
│   ├── 03_preprocessing.Rmd          # feature table filtering and identification of contaminants    
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
* [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) is located in the HOME directory.
* [grabseqs (0.7.0)](https://github.com/louiejtaylor/grabseqs) should be installed via conda.
* [QIIME2 (2020.2)](https://docs.qiime2.org/2020.2/) should be installed within a conda environment (qiime2-2020.2).
* R (3.6.0) should be located in the user's PATH.
* R packages: refer to the `SessionInfo` in the html files rendered from the RMarkdown files.
  
#### Running analysis
Download or clone the github repository to a local working directory.
```
git clone https://github.com/yanxianl/AquaFly-SeawaterGutHealth-Aquaculture-2019.git
```
Run codes in the `code/` folder from the project's root directory. The codes are numbered by the order of execution.

### To-do list
* Add a driver script to automate all the analysis, e.g., `make ` or `snakemake`.

### Acknowledgements
The initial file and directory structure of this project is based on the [template](https://github.com/SchlossLab/new_project/releases/latest) shared by [Dr. Pat Schloss](http://www.schlosslab.org/) to improve the reproducibility of microbiome data analysis. For trainings and tutorials on reproducible data analysis in microbiome research, check the [*Riffomonas*](http://www.riffomonas.org/) project.
