<!-- badges: start -->
  [![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/yanxianl/Li_AqFl2-Microbiota_ASM_2020/master?urlpath=rstudio)
<!-- badges: end -->

## Differential response of digesta- and mucosa-associated intestinal microbiota to dietary black soldier fly (*Hermetia illucens*) larvae meal in seawater phase Atlantic salmon (*Salmo salar*)

Knowledge on the nutritional value of black soldier fly (Hermetia illucens) as a feed ingredient for fish is accumulating, yet less is known regarding how it may modulate fish intestinal microbiota. In a 16-week seawater feeding trial, Atlantic salmon (Salmo salar) were fed either a commercially-relevant reference diet, or an insect meal diet wherein black soldier fly larvae meal was included to replace all the fish meal and most of the pea protein concentrate in the reference diet. The digesta- and mucosa-associated distal intestinal microbiota were profiled by sequencing V1-2 regions of the 16S rRNA gene. Regardless of diet, we observed substantial differences between digesta- and mucosa-associated intestinal microbiota. Microbial richness and diversity were much higher in the digesta than the mucosa. The insect meal diet altered the distal intestinal microbiota assemblage resulting in higher microbial richness and diversity. The diet effect was dependent on the sample origin, with digesta-associated intestinal microbiota showing more pronounced changes than the mucosa-associated microbiota. Multivariate association analysis identified two mucosa-enriched taxa, Brevinema andersonii and unclassified Spirochaetaceae, associated with the expression of genes related to immune responses and barrier function in the distal intestine, respectively. The insect meal diet increased the relative abundance of several taxa that have been reported in fish, including Actinomyces, Bacillus, Brevibacterium, Corynebacterium 1 and Enterococcus. In conclusion, the response of digesta- and mucosa-associated intestinal microbiota to dietary inclusion of insect meal was different, with the latter being more resilient to the dietary change.

### Overview

```
root
├── code
│   ├── 00_setup.bash
│   ├── 01_dada2.html
│   ├── 01_dada2.Rmd
│   ├── 02_qiime2_part1.bash
│   ├── 03_preprocessing.html
│   ├── 03_preprocessing.Rmd
│   ├── 04_qiime2_part2.bash
│   ├── 05_qiime2R.html
│   ├── 05_qiime2R.Rmd
│   ├── 06_taxonomy.html
│   ├── 06_taxonomy.Rmd
│   ├── 07_alpha-diversity.html
│   ├── 07_alpha-diversity.Rmd
│   ├── 08_beta-diversity.html
│   ├── 08_beta-diversity.Rmd
│   ├── 09_metadata_association.html
│   ├── 09_metadata_association.Rmd
│   ├── functions
│   ├── utilities
│   └── README.md
├── data
│   ├── metadata.tsv
│   ├── raw 
│   ├── reference 
│   ├── qPCR 
│   ├── dada2 
│   ├── qiime2 
│   ├── preprocessing
│   ├── qiime2R
│   ├── permanova
│   └── maaslin2
├── image
│   └── SBMIE.png
├── result
│   ├── figures
│   ├── README.md
│   └── tables
├── LICENSE.md
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
Download or clone the github repository to local working directory.
```
git clone https://github.com/yanxianl/AquaFly-SeawaterGutHealth-Aquaculture-2019.git
```
Run all the codes from the project's root directory. The codes are numbered by the order of execution

### To-do list
* Add a driver script to automate all the analysis, e.g., `make `or `snakemake`.

### Acknowledgements
The initial file and directory structure of this project was developed by a group of participants in the Reproducible Science Curriculum Workshop, held at [NESCent] in December 2014 ([rr-init repository]). The structure is based on, and heavily follows the one proposed by [Noble 2009], with a few but small modifications. All copyright and related and neighboring rights to the original template were dedicated to the public domain worldwide under the [CC0 Public Domain Dedication]. The template and its derivatives are distributed without any warranty. It has been further modified by Pat Schloss to fit the needs of his research group.
