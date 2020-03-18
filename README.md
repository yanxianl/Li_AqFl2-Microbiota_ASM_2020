## Differential response of digesta- and mucosa-associated intestinal microbiota to dietary black soldier fly (*Hermetia illucens*) larvae meal in seawater phase Atlantic salmon (*Salmo salar*)

Knowledge on the nutritional value of black soldier fly (Hermetia illucens) as a feed ingredient for fish is accumulating, yet less is known regarding how it may modulate fish intestinal microbiota. In a 16-week seawater feeding trial, Atlantic salmon (Salmo salar) were fed either a commercially-relevant reference diet, or an insect meal diet wherein black soldier fly larvae meal was included to replace all the fish meal and most of the pea protein concentrate in the reference diet. The digesta- and mucosa-associated distal intestinal microbiota were profiled by sequencing V1-2 regions of the 16S rRNA gene. Regardless of diet, we observed substantial differences between digesta- and mucosa-associated intestinal microbiota. Microbial richness and diversity were much higher in the digesta than the mucosa. The insect meal diet altered the distal intestinal microbiota assemblage resulting in higher microbial richness and diversity. The diet effect was dependent on the sample origin, with digesta-associated intestinal microbiota showing more pronounced changes than the mucosa-associated microbiota. Multivariate association analysis identified two mucosa-enriched taxa, Brevinema andersonii and unclassified Spirochaetaceae, associated with the expression of genes related to immune responses and barrier function in the distal intestine, respectively. The insect meal diet increased the relative abundance of several taxa that have been reported in fish, including Actinomyces, Bacillus, Brevibacterium, Corynebacterium 1 and Enterococcus. In conclusion, the response of digesta- and mucosa-associated intestinal microbiota to dietary inclusion of insect meal was different, with the latter being more resilient to the dietary change.

### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) is located in the HOME directory.
* [grabseqs (0.7.0)](https://github.com/louiejtaylor/grabseqs) should be installed via conda.
* [QIIME2 (2020.2)](https://docs.qiime2.org/2020.2/) should be installed within a conda environment (qiime2-2020.2).
* R (3.6.3) should be located in the user's PATH.
* R packages: refer to the `SessionInfo`contained in the html files rendered from the RMarkdown files.
  
#### Running analysis
All analysis should be run from the project's root directory.

To begin with, run the following bash command to download


```
git clone https://github.com/SchlossLab/LastName_BriefDescription_Journal_Year.git
make write.paper
```
### To-do list
* Add a driver script to automate all the analysis, e.g., `make `or `snakemake`.
* 
### Acknowledgements
The initial file and directory structure of this project was developed by a group of participants in the Reproducible Science Curriculum Workshop, held at [NESCent] in December 2014 ([rr-init repository]). The structure is based on, and heavily follows the one proposed by [Noble 2009], with a few but small modifications. All copyright and related and neighboring rights to the original template were dedicated to the public domain worldwide under the [CC0 Public Domain Dedication]. The template and its derivatives are distributed without any warranty. It has been further modified by Pat Schloss to fit the needs of his research group.
