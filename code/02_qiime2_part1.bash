#!/bin/bash

#======================================================================================================================
# Processing marker-gene data in QIIME2, part1
#======================================================================================================================
# Activate QIIME2-2020.2 
source $HOME/miniconda3/bin/activate
conda activate qiime2-2020.2

#======================================================================================================================
# Import feature table and representative sequences from dada2
#======================================================================================================================
# Import feature table
qiime tools import \
  --input-path data/dada2/table.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path data/qiime2/table.qza

# Import representative sequences
qiime tools import \
  --input-path data/dada2/rep-seqs.fna \
  --type 'FeatureData[Sequence]' \
  --output-path data/qiime2/rep-seqs.qza

# Check the imported feature table and representative sequences
qiime feature-table summarize \
  --i-table data/qiime2/table.qza \
  --m-sample-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/table.qzv 

qiime feature-table tabulate-seqs \
  --i-data data/qiime2/rep-seqs.qza \
  --o-visualization data/qiime2/rep-seqs.qzv

#======================================================================================================================
# Taxanomic assignment
#======================================================================================================================
# Import reference sequence and taxonomy to train the feature-classifier
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path data/reference/silva_132_99_16S.fna \
  --output-path data/qiime2/99-otus-silva132.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path data/reference/silva_132_consensus_taxonomy_l7.txt \
  --output-path data/qiime2/ref-taxonomy-silva132.qza

qiime feature-classifier extract-reads \
  --i-sequences data/qiime2/99-otus-silva132.qza \
  --p-f-primer AGAGTTTGATCMTGGCTCAG \
  --p-r-primer GCWGCCWCCCGTAGGWGT \
  --p-n-jobs 16 \
  --o-reads data/qiime2/ref-seqs-silva132.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads data/qiime2/ref-seqs-silva132.qza \
  --i-reference-taxonomy data/qiime2/ref-taxonomy-silva132.qza \
  --o-classifier data/qiime2/silva132-99otu-27-338-classifier.qza

# Assign taxanomy using the trained featureClassifier
qiime feature-classifier classify-sklearn \
  --i-classifier data/qiime2/silva132-99otu-27-338-classifier.qza \
  --i-reads data/qiime2/rep-seqs.qza \
  --o-classification data/qiime2/taxonomy-silva132.qza

qiime metadata tabulate \
  --m-input-file data/qiime2/taxonomy-silva132.qza \
  --o-visualization data/qiime2/taxonomy-silva132.qzv

# Taxonomic barplot 
qiime taxa barplot \
  --i-table data/qiime2/table.qza \
  --i-taxonomy data/qiime2/taxonomy-silva132.qza \
  --m-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/taxa-bar-plots.qzv

# exit qiime2-2020.2 
conda deactivate
