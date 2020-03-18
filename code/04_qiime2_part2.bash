#!/bin/bash

#======================================================================================================================
# Processing marker-gene data in QIIME2, part2
#======================================================================================================================
# Activate QIIME2-2020.2 
source $HOME/miniconda3/bin/activate
conda activate qiime2-2020.2

# Import the filtered feature table
qiime tools import \
  --input-path data/preprocessing/table-filtered.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path data/qiime2/table-filtered.qza

# Feature table summary 
qiime feature-table summarize \
  --i-table data/qiime2/table-filtered.qza \
  --m-sample-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/table-filtered.qzv 

# Taxonomic barplot 
qiime taxa barplot \
  --i-table data/qiime2/table-filtered.qza \
  --i-taxonomy data/qiime2/taxonomy-silva132.qza \
  --m-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/taxa-bar-plots-filtered.qzv

#======================================================================================================================
# Phylogeny
#======================================================================================================================
# Filter rep-seqs based on the filtered feature table
qiime feature-table filter-seqs \
  --i-data data/qiime2/rep-seqs.qza \
  --i-table data/qiime2/table-filtered.qza \
  --p-no-exclude-ids \
  --o-filtered-data data/qiime2/rep-seqs-filtered.qza

# Reference-based fragment insertion with SEPP
qiime fragment-insertion sepp \
  --i-representative-sequences data/qiime2/rep-seqs-filtered.qza \
  --i-reference-database data/reference/sepp-refs-silva-128.qza \
  --o-tree data/qiime2/insertion-tree.qza \
  --o-placements data/qiime2/tree-placements.qza \
  --p-threads 16 \
  --p-debug

# Filter the uninserted sequences from the feature table
qiime fragment-insertion filter-features \
  --i-table data/qiime2/table-filtered.qza \
  --i-tree data/qiime2/insertion-tree.qza \
  --o-filtered-table data/qiime2/table-filtered-sepp-inserted.qza \
  --o-removed-table data/qiime2/table-filtered-sepp-uninserted.qza \
  --verbose

# Feature table summary 
qiime feature-table summarize \
  --i-table data/qiime2/table-filtered-sepp-inserted.qza \
  --m-sample-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/table-filtered-sepp-inserted.qzv 

qiime feature-table summarize \
  --i-table data/qiime2/table-filtered-sepp-uninserted.qza \
  --m-sample-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/table-filtered-sepp-uninserted.qzv 

#======================================================================================================================
# Quality control: taxonomic composition of the mock
#======================================================================================================================
# Import the expected taxonomic composition of the mock ===============================================================
biom convert \
  -i data/reference/mock_expected.tsv \
  -o data/reference/mock-expected.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path data/reference/mock-expected.biom \
  --type 'FeatureTable[RelativeFrequency]' \
  --input-format BIOMV210Format \
  --output-path data/qiime2/quality-control/mock-expected.qza

# Get the observed taxonomic composition of the mock samples ==========================================================
# Filter feature table to contain the mock samples only 
qiime feature-table filter-samples \
  --i-table data/qiime2/table-filtered.qza \
  --m-metadata-file data/metadata.tsv \
  --p-where "SampleType='Mock'" \
  --p-no-exclude-ids \
  --o-filtered-table data/qiime2/quality-control/mock-observed.qza

# Inspect the taxonomic composition of mock samples
qiime taxa barplot \
  --i-table data/qiime2/quality-control/mock-observed.qza \
  --i-taxonomy data/qiime2/taxonomy-silva132.qza \
  --m-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/quality-control/mock-observed.qzv

# Agglomerate taxa at species level
qiime taxa collapse \
  --i-table data/qiime2/quality-control/mock-observed.qza \
  --i-taxonomy data/qiime2/taxonomy-silva132.qza \
  --p-level 7 \
  --o-collapsed-table data/qiime2/quality-control/mock-observed-l7.qza

# Convert sequence counts into relative abundances
qiime feature-table relative-frequency \
  --i-table data/qiime2/quality-control/mock-observed-l7.qza \
  --o-relative-frequency-table data/qiime2/quality-control/mock-observed-l7-rel.qza

# Compare observed taxonomic composition to the expected values =======================================================
qiime quality-control evaluate-composition \
  --i-expected-features data/qiime2/quality-control/mock-expected.qza \
  --i-observed-features data/qiime2/quality-control/mock-observed-l7-rel.qza \
  --o-visualization data/qiime2/quality-control/mock-comparison.qzv

#======================================================================================================================
# Alpha and beta diversity analysis 
#======================================================================================================================
# Exclude mock and negative controls from downstream analysis
qiime feature-table filter-samples \
  --i-table data/qiime2/table-filtered-sepp-inserted.qza \
  --m-metadata-file data/metadata.tsv \
  --p-where "SampleType IN ('Mock', 'Blank-extraction', 'Blank-library')" \
  --p-exclude-ids \
  --o-filtered-table data/qiime2/table-filtered-sepp-inserted-no-control.qza

# Rarefaction analysis 
qiime diversity alpha-rarefaction \
  --i-table data/qiime2/table-filtered-sepp-inserted-no-control.qza \
  --i-phylogeny data/qiime2/insertion-tree.qza \
  --p-max-depth 10000 \
  --p-steps 20 \
  --m-metadata-file data/metadata.tsv \
  --o-visualization data/qiime2/alpha-rarefaction.qzv

# Generating core metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny data/qiime2/insertion-tree.qza \
  --i-table data/qiime2/table-filtered-sepp-inserted-no-control.qza \
  --m-metadata-file data/metadata.tsv \
  --p-sampling-depth 2047 \
  --output-dir data/qiime2/core-metrics-results

# Comparing beta-diversity using robust Aitchison PCA 
qiime deicode rpca \
  --i-table data/qiime2/table-filtered-sepp-inserted-no-control.qza \
  --p-min-feature-count 10 \
  --p-min-sample-count 500 \
  --output-dir data/qiime2/robust-Aitchison-pca

qiime emperor biplot \
  --i-biplot data/qiime2/robust-Aitchison-pca/biplot.qza \
  --m-sample-metadata-file data/metadata.tsv \
  --m-feature-metadata-file data/qiime2/taxonomy-silva132.qza \
  --o-visualization data/qiime2/robust-Aitchison-pca/biplot.qzv \
  --p-number-of-features 8

# exit qiime2-2020.2 
conda deactivate
