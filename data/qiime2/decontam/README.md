This folder holds data and visualization files relevant to the screening of contaminating features in the data.

1. Data output files
    *`table-no-chlo-mito-lowPre-with-phyla-tax-rel.tsv`: total sum scaled feature table for the contaminating features screening, which has been filtered to exclude 1)chloroplast/mitochondria sequences and those without a phylum-level taxonomic annotation; 2)low-prevalence features that only present in one sample.
    *`contaminants-sample.tsv`: feature ID and taxonomy of contaminating sequences to be removed from the biological samples.
    *`contaminants-mock.tsv`: feature ID and taxonomy of contaminating sequences (including cross-contaminants) to be removed from the mock samples.

2. Visualization files
    *`feature_prevalence.pdf`: barplots showing the prevalence and abundance of features found in the negative controls.
    *`correlation.pdf`: correlation plots showing the correlation between sample bacterial DNA concentration and the abundance of features found in the negative controls.