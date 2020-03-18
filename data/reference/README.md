This folder holds all the reference data used for the analyses.

Run the following codes to download the reference sequence and taxonomy of SILVA132.
```bash
# Download the reference database, SILVA132
wget -P data/reference https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip

# Decompress and delete the downloaded zip file 
unzip data/reference/Silva_132_release.zip -d data/reference/silva_132 && rm -f data/reference/Silva_132_release.zip

# Copy and rename the reference sequence and taxonomy 
cp data/reference/silva_132/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna data/reference
cp data/reference/silva_132/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt data/reference
mv data/reference/consensus_taxonomy_7_levels.txt data/reference/silva_132_consensus_taxonomy_l7.txt

# Delete the folder
rm -rf data/reference/silva_132
```

Run the following code to download the SILVA128 reference tree for phylogenetic placement of amplicon sequence variants (ASVs).
```bash
# Download the reference database, SILVA_132
wget -P data/reference https://data.qiime2.org/2020.2/common/sepp-refs-silva-128.qza
```