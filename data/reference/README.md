This folder holds all the reference data used for the analyses.

Run the following codes to download the reference sequence and taxonomy of SILVA132.
```bash
# Download the reference database, SILVA_132
wget https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip

# Decompress and delete the downloaded zip file 
unzip Silva_132_release.zip -d silva_132 && rm -f Silva_132_release.zip

# Copy and rename the reference sequence and taxonomy 
cp silva_132/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna .
cp silva_132/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt .
mv consensus_taxonomy_7_levels.txt silva_132_consensus_taxonomy_l7.txt

# Delete the folder
rm -rf silva_132
```
