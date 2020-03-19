#!/bin/bash

#==================================================================================================
# Download the raw sequence data and the associated metadata from NCBI SRA database               
#==================================================================================================
# Activate miniconda 
source $HOME/miniconda3/bin/activate
# Download raw sequence data using grabseqs tool
grabseqs sra -m metadata_sra.csv -o data/raw/casava-18-paired-end-demultiplexed/ -r 3 PRJNA555355
# Rename the downloaded fastq files
Rscript code/utilities/rename_sra_fastq.R
#==================================================================================================
# Download reference sequences and taxonomy               
#==================================================================================================
# Download the reference database, SILVA132
wget -P data/reference https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip
# Decompress and delete the downloaded zip file 
unzip data/reference/Silva_132_release.zip -d data/reference/silva_132 && rm -f data/reference/Silva_132_release.zip
# Copy and rename the reference sequence and taxonomy file
cp data/reference/silva_132/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna data/reference
cp data/reference/silva_132/SILVA_132_QIIME_release/taxonomy/16S_only/99/consensus_taxonomy_7_levels.txt data/reference
mv data/reference/consensus_taxonomy_7_levels.txt data/reference/silva_132_consensus_taxonomy_l7.txt
# Delete data to free up disk space
rm -rf data/reference/silva_132
#==================================================================================================
# Download the SILVA128 reference tree for phylogenetic placement of ASVs           
#==================================================================================================
wget -P data/reference https://data.qiime2.org/2020.2/common/sepp-refs-silva-128.qza
