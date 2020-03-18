This folder holds the raw data to reproduce the results presented in the paper. 

The qPCR assays files can be viewed and edited by the [Roche LightCycler 96 System](https://lifescience.roche.com/en_no/brands/realtime-pcr-overview.html#software).

Run the following code to download the raw sequence data deposited at the [NCBI SRA database](https://www.ncbi.nlm.nih.gov/bioproject/) under the BioProject PRJNA555355.  
```bash
# download sequence and metadata from SRA
grabseqs sra -m metadata_sra.csv -o data/raw/casava-18-paired-end-demultiplexed/ -r 3 PRJNA555355

# rename files

```