library(here) # [CRAN] A replacement for 'file.path()', locating the files relative to the project root
library(tidyverse) # [CRAN] # Easily Install and Load the 'Tidyverse'

# get the absolute file paths
files <- list.files(here("data", "raw", "casava-18-paired-end-demultiplexed"), full.names = TRUE)

# read in sra metadata
mdat <- read_csv(here("data", "raw", "casava-18-paired-end-demultiplexed", "metadata_sra.csv"))

# make a lookup table for renaming fastq files
filenames <- mdat %>%
  select(Run, SampleName) %>%
  uncount(2)  %>%
  mutate(index = rep(c(1, 2), nrow(.)/2),
         Run = paste0(dirname(files[1]), "/", Run, "_", index, ".fastq.gz"),
         SampleName = paste0(dirname(files[1]), "/", SampleName, "_R", index, "_001.fastq.gz"))

# rename fastq files
file.rename(from = filenames[["Run"]], to = filenames[["SampleName"]])

