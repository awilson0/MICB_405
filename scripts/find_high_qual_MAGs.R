library(tidyverse)

#load Medium and High quality MAGs data
raw_dat <- read_tsv(file = 'MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv')

# Find High Quality MAGs (Completeness > 90% and Contamination < 5%)
dat <- filter(raw_dat, Completeness > 90, Contamination < 5)