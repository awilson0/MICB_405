library(tidyverse)

raw_dat <- read_tsv("MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv")

dat <- raw_dat

ggplot(dat) 

dat %>%
  ggplot(aes(x=Contamination, y=Completeness)) + geom_point()