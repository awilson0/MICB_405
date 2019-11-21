## Scripts for MICB 405 Project 2
## JJ Hum

#library(tidyverse)
#library(pathview)

checkm <- read_tsv("MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv")
arc <- read_tsv("gtdbtk.ar122.classification_pplacer.tsv", col_names = FALSE)
bac <- read_tsv("gtdbtk.bac120.classification_pplacer.tsv", col_names = FALSE)
rpkm <- read_csv("SaanichInlet_10m_binned.rpkm.csv")
rpkm2 <- read_csv("SI072_10m_MAG_RPKM.csv")

# Separate out taxonomy info, combine into 1 dataframe
a1 <- arc$X2
a1 <- str_remove_all(a1, ".__")
arc$X3 <- a1
arc1 <- arc %>% 
  separate(X3, into = c("Domain", "Phylum", "Class", 
                        "Order", "Family", "Genus", "Species"), sep=";")
b1 <- bac$X2
b1 <- str_remove_all(b1, ".__")
bac$X3 <- b1
bac1 <- bac %>% 
  separate(X3, into = c("Domain", "Phylum", "Class", 
                        "Order", "Family", "Genus", "Species"), sep=";")

tax <- bind_rows(arc1, bac1) %>% 
  select(-X2) %>% 
  rename(`Bin ID` = "X1") %>% 
  mutate_if(is_character, funs(na_if(.,"")))

# Sum rpkms and add column matching `Bin id`
r2 <- rpkm %>% separate(Sequence_name, into = c("Bin", "num"), sep="(?<=SaanichInlet_10m_...)") %>% 
  group_by(Bin) %>% 
  summarize(RPKM =  sum(RPKM)) 


