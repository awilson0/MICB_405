## Scripts for MICB 405 Project 2
## JJ Hum

library(tidyverse)
library(pathview)
setwd("~/micb_405/MICB_405/")
checkm <- read_tsv("data/MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv")
arc <- read_tsv("data/gtdbtk.ar122.classification_pplacer.tsv", col_names = FALSE)
bac <- read_tsv("data/gtdbtk.bac120.classification_pplacer.tsv", col_names = FALSE)
rpkm <- read_csv("data/SaanichInlet_10m_binned.rpkm.csv")
rpkm2 <- read_csv("data/SI072_10m_MAG_RPKM.csv")

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
  dplyr::select(-X2) %>% 
  dplyr::rename(`Bin Id` = "X1") %>% 
  mutate_if(is_character, funs(na_if(.,"")))

# Sum rpkms and add column matching `Bin id`
r2 <- rpkm2 %>% separate(Sequence_name, into = c("Bin", "num"), sep="(?<=SaanichInlet_10m_...)") %>% 
  group_by(Bin) %>% 
  summarize(RPKM =  sum(RPKM)) %>% 
  mutate(Bin = gsub('m_', 'm.', Bin)) %>% 
  dplyr::rename(`Bin Id` = Bin)

# Combine all data frames for plotting
dat <- full_join(r2, tax, by = "Bin Id") %>% full_join(checkm, by = "Bin Id") %>% 
  filter(!is.na(Domain))

# Find most represented phyla
top_phyla <- as.vector(table(tax$Phylum))
names(top_phyla) <- names(table(tax$Phylum))
top_phyla <- sort(top_phyla, decreasing = TRUE)
top_phyla <- names(top_phyla[1:10])
top_phyla

dat$Phylum <- sapply(dat$Phylum, FUN = function(p){
  if(p %in% top_phyla){
    p
  } else {
    "Other"
  }
})

dat$Phylum <- factor(dat$Phylum, levels = c("Actinobacteriota", "Bacteroidota", "Chloroflexota",
                                            "Dependentiae", "Marinisomatota", "Nanoarchaeota",
                                            "Patescibacteria", "Proteobacteria", "Thermoplasmatota",
                                            "Verrucomicrobiota", "Other"))


contamination_plot <- ggplot(dat, aes(x = Completeness, y = Contamination)) +
  geom_point(aes(colour = Phylum, size = RPKM, shape = Domain))+
  scale_shape_manual(values = c(17, 16)) + 
  geom_segment(aes(x=90,y=1.25,xend=100,yend=1.25)) +
  geom_segment(aes(x=90,y=0,xend=90,yend=1.25))+
  geom_segment(aes(x=100, y=0, xend = 100, yend = 1.25)) +
  geom_segment(aes(x = 90, y=0, xend = 100, yend = 0))
contamination_plot

