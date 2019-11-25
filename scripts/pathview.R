library(tidyr)
library(dplyr)
library(pathview)
library(RColorBrewer)
library(knitr)

data_path = "../data/files_for_pathview/"

# metabolic reconstruction
ko <- read.table(paste(data_path,"SaanichInlet_MAGx_ORFs_ko.cleaned.txt",sep=""), sep="\t") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)
# combine rpkm result for multiple cruise, result from map bwa -p??
metat_rpkm <- read.table(paste(data_path,"SI072_SaanichInlet_MAG_ORFs.sam_RPKM.csv",sep=""), sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm = V2)

prokka_mag_map <- read.table(paste(data_path,"Prokka_MAG_map.csv",sep=""), header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)
# todo inconsistent mag name, change . to _ use gsub()?
arc_class <- read.table(paste(data_path,"gtdbtk.ar122.classification_pplacer.tsv",sep=""), sep="\t")
bac_class <- read.table(paste(data_path,"gtdbtk.bac120.classification_pplacer.tsv",sep=""), sep="\t")
gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  mutate(mag = gsub('m.', 'm_', mag))

checkm_dat <- read.table(paste(data_path,"MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv",sep=""),
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination) %>% 
  mutate(mag = gsub('m.', 'm_', mag))

# change header name!!
# Due to a bug in the renaming script we have to rename the bins. Its a bit hacky but works using tidyverse functions
metag_rpkm <- read.table(paste(data_path,"SaanichInlet_10m_binned.rpkm.csv",sep=""), header=T, sep=',') %>% 
  #mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  #mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence_name, into=c("mag1","mag2","mag3","contig"), sep='_', extra="merge") %>% 
  unite('mag',c('mag1','mag2','mag3')) %>% 
  group_by(Sample_name, mag) %>% 
  summarise(g_rpkm = sum(RPKM))
  #mutate(mag = gsub('Inlet.', 'Inlet_', mag))

# determine the number of Phyla present in our bins, and how many MAGs are representing each
gtdb_dat %>% 
  group_by(Phylum) %>% 
  summarise(count = n_distinct(mag)) %>% 
  kable()
##?
gtdb_dat <- dplyr::select(gtdb_dat, mag, Kingdom, Phylum, Class, Order, Family)

rpkm_dat <- left_join(ko, metat_rpkm, by="orf") %>%
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  left_join(checkm_dat, by="mag")

# If you also wanted to add the RPKM abundance values from the metagenome:
# left_join(metag_rpkm, by="mag")

head(rpkm_dat) %>% kable()

# todo
# aggregation

# # Subset by taxon
# ko_rpkm <- rpkm_dat %>%
#   filter(Phylum %in% c("p__Proteobacteria", "p__Nanoarchaeota", "p__Thermoplasmatota")) %>%
#   group_by(mag, ko) %>% 
#   summarise(t_rpkm = sum(rpkm)) %>% 
#   spread(key = mag, value = t_rpkm)
# 
# Subset by completeness and contamination
ko_rpkm <- rpkm_dat %>%
  filter(Completeness >= 90 & Contamination < 5) %>%
  group_by(mag, ko) %>%
  summarise(t_rpkm = sum(rpkm)) %>%
  spread(key = mag, value = t_rpkm)

# # Aggregate by a taxonomy, still summing RPKM of each KO number. You could use mean() instead.
# ko_rpkm <- rpkm_dat %>%
#   group_by(Class, ko) %>% 
#   summarise(t_rpkm = sum(rpkm)) %>% 
#   spread(key = Class, value = t_rpkm)
# 
pv_mat <- dplyr::select(ko_rpkm, -ko)
rownames(pv_mat) <- ko_rpkm$ko

# Nitrogen metabolism
# todo
# need to change kegg.dir path
pv.out <- pathview(gene.data = pv_mat,
                   limit = list(gene = c(0,10)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"),
                   species = "ko",
                   pathway.id="01200",
                   kegg.dir = "../data")

