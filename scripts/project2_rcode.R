## Scripts for MICB 405 Project 2
## JJ Hum, Andrew Wilson

library(tidyverse)
library(pathview)
library(cowplot)
library(KEGGREST)
library(knitr)

setwd("~/micb_405/MICB_405/")
checkm <- read_tsv("data/MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv")
arc <- read_tsv("data/gtdbtk.ar122.classification_pplacer.tsv", col_names = FALSE)
bac <- read_tsv("data/gtdbtk.bac120.classification_pplacer.tsv", col_names = FALSE)
rpkm_g <- read_csv("data/SaanichInlet_10m_binned.rpkm.csv")
rpkm_t <- read_csv("data/files_for_pathview/SI072_SaanichInlet_MAG_ORFs.sam_RPKM.csv", col_names = FALSE)
prk <- read_csv("data/Prokka_MAG_map.csv", col_names = FALSE)
ko <- read_tsv("data/files_for_pathview/SaanichInlet_MAGx_ORFs_ko.cleaned.txt", col_names = FALSE)

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
r2 <- rpkm_g %>% separate(Sequence_name, into = c("Bin", "num"), sep="(?<=SaanichInlet_10m_...)") %>% 
  group_by(Bin) %>% 
  summarize(RPKM =  sum(RPKM)) %>% 
  mutate(Bin = gsub('m_', 'm.', Bin)) %>% 
  dplyr::rename(`Bin Id` = Bin)

# Combine all data frames for plotting
dat <- full_join(r2, tax, by = "Bin Id") %>% full_join(checkm, by = "Bin Id") %>% 
  filter(!is.na(Domain))

## Andrew's code
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
# Assign nice colours to each phylum
phylum_colours <- c("#9c47d5","#a75e3d","#697d36","#c08426","#725ec2","#d6482e","#6778b6","#be4964","#5ea536","#b44d95","#3c9472")

color1 <- c("#00cf64","#e60096","#9688ff","#ff8f0c","#7ad0ff","#a46c00","#46044d","#76d9a8","#650017","#ffab7b","#ff949b")

# Plot Contamination vs Completion
contamination_plot <- ggplot(dat, aes(x = Completeness, y = Contamination)) +
  geom_point(aes(colour = Phylum, size = RPKM, shape = Domain))+
  scale_shape_manual(values = c(17, 16)) + 
  scale_colour_manual(values = color1) +
  geom_hline(yintercept = 0) + geom_hline(yintercept = 5) +
  geom_vline(xintercept = 90) + geom_vline(xintercept = 100)+
  geom_segment(aes(x=90,y=5,xend=100,yend=5), color = "blue") +
  geom_segment(aes(x=90,y=0,xend=90,yend=5), color = "blue")+
  geom_segment(aes(x=100, y=0, xend = 100, yend = 5), color = "blue") +
  geom_segment(aes(x = 90, y=0, xend = 100, yend = 0), color = "blue") +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 4)))

contamination_plot
## End Andrew's code

# High quality contamination vs. completion plot
dat1 <- dat %>% filter(Completeness > 90) %>% filter(Contamination < 5) %>% 
  # Put Phylum name back after it was overwritten with other for the larger plot
  mutate(Phylum = gsub("Other", "Desulfobacterota", Phylum))

zoomplot <- ggplot(dat1, aes(x = Completeness, y = Contamination)) +
  geom_point(aes(colour = Phylum, size = RPKM))+
  scale_shape_manual(values = c(17, 16)) + 
  scale_colour_manual(values = color1) +
  theme_cowplot() +
  guides(color = guide_legend(override.aes = list(size = 4)))
zoomplot

# Counting total MAGs with RPKM + taxonomy data
dat2 <- dat %>% filter(!is.na(RPKM))
table(dat1$Phylum)
table(dat1$Class)

## Prokka data wrangling
ko1 <- ko %>% dplyr::rename("ORF" = "X1", "Kegg_ID" = "X2")
prk1 <- prk %>% mutate(X2 = gsub('m_', 'm.', X2)) %>% dplyr::rename("Prokka_Id" = "X1", `Bin Id` = "X2")
rpkm_t2 <- rpkm_t %>% dplyr::rename("ORF" = "X1", "RPKM" = "X2")

# Join all data frames
all <- left_join(ko1, rpkm_t2, by = "ORF") %>% 
  separate(ORF, into = c("Prokka_Id", "ORF")) %>% 
  left_join(prk1, by = "Prokka_Id") %>% 
  left_join(tax, by = "Bin Id") %>% 
  left_join(checkm, by = "Bin Id") %>% 
  dplyr::select(-"Marker lineage", -"# genomes", -"# markers", -"# marker sets", 
         -"0", -"1", -"2", -"3", -"4", -"5+", -"Strain heterogeneity")

# Get KEGGREST lookups and reformat into tables with proper headings
zyms <-keggLink("enzyme","ko")
zym1 <- keggList("enzyme")
paths <-keggLink("enzyme", "pathway")

enzymes <- enframe(zyms) %>% 
  dplyr::rename("Kegg_ID" = "name", "num" = "value") %>% 
  mutate(Kegg_ID = gsub('ko:K', 'K', Kegg_ID)) 
zym2 <- enframe(zym1) %>% 
  dplyr::rename("num" = "name", "Enzyme" = "value") %>% 
  separate(Enzyme, sep = ";", into = c("Enzyme", "Enzyme2"), extra = "merge")
pathways <- enframe(paths) %>% dplyr::rename("Path" = "name", "num" = "value")

# Get gene names by uploading KO list to KEGG mapper
# Copy pasted list of genes "names.csv"
names <- read_csv("data/names.csv", col_names = FALSE) %>% 
  mutate(X1 = gsub("ko:", "", X1)) %>% 
  separate(X1, into = c("Kegg_ID", "Gene"), sep = "(?<=K.....)") %>% 
  separate(Gene, into = c("Gene", "Enzyme2"), sep = ";", extra = "merge")

# Combine tables
alldat <- left_join(enzymes, zym2, by = "num") %>% 
  left_join(pathways, by = "num") %>% 
  left_join(all, by = "Kegg_ID") %>% 
  left_join(names, by = "Kegg_ID") %>% 
  filter(RPKM != 0)

# Most abundant pathways?
toppath <- alldat %>% 
  group_by(Path) %>% 
  summarize(RPKM =  sum(RPKM)) %>% 
  arrange(desc(RPKM)) 

# Looking at C fixation pathway
cfix <- alldat %>% filter(Path == "path:map00720")

cfplot <- 
  ggplot(cfix, aes(y = Gene, x = Class, color = Phylum)) +
  geom_point(aes(size = RPKM)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.major.y = element_line(color = "gray")) +
  guides(color = guide_legend(override.aes = list(size = 4)))
cfplot

topcf <- cfix %>% 
  group_by(Class) %>% 
  summarize(RPKM =  sum(RPKM))

# Looking at Calvin cycle pathway
cph <- alldat %>% filter(Path == "path:map00710")

cphplot <- 
  ggplot(cph, aes(y = Gene, x = Class, color = Phylum)) +
  geom_point(aes(size = RPKM)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        panel.grid.major.y = element_line(color = "gray")) +
  guides(color = guide_legend(override.aes = list(size = 4)))

cphplot

topcph <- cph %>% 
  group_by(Class) %>% 
  summarize(RPKM =  sum(RPKM))
  

## Pathview visualization

ko_rpkm <- alldat %>%
  group_by(Class, Kegg_ID) %>% 
   summarise(t_rpkm = sum(RPKM)) %>% 
  spread(key = Class, value = t_rpkm)

pv_mat <- dplyr::select(ko_rpkm, -Kegg_ID)
rownames(pv_mat) <- ko_rpkm$Kegg_ID

# Carbon fixation by prokaryotes pathview
pv.cfix <- pathview(gene.data = pv_mat,
                   limit = list(gene = c(0,10)),
                   low = list(gene = "#91bfdb"),
                   mid = list(gene = "#ffffbf"),
                   high = list(gene = "#fc8d59"),
                   species = "ko",
                   pathway.id="00720")

# Looking at interesting spots
cf1 <- alldat %>% filter(num == "ec:1.2.7.4") %>% distinct(Class, .keep_all = TRUE)
cf2 <- alldat %>% filter(num == "ec:2.3.1.169") %>% distinct(Class, .keep_all = TRUE)

# Looking for epsilon-proteobacteria
table(dat$Class)

# Carbon fixation by photosynthetes pathview
pv.cpho <- pathview(gene.data = pv_mat,
                    limit = list(gene = c(0,10)),
                    low = list(gene = "#91bfdb"),
                    mid = list(gene = "#ffffbf"),
                    high = list(gene = "#fc8d59"),
                    species = "ko",
                    pathway.id="00710")

# Interesting enzymes
cp1 <- alldat %>% filter(num == "ec:2.7.1.19") %>% distinct(Class, .keep_all = TRUE)
cp2 <- alldat %>% filter(num == "ec:4.1.1.39") %>% distinct(Class, .keep_all = TRUE)

## Plotting geochemical data
geodata <- read_csv("data/Saanich_TimeSeries_Chemical_Data.txt")

g2 <- geodata %>% 
  select(Cruise, Date, Depth, Mean_O2, Mean_H2S, 
         Mean_co2, Mean_CH4, Mean_N2, Mean_NH4, Mean_NO2) %>% 
  rename("Mean CH4" = Mean_CH4, "Mean CO2" = Mean_co2, "Mean H2S" = Mean_H2S,
         "Mean N2" = Mean_N2, "Mean NH4" = Mean_NH4,
         "Mean NO2" = Mean_NO2, "Mean O2" = Mean_O2) %>% 
  gather(key="Chemical", value="Concentration", -Cruise, -Date, -Depth) 

geoplot <-
ggplot(g2, aes(x=Concentration, y=Depth)) +
  geom_point(aes(colour=Chemical)) +
  scale_y_reverse(limits=c(200,0)) +
  facet_wrap(~Chemical, scales="free_x", nrow = 2) +
  xlab("Concentration, uM") + ylab("Depth, m") +
  theme_cowplot() +
  theme(legend.position = "none", 
        panel.grid.major = element_line(color = "gray"),
        panel.spacing = unit(1.5, "lines"))
geoplot