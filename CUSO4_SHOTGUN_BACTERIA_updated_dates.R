# setwd,load libraries, source functions ####
setwd('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Analysis/Bacteria_archaea/Figures')

# install.packages("devtools")
# devtools::install_github("vmikk/metagMisc")
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install(version = "3.23")
# # BiocManager::install("phyloseq")
#BiocManager::install("metagenomeSeq")
# BiocManager::install("ANCOMBC")
# BiocManager::install("maaslin3")
# BiocManager::install("microbiome")
# BiocManager::install("MicrobiotaProcess")
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install.packages("svglite")
#install.packages("writexl")
#install.packages("pals")
#install.packages("changepoint")

library(phyloseq); library (tidyverse); library(ggplot2);  library(stringr); 
library(dplyr);library(metagMisc); library(metagenomeSeq); library(vegan); library(cowplot);
library(ggdendro); library(pairwiseAdonis); library(randomcoloR); library(ggpubr); library(ppcor)
library(ggsignif); library (ANCOMBC);library(maaslin3); library (UpSetR); library(MicrobiotaProcess); library(microbiome)
library(ggtext); library(ggnewscale); library(rstatix); library(ggrepel); library(ggh4x); library(svglite);
library(lmerTest); library(mgcv); library(rmcorr); library(patchwork); library(colorspace)
library(writexl)
library(pals); library(changepoint); library(paletteer);library(RColorBrewer)


##Source functions
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbundanceOthersPercentage.R')
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/top_taxa_legend_updated.R")
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/fill_taxonomy_updated.R')
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbun_group_microbiome.R')
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/R_functions/MergeLowAbun_group_ARG.R")


#Importing data from kraken output nt_core - counts will be classified reads#### 
counts <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/Paired_end_mode_GTDB_updated_20260602/Conf_005/kraken_analytic_matrix.conf_005.csv')
#There are some extra samples that will not be used for this project. Filtering those.
dropping_samples <- c("H21_0912", "H21_1005", "P1_0420", "P1_0427", 
                      "P1_0504", "H21_1202b", "H21_1021a", "H21_1021b")
#Dropping them from the count matrix
counts <- counts %>%
  select(-all_of(dropping_samples)) %>% 
  filter(rowSums(select(., -taxa)) > 0)
ncol(counts) #235 = 219 samples + 3 mocks + 12 negative controls + "taxa" column

##Separating into taxonomy levels
counts_separated_tax <- counts %>%
  separate(taxa, 
           into = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|",#splits the strings by the "|" symbol.
           fill = "right") #fill = "right", missing components are added as "NA" to the right (last columns) instead of to the left 

##Extracting just taxonomy  (columns 1:8 are taxonomy, the rest are counts)
tax.table<- counts_separated_tax %>%
  dplyr::select(3:8) #GTDB does not give Domain or Kingdom
tax.table

##Filling up actual NAs and string "NA"s in the taxonomy table
filled_taxonomy <- fill_taxonomy(tax.table) ##apply the function to the taxonomy table
anyNA(filled_taxonomy) ##OK, no NAs now 
grep("^NA$", filled_taxonomy, value = T) ##OK, no "NA" strings now

##Now, to add the row names as "OTU1, OTU2, etc..." for phyloseq later on
filled_taxonomy_2<- filled_taxonomy %>%
  mutate(OTU = paste0("OTU", 1:nrow(filled_taxonomy))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##convert into matrix for phyloseq
filled_taxonomy_2

#Make a csv file for the kraken taxonomy table
write.csv(filled_taxonomy_2,
          "kraken_taxonomy.csv",
          row.names = F)

###OTU table #####
otu_table <- counts[, -1]%>% #Excludes the first column (taxonomy)
  mutate(OTU = paste0("OTU", 1:nrow(counts))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##make into matrix so it is compatible with otu_table function from phyloseq
otu_table

#IMPORT METADATA####
#This comes from an already clean metadata file with data for both systems, as well as positive and negative controls from the "Metadata_cleaning.R" script
metadata <- read_csv("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/metadata_all_systems_phyloseq.csv")

#Factor ordering
metadata <- metadata %>%
  # Convert Collection_Date to Date if it isn’t already
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Date = as.Date(Collection_Date))%>%
  mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
                                       levels = c("Low Copper Levels (Day 1-27)", 
                                                  "Transition Period 1 (Day 29-38)", 
                                                  "First Tx Copper Exposure (Day 39-51)", 
                                                  "Transition Period 2 (Day 52-81)", 
                                                  "Transition Period 3 (Day 86-108)", 
                                                  "Second Tx Copper Exposure (Day 109-135)", 
                                                  "Post-Treatment (Day 136-146)")))%>%
  mutate(Date_num_phase_naive_abbrv = factor(Date_num_phase_naive_abbrv, 
                                             levels = c("L",
                                                        "T1",
                                                        "E1",
                                                        "T2", 
                                                        "T3",
                                                        "E2", 
                                                        "P")))%>%
  mutate(Date_num_phase_established = factor(Date_num_phase_established, 
                                             levels = c("Low Copper Levels (Day 1-53)",
                                                        "Transition Period 1 (Day 54-65)",
                                                        "Tx Copper Exposure (Day 66-104)",
                                                        "Post-Treatment (Day 105-169)")))%>%
  mutate(Date_num_phase_established_abbrv = factor(Date_num_phase_established_abbrv, 
                                                   levels = c("L", 
                                                              "T1", 
                                                              "E", 
                                                              "P")))%>%
  mutate(Date_num_phase = factor(Date_num_phase, 
                                 levels = c(
                                   #Naive
                                   "Low Copper Levels (Day 1-27)", 
                                   "Transition Period 1 (Day 29-38)", 
                                   "First Tx Copper Exposure (Day 39-51)",
                                   "Transition Period 2 (Day 52-81)", 
                                   "Transition Period 3 (Day 86-108)", 
                                   "Second Tx Copper Exposure (Day 109-135)", 
                                   "Post-Treatment (Day 136-146)", 
                                   #Established
                                   "Low Copper Levels (Day 1-53)",
                                   "Transition Period 1 (Day 54-65)",
                                   "Tx Copper Exposure (Day 66-104)",
                                   "Post-Treatment (Day 105-169)")))%>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = c("L", "T1", "E", 
                                                  "E1", "T2", "T3", "E2", 
                                                  "P")))

##Making into phyloseq-compatible object
sampledata_phyloseq <- metadata %>%
  mutate(rows = SampleID)%>%
  column_to_rownames(var= "rows") %>%##Make sampleID column into row names, so they match sample_names() with OTU and TAX
  sample_data(metadata) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

#PHYLOSEQ####
#Make phyloseq object
OTU <-phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
TAX <-phyloseq::tax_table(filled_taxonomy_2)
phyloseq <- phyloseq(OTU, TAX, sampledata_phyloseq)

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU), metadata$SampleID) #No

#Are there samples in metadata that don't have sequencing data?
setdiff(metadata$SampleID, sample_names(OTU)) #Yes, "P1_1126", "P1_1203", 
#"P1_1216", "P1_1225", "P1_1228", "P1_0104", "P1_0112", "P1_0115", 
#"P1_0205", "P1_0212", "P1_0218", "H21_1021", "H21_1122b"

#COLOR PALETTES#####
enclosure.palette <- c("H21" = "#fc8d62",  
                       "P1"  = "#8da0cb" )

reads.palette <- c("Raw" ="#E69F00", 
                   "QC - Trimmomatic" = "#0072B2")


naive_phase_palette <- c("L" = "#A3A33D",
                         "T1" = "#0CE3A3",
                         "E1" = "#7EDC8A",
                         "T2" = "#066E10",
                         "T3" = "#9E2307B4",
                         # "T3" = "#E2762C",
                         "E2" = "#DE1B1BE6",
                         # "E2" = "#DE4A1B", 
                         "P" = "#88CCEE")

established_phase_palette <- c("L" = "#7F7F2A",
                               "T1" = "#B85114",
                               "E" = "#9E2007", 
                               "P" = "#0C8AE3")

phases_naive_established_palette <- c("L_Naive" = "#A3A33D",
                                      "T1_Naive" = "#0CE3A3",
                                      "E1_Naive" = "#7EDC8A",
                                      "T2_Naive" = "#066E10",
                                      "T3_Naive" = "#9E2307B4",
                                      "E2_Naive" = "#DE1B1BE6", 
                                      "P_Naive" = "#88CCEE",
                                      "L_Established" = "#7F7F2A",
                                      "T1_Established" = "#B85114",
                                      "E_Established" = "#9E2007", 
                                      "P_Established" = "#0C8AE3")
#PREPROCESSING ####
phyloseq #160,115 taxa and 234 samples 
      
##Selecting only Bacteria/Archaea - GTDB only has bacteria and archaea####
phyloseq.bacteria <- phyloseq

#WORKING ON BACTERIA/ARCHAEA ONLY####
# some QC checks of the "classified" reads per samples
min(sample_sums(phyloseq.bacteria)) # 1 (P1_0308)
max(sample_sums(phyloseq.bacteria)) # 109,875,085  (H21_0119) 
mean(sample_sums(phyloseq.bacteria)) #27,862,505
median(sample_sums(phyloseq.bacteria)) # 23,288,280
sort(sample_sums(phyloseq.bacteria))

##ZYMO MOCK COMMUNITIES AND NEGATIVE CONTROLS####
### Getting samples from ZYMOs and EB, NTC
phyloseq.bacteria.controls <- subset_samples(phyloseq.bacteria, 
  grepl("NTC|EB|Zymo", sample_names(phyloseq.bacteria)))
phyloseq.bacteria.controls <- prune_taxa(taxa_sums(phyloseq.bacteria.controls) > 0, phyloseq.bacteria.controls) 
phyloseq.bacteria.controls #20616 taxa, 15 samples(NTC, EB and Zymos)


###Zymo Mock communities compositions#####
phyloseq.bacteria.controls.zymo <- subset_samples(phyloseq.bacteria.controls, 
                                                  grepl("Zymo", sample_names(phyloseq.bacteria.controls)))
phyloseq.bacteria.controls.zymo <- prune_taxa(taxa_sums(phyloseq.bacteria.controls.zymo) > 0, 
                                              phyloseq.bacteria.controls.zymo) 
phyloseq.bacteria.controls.zymo #20473 taxa, 3 samples

#Relative abundance
phyloseq.bacteria.controls.zymo.ra <- transform_sample_counts(phyloseq.bacteria.controls.zymo, 
                                                        function(x) x/sum(x)*100) ##Relative abundance 
#GENUS LEVEL 
phyloseq.bacteria.controls.zymo.ra.genus <- tax_glom(phyloseq.bacteria.controls.zymo.ra, taxrank = "Genus", NArm = F)

#Melt to long format at genus level
phyloseq.bacteria.controls.zymo.ra.genus.melt <- psmelt(phyloseq.bacteria.controls.zymo.ra.genus)

#What are the top genera
top_genera_zymo <- phyloseq.bacteria.controls.zymo.ra.genus.melt %>%
  group_by(Genus) %>%
  summarise(`Mean Relative Abundance (%)` = mean(Abundance, na.rm = TRUE),
            `Standard Deviation` = sd(Abundance, na.rm = TRUE)) %>%
  arrange(desc(`Mean Relative Abundance (%)`))%>%
  head(n = 30)
top_genera_zymo

#SPECIES LEVEL
phyloseq.bacteria.controls.zymo.ra.melt <- psmelt(phyloseq.bacteria.controls.zymo.ra)

#What are the top species?
top_species_zymo <- phyloseq.bacteria.controls.zymo.ra.melt %>%
  group_by(Family, Species) %>%
  summarise(
    `Mean Relative Abundance (%)` = mean(Abundance, na.rm = TRUE),
    `Standard Deviation` = sd(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(`Mean Relative Abundance (%)`)) %>%
  slice_head(n = 30)%>%
  group_by(Family)%>%
  arrange(Family)
top_species_zymo

#####SUPPLEMENTARY TABLE 4 - MOCK COMMUNITIES SPECIES####
write_xlsx(top_species_zymo, 
           "SupplementaryTable4.xlsx")

#Plot - GENUS 
phyloseq.bacteria.controls.zymo.ra.genus.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.controls.zymo.ra.genus, 
                                                                          variable = "Enclosure",
                                                                          level = "Genus", 
                                                                          threshold = 0.5)
phyloseq.bacteria.controls.zymo.ra.genus.filt #10 genera over 0.5% mean RA

phyloseq.bacteria.controls.zymo.ra.genus.filt.melt <- psmelt(phyloseq.bacteria.controls.zymo.ra.genus.filt)%>%
  mutate(Genus = factor(Genus, 
                        levels = c(setdiff(Genus, 
                                           unique(grep("Others", Genus, value = TRUE))), 
                                   unique(grep("Others", Genus, value = TRUE)))))##Factoring the Phylum column so that "Others.." is the last category
levels(phyloseq.bacteria.controls.zymo.ra.genus.filt.melt$Genus) ##ok

#Palette
library(Polychrome)
#zymo.genus.palette <- as.character(brewer.pal(n = 10, name = "Accent")) #10 colors
zymo.genus.palette <- as.character(palette36.colors(10))
names(zymo.genus.palette) <- unique(phyloseq.bacteria.controls.zymo.ra.genus.filt.melt$Genus)
zymo.genus.palette$'Listeria' <- "#16FF32"
zymo.genus.palette$'Limosilactobacillus' <- "#5A5156"
zymo.genus.palette$'Staphylococcus' <- "darkblue" 
zymo.genus.palette$'Others <0.5% RA' <- 'grey90'

phyloseq.bacteria.controls.zymo.ra.genus.plot <- ggplot(phyloseq.bacteria.controls.zymo.ra.genus.filt.melt, aes(x=Sample, y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Zymo Mock Communities") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values =zymo.genus.palette) +
  guides(fill=guide_legend(title.position="top", nrow = 7))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 30, colour = "black", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
phyloseq.bacteria.controls.zymo.ra.genus.plot


#Plot - SPECIES 
phyloseq.bacteria.controls.zymo.ra.species.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.controls.zymo.ra, 
                                                                                variable = "Enclosure",
                                                                                level = "Species", 
                                                                                threshold = 0.5)
phyloseq.bacteria.controls.zymo.ra.species.filt #17 species over 0.5% mean RA

phyloseq.bacteria.controls.zymo.ra.species.filt.melt <- psmelt(phyloseq.bacteria.controls.zymo.ra.species.filt)
phyloseq.bacteria.controls.zymo.ra.species.filt.melt <- phyloseq.bacteria.controls.zymo.ra.species.filt.melt%>%
  mutate(Species = factor(Species, levels = c(
    "Limosilactobacillus fermentum",
    "Listeria monocytogenes_B", 
    "unclassified Listeria", 
    "Pseudomonas aeruginosa", 
    "Bacillus spizizenii", 
    "unclassified Bacillus",
    "Escherichia coli", 
    "Escherichia coli_F", 
    "Escherichia sp004211955",
    "unclassified Escherichia", 
    "Salmonella enterica",
    "unclassified Salmonella", 
    "Lactobacillus fermentum", 
    "Enterococcus faecalis", 
    "Staphylococcus aureus",
    "unclassified Staphylococcus", 
    "unclassified Enterobacteriaceae",
    "Others <0.5% RA"
    )))
  # mutate(Species = factor(Species, 
  #                       levels = c(setdiff(Species, 
  #                                          unique(grep("Others", Species, value = TRUE))), 
  #                                  unique(grep("Others", Species, value = TRUE)))))##Factoring the Species column so that "Others.." is the last category
levels(phyloseq.bacteria.controls.zymo.ra.species.filt.melt$Species) ##ok


#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
zymo_genus_base_colors <- zymo.genus.palette
#Make hues based on families within each ammonia-nitrite oxidizing group
palette_zymo_genus_df <- phyloseq.bacteria.controls.zymo.ra.species.filt.melt %>% 
  distinct(Genus, Species) %>%
  group_by(Genus) %>%  
  arrange(Species) %>%   
  mutate(
    base_color = zymo_genus_base_colors[Genus],
    shade = seq(-0.2, 0.2, length.out = n()),  # slightly wider range helps
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_zymo_genus_df

#Set up final palette
palette_zymo_genus <- setNames(
  palette_zymo_genus_df$color,
  palette_zymo_genus_df$Species)
palette_zymo_genus
palette_zymo_genus$'Others <0.5% RA' <- 'grey90'


#PLOT
phyloseq.bacteria.controls.zymo.ra.species.plot <-
  ggplot(phyloseq.bacteria.controls.zymo.ra.species.filt.melt, aes(x=Sample, y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Zymo Mock Communities") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = palette_zymo_genus) +
  guides(fill=guide_legend(title.position="top", ncol = 2))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 30, colour = "black", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) 
phyloseq.bacteria.controls.zymo.ra.species.plot

#####SUPPLEMENTARY FIGURE S2 #####
sfigure2 <- phyloseq.bacteria.controls.zymo.ra.species.plot
ggsave("SupplementaryFigure2.png", 
       sfigure2, 
       device = "png", 
       width = 10, height =10, 
       dpi = 500)

##SAMPLES#####
##New phyloseq of just samples
phyloseq.bacteria.samples <- subset_samples(phyloseq.bacteria, 
                                             !grepl("NTC|EB|Zymo", sample_names(phyloseq.bacteria)))
phyloseq.bacteria.samples #160115  taxa and 219 samples
sort(sample_sums(phyloseq.bacteria.samples))
#Taking out those with low counts
phyloseq.bacteria.samples <- prune_samples(sample_sums(phyloseq.bacteria.samples) > 300000, phyloseq.bacteria.samples) 
phyloseq.bacteria.samples <- prune_taxa(taxa_sums(phyloseq.bacteria.samples) > 0, phyloseq.bacteria.samples) 
phyloseq.bacteria.samples #160,112 taxa and 216 samples  (dropped P1_0308, H21_0109, and H21_0120)
sort(sample_sums(phyloseq.bacteria.samples)) #OK

###COMPARING SEQUENCING DEPTHS#######
cuso4_raw_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Raw/CuSO4_clean_raw_read_counts.csv')
str(cuso4_raw_read_counts)

#Keeping just those samples that I'm analyzing in phyloseq.bacteria.samples
sampleIDs_phyloseq.bacteria.samples <- phyloseq.bacteria.samples@sam_data$SampleID 
length(sampleIDs_phyloseq.bacteria.samples)#216 samples

#Filtering just phyloseq.bacteria.samples SampleIDs, adding metadata
cuso4_raw_read_counts_samples_metadata <- cuso4_raw_read_counts %>%
  filter(SampleID %in% sampleIDs_phyloseq.bacteria.samples)%>%
  dplyr::left_join(metadata, by = "SampleID")
nrow(cuso4_raw_read_counts_samples_metadata) #Ok, 216 samples

#Summary of raw reads
summary(cuso4_raw_read_counts_samples_metadata$Num_Reads_Forward_Raw)
sd(cuso4_raw_read_counts_samples_metadata$Num_Reads_Forward_Raw)
####Established vs Naive####
sequencing_depth_P1vsH21<- ggplot(cuso4_raw_read_counts_samples_metadata, 
                                  aes(x = Enclosure, y= Num_Reads_Forward_Raw, 
                                      color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Paired-end Reads", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  # geom_text_repel(aes(label = SampleID),   
  #                 size = 3)+
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "none",
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T) 
sequencing_depth_P1vsH21
wilcox_test(cuso4_raw_read_counts_samples_metadata, Num_Reads_Forward_Raw ~ Enclosure)

####Established and Naive over time####
sequencing_depth_P1andH21_overtime<-  ggplot(cuso4_raw_read_counts_samples_metadata, 
                                             aes(x = factor(Date_num), 
                                                 y= Num_Reads_Forward_Raw, 
                                                 color = Enclosure, 
                                                 fill = Enclosure)) +
  theme_bw()+
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Paired-End Reads", color = "System",
       x = "Day") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_color_manual(values=enclosure.palette)+
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "none",
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(size = 28, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text= element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
sequencing_depth_P1andH21_overtime

###COMPARING TRIMMED READS#######
cuso4_trimmed_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Trimmed/CuSO4_clean_trimmed_read_counts.csv')
str(cuso4_trimmed_read_counts)

#Filtering just phyloseq.bacteria.samples SampleIDs, adding metadata
cuso4_trimmed_read_counts_samples_metadata <- cuso4_trimmed_read_counts %>%
  filter(SampleID %in% sampleIDs_phyloseq.bacteria.samples)%>%
  dplyr::left_join(cuso4_raw_read_counts_samples_metadata, by = "SampleID")
nrow(cuso4_trimmed_read_counts_samples_metadata) #Ok, 216 samples

####Established vs Naive####
trimmed_reads_P1vsH21<- ggplot(cuso4_trimmed_read_counts_samples_metadata, 
                               aes(x = Enclosure, y= Num_Reads_Forward_Trimmed_Paired, 
                                   color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Paired-end Reads", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "none",
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
trimmed_reads_P1vsH21

####Established and Naive, Raw and Trimmed over time####
cuso4_trimmed_read_counts_samples_metadata_long <- 
  cuso4_trimmed_read_counts_samples_metadata %>%
  pivot_longer(cols = c(Num_Reads_Forward_Trimmed_Paired,
                        Num_Reads_Forward_Raw),
               names_to = "Read_Status",
               values_to = "Num_Paired_Reads") %>%
  mutate(
    Read_Status = dplyr::recode(Read_Status,
                                "Num_Reads_Forward_Trimmed_Paired" = "QC - Trimmomatic",
                                "Num_Reads_Forward_Raw" = "Raw"))%>%
  mutate(
    Read_Status = factor(Read_Status, levels = c("Raw",
                                                 "QC - Trimmomatic"))
  )
#Now, plot
trimmed_and_raw_reads_time<- ggplot(cuso4_trimmed_read_counts_samples_metadata_long, 
                                    aes(x = factor(Date_num), 
                                        y= Num_Paired_Reads, 
                                        color = Read_Status)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Paired-End Reads", color = "Read Status",
       x = "Day") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_color_manual(values=reads.palette)+
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(size = 28, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text= element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
trimmed_and_raw_reads_time

####Established and Naive, Raw vs Trimmed####
trimmed_and_raw_reads_P1vsH21 <- ggplot(cuso4_trimmed_read_counts_samples_metadata_long, 
                                        aes(x = Enclosure, 
                                            y= Num_Paired_Reads, 
                                            color = Read_Status)) +
  theme_bw() +
  labs(y= "Paired-End Reads", color = "Read Status",
       x = "System") +
  geom_boxplot(position = position_dodge(width = 0.8), 
               alpha = 0.3) +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.4, 
    dodge.width = 0.8),
    size = 3, shape = 18, alpha = 0.8) +
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  scale_color_manual(values=reads.palette)+
  scale_x_discrete(labels= c("H21" = "Naive", 
                             "P1" = "Established"))+
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(size = 28, colour = "black", face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.y= element_text(colour = "black", size = 20),
        axis.text.x= element_text(colour = "black", size = 25),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
trimmed_and_raw_reads_P1vsH21

###COMPARING CLASSIFIED READS BY KRAKEN#######
kraken_unclassified_reads <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/Paired_end_mode_GTDB_updated_20260602/Conf_005/unclassifieds_kraken_analytic_matrix.conf_005.csv')

#Mock communities 
kraken_unclassified_reads_pos_control_metadata <- kraken_unclassified_reads %>%
  filter(grepl("Zymo", SampleID))%>%
  dplyr::left_join(cuso4_trimmed_read_counts_samples_metadata, by = "SampleID")%>%
  rename(Kraken2_Input_PairedEnd_Reads = Total, 
         Kraken2_Unclassified_PairedEnd_Reads = NumberUnclassified, 
         Kraken2_Unclassified_Percentage_Reads = PercentUnclassified)%>%
  mutate(Kraken2_Classified_Percentage_Reads = (100 - Kraken2_Unclassified_Percentage_Reads))
nrow(kraken_unclassified_reads_pos_control_metadata) #Ok, 3 mock communities
kraken_unclassified_reads_pos_control_metadata$Kraken2_Classified_Percentage_Reads
#93.92, 93.97, 93.94

#Filtering just samples included in phyloseq.bacteria.samples, adding metadata, calculating percentage classified
kraken_unclassified_reads_samples_metadata <- kraken_unclassified_reads %>%
  filter(SampleID %in% sampleIDs_phyloseq.bacteria.samples)%>%
  dplyr::left_join(cuso4_trimmed_read_counts_samples_metadata, by = "SampleID")%>%
  rename(Kraken2_Input_PairedEnd_Reads = Total, 
         Kraken2_Unclassified_PairedEnd_Reads = NumberUnclassified, 
         Kraken2_Unclassified_Percentage_Reads = PercentUnclassified)%>%
  mutate(Kraken2_Classified_Percentage_Reads = (100 - Kraken2_Unclassified_Percentage_Reads))
nrow(kraken_unclassified_reads_samples_metadata) #Ok, 216 samples

#Descriptive stats
summary(kraken_unclassified_reads_samples_metadata$Kraken2_Classified_Percentage_Reads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 29.52   44.30   62.23   58.87   74.48   82.88 
sd(kraken_unclassified_reads_samples_metadata$Kraken2_Classified_Percentage_Reads) #16.24
#Descriptive stats per group
kraken_unclassified_reads_samples_metadata %>%
  group_by(Enclosure)%>%
  summarise(mean_percentage_classified_reads = mean(Kraken2_Classified_Percentage_Reads),
            sd_percentage_classified_reads = sd(Kraken2_Classified_Percentage_Reads), 
            min_percentage_classified_reads = min(Kraken2_Classified_Percentage_Reads), 
            max_percentage_classified_reads = max(Kraken2_Classified_Percentage_Reads))
# Enclosure mean_percentage_classified_reads sd_percentage_classified_reads min_percentage_classifi…¹ max_percentage_class…²
# H21                                 45.5                           9.73                      29.5                   80.2
# P1                                  68.8                          13.0                       30.4                   82.9


####Kraken2 Classified Percentages Established vs Naive####
kraken2_classified_read_percentages_P1vsH21<- ggplot(kraken_unclassified_reads_samples_metadata, 
                                                     aes(x = Enclosure, 
                                                         y= Kraken2_Classified_Percentage_Reads, 
                                                         color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Percentage (%) Classified Reads", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "none",
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 22, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
kraken2_classified_read_percentages_P1vsH21
wilcox_test(kraken_unclassified_reads_samples_metadata, Kraken2_Classified_Percentage_Reads ~ Enclosure) #S., 9.58e-23

####Kraken2 Classified Percentages Established and Naive over time#### 
kraken2_classified_read_percentages_P1andH21_overtime <- ggplot(kraken_unclassified_reads_samples_metadata,
                                                                aes(x = factor(Date_num), 
                                                                    y= Kraken2_Classified_Percentage_Reads, 
                                                                    color = Enclosure)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Percentage (%) Classified Reads",
       x = "Day") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_color_manual(values=enclosure.palette)+
  theme(
    plot.title = element_text(colour = "black", size = 32, face = "bold"),
    legend.position = "none",
    # legend.text = element_text(size = 20),
    # legend.title = element_text(size = 22, face = "bold"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    strip.background = element_rect(fill = "black"),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text = element_text(colour = "white", size = 28, face = "bold"),
    axis.title = element_text(size = 28, colour = "black"),
    axis.ticks.x = element_blank(),
    axis.text= element_text(colour = "black", size = 20),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
kraken2_classified_read_percentages_P1andH21_overtime


###COMPARING SAMPLE SUMS OTUs (CLASSIFIED READS FROM KRAKEN)#######
##ALL TAXA#####
sample.sums <- sample_sums(phyloseq.bacteria.samples) #making a sample sums object
phyloseq.bacteria.samples.samplessums.df <- cbind(phyloseq.bacteria.samples@sam_data, 
                                                        sample.sums) #combining sample sums with metaphyloseq
nrow(phyloseq.bacteria.samples.samplessums.df) #216 samples, OK 
phyloseq.bacteria.samples.samplessums.df$sampleID <- rownames(phyloseq.bacteria.samples.samplessums.df) ##making a sampleID column

#Descriptive stats
summary(phyloseq.bacteria.samples.samplessums.df$sample.sums)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 6510692  14633429  23745867  29392014  40452430 109875085 
sd(phyloseq.bacteria.samples.samplessums.df$sample.sums) #19,009,601
#Descriptive stats per group
phyloseq.bacteria.samples.samplessums.df %>%
  group_by(Enclosure)%>%
  summarise(mean_OTU_counts = mean(sample.sums),
            sd_OTU_counts = sd(sample.sums), 
            min_OTU_counts = min(sample.sums), 
            max_OTU_counts = max(sample.sums))
# Enclosure mean_OTU_counts sd_OTU_counts min_OTU_counts max_OTU_counts
# H21             18687286.     13700008.        6510692      109875085
# P1              37334231.     18533532.        7676690      101124237

#Samplesums of Zymo mock communities
sample_sums(phyloseq.bacteria.controls)

###Established vs Naive####
bacteria_archaea_samplesums_P1vsH21<- ggplot(phyloseq.bacteria.samples.samplessums.df, 
                                             aes(x = Enclosure, y= sample.sums, 
                                                 color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "OTUs", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
bacteria_archaea_samplesums_P1vsH21

##Stats 
wilcox_test(phyloseq.bacteria.samples.samplessums.df, sample.sums~Enclosure) #S. p = 2.83e-18

###Established and Naive over time####
bacteria_archaea_samplesums_P1andH21_overtime<- ggplot(phyloseq.bacteria.samples.samplessums.df, 
                                                       aes(x = factor(Date_num), 
                                                           y= sample.sums, 
                                                           color = Enclosure)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "OTUs", x = "Day", color = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_color_manual(values=enclosure.palette)+
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "none",
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(size = 28, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text= element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
bacteria_archaea_samplesums_P1andH21_overtime

###H21####
phyloseq.bacteria.samples_H21 <- subset_samples(phyloseq.bacteria.samples, Enclosure == "H21")
phyloseq.bacteria.samples_H21 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples_H21) > 0, 
                                                  phyloseq.bacteria.samples_H21)
phyloseq.bacteria.samples_H21 #159,004 taxa and 92 samples
range(phyloseq.bacteria.samples_H21@sam_data$Collection_Date)#OK, "2023-10-09" through "2024-03-02"

###P1####
phyloseq.bacteria.samples_P1 <- subset_samples(phyloseq.bacteria.samples, Enclosure == "P1")
phyloseq.bacteria.samples_P1 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples_P1) > 0, 
                                                 phyloseq.bacteria.samples_P1)
phyloseq.bacteria.samples_P1 #158,612 taxa and 124 samples
range(phyloseq.bacteria.samples_P1@sam_data$Collection_Date)#OK, "2023-11-14" through "2024-04-30"


##NITRIFYING TAXA####
nitrifiers_all <- subset_taxa(phyloseq.bacteria.samples, 
                              Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                            Family == "Chromatiaceae" | # AOB (Nitrosococcus genus)
                            Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales
                            Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales
                            Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales
                            Genus == "JBMCYV01" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is JBMCYV01 sp04919959
                            Order == "Nitrosomirales" | #AOA
                            Order == "Candidatus Nitrosocaldales" | #AOA
                            Family == "Nitrospiraceae" | # NOB/Commamox; some!
                            Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                            Family == "Nitrobacteraceae" | #NOB
                            Family == "Gallionellaceae" | #NOB none
                            Family == "Nitrospinaceae") # NOB; some, plus a new one!

nitrifiers_all #1437 taxa and 216 samples
nitrifiers <- subset_samples(nitrifiers_all, sample_sums(nitrifiers_all) > 0)
nitrifiers #1437 taxa and 216 samples 

##QC checks again
min(sample_sums(phyloseq.bacteria.samples)) #6,510,692 (H21_1101)
max(sample_sums(phyloseq.bacteria.samples)) #109,875,085(H21_0119) 
mean(sample_sums(phyloseq.bacteria.samples)) #29,392,014
median(sample_sums(phyloseq.bacteria.samples)) #23,745,867
sort(sample_sums(phyloseq.bacteria.samples)) 


##NITRIFIERS#####
sample.sums.nit <- sample_sums(nitrifiers) #making a sample sums object
nitrifiers.samplesums.df <- cbind(nitrifiers@sam_data, 
                                      sample.sums.nit) #combining sample sums with metadata
nitrifiers.samplesums.df
nitrifiers.samplesums.df$SampleID <- rownames(nitrifiers.samplesums.df) ##making a sampleID column

###Established vs Naive####
nitrifier_bacteria_archaea_samplesums_P1vsH21<- ggplot(nitrifiers.samplesums.df, 
                                  aes(x = Enclosure, y= sample.sums.nit, 
                                      color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "OTUs", color = "System", fill = "System",) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
nitrifier_bacteria_archaea_samplesums_P1vsH21
##Stats 
wilcox_test(nitrifiers.samplesums.df, sample.sums.nit~Enclosure) #S. p = 0.000309

#RELATIVE ABUNDANCE####
any(sample_sums(phyloseq.bacteria.samples)== 0) ## no samples with 0 OTUs
phyloseq.bacteria.samples.ra <- transform_sample_counts(phyloseq.bacteria.samples, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data
##CLASSIFICATION PERCENTAGES AT DIFFERENT TAXONOMIC LEVELS####
###PHYLUM######
phyloseq.bacteria.samples_phylum.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Phylum", NArm = F) 
phyloseq.bacteria.samples_phylum.ra #204 phyla and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"])) #204 phyla (so No duplicates)

Unknown_phylum_abundance <- phyloseq.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Phylum, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_phylum_abundance ##0% abundance by Unknown Phyla

Unclassified_phylum_abundance <- phyloseq.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##0% abundance by Unclassified Phyla

Classified_phylum_abundance <- phyloseq.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##100% abundance by Classified Phyla

##Checking on excel
write.csv(phyloseq.bacteria.samples_phylum.ra@otu_table, "phylum_otus.csv")
write.csv(phyloseq.bacteria.samples_phylum.ra@tax_table, "phylum_taxa.csv")  

# #How many unclassified? -None
# phyloseq.bacteria.samples_phylum.unclassified.ra <- prune_taxa(
#   grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"]),
#   phyloseq.bacteria.samples_phylum.ra)
phyloseq.bacteria.samples_phylum.unclassified.ra  <- 0

# #How many unknown? - none 
# phyloseq.bacteria.samples_phylum.unknown.ra <- prune_taxa(
#   grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"]),
#   phyloseq.bacteria.samples_phylum.ra)
# phyloseq.bacteria.samples_phylum.unknown.ra

#Keep just classified Phyla - ALL
phyloseq.bacteria.samples_phylum.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples_phylum.ra)
phyloseq.bacteria.samples_phylum.classified.ra ##204 classified (not unknown or unclassified) Phyla

###CLASS#####
phyloseq.bacteria.samples_class.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Class", NArm = F) 
phyloseq.bacteria.samples_class.ra #687 classes and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"])) #687 classes (so No duplicates)

Unknown_class_abundance <- phyloseq.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Class, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_class_abundance #0% Abundance by Unknown classes

Unclassified_class_abundance <- phyloseq.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##3.2% Abundance by Unclassified Classes

Classified_class_abundance <- phyloseq.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##96.8% Abundance by Classified classes

##Checking on excel
write.csv(phyloseq.bacteria.samples_class.ra@otu_table, "class_otus.csv")
write.csv(phyloseq.bacteria.samples_class.ra@tax_table, "class_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"]),
  phyloseq.bacteria.samples_class.ra)
phyloseq.bacteria.samples_class.unclassified.ra #49 unclassified classes

# #How many unknown? - NONE
# phyloseq.bacteria.samples_class.unknown.ra <- prune_taxa(
#   grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"]),
#   phyloseq.bacteria.samples_class.ra)
# phyloseq.bacteria.samples_class.unknown.ra 

#Keep just classified Classes
phyloseq.bacteria.samples_class.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"]),
  phyloseq.bacteria.samples_class.ra)
phyloseq.bacteria.samples_class.classified.ra #638 classified classes

###ORDER######
phyloseq.bacteria.samples_order.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Order", NArm = F) 
phyloseq.bacteria.samples_order.ra #2405 orders

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"])) #2404 orders (1 duplicate)
order_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"])
unique(order_taxa_vec[duplicated(order_taxa_vec)])
#"unclassified WOR-3" is duplicated

Unknown_order_abundance <- phyloseq.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Order, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_order_abundance ##0% abundance by Unknown Orders

Unclassified_order_abundance <- phyloseq.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##11.3% abundance by Unclassified Orders

Classified_order_abundance <- phyloseq.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##88.7% abundance by Classified orders

#Checking on excel
write.csv(phyloseq.bacteria.samples_order.ra@otu_table, "order_otus.csv")
write.csv(phyloseq.bacteria.samples_order.ra@tax_table, "order_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"]),
  phyloseq.bacteria.samples_order.ra)
phyloseq.bacteria.samples_order.unclassified.ra #241 unclassified orders

# #How many unknown? -NONE
# phyloseq.bacteria.samples_order.unknown.ra <- prune_taxa(
#   grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"]),
#   phyloseq.bacteria.samples_order.ra)
# phyloseq.bacteria.samples_order.unknown.ra 

#Keep just classified Orders
phyloseq.bacteria.samples_order.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"]),
  phyloseq.bacteria.samples_order.ra)
phyloseq.bacteria.samples_order.classified.ra #2164 classified orders
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_order.classified.ra)[, "Order"])) ##2164 classified orders (unique - without duplicates)

###FAMILY######
phyloseq.bacteria.samples_family.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Family", NArm = F) 
phyloseq.bacteria.samples_family.ra #6755 families
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"])) #6728 taxa (there are duplicates)
family_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"])
unique(family_taxa_vec[duplicated(family_taxa_vec)]) #Duplicates:
# "unclassified UBA1135"         "unclassified J058"            "unclassified UBA8108"         "unclassified SZUA-567"        "unclassified UBA796"         
# "unclassified UBA9160"         "unclassified UBA2968"         "unclassified UBA2361"         "unclassified RBG-16-71-46"    "unclassified WGA-4E"         
# "unclassified GWC2-55-46"      "unclassified UBA1018"         "unclassified UBA10199"        "unclassified JALEGL01"        "unclassified CADDZG01"       
# "unclassified XYA12-FULL-58-9" "unclassified UBA11872"        "unclassified UBA5829"         "unclassified MSB-5A5"         "unclassified SAR324"         
# "unclassified 4484-113"        "unclassified WOR-3"           "unclassified UBA7883"         "unclassified DTGP01"          "unclassified UBA6099"        
# "unclassified UBA1177"         "unclassified AC1"   

Unknown_family_abundance <- phyloseq.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_family_abundance #0% abundance by Unknown Families

Unclassified_family_abundance <- phyloseq.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##13.7% abundance by Unclassified Families

Classified_family_abundance <- phyloseq.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##86.3% abundance by Classified Families

#Checking on excel
write.csv(phyloseq.bacteria.samples_family.ra@otu_table, "family_otus.csv")
write.csv(phyloseq.bacteria.samples_family.ra@tax_table, "family_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"]),
  phyloseq.bacteria.samples_family.ra)
phyloseq.bacteria.samples_family.unclassified.ra #823 unclassified families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.unclassified.ra)[, "Family"])) ##796 classified families (unique - without duplicates)

# #How many unknown? - NONE
# phyloseq.bacteria.samples_family.unknown.ra <- prune_taxa(
#   grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"]),
#   phyloseq.bacteria.samples_family.ra)
# phyloseq.bacteria.samples_family.unknown.ra #3298 "unknown" families
# length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.unknown.ra)[, "Family"]))

#Keep just classified Families
phyloseq.bacteria.samples_family.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"]),
  phyloseq.bacteria.samples_family.ra)
phyloseq.bacteria.samples_family.classified.ra #5932 classified families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.classified.ra)[, "Family"]))#5932 classified families (unique - without duplicates)

###GENUS ######
phyloseq.bacteria.samples_genus.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Genus", NArm = F) 
phyloseq.bacteria.samples_genus.ra #32836 genera
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"])) #32604 taxa (there are duplicates)
genus_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"])
unique(genus_taxa_vec[duplicated(genus_taxa_vec)]) #216 duplicated unique ones

Unknown_genus_abundance <- phyloseq.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_genus_abundance ##0%  abundance by unknown genera

Unclassified_genus_abundance <- phyloseq.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance)) #Sum across OTUs
Unclassified_genus_abundance ##26.3% abundance by unclassified genera

Classified_genus_abundance <- phyloseq.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##73.7% abundance by Classified Genera

#Checking on excel
write.csv(phyloseq.bacteria.samples_genus.ra@otu_table, "genus_otus.csv")
write.csv(phyloseq.bacteria.samples_genus.ra@tax_table, "genus_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples_genus.ra)
phyloseq.bacteria.samples_genus.unclassified.ra #3431 unclassified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.unclassified.ra)[, "Genus"])) ##3199 unclassified genera (unique - without duplicates)

# #How many unknown? -NONE
# phyloseq.bacteria.samples_genus.unknown.ra <- prune_taxa(
#   grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"]),
#   phyloseq.bacteria.samples_genus.ra)
# phyloseq.bacteria.samples_genus.unknown.ra #3755 "unknown" genera
# length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.unknown.ra)[, "Genus"])) ##3755 unknown genera (unique - without duplicates)


#Keep just classified Genera
phyloseq.bacteria.samples_genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples_genus.ra)
phyloseq.bacteria.samples_genus.classified.ra #29405 classified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.classified.ra)[, "Genus"])) ##29405 classified genera (unique - without duplicates)


###SPECIES######
phyloseq.bacteria.samples.ra ##160,112  Species- OTUs
phyloseq.bacteria.samples_species.ra <- phyloseq.bacteria.samples.ra
phyloseq.bacteria.samples_species.ra #160112 Species

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_species.ra)[, "Species"])) #158975 species (there are duplicates)
species_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_species.ra)[, "Species"])
unique(species_taxa_vec[duplicated(species_taxa_vec)]) #duplicated unique ones

#Extract OTU and tax table
species_level_otu <- as(phyloseq::otu_table(phyloseq.bacteria.samples_species.ra), "matrix")
species_level_tax <- as(phyloseq::tax_table(phyloseq.bacteria.samples_species.ra), "matrix")

# Identify unclassified and classifiedspecies
unclassified_idx_species <- grepl("unclassified", species_level_tax[, "Species"], ignore.case = TRUE)
classified_idx_species <- !grepl("unclassified|unknown", species_level_tax[, "Species"], ignore.case = TRUE)

# Subset OTUs
otu_unclassified_species <- species_level_otu[unclassified_idx_species, , drop = FALSE]
otu_classified_species <- species_level_otu[classified_idx_species, , drop = FALSE]

# Mean abundance per OTU (row means)
otu_unclassified_species_means <- rowMeans(otu_unclassified_species)
otu_classified_species_means <- rowMeans(otu_classified_species)

# Sum across OTUs
unclassified_species_sum <- sum(otu_unclassified_species_means) #30.91% abundance by unclassified species
classified_species_sum <- sum(otu_classified_species_means) #69.1% abundance by classified species

#Checking on excel
write.csv(phyloseq.bacteria.samples_species.ra@otu_table, "species_otus.csv")
write.csv(phyloseq.bacteria.samples_species.ra@tax_table, "species_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_species.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_species.ra)[, "Species"]),
  phyloseq.bacteria.samples_species.ra)
phyloseq.bacteria.samples_species.unclassified.ra #16,655 unclassified species
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_species.unclassified.ra)[, "Species"])) ##15,518 unclassified species (unique - without duplicates)

# #How many unknown? -NONE
# phyloseq.bacteria.samples_species.unknown.ra <- prune_taxa(
#   grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_species.ra)[, "Species"]),
#   phyloseq.bacteria.samples_species.ra)
# phyloseq.bacteria.samples_species.unknown.ra #4 "unknown" species
# length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_species.unknown.ra)[, "Species"])) 

#Keep just classified Genera
phyloseq.bacteria.samples_species.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_species.ra)[, "Species"]),
  phyloseq.bacteria.samples_species.ra)
phyloseq.bacteria.samples_species.classified.ra #143,457 classified species
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_species.classified.ra)[, "Species"])) ##143,457 classified species (unique - without duplicates)


###SUPPLEMENTARY TABLE 3####
stable3 <- data.frame(
  "Taxonomic level" = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  # Total taxa at each level
  "Number of Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples_phylum.ra)),
    length(taxa_names(phyloseq.bacteria.samples_class.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_order.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_family.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_genus.ra)),
    length(taxa_names(phyloseq.bacteria.samples_species.ra))),

  
  # Unclassified taxa
  "Number of Unclassified Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples_phylum.unclassified.ra)),
    length(taxa_names(phyloseq.bacteria.samples_class.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_order.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_family.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_genus.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_species.unclassified.ra))),
  
  # Classified taxa
  "Number of Classified Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples_phylum.classified.ra)),
    length(taxa_names(phyloseq.bacteria.samples_class.classified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_order.classified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_family.classified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples_genus.classified.ra)),
    length(taxa_names(phyloseq.bacteria.samples_species.classified.ra))),
  
  # Relative abundances
  "Mean Relative Abundance (%) of Unclassified Taxa" = c(
    Unclassified_phylum_abundance$Unclassified_sum, 
    Unclassified_class_abundance$Unclassified_sum,
    Unclassified_order_abundance$Unclassified_sum,
    Unclassified_family_abundance$Unclassified_sum,
    Unclassified_genus_abundance$Unclassified_sum,
    unclassified_species_sum
  ),
  
  "Mean Relative Abundance (%) of Classified Taxa" = c(
    Classified_phylum_abundance$Classified_sum,
    Classified_class_abundance$Classified_sum, 
    Classified_order_abundance$Classified_sum,
    Classified_family_abundance$Classified_sum,
    Classified_genus_abundance$Classified_sum,
    classified_species_sum),
  check.names = FALSE
) %>%
  mutate(`Percentage of Classified Taxa` = 
           (`Number of Classified Taxa` / `Number of Taxa`) * 100)
stable3
#Make into excel file
write_xlsx(stable3, 
          "SupplementaryTable3.xlsx")



#ALPHA DIVERSITY ######
## ALL COMMUNITIES#####
alpha_div1 <- phyloseq::estimate_richness(phyloseq.bacteria.samples, 
                                          measures = c("Observed", "Shannon")) # richness, diversity
alpha_div2 <- microbiome::evenness(phyloseq.bacteria.samples, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness

#Combine alpha div metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div

#Metadata and div metrics
alpha_div_meta <- cbind(phyloseq.bacteria.samples@sam_data, 
                        alpha_div) 
alpha_div_meta 

#Pivot to long format 
alpha_div_meta_long <- 
  alpha_div_meta %>%
  pivot_longer(cols = c(Observed, Shannon, pielou),  
               names_to = "alpha_div_metric", 
               values_to = "alpha_div_value") 
alpha_div_meta_long

##Factoring alpha div metrics
alpha_div_meta_long$alpha_div_metric <- factor(alpha_div_meta_long$alpha_div_metric, levels = c("Observed","pielou", "Shannon"))

###P1 vs H21#####
alpha_div_P1vsH21 <- ggplot(alpha_div_meta_long, 
                            aes(x = Enclosure, y= alpha_div_value, 
                                fill= Enclosure, colour = Enclosure)) +
  theme_bw() +
  labs(color = "System", fill = "System") +
  facet_wrap(~alpha_div_metric, 
             scales = "free",
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS"))) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = enclosure.palette,labels = c("H21" = "Naive", "P1" = "Established")) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))+
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_P1vsH21

####Alpha diversity indexes and Water quality over time#####
#Editing alpha diversity and metadata dataframe
alpha_div_wq_time_long <- alpha_div_meta %>%
  mutate(Copper_keep = Copper_level_mg_L) %>%   #keep a column as copper to be able to color the plot 
  pivot_longer(cols = c("Copper_level_mg_L",
                        "Temperature_F",
                        "Chlorine_mg_L", 
                        "pH_spu", 
                        "Ammonia_mg_L",
                        "Nitrite_mg_L",
                        "Nitrate_UV_mg_L", 
                        "Salinity_ppt",
                        "Alkalinity_mg_L", 
                        "Shannon",
                        "Observed",
                        "pielou"),
               names_to = "Index",
               values_to = "Index_value")%>%
  mutate(Index = factor(Index, levels = c(
    "Observed",
    "Shannon",
    "pielou",
    "Copper_level_mg_L",
    "Ammonia_mg_L",
    "Nitrite_mg_L",
    "Nitrate_UV_mg_L",
    "pH_spu", 
    "Salinity_ppt", 
    "Temperature_F",
    "Chlorine_mg_L", 
    "Alkalinity_mg_L"
  )))

#####Adding vertical dashed line for when "phase" changes#######
metadata_phyloseq.bacteria.samples <- data.frame(phyloseq.bacteria.samples@sam_data)
#When do phases switch?
#Naive
metadata_phyloseq.bacteria.samples %>%
  filter(Enclosure == "H21") %>%
  group_by(Date_num_phase_naive) %>%
  summarise(
    min_naive = min(Date_num_naive, na.rm = TRUE),
    max_naive = max(Date_num_naive, na.rm = TRUE)
  )

#Established
metadata_phyloseq.bacteria.samples %>%
  filter(Enclosure == "P1") %>%
  group_by(Date_num_phase_established) %>%
  summarise(
    min_naive = min(Date_num_established, na.rm = TRUE),
    max_naive = max(Date_num_established, na.rm = TRUE)
  )


#Vertical lines to be places
line_breaks_phases <- data.frame(
  Enclosure = c("H21", "H21", "H21", "H21", "H21", "H21","H21",
                "H21", "H21", "H21", "H21", "H21", "H21","H21",
                "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1"), 
  Date_num = c("1", "27", "29", "38", "39", "51", "52", "81", "86", "108", "109", "135", "136", "146",
               "1", "53", "54", "65", "66", "104", "105", "169"))
line_breaks_phases

#####Plot - Just  Shannon, Copper and Ammonia levels - Date number since start of sampling as as.factor####
alpha_div_wq_date_num_factor <- ggplot(alpha_div_wq_time_long%>%
                                filter(Index %in% c("Copper_level_mg_L",
                                                    "Shannon",
                                                    "Ammonia_mg_L")),
                              aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "MICROBIOME\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)")))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme(
    legend.position = c(0.85, 1.7),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.y  = element_text(colour = "white", size = 30, 
                                 face = "bold", angle = 0),
    strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_wq_date_num_factor
ggsave("alpha_div_wq_date_num_factor.png",
       alpha_div_wq_date_num_factor, 
       device = "png", 
       dpi = 600, 
       height = 12, 
       width = 25)

#####Plot - Other Water quality measures- Date number since start of sampling as as.factor#####
alpha_div_wq_date_num_factor_other_metadata <- ggplot(alpha_div_wq_time_long%>%
                                                        filter(Index %in% c("Copper_level_mg_L",
                                                                            "Shannon",
                                                                            "Ammonia_mg_L",
                                                                            #"Chlorine_mg_L", 
                                                                            #"Alkalinity_mg_L",
                                                                            "Temperature_F",
                                                                            "pH_spu",
                                                                            "Salinity_ppt"
                                                                            )),
                              aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black", 
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  # 
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  # 
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  # 
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  # 
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "MICROBIOME\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)",
                                      "Temperature_F"= "Temperature (F)",
                                      "Salinity_ppt" = "Salinity (ppt)", 
                                      "pH_spu" = "pH (spu)"
                                      #"Chlorine_mg_L" = "Chlorine (mg/L)",
                                      #"Alkalinity_mg_L" = "Alkalinity (mg/L)"
                                      )))+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme(
    legend.position = c(0.95, 1.3),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.y  = element_text(colour = "white", size = 30, 
                                 face = "bold", angle = 0),
    strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_wq_date_num_factor_other_metadata

###Copper vs Shannon levels #####
##Time series
copper_shannon <- ggplot(alpha_div_meta,
                      aes(x = Copper_level_mg_L, 
                          y = Shannon,
                          color = Collection_Month)) +
  geom_point(size = 3, shape = 18) +
  labs(title= "Shannon's Diversity vs Copper levels (mg/L) \nBacterial - Archaeal Communities", 
       y = "SHANNON'S INDEX",
       x = "Copper levels (mg/L)",
       color = "Month") +
  facet_grid(~Enclosure,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  theme_bw() +
  #geom_smooth(method="loess", se=TRUE) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12, angle = 45, vjust = 0.5),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(colour = "black", size = 20),
        axis.text.x = element_text(colour = "black", size = 20, angle = 45,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
copper_shannon

####H21########
alpha_div_meta_overall_correlation_H21_df <- alpha_div_meta%>%
  filter(Enclosure == "H21")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month),
         Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%  # convert to factor 
  filter(!Collection_Month %in% c("2023-09", "2024-03"))%>%
  filter(!is.na(Copper_level_mg_L),
         !is.na(Collection_Date))%>%
  mutate(Attempt = case_when(is.na(Attempt) ~ "No CuSO4 Addition",
                             Attempt == "1" ~ "PHASE 1",
                             Attempt == "2" ~ "PHASE 2",
                             Attempt == "3" ~ "PHASE 3",
                             TRUE ~ as.character(Attempt)))%>%
  mutate(Attempt = factor(Attempt, 
                          levels = c("PHASE 1", "PHASE 2", "PHASE 3", 
                                     "No CuSO4 Addition")))%>%
  mutate(Collection_Month = case_when(Collection_Month == "2023-10"~ "Oct-2023", 
                                      Collection_Month == "2023-11"~ "Nov-2023", 
                                      Collection_Month == "2023-12"~ "Dec-2023", 
                                      Collection_Month == "2024-01"~ "Jan-2024", 
                                      Collection_Month == "2024-02"~ "Feb-2024", 
                                      Collection_Month == "2024-03"~ "Mar-2024",
                                      TRUE ~ as.character(Collection_Month)))%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("Oct-2023",
                                                                "Nov-2023", 
                                                                "Dec-2023", 
                                                                "Jan-2024", 
                                                                "Feb-2024", 
                                                                "Mar-2024")))

#Plot
alpha_div_meta_overall_correlation_H21_plot <- ggplot(alpha_div_meta_overall_correlation_H21_df,
                                           aes(x = Copper_level_mg_L, 
                                               y = Shannon)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_color_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Shannon",
       x = "Copper levels (mg/L)", 
       title = "NAIVE", 
       color = "Collection Month") +
  theme_bw() +
  geom_smooth(method="loess", 
              se=TRUE,
              color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(method = "spearman",
           label.x.npc = 0.65,
           label.y.npc = 1,
           size = 8) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
alpha_div_meta_overall_correlation_H21_plot
ggsave("alpha_div_meta_overall_correlation_H21_plot.png", 
       alpha_div_meta_overall_correlation_H21_plot, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)

####P1########
alpha_div_meta_overall_correlation_P1_df <- alpha_div_meta%>%
  filter(Enclosure == "P1")%>%
  filter(Collection_Date > "2023-06-01")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month))%>%  # convert to factor 
  mutate(Attempt = case_when(is.na(Attempt) ~ "No CuSO4 Addition",
                             Attempt == "1" ~ "PHASE 1",
                             Attempt == "2" ~ "PHASE 2",
                             Attempt == "3" ~ "PHASE 3",
                             TRUE ~ as.character(Attempt)))%>%
  mutate(Attempt = factor(Attempt, 
                          levels = c("PHASE 1", "PHASE 2", "PHASE 3", 
                                     "No CuSO4 Addition")))%>%
  mutate(Collection_Month = case_when(Collection_Month == "2023-11"~ "Nov-2023", 
                                      Collection_Month == "2023-12"~ "Dec-2023", 
                                      Collection_Month == "2024-01"~ "Jan-2024", 
                                      Collection_Month == "2024-02"~ "Feb-2024", 
                                      Collection_Month == "2024-03"~ "Mar-2024",
                                      Collection_Month == "2024-04"~ "Apr-2024",
                                      TRUE ~ as.character(Collection_Month)))%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("Nov-2023", 
                                                                "Dec-2023", 
                                                                "Jan-2024", 
                                                                "Feb-2024", 
                                                                "Mar-2024", 
                                                                "Apr-2024")))

#Plot
alpha_div_meta_overall_correlation_P1_plot <- ggplot(alpha_div_meta_overall_correlation_P1_df,
                                                      aes(x = Copper_level_mg_L, 
                                                          y = Shannon)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_color_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Shannon",
       x = "Copper levels (mg/L)", 
       title = "ESTABLISHED", 
       color = "Collection Month") +
  theme_bw() +
  geom_smooth(method="loess", 
              se=TRUE,
              color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(method = "spearman",
           label.x.npc = 0,
           label.y.npc = 1,
           size = 8) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 23),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
alpha_div_meta_overall_correlation_P1_plot
ggsave("alpha_div_meta_overall_correlation_P1_plot.png", 
       alpha_div_meta_overall_correlation_P1_plot, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)

##Linear models ####
alpha_div_meta_clean <- alpha_div_meta %>%
    filter(!is.na(Shannon),
           !is.na(Copper_level_mg_L),
           !is.na(Collection_Date))%>%
  arrange(Enclosure, Collection_Date)%>%
  mutate(Enclosure = factor(Enclosure, levels = c("P1", "H21")),
         Collection_Date = factor(Collection_Date),
         Collection_Month = factor(Collection_Month)) #Make collection month and date a factor for models

#Per enclosure 
alpha_div_meta_clean_H21 <- alpha_div_meta_clean%>%
  filter(Enclosure =="H21")
alpha_div_meta_clean_P1 <- alpha_div_meta_clean%>%
  filter(Enclosure =="P1")

#Overall linear model
lm_model_copper_shannon <- lm(Shannon ~ Copper_level_mg_L * Enclosure, 
                                  data = alpha_div_meta_clean)
summary(lm_model_copper_shannon) ##Marginal effect of enclosure. Copper and Enclosure Interaction not significant 

#Linear model H21
lm_model_copper_shannon_H21 <- lm(Shannon ~ Copper_level_mg_L, 
                                  data = alpha_div_meta_clean_H21)
summary(lm_model_copper_shannon_H21) #NS 
summary(lm_model_copper_shannon_H21)$r.squared #0.011

##LOESS model H21
loess_model_copper_shannon_H21 <- loess(Shannon ~ Copper_level_mg_L, 
                     data = alpha_div_meta_clean_H21)
pred_copper_shannon_H21 <- predict(loess_model_copper_shannon_H21)
R2_copper_shannon_H21 <- 1 - sum((alpha_div_meta_clean_H21$Shannon - pred_copper_shannon_H21)^2) / 
  sum((alpha_div_meta_clean_H21$Shannon - mean(alpha_div_meta_clean_H21$Shannon))^2)
R2_copper_shannon_H21 ##0.08

#Linear model P1
lm_model_copper_shannon_P1 <- lm(Shannon ~ Copper_level_mg_L, 
                                  data = alpha_div_meta_clean_P1)
summary(lm_model_copper_shannon_P1) #NS
summary(lm_model_copper_shannon_P1)$r.squared #0.018

##LOESS model P1
loess_model_copper_shannon_P1 <- loess(Shannon ~ Copper_level_mg_L, 
                                        data = alpha_div_meta_clean_P1)
pred_copper_shannon_P1 <- predict(loess_model_copper_shannon_P1)
R2_copper_shannon_P1 <- 1 - sum((alpha_div_meta_clean_P1$Shannon - pred_copper_shannon_P1)^2) / 
  sum((alpha_div_meta_clean_P1$Shannon - mean(alpha_div_meta_clean_P1$Shannon))^2)
R2_copper_shannon_P1 ##0.04


###GAM Generalized Additive Model####
#The “additive” part means the effects of predictors are added together, but each effect can be nonlinear#
##With smooth: 
#Which distribution?
ggplot(alpha_div_meta_clean, aes(x= Shannon))+
  geom_histogram()  #Continuous, strictly positive, left-skewed responses

#Reflect to make right skewed?
# Choose a constant greater than the maximum value, add 1 or any other nuumber. Then from that, substract values
alpha_div_meta_clean_shannon_ref <- alpha_div_meta_clean %>%
  mutate(Shannon_reflected = (max(Shannon) + 1) - Shannon)

ggplot(alpha_div_meta_clean_shannon_ref, aes(x= Shannon_reflected))+
  geom_histogram() #Continuous, strictly positive, right-skewed responses: GAMMA


###
# Creates separate smooths for Copper_level_mg_L for each level of Enclosure, model will fit two different smooth
# curves for Shannon vs Copper_level_mg_L, one for each enclosure. Also, + Enlcosure tests for the main effect of 
# enclosure independent of Copper_level_mg_L (baseline differences between enclosures)
# Gamma-distributed response with an identity link, so the expected value of Shannon is 
# modeled directly as a sum of smooths and parametric terms. 
# Using a Gamma distribution (positive, right-skewed data)
# The identity link means the model predicts the response directly (no log or inverse transformation).
# Default method for estimating smooth is GCV (Generalized Cross- Validation)
###
gam_model <- gam(Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
                 family = Gamma(link="identity"),
                 data = alpha_div_meta_clean_shannon_ref)
summary(gam_model)

###
#Method REML uses restricted maximum likelihood to estimate the smoothness of each term. 
#Tends to avoid overfitting better than GCV 
####Overall model- with interaction #####
gam_model_interaction <- gam(Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
                 data = alpha_div_meta_clean_shannon_ref,
                 family = Gamma(link="identity"),
                 method = "REML")
summary(gam_model_interaction)
gam.check(gam_model_interaction)

####Overall model- without interaction #####
gam_model_no_interaction <- gam(Shannon_reflected ~ s(Copper_level_mg_L) + Enclosure,
                 data = alpha_div_meta_clean_shannon_ref,
                 family = Gamma(link="identity"),
                 method = "REML")
summary(gam_model_no_interaction)
gam.check(gam_model_no_interaction)

anova(gam_model_no_interaction, gam_model_interaction) ##No difference between the one with the interaction and te one without. So no difference 
##in effect of Copper on Shannons between enclosures (the curves don’t need to be different)

####P1 ######
gam_model_P1 <- gam(Shannon_reflected ~ s(Copper_level_mg_L),
                 data = alpha_div_meta_clean_shannon_ref%>%filter(Enclosure == "P1"),
                 family = Gamma(link="identity"),
                 method = "REML")
summary(gam_model_P1)

####H21 #####
gam_model_H21 <- gam(Shannon_reflected ~ s(Copper_level_mg_L),
                 data = alpha_div_meta_clean_shannon_ref%>%filter(Enclosure == "H21"),
                 family = Gamma(link="identity"),
                 method = "REML")
summary(gam_model_H21)

####Gam - Collection_Date as random#####
#Including interaction
gamm_model_interaction_random <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
  random = list(Collection_Month = ~1),
  family = Gamma(link="identity"),
  data = alpha_div_meta_clean_shannon_ref
)
summary(gamm_model_interaction_random$gam) #NS
summary(gamm_model_interaction_random$lme)#NS

#Not including interaction
gamm_model_random <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L) + Enclosure,
  random = list(Collection_Month = ~1),
  family = Gamma(link="identity"),
  data = alpha_div_meta_clean_shannon_ref
)
summary(gamm_model_random$gam)#NS
summary(gamm_model_random$lme)

anova(gamm_model_interaction_random$lme, gamm_model_random$lme) #no difference
anova(gamm_model_interaction_random$gam, gamm_model_random$gam) #adding separate smooths for each enclosure didn’t improve the model fit enough to be significant.

#####P1######
gamm_model_P1_random <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L),
  random = list(Collection_Month = ~1),
  family = Gamma(link="identity"),
  data = alpha_div_meta_clean_shannon_ref%>%filter(Enclosure == "P1")
)
summary(gamm_model_P1_random$gam)

#####H21#####
gamm_model_H21_random <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L),
  random = list(Collection_Month = ~1),
  family = Gamma(link="identity"),
  data = alpha_div_meta_clean_shannon_ref%>%filter(Enclosure == "H21")
)
summary(gamm_model_H21_random$gam)


####Gam -  Collection_Date as fixed#####
gamm_model_date <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + s(as.numeric(Collection_Month)) + Enclosure,
  family = Gamma(link = "identity"),  # log link is more stable for Gamma
  data = alpha_div_meta_clean_shannon_ref)
  
summary(gamm_model_date$gam)

###GEE MODEL #####
#install.packages("geepack")
library(geepack)

#Arrange by date order
alpha_div_meta_clean_2 <- alpha_div_meta_clean %>%
  mutate(Enclosure = factor(Enclosure, levels = c("H21", "P1")),
         Collection_Date = as.Date(Collection_Date))%>% #Enclosure needs to be a factor for GEE 
  arrange(Enclosure, Collection_Date) #row 1 = earliest date, row 2 = next date, etc

alpha_div_meta_clean_shannon_ref_2 <- alpha_div_meta_clean_shannon_ref %>%
  mutate(Enclosure = factor(Enclosure, levels = c("H21", "P1")),
                            Collection_Date = as.Date(Collection_Date))%>% #Enclosure needs to be a factor for GEE 
  arrange(Enclosure, Collection_Date) 

##GEE will then implicitly use row order as the time order for AR1
gee_model <- geeglm(
  Shannon_reflected ~ Copper_level_mg_L * Enclosure,
  data = alpha_div_meta_clean_shannon_ref_2,
  id = Enclosure,                  # repeated measures within enclosure
  family = Gamma(link = "log"),    # ensures positive fitted values
  corstr = "ar1"                   # AR1 correlation for time-ordered data
)
summary(gee_model) #Small alpha indicates that repeated measurements are almost independent, 
#so AR1 doesn’t have much impact on standard errors here

#Shannon as left skewed
gee_model <- geeglm(
  Shannon~ Copper_level_mg_L * Enclosure,
  data = alpha_div_meta_clean_2,
  id = Enclosure,                  # repeated measures within enclosure
  family = Gamma(link = "log"),    # ensures positive fitted values
  corstr = "ar1"                   # AR1 correlation for time-ordered data
)
summary(gee_model)


###LME#####
##Collection_month: trying to accounts for repeated measurements or clustering by month: each month may have its own baseline Shannon diversity.
#Tells me the average effect of copper, enclosure, and their interaction on Shannon diversity across all months.
#This assumes a relationship that is aprox linear!
model_lme <- lmer(Shannon ~ Copper_level_mg_L * Enclosure + (1 | Collection_Date),
                  data = alpha_div_meta_clean)
summary(model_lme) #Only marginal effect of exposure, not of copper or interaction

##Rm correlation#####
rmcorr(Enclosure, Copper_level_mg_L, Shannon, data = alpha_div_meta_clean) #Slight trend of Shannon decreasing with Higher copper levels 
#The grouping factor (Enclosure) is not being “tested”.
#It’s used to control for non-independence in repeated measures.
rmcorr(Enclosure, A, Shannon, data = alpha_div_meta_clean) 

##Spearman correlation#####
cor.test(x = alpha_div_meta_clean_H21$Shannon, 
         y = alpha_div_meta_clean_H21$Copper_level_mg_L, 
         method = 'spearman')

cor.test(x = alpha_div_meta_clean_P1$Shannon, 
         y = alpha_div_meta_clean_P1$Copper_level_mg_L, 
         method = 'spearman')

## Partial correlation of Shannon vs Copper, controlling for Date#####
###Turning date into numeric to add to partial correlation 
alpha_div_meta_clean_H21 <- alpha_div_meta_clean_H21 %>%
  mutate(
    Collection_Date = as.Date(Collection_Date),  # make sure it's a Date
    Date_num = as.numeric(Collection_Date - min(Collection_Date))) #Calculating dates from date #1

alpha_div_meta_clean_P1 <- alpha_div_meta_clean_P1 %>%
  mutate(
    Collection_Date = as.Date(Collection_Date),  # make sure it's a Date
    Date_num = as.numeric(Collection_Date - min(Collection_Date))) #Calculating dates from date #1
# alpha_div_meta_clean_P1$Date_num <- as.numeric(as.Date(alpha_div_meta_clean_P1$Collection_Date))
# alpha_div_meta_clean_H21$Date_num <- as.numeric(as.Date(alpha_div_meta_clean_H21$Collection_Date))

#H21
pcor.test(x = alpha_div_meta_clean_H21$Shannon,
          y = alpha_div_meta_clean_H21$Copper_level_mg_L,
          z = alpha_div_meta_clean_H21["Date_num"],
          method = "spearman")


#P1
pcor.test(alpha_div_meta_clean_P1$Shannon,
          alpha_div_meta_clean_P1$Copper_level_mg_L,
          alpha_div_meta_clean_P1["Date_num"],
          method = "spearman")


##MODEL
# #Collection_date as factor to include as random effect
# alpha_div_meta2 <- alpha_div_meta %>%
#   #mutate(Collection_Date=as.factor(Collection_Date))
# alpha_div_meta_P1 <- alpha_div_meta %>%
#   #mutate(Collection_Date=as.factor(Collection_Date))%>%
#   filter(Enclosure == "P1")
# alpha_div_meta_H21 <- alpha_div_meta %>%
#   #mutate(Collection_Date=as.factor(Collection_Date))%>%
#   filter(Enclosure == "H21")
# 
# library(nlme)
# 
# alpha_div_meta_clean <- alpha_div_meta %>%
#   filter(!is.na(Shannon), 
#          !is.na(Copper_level_mg_L), 
#          !is.na(Collection_Date))
# 
# alpha_div_meta_mean <- alpha_div_meta_clean %>%
#   group_by(Collection_Date, Enclosure) %>%
#   summarise(Shannon = mean(Shannon, na.rm = TRUE),
#             Copper_level_mg_L = mean(Copper_level_mg_L, na.rm = TRUE))
# 
# 
# model <- lme(Shannon ~ Copper_level_mg_L, 
#              random = ~ 1 | Enclosure, 
#              correlation = corAR1(form = ~ Collection_Date),
#              data = alpha_div_meta_mean)
# 
# cor.test(alpha_div_meta_P1$Shannon, alpha_div_meta_P1$Copper_level_mg_L, method = "spearman")
# 
# 
# model <- lmer(Shannon ~ Copper_level_mg_L + (1 | Collection_Date), data = alpha_div_meta_P1)
# summary(model)
# 
# setdiff(alpha_div_meta$Copper_level_mg_L,alpha_div_meta$Copper_mg_L)
# 
# library(mgcv)   # for GAMs

# 
# H21_data_alpha_div <- alpha_div_meta_clean[alpha_div_meta_clean$Enclosure == "H21", ]
# P1_data_alpha_div  <- alpha_div_meta_clean[alpha_div_meta_clean$Enclosure == "P1", ]
# 
# 
# 
# lm_log <- lm(log(Shannon + 0.01) ~ Copper_level_mg_L, data = H21_data_alpha_div)
# summary(lm_log)
# qqnorm(residuals(lm_log))
# qqline(residuals(lm_log))
# shapiro.test(residuals(lm_log))
# 
# 
# 
# lm_model <- lm(Shannon ~ Copper_level_mg_L, data = P1_data_alpha_div)
# summary(lm_model)
# 
# 
# lm_model <- lm(Shannon ~ Copper_level_mg_L, data = alpha_div_meta_clean)
# shapiro.test(residuals(lm_model))
# 
# 
# # Gaussian
# glm_gaussian <- glm(Shannon ~ Copper_level_mg_L * Enclosure, family=gaussian, data=alpha_div_meta_clean)
# 
# # Gamma with log link
# glm_gamma <- glm(Shannon ~ Copper_level_mg_L * Enclosure, family=Gamma(link="log"), data=alpha_div_meta_clean)
# 
# # Quasi
# glm_quasi <- glm(Shannon ~ Copper_level_mg_L * Enclosure, family=quasi(link="log"), data=alpha_div_meta_clean)
# 
# par(mfrow=c(1,3))
# plot(residuals(glm_gaussian, type="deviance"), main="Gaussian Residuals")
# plot(residuals(glm_gamma, type="deviance"), main="Gamma Residuals")
# plot(residuals(glm_quasi, type="deviance"), main="Quasi Residuals")
# 
# AIC(glm_gaussian, glm_gamma)   # lower is better
# 
# 
# plot(alpha_div_meta_clean$Shannon, predict(glm_quasi, type="response"),
#      xlab="Observed Shannon", ylab="Predicted Shannon")
# abline(a=0, b=1, col="red")
# 
# plot(residuals(glm_gamma, type="deviance"))
# 
# 
# gam_model <- gam(Shannon ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
#                  family = Gamma(link="log"),
#                  data = alpha_div_meta_clean)
# 
# plot(gam_model, shade=TRUE, pages=1)




## NITRIFIERS #####
nitrifiers #1437 taxa and 216 samples
alpha_div1_nit <- phyloseq::estimate_richness(nitrifiers, 
                                              measures = c("Observed", "Shannon")) # richness, diversity
#Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
alpha_div2_nit <- microbiome::evenness(nitrifiers, index = "pielou", 
                                   zeroes = TRUE, 
                                   detection = 0) ##evenness

# combine metrics with metadata
alpha_div_nit <- cbind(alpha_div1_nit, alpha_div2_nit)
alpha_div_nit

##Add metadata
alpha_div_nit_meta <- cbind(nitrifiers@sam_data, 
                        alpha_div_nit) 
alpha_div_nit_meta # metadata and div metrics


#Pivot to long format 
alpha_div_nit_meta_long <- 
  alpha_div_nit_meta %>%
  pivot_longer(cols = c(Observed, Shannon, pielou),  
               names_to = "alpha_div_metric", 
               values_to = "alpha_div_value") 
alpha_div_nit_meta_long

##Factoring alpha div metrics
alpha_div_nit_meta_long$alpha_div_metric<- factor(alpha_div_nit_meta_long$alpha_div_metric, levels = c("Observed","pielou", "Shannon"))

###P1 vs H21#####
alpha_div_nit_P1vsH21 <- ggplot(alpha_div_nit_meta_long, 
                            aes(x = Enclosure, y= alpha_div_value, 
                                fill= Enclosure, colour = Enclosure)) +
  theme_bw() +
  labs(color = "System", fill = "System") +
  facet_wrap(~alpha_div_metric, 
             scales = "free",
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS"))) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.1) +
  scale_color_manual(values = enclosure.palette,labels = c("H21" = "Naive", "P1" = "Established")) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established")) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))+
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_nit_P1vsH21

###Time series#####
####Alpha div indexes####
alpha_div_nit_time <- ggplot(alpha_div_nit_meta_long, 
                         aes(x = Collection_Date, y= alpha_div_value)) +
  theme_bw() +
  labs(title= "ALPHA DIVERSITY") +
  facet_nested(alpha_div_metric ~Enclosure,
               scales = "free",
               #switch = "y", 
               labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                        "Shannon" = "DIVERSITY\n(SHANNON)",
                                        "pielou" = "EVENNESS",
                                        "P1" = "Established",
                                        "H21" = "Naive"))) +
  #geom_boxplot(alpha = 0.1) +
  geom_smooth(method="loess", se=TRUE) +
  geom_point(size = 3, shape = 18)+
  #geom_text(aes(label = SampleID), vjust = -0.5, size = 3, angle = 90)+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 8, angle = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_nit_time


####Alpha diversity indexes and Water quality over time#####
#Editing alpha diversity and metadata dataframe
alpha_div_nit_wq_time_long <- alpha_div_nit_meta %>%
  mutate(Copper_keep = Copper_level_mg_L) %>%   #keep a column as copper to be able to color the plot 
  pivot_longer(cols = c("Copper_level_mg_L",
                        "Temperature_F",
                        "Chlorine_mg_L", 
                        "pH_spu", 
                        "Ammonia_mg_L",
                        "Nitrite_mg_L",
                        "Nitrate_UV_mg_L", 
                        "Salinity_ppt",
                        "Alkalinity_mg_L", 
                        "Shannon",
                        "Observed",
                        "pielou"),
               names_to = "Index",
               values_to = "Index_value")%>%
  mutate(Index = factor(Index, levels = c(
    "Observed",
    "pielou",
    "Shannon",
    "Copper_level_mg_L",
    "Ammonia_mg_L",
    "Nitrite_mg_L",
    "Nitrate_UV_mg_L",
    "pH_spu", 
    "Salinity_ppt", 
    "Temperature_F",
    "Chlorine_mg_L", 
    "Alkalinity_mg_L"
  )))

####Adding vertical dashed line for when "phase" changes#####
metadata_phyloseq.bacteria.samples <- data.frame(phyloseq.bacteria.samples@sam_data)
#When do phases switch?
#Naive
metadata_phyloseq.bacteria.samples %>%
  filter(Enclosure == "H21") %>%
  group_by(Date_num_phase_naive) %>%
  summarise(
    min_naive = min(Date_num_naive, na.rm = TRUE),
    max_naive = max(Date_num_naive, na.rm = TRUE)
  )

#Established
metadata_phyloseq.bacteria.samples %>%
  filter(Enclosure == "P1") %>%
  group_by(Date_num_phase_established) %>%
  summarise(
    min_naive = min(Date_num_established, na.rm = TRUE),
    max_naive = max(Date_num_established, na.rm = TRUE)
  )

#Vertical lines to be places
line_breaks_phases <- data.frame(
  Enclosure = c("H21", "H21", "H21", "H21", "H21", "H21","H21",
                "H21", "H21", "H21", "H21", "H21", "H21","H21",
                "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1"), 
  Date_num = c("1", "27", "29", "38", "39", "51", "52", "81", "86", "108", "109", "135", "136", "146",
               "1", "53", "54", "65", "66", "104", "105", "169"))

#####Plot - Just  Shannon, Copper and Ammonia levels - Date number since start of sampling as as.factor#####
#Plot
alpha_div_nit_wq_date_num_factor <- ggplot(alpha_div_nit_wq_time_long%>%
                                             filter(Index %in% c("Copper_level_mg_L",
                                                                 "Shannon",
                                                                 "Ammonia_mg_L"
                                             )),
                                           aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black", 
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "NITRIFYING TAXA\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)")))+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme(
        legend.position = c(0.85, 1.6),
        legend.justification = c(0, 1),
        legend.title.position = "left",
        legend.direction = "horizontal",
        legend.text = element_text(size = 15, angle = 45,
                                   hjust = 0.5, vjust = 0.5),
        legend.title = element_text(size = 18, face = "bold"),
        strip.background = element_rect(fill = "black"),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text.y  = element_text(colour = "white", size = 30, 
                                     face = "bold", angle = 0),
        strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 18),
        axis.ticks.x = element_line(colour = "black", linewidth = 1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_nit_wq_date_num_factor


#####Plot - Shannon with Other Water quality measures- Date number since start of sampling as as.factor#####
alpha_div_nit_wq_date_num_factor_other_metadata <- ggplot(alpha_div_nit_wq_time_long%>%
                                             filter(Index %in% c("Copper_level_mg_L",
                                                                 "Shannon",
                                                                 "Ammonia_mg_L",
                                                                 #"Chlorine_mg_L", 
                                                                 #"Alkalinity_mg_L",
                                                                 "Temperature_F",
                                                                 "pH_spu",
                                                                 "Salinity_ppt"
                                             )),
                                           aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black", 
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "NITRIFYING TAXA\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)",
                                      "Temperature_F"= "Temperature (F)",
                                      "Salinity_ppt" = "Salinity (ppt)", 
                                      "pH_spu" = "pH (spu)"
                                      #"Chlorine_mg_L" = "Chlorine (mg/L)",
                                      #"Alkalinity_mg_L" = "Alkalinity (mg/L)"
             )))+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme(
    legend.position = c(0.85, 1.6),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.y  = element_text(colour = "white", size = 30, 
                                 face = "bold", angle = 0),
    strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.y = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_nit_wq_date_num_factor_other_metadata


#####Plot - Shannon with Other Water quality measures- Date number since start of sampling as as.factor#####
alpha_div_all_nit_wq_date_num_factor_other_metadata <- ggplot(alpha_div_nit_wq_time_long%>%
                                                            filter(Index %in% c("Copper_level_mg_L",
                                                                                "Shannon",
                                                                                "pielou", 
                                                                                "Observed",
                                                                                "Ammonia_mg_L",
                                                                                #"Chlorine_mg_L", 
                                                                                #"Alkalinity_mg_L",
                                                                                "Temperature_F",
                                                                                "pH_spu",
                                                                                "Salinity_ppt"
                                                            )),
                                                          aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black", 
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "NITRIFYING TAXA\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Observed" = "Richness", 
                                      "pielou" = "Evenness",
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)",
                                      "Temperature_F"= "Temperature (F)",
                                      "Salinity_ppt" = "Salinity (ppt)", 
                                      "pH_spu" = "pH (spu)"
                                      #"Chlorine_mg_L" = "Chlorine (mg/L)",
                                      #"Alkalinity_mg_L" = "Alkalinity (mg/L)"
             )))+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme(
    legend.position = c(0.85, 1.6),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.y  = element_text(colour = "white", size = 30, 
                                 face = "bold", angle = 0),
    strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.y = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_all_nit_wq_date_num_factor_other_metadata

###Linear models ####
#Have to remove missing values
alpha_div_nit_meta_clean <- alpha_div_nit_meta %>%
  filter(!is.na(Shannon),
         !is.na(Copper_level_mg_L),
         !is.na(Collection_Date))%>%
  arrange(Enclosure, Collection_Date)%>%
  mutate(Enclosure = factor(Enclosure, levels = c("P1", "H21")))
alpha_div_nit_meta_clean

#Per enclosure
alpha_div_nit_meta_clean_H21 <- alpha_div_nit_meta_clean%>%
  filter(Enclosure =="H21")%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("2023-10", 
                                                                "2023-11", "2023-12", 
                                                              "2024-01", "2024-02", 
                                                              "2024-03")))

alpha_div_nit_meta_clean_P1 <- alpha_div_nit_meta_clean%>%
  filter(Enclosure =="P1")%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("2023-11", "2023-12", 
                                                              "2024-01", "2024-02", 
                                                              "2024-03", "2024-04",
                                                              "2024-05")))

####Overall linear model#####
lm_model_copper_shannon_nit <- lm(Shannon ~ Copper_level_mg_L * Enclosure, 
                                      data = alpha_div_nit_meta_clean)
summary(lm_model_copper_shannon_nit) ##Effect of copper, enclosure, and their interaction

#Linear model H21
lm_model_copper_shannon_H21_nit <- lm(Shannon ~ Copper_level_mg_L, 
                                  data = alpha_div_nit_meta_clean_H21)
summary(lm_model_copper_shannon_H21_nit) #S (For every 1 mg/L increase in copper, the Shannon diversity decreases by around 3.75 units)
summary(lm_model_copper_shannon_H21_nit)$r.squared #0.0986
qqnorm(residuals(lm_model_copper_shannon_H21_nit))
qqline(residuals(lm_model_copper_shannon_H21_nit))

##LOESS model H21
loess_model_copper_shannon_H21_nit <- loess(Shannon ~ Copper_level_mg_L, 
                                        data = alpha_div_nit_meta_clean_H21)
summary(loess_model_copper_shannon_H21_nit)
pred_copper_shannon_H21_nit <- predict(loess_model_copper_shannon_H21_nit)
R2_copper_shannon_H21_nit <- 1 - sum((alpha_div_nit_meta_clean_H21$Shannon - pred_copper_shannon_H21_nit)^2) / 
  sum((alpha_div_nit_meta_clean_H21$Shannon - mean(alpha_div_nit_meta_clean_H21$Shannon))^2)
R2_copper_shannon_H21_nit ##0.241

#Linear model P1
lm_model_copper_shannon_P1_nit <- lm(Shannon ~ Copper_level_mg_L, 
                                 data = alpha_div_nit_meta_clean_P1)
summary(lm_model_copper_shannon_P1_nit) #S
summary(lm_model_copper_shannon_P1_nit)$r.squared #0.107
qqnorm(residuals(lm_model_copper_shannon_P1_nit))
qqline(residuals(lm_model_copper_shannon_P1_nit))

##LOESS model P1
loess_model_copper_shannon_P1_nit <- loess(Shannon ~ Copper_level_mg_L, 
                                       data = alpha_div_nit_meta_clean_P1)
pred_copper_shannon_P1_nit <- predict(loess_model_copper_shannon_P1_nit)
R2_copper_shannon_P1_nit <- 1 - sum((alpha_div_nit_meta_clean_P1$Shannon - pred_copper_shannon_P1_nit)^2) / 
  sum((alpha_div_nit_meta_clean_P1$Shannon - mean(alpha_div_nit_meta_clean_P1$Shannon))^2)
R2_copper_shannon_P1_nit ##0.116

###LME#####
##Collection_month: trying to accounts for repeated measurements or clustering by month: each month may have its own baseline Shannon diversity.
#Tells me the average effect of copper, enclosure, and their interaction on Shannon diversity across all months.
#This assumes a relationship that is aprox linear!
model_lme_nit <- lmer(Shannon ~ Copper_level_mg_L * Enclosure + (1 | Collection_Month),
                  data = alpha_div_nit_meta_clean)
summary(model_lme_nit) ##Effect of copper, enclosure, and their interaction

lmer(
  Shannon ~ Copper_level_mg_L + scale(Date_num) + (1 | Enclosure),
  data = alpha_div_nit_meta_clean
)

####NAIVE #######
model_lme_nit_H21 <- lmer(Shannon ~ Copper_level_mg_L + (1 | Collection_Month),
                      data = alpha_div_nit_meta_clean_H21)
summary(model_lme_nit_H21) ##Effect of copper 

alpha_div_nit_meta_clean_H21_filt <- alpha_div_nit_meta_clean_H21%>%
  filter(!Collection_Month %in% c("2024-03"))

model_lm_nit_H21 <- lm(Shannon ~ poly(Copper_level_mg_L,2, raw = TRUE) * Collection_Month,
                       alpha_div_nit_meta_clean_H21_filt)
summary(model_lm_nit_H21)

model_lm_nit_H21 <- lm(Shannon ~ Copper_level_mg_L + I(Copper_level_mg_L^2) +
    Collection_Month +
    Copper_level_mg_L:Collection_Month +
    I(Copper_level_mg_L^2):Collection_Month,
  data = alpha_div_nit_meta_clean_H21_filt)
summary(model_lm_nit_H21)

#Confidence Intervals
confint(model_lm_nit_H21)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole.
Anova(model_lm_nit_H21, type = "III") 

# Add fitted values to your data
alpha_div_nit_meta_clean_H21_filt$fitted <- fitted(model_lm_nit_H21)

# Plot actual vs fitted
ggplot(alpha_div_nit_meta_clean_H21_filt, aes(x = fitted, y = Shannon)) +
  geom_point(color = "steelblue", size = 2) +     # points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # 1:1 line
  theme_minimal() +
  labs(
    x = "Fitted Shannon",
    y = "Observed Shannon",
    title = "Fitted vs Observed Shannon Diversity")

gam_model_nit_H21 <- gam(Shannon ~ s(Copper_level_mg_L, by = Collection_Month) + Collection_Month,
                     data = alpha_div_nit_meta_clean_H21)
summary(gam_model_nit_H21) ##No effect of enclosure, but copper effect did vary between enclosures 
plot(gam_model_nit_H21, pages = 1, shade = TRUE)

##Trying to plot? Need copper levels to predict at 
# copper_seq <- seq(
#   min(alpha_div_nit_meta_clean_H21_filt$Copper_level_mg_L, na.rm = TRUE),
#   max(alpha_div_nit_meta_clean_H21_filt$Copper_level_mg_L, na.rm = TRUE),
#   length.out = 20)

emmeans_h21 <- emmeans(model_lm_nit_H21, ~ Copper_level_mg_L + I(Copper_level_mg_L^2) | Collection_Month,
                       #at = list(Copper_level_mg_L = c(0, 0.05, 0.1, 0.15, 0.2)))%>%  # values to predict
                       at = list(Copper_level_mg_L = unique(alpha_div_nit_meta_clean_H21_filt$Copper_level_mg_L)))
                       #at = list(Copper_level_mg_L = copper_seq))%>%
emmeans_h21


emmeans_h21_filt <- emmeans_h21 %>%
  data.frame()%>%
  filter(
    !(Collection_Month == "2023-10" & Copper_level_mg_L > 0.02),
    !(Collection_Month == "2023-11" & Copper_level_mg_L > 0.17),
    !(Collection_Month == "2023-12" & Copper_level_mg_L > 0.10),
    !(Collection_Month == "2024-01" & Copper_level_mg_L > 0.22),
    !(Collection_Month == "2024-02" & (Copper_level_mg_L > 0.22 | Copper_level_mg_L < 0.06)))

  
emmeans_h21_plot <- emmeans_h21_filt %>%
  ggplot(aes(x = Copper_level_mg_L, y = emmean, color = Collection_Month, fill = Collection_Month)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.2, 
              color = NA) +
  geom_point(aes(x = Copper_level_mg_L,
                  y = Shannon),
              alpha = 0.35,
              size = 3,
              data = alpha_div_nit_meta_clean_H21_filt) +#raw data
  #geom_point(size = 4, shape = 20) + ##emmean
  # geom_errorbar(aes(ymin = `lower.CL`, ymax = `upper.CL`), 
  #               #position = position_dodge(width = 0.5), 
  #               width = 0.03,
  #               linewidth = 0.3) + #error bars for confidence intervals
  facet_grid(~ Collection_Month,
             scales = "free", labeller = as_labeller(c( "2023-10" = "Oct 2023",
                                                      "2023-11" = "Nov 2023", 
                                                      "2023-12" = "Dec 2023",
                                                      "2024-01" = "Jan 2024",
                                                      "2024-02" = "Feb 2024")))+
  theme_bw() +
  labs(title= "NITRIFYING TAXA",
       y = "Shannon's Diversity", 
       x = "Copper Levels (mg/L)") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text = element_text(colour = "white", size = 28, face = "bold"),
    axis.title = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.text.x = element_text(colour = "black", size = 20, 
                               angle = 45, 
                               vjust = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 30, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) 
emmeans_h21_plot

#Predict values from lme model
predict_values_lme_H21 <- ggpredict(model_lm_nit_H21, terms = "Copper_level_mg_L")


#Plot
ggplot(predict_values_lme_H21, aes(x = x, y = predicted)) +
  geom_point(size = 3, shape = 18, color = "dodgerblue") + #Predicted values
  geom_point(data = alpha_div_nit_meta_clean_H21, 
             aes(x = Copper_level_mg_L, y = Shannon), size = 3, shape = 18,
             alpha = 0.5) + #Raw data
  geom_line(linewidth = 1, color = "dodgerblue",
            alpha = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "dodgerblue",
              alpha = 0.2) +
  theme_bw() +
  labs(title= "NAIVE SYSTEM\nNitrifiers", 
       y = "Predicted Shannon's Diversity",
       x = "Copper level (mg/L)") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12, angle = 45, vjust = 0.5),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(colour = "black", size = 20),
        axis.text.x = element_text(colour = "black", size = 20, angle = 45,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))
  
####ESTABLISHED #######
model_lme_nit_P1 <- lmer(Shannon ~ Copper_level_mg_L + (1 | Collection_Month),
                          data = alpha_div_nit_meta_clean_P1)
summary(model_lme_nit_P1) ##Effect of copper 

lm(Shannon ~ Copper_level_mg_L + Date_num,
   data = alpha_div_nit_meta_clean_P1)

#Predict values from lme model
predict_values_lme_P1 <- ggpredict(model_lme_nit_P1, terms = "Copper_level_mg_L")

#Plot
ggplot(predict_values_lme_P1, aes(x = x, y = predicted)) +
  geom_point(size = 3, shape = 18, color = "dodgerblue") + #Predicted values
  geom_point(data = alpha_div_nit_meta_clean_P1, 
             aes(x = Copper_level_mg_L, y = Shannon), size = 3, shape = 18,
             alpha = 0.5) + #Raw data
  geom_line(linewidth = 1, color = "dodgerblue",
            alpha = 0.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "dodgerblue",
              alpha = 0.2) +
  theme_bw() +
  labs(title= "ESTABLISHED SYSTEM\nNitrifiers", 
       y = "Predicted Shannon's Diversity",
       x = "Copper level (mg/L)") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12, angle = 45, vjust = 0.5),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(colour = "black", size = 20),
        axis.text.x = element_text(colour = "black", size = 20, angle = 45,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))



# #This assumes a relationship that is aprox linear!
# model_lme_nit_H21 <- lmer(Shannon ~ Copper_level_mg_L + (1 | Collection_Date),
#                       data = alpha_div_nit_meta_clean_H21)
# summary(model_lme_nit_H21) ##Effect of copper 
# 
# model_lme_nit_P1 <- lmer(Shannon ~ Copper_level_mg_L + (1 | Collection_Date),
#                           data = alpha_div_nit_meta_clean_P1)
# summary(model_lme_nit_P1) ##Effect of copper 

###GAM Generalized Additive Model####
#The “additive” part means the effects of predictors are added together, but each effect can be nonlinear#
##With smooth: 
#Which distribution?
ggplot(alpha_div_nit_meta_clean, aes(x= Shannon))+
  geom_histogram()  #Continuous, strictly positive, left-skewed responses

#Reflect to make right skewed?
# Choose a constant greater than the maximum value, add 1 or any other nuumber. Then from that, substract values
alpha_div_nit_meta_clean_shannon_ref <- alpha_div_nit_meta_clean %>%
  mutate(Shannon_reflected = (max(Shannon) + 1) - Shannon)

ggplot(alpha_div_nit_meta_clean_shannon_ref, aes(x= Shannon_reflected))+
  geom_histogram() #Continuous, strictly positive, right-skewed responses: GAMMA


###
# Creates separate smooths for Copper_level_mg_L for each level of Enclosure, model will fit two different smooth
# curves for Shannon vs Copper_level_mg_L, one for each enclosure. Also, + Enlcosure tests for the main effect of 
# enclosure independent of Copper_level_mg_L (baseline differences between enclosures)
# Gamma-distributed response with an identity link, so the expected value of Shannon is 
# modeled directly as a sum of smooths and parametric terms. 
# Using a Gamma distribution (positive, right-skewed data)
# The identity link means the model predicts the response directly (no log or inverse transformation).
# Default method for estimating smooth is GCV (Generalized Cross- Validation)
###
gam_model_nit <- gam(Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
                 family = Gamma(link="identity"),
                 data = alpha_div_nit_meta_clean_shannon_ref)
summary(gam_model_nit) ##No effect of enclosure, but copper effect did vary between enclosures 

###
#Method REML uses restricted maximum likelihood to estimate the smoothness of each term. 
#Tends to avoid overfitting better than GCV 
####Overall model- with interaction #####
gam_model_interaction_nit <- gam(Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
                             data = alpha_div_nit_meta_clean_shannon_ref,
                             family = Gamma(link="identity"),
                             method = "REML")
summary(gam_model_interaction_nit) #NS effect of enclosure, but copper effect did vary between enclosures 

####Overall model- without interaction #####
gam_model_no_interaction_nit <- gam(Shannon_reflected ~ s(Copper_level_mg_L) + Enclosure,
                                data = alpha_div_nit_meta_clean_shannon_ref,
                                family = Gamma(link="identity"),
                                method = "REML")
summary(gam_model_no_interaction_nit)
gam.check(gam_model_no_interaction_nit)

anova(gam_model_no_interaction_nit, gam_model_interaction_nit) ##Difference between the one with the interaction and the one without.
##Different effect of Copper on Shannons between enclosures (the curves are different)

####P1 ######
gam_model_P1_nit <- gam(Shannon_reflected ~ s(Copper_level_mg_L),
                    data = alpha_div_nit_meta_clean_shannon_ref%>%filter(Enclosure == "P1"),
                    family = Gamma(link="identity"),
                    method = "REML")
summary(gam_model_P1_nit) #Effect of copper on shannons in P1

####H21 #####
gam_model_H21 <- gam(Shannon_reflected ~ s(Copper_level_mg_L),
                     data = alpha_div_nit_meta_clean_shannon_ref%>%filter(Enclosure == "H21"),
                     family = Gamma(link="identity"),
                     method = "REML")
summary(gam_model_H21) #Effect of copper on shannons in H21

####Gam - Collection_Date as random#####
#Including interaction
gamm_model_interaction_random_nit <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L, by = Enclosure) + Enclosure,
  random = list(Collection_Month = ~1),
  family = Gamma(link="identity"),
  data = alpha_div_nit_meta_clean_shannon_ref
)
summary(gamm_model_interaction_random_nit$gam) #NS
summary(gamm_model_interaction_random_nit$lme)#NS

#Not including interaction
gamm_model_random_nit <- gamm(
  Shannon_reflected ~ s(Copper_level_mg_L) + Enclosure,
  random = list(Collection_Month = ~1),
  family = Gamma(link="identity"),
  data = alpha_div_nit_meta_clean_shannon_ref
)
summary(gamm_model_random_nit$gam)#NS
summary(gamm_model_random_nit$lme)

anova(gamm_model_interaction_random_nit$lme, gamm_model_random_nit$lme) #S difference
anova(gamm_model_interaction_random_nit$gam, gamm_model_random_nit$gam) #S difference. 



###GEE MODEL #####
#install.packages("geepack")
library(geepack)

#Arrange by date order
alpha_div_nit_meta_clean_2 <- alpha_div_nit_meta_clean %>%
  mutate(Enclosure = factor(Enclosure, levels = c("H21", "P1")),
         Collection_Date = as.Date(Collection_Date))%>% #Enclosure needs to be a factor for GEE 
  arrange(Enclosure, Collection_Date) #row 1 = earliest date, row 2 = next date, etc

alpha_div_nit_meta_clean_shannon_ref_2 <- alpha_div_nit_meta_clean_shannon_ref %>%
  mutate(Enclosure = factor(Enclosure, levels = c("H21", "P1")),
         Collection_Date = as.Date(Collection_Date))%>% #Enclosure needs to be a factor for GEE 
  arrange(Enclosure, Collection_Date) 

##GEE will then implicitly use row order as the time order for AR1
gee_model_nit_refl <- geeglm(
  Shannon_reflected ~ Copper_level_mg_L * Enclosure,
  data = alpha_div_nit_meta_clean_shannon_ref_2,
  id = Enclosure,                  # repeated measures within enclosure
  family = Gamma(link = "log"),    # ensures positive fitted values
  corstr = "ar1"                   # AR1 correlation for time-ordered data
)
summary(gee_model_nit_refl) #Small alpha indicates that repeated measurements are almost independent, 
#so AR1 doesn’t have much impact on standard errors here

#Shannon as left skewed
gee_model_nit_link <- geeglm(
  Shannon~ Copper_level_mg_L * Enclosure,
  data = alpha_div_nit_meta_clean_2,
  id = Enclosure,                  # repeated measures within enclosure
  family = Gamma(link = "log"),    # ensures positive fitted values
  corstr = "ar1"                   # AR1 correlation for time-ordered data
)
summary(gee_model_nit_link)


###Rm correlation#####
rmcorr(Enclosure, Copper_level_mg_L, Shannon, data = alpha_div_nit_meta_clean) 
#No trend. 

###Spearman correlation#####
cor.test(x = alpha_div_nit_meta_clean_H21$Shannon, 
         y = alpha_div_nit_meta_clean_H21$Copper_level_mg_L, 
         method = 'spearman') #Significant for naive one

cor.test(x = alpha_div_nit_meta_clean_P1$Shannon, 
         y = alpha_div_nit_meta_clean_P1$Copper_level_mg_L, 
         method = 'spearman') #Not significant for established 

### Partial correlation of Shannon vs Copper, controlling for Date#####
##Turning date into numeric to add to partial correlation 
alpha_div_nit_meta_clean_H21 <- alpha_div_nit_meta_clean_H21 %>%
  mutate(
    Collection_Date = as.Date(Collection_Date),  # make sure it's a Date
    Date_num = as.numeric(Collection_Date - min(Collection_Date))) #Calculating dates from date #1

alpha_div_nit_meta_clean_P1 <- alpha_div_nit_meta_clean_P1 %>%
  mutate(
    Collection_Date = as.Date(Collection_Date),  # make sure it's a Date
    Date_num = as.numeric(Collection_Date - min(Collection_Date))) #Calculating dates from date #1

####H21#######
H21_pcor <- pcor.test(x = alpha_div_nit_meta_clean_H21$Shannon,
          y = alpha_div_nit_meta_clean_H21$Copper_level_mg_L,
          z = alpha_div_nit_meta_clean_H21$Date_num,
          method = "pearson")

# residuals
res_shannon_H21 <- resid(lm(Shannon ~ Date_num, data = alpha_div_nit_meta_clean_H21))
res_copper_H21 <- resid(lm(Copper_level_mg_L ~ Date_num, data = alpha_div_nit_meta_clean_H21))


plot(res_copper_H21, res_shannon_H21,
     xlab = "Copper level (residuals)",
     ylab = "Shannon diversity (residuals)",
     main = paste0("Partial correlation: r = ", round(H21_pcor$estimate, 2)))
abline(lm(res_shannon_H21 ~ res_copper_H21), col = "red", lwd = 2)
summary(lm(res_shannon_H21 ~ res_copper_H21))

plot(rank(res_copper_H21), rank(res_shannon_H21),
     xlab = "Copper level (residuals)",
     ylab = "Shannon diversity (residuals)",
     main = paste0("Partial correlation: r = ", round(H21_pcor$estimate, 2)))
abline(lm(rank(res_shannon_H21) ~ rank(res_copper_H21)), col = "red")
summary(lm(rank(res_shannon_H21) ~ rank(res_copper_H21)))


####P1#######
P1_pcor <- pcor.test(x = alpha_div_nit_meta_clean_P1$Shannon,
                      y = alpha_div_nit_meta_clean_P1$Copper_level_mg_L,
                      z = alpha_div_nit_meta_clean_P1$Date_num,
                      method = "pearson")

# residuals
res_shannon_P1 <- resid(lm(Shannon ~ Date_num, data = alpha_div_nit_meta_clean_P1))
res_copper_P1 <- resid(lm(Copper_level_mg_L ~ Date_num, data = alpha_div_nit_meta_clean_P1))


plot(res_copper_P1, res_shannon_P1,
     xlab = "Copper level (residuals)",
     ylab = "Shannon diversity (residuals)",
     main = paste0("Partial correlation: r = ", round(P1_pcor$estimate, 2)))
abline(lm(res_shannon_P1 ~ res_copper_P1), col = "red", lwd = 2)
summary(lm(res_shannon_P1 ~ res_copper_P1))

plot(rank(res_copper_P1), rank(res_shannon_P1),
     xlab = "Copper level (residuals)",
     ylab = "Shannon diversity (residuals)",
     main = paste0("Partial correlation: r = ", round(P1_pcor$estimate, 2)))
abline(lm(rank(res_shannon_P1) ~ rank(res_copper_P1)), col = "red")
summary(lm(rank(res_shannon_P1) ~ rank(res_copper_P1)))



#RELATIVE ABUNDANCE######
#ALL TAXA######
## ORDER #####
phyloseq.bacteria.samples_order.ra #2405 taxa and 216 samples 

#Grouping the low abundance orders into one category
phyloseq.bacteria.samples.order.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.samples_order.ra, 
                                                                        "Enclosure", 
                                                                        level = "Order", 
                                                                        threshold = 0.5)
phyloseq.bacteria.samples.order.filt #24 orders over 0.5% mean RA
phyloseq.bacteria.samples.order.filt.melt <- psmelt(phyloseq.bacteria.samples.order.filt)%>%
  mutate(Order = factor(Order, 
                         levels = c(setdiff(Order, 
                                            unique(grep("Others", Order, value = TRUE))), 
                                    unique(grep("Others", Order, value = TRUE)))))##Factoring the Order column so that "Others.." is the last category
levels(phyloseq.bacteria.samples.order.filt.melt$Order) ##ok

##Create color palette
#order.filt.palette <- distinctColorPalette(length(unique(phyloseq.bacteria.samples.order.filt.melt$Order)))
order.filt.palette <- unname(alphabet2())
order.filt.palette <- unname(polychrome())

order.filt.palette <- c(
  "#6F00FF", "#BFFF00", "#FF6F61", "#014D4E", "#E97451",
  "#191970", "#DAA520", "#FF00AF", "#046307", "#71797E",
  "#FF7518", "#7851A9", "#00CED1", "#990000", "#B497BD",
  "#DFFF00", "#6A5ACD", "#4B2E2B", "#C2B280", "#00FFFF",
  "#6B8E23", "#FFB07C", "#0F52BA", "#FF0800", "#8A9A5B",
  "#7B68EE", "#FFBF00", "#36454F", "#2E8B8B", "#FF00FF"
)

order.filt.palette <- c(
  "#E6194B", # red
  "#3CB44B", # green
  "#FFE119", # yellow
  "#4363D8", # blue
  "#F58231", # orange
  "#911EB4", # purple
  "#46F0F0", # cyan
  "#F032E6", # magenta
  "#BCF60C", # lime
  "#FABEBE", # pink
  "#008080", # teal
  "#E6BEFF", # lavender
  "#9A6324", # brown
  #"#FFFAC8", # beige
  #"#800000", # maroon
  "#C04060",
  "#AAFFC3", # mint
  "#808000", # olive
  #"#FFD8B1", # apricot
  "#000075", # navy
  "#808080", # gray
  "#FF4500", # orange-red
  #"#3CB371", #sea green
  "#2A6F62",
  #"#00FF7F", # spring green
  "#00BFA5",
  #"#8B0000", # dark red
  "#FFD8B1", # apricot, 
  "#00BFFF", # deep sky blue
  "#FFD700", # gold
  "#4B0082", # indigo
  "#FFFAC8", # beige
  "#FF1493", # deep pink
  "#228B22", # forest green
  "#A52A2A", # brown-red
  "#00CED1"  # dark turquoise
)



order_filt_names <- unique(phyloseq.bacteria.samples.order.filt.melt$Order)# Create a named vector for the palette, where the names correspond to phlyum names
order_named_palette <- setNames((order.filt.palette)[1:length(order_filt_names)], order_filt_names)
order_named_palette$'Others <0.5% RA' <- "grey95"
order_named_palette$'Flavobacteriales' <-  "#63A184"
order_named_palette$'Rhodobacterales' <- "#E3B199"
order_named_palette$'Nitrososphaerales' <- "#951942"
# order_named_palette$'unclassified Bacteria' <- "darkred"
order_named_palette$'unclassified Alphaproteobacteria' <- "dodgerblue"
# order_named_palette$'Mycobacteriales' <- "#8B005D"

# ##Apply the function to obtain top orders (n=15)
# top_orders <- top_taxa_legend(phyloseq.bacteria.samples.order.filt.melt, 
#                               taxlevel = "Order", n = 15)
# top_orders

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_order_enclosures_overall_plot_datenum <- ggplot(phyloseq.bacteria.samples.order.filt.melt,
                                     aes(x=factor(Date_num), y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black", 
  #            alpha = 0.8) +
  # #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = order_named_palette,
                    #breaks = top_orders,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
        legend.position = "right",
        #legend.position = c(1.09, 0.5),  # x, y inside plot
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(0.6, "cm"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
        strip.text.x  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 20,
                                   vjust = 0.5, hjust = 0.5),
        axis.title = element_text(colour = "black", size = 22),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_order_enclosures_overall_plot_datenum

## FAMILY #####
phyloseq.bacteria.samples_family.ra #6755 families and 216 samples 

phyloseq.bacteria.samples.family.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.samples_family.ra, 
                                                                             "Enclosure", 
                                                                             level = "Family", 
                                                                             threshold = 0.5)
phyloseq.bacteria.samples.family.filt #27 families over 0.5% mean RA
phyloseq.bacteria.samples.family.filt.melt <- psmelt(phyloseq.bacteria.samples.family.filt)%>%
  mutate(Family = factor(Family, 
                        levels = c(setdiff(Family, 
                                           unique(grep("Others", Family, value = TRUE))), 
                                   unique(grep("Others", Family, value = TRUE)))))##Factoring the Family column so that "Others.." is the last category
levels(phyloseq.bacteria.samples.family.filt.melt$Family) ##ok

####SUPPLEMENTARY TABLE 5 #####
stable5 <-  phyloseq.bacteria.samples.family.filt.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Family, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  rename(Phase = Date_num_phase_abbrv, 
         `Mean Relative Abundance (%) Within Overall Microbial Community` = mean_abun,  
         `Standard Deviation (%) Within Overall Microbial Community` = sd_abun, 
         `Min Relative Abundance (%) Within Overall Microbial Community` = min_abun, 
         `Max Relative Abundance (%) Within Overall Microbial Community` = max_abun)
stable5

write_xlsx(stable5, 
           "SupplementaryTable5.xlsx")


##Create color palette - based on families within the same order
palette_family_level_df <- phyloseq.bacteria.samples.family.filt.melt %>% 
  arrange(Order, Family) %>%   # ensure consistent shading order
  group_by(Order) %>%
  mutate(
    base_color = order_named_palette[Order],
    shade = seq(0.02, 0.5, length.out = n()),  # avoid extremes
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_family_level_df

#Set up final palette
family_named_palette <- setNames(
  palette_family_level_df$color,
  palette_family_level_df$Family)
family_named_palette
family_named_palette$'Nitrosopumilaceae' <- "#951942"
family_named_palette$'Others <0.5% RA' <- "grey95"

##Apply the function to obtain top familys (n=15)
top_families <- top_taxa_legend(phyloseq.bacteria.samples.family.filt.melt, 
                              taxlevel = "Family", n = 20)
top_families

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_family_enclosures_overall_plot_datenum <- ggplot(phyloseq.bacteria.samples.family.filt.melt,
                                             aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = family_named_palette,
                    breaks = top_families,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.08, 0.5),  # x, y inside plot
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 17, face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_enclosures_overall_plot_datenum


###Plot RA of just Rhodobacteraceae between phases ##############
phyloseq.bacteria.samples_family.ra #6755 families, 216 samples

#Out of this overall communities object, select only Nitrosopumilaceae
phyloseq.bacteria.samples_family.ra.rhodobacteraceae <- subset_taxa(phyloseq.bacteria.samples_family.ra, 
                                                                     Family == "Rhodobacteraceae") 
phyloseq.bacteria.samples_family.ra.rhodobacteraceae <- subset_samples(phyloseq.bacteria.samples_family.ra.rhodobacteraceae, 
                                                                        sample_sums(phyloseq.bacteria.samples_family.ra.rhodobacteraceae) > 0)
phyloseq.bacteria.samples_family.ra.rhodobacteraceae #Rhodobacteraceae (1 taxa) in 216 samples 

#Melt to long format
phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt <- psmelt(phyloseq.bacteria.samples_family.ra.rhodobacteraceae)

#Going to make a single column for abbreviated phases, so I can color by this column 
phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt <-  phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt%>%
  mutate(System = ifelse(Enclosure == "H21", "Naive", "Established"), 
         System = factor(System, levels = c("Naive", "Established")))%>%
  mutate(Date_num_colors_systems = case_when(
    System == "Naive" & !is.na(Date_num_phase_naive_abbrv) ~
      paste0(Date_num_phase_naive_abbrv, "_", System),
    
    System == "Established" & !is.na(Date_num_phase_established_abbrv) ~
      paste0(Date_num_phase_established_abbrv, "_", System),
    
    TRUE ~ NA_character_
  ))%>%
  mutate(Date_num_colors_systems = factor(Date_num_colors_systems, 
                                          levels = c("L_Naive", 
                                                     "T1_Naive",
                                                     "E1_Naive", 
                                                     "T2_Naive", 
                                                     "T3_Naive", 
                                                     "E2_Naive",
                                                     "P_Naive", 
                                                     "L_Established", 
                                                     "T1_Established", 
                                                     "E_Established", 
                                                     "P_Established"
                                          )))

labels_rhodobacteraceae <- phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt %>%
  distinct(Date_num_colors_systems, Date_num_phase_abbrv) %>%
  mutate(Date_num_phase_abbrv = as.character(Date_num_phase_abbrv)) %>%
  tibble::deframe()

#Get wilcox stats
#Have to add another variable for each system to get their respective stats 
phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt_2 <- phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt%>%
  mutate(Date_num_phase_naive_abbrv_stats = paste0(Date_num_phase_naive_abbrv, "_Naive"), 
         Date_num_phase_naive_abbrv_stats = factor(Date_num_phase_naive_abbrv_stats, 
                                                   levels = c(
                                                     "L_Naive", 
                                                     "T1_Naive",
                                                     "E1_Naive", 
                                                     "T2_Naive", 
                                                     "T3_Naive", 
                                                     "E2_Naive",
                                                     "P_Naive"))
  )%>%
  mutate(Date_num_phase_established_abbrv_stats = paste0(Date_num_phase_established_abbrv, "_Established"), 
         Date_num_phase_established_abbrv_stats = factor(Date_num_phase_established_abbrv_stats, 
                                                         levels = c(
                                                           "L_Established", 
                                                           "T1_Established", 
                                                           "E_Established", 
                                                           "P_Established")))

#Get stats for naive system - only comparisons in order of time (L vs T1, T1 vs E1, E1 vs T2, T2 vs T3, T3 vs E2, E2 vs P)
stat_wilcox_test_rhodobacteraceae_H21 <- phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt_2 %>%
  filter(Enclosure == "H21")%>%
  wilcox_test(
    Abundance ~ Date_num_phase_naive_abbrv_stats,
    comparisons = list(
      c("L_Naive", "T1_Naive"),
      c("T1_Naive", "E1_Naive"),
      c("E1_Naive", "T2_Naive"),
      c("T2_Naive", "T3_Naive"),
      c("T3_Naive", "E2_Naive"),
      c("E2_Naive", "P_Naive")
    ),
    p.adjust.method = "BH"
  ) %>%
  filter(p.adj <= 0.05) 

stat_wilcox_test_rhodobacteraceae_df_H21 <- stat_wilcox_test_rhodobacteraceae_H21 %>%
  add_xy_position(x = "Date_num_phase_naive_abbrv_stats", 
                  step.increase = 0.10, 
                  scales = "free"
  )%>%
  mutate(System = "Naive")

#Get stats for established system - only comparisons in order of time (L vs T1, T1 vs E, E vs P)
stat_wilcox_test_rhodobacteraceae_P1 <- 
  phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt_2 %>%
  filter(Enclosure == "P1") %>%
  wilcox_test(
    Abundance ~ Date_num_phase_established_abbrv_stats,
    comparisons = list(
      c("L_Established", "T1_Established"),
      c("T1_Established", "E_Established"),
      c("E_Established", "P_Established")
    ),
    p.adjust.method = "BH"
  ) %>%
  filter(p.adj <= 0.05)

stat_wilcox_test_rhodobacteraceae_df_P1 <- stat_wilcox_test_rhodobacteraceae_P1 %>%
  add_xy_position(x = "Date_num_phase_established_abbrv_stats",
                  step.increase = 0.30)%>%
  mutate(System = "Established")

stat_wilcox_test_rhodobacteraceae_df_both_H21_P1 <- bind_rows(stat_wilcox_test_rhodobacteraceae_df_P1, 
                                                               stat_wilcox_test_rhodobacteraceae_df_H21)%>%
  mutate(System = factor(System, levels = c("Naive", "Established")))

#Plot
rhodobacteraceae.RA.phases.boxplot <- ggplot(phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt, 
                                              aes(x = Date_num_colors_systems, 
                                                  y= Abundance, 
                                                  color = Date_num_colors_systems, 
                                                  fill = Date_num_colors_systems)) +
  theme_bw() +
  facet_wrap(~System, 
             scales = "free")+
  labs(x = "Phase", 
       y = expression(italic(Rhodobacteraceae) ~ "RA (%)"),
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_x_discrete(labels = labels_rhodobacteraceae)+
  scale_fill_manual(values = phases_naive_established_palette)+
  scale_color_manual(values = phases_naive_established_palette)+
  scale_y_continuous(expand= c(0.03,0,0.1,0)) +
  theme(
    legend.position = "none",
    # legend.text = element_text(size = 20),
    # legend.title = element_text(size = 22, face = "bold"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 32),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    axis.title = element_text(size = 24, colour = "black", face = "bold"),
    axis.text = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.5))+
  stat_pvalue_manual(stat_wilcox_test_rhodobacteraceae_df_both_H21_P1,
                   label = "p.adj.signif",
                   hide.ns = TRUE,
                   tip.length = 0.01)
rhodobacteraceae.RA.phases.boxplot

ggsave("rhodobacteraceae_RA_between_phases.png", 
       rhodobacteraceae.RA.phases.boxplot, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 16)


#NITRIFIERS WITHIN THE OVERALL COMMUNITY#####
## FAMILY #######
phyloseq.bacteria.samples_family.ra #6755 families

##Which families are nitrifiers? 
nitrifiers.melt <- psmelt(nitrifiers)
unique(nitrifiers.melt$Family) #"Nitrosopumilaceae" , "Nitrospinaceae" , "Nitrospiraceae", "Nitrosomonadaceae",  "Chromatiaceae" , "Gallionellaceae", "Ectothiorhodospiraceae"
#"Nitrososphaeraceae","Nitrosocaldaceae"

#Out of this overall communities object, select only nitrifiers 
phyloseq.bacteria.samples_family.ra.nitrifiers <- subset_taxa(phyloseq.bacteria.samples_family.ra, 
                                                              Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                                                                Family == "Chromatiaceae" | # AOB (Nitrosococcus genus)
                                                                Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales
                                                                Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales
                                                                Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales
                                                                Genus == "JBMCYV01" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is JBMCYV01 sp04919959
                                                                Order == "Nitrosomirales" | #AOA
                                                                Order == "Candidatus Nitrosocaldales" | #AOA
                                                                Family == "Nitrospiraceae" | # NOB/Commamox; some!
                                                                Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                                                                Family == "Nitrobacteraceae" | #NOB
                                                                Family == "Gallionellaceae" | #NOB none
                                                                Family == "Nitrospinaceae") # NOB; some, plus a new one!
phyloseq.bacteria.samples_family.ra.nitrifiers <- subset_samples(phyloseq.bacteria.samples_family.ra.nitrifiers, 
                                                                 sample_sums(phyloseq.bacteria.samples_family.ra.nitrifiers) > 0)
phyloseq.bacteria.samples_family.ra.nitrifiers #9 nitrifying families in 216 samples 


#Melt to plot 
phyloseq.bacteria.samples_family.ra.nitrifiers.melt <- psmelt(phyloseq.bacteria.samples_family.ra.nitrifiers)

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.bacteria.samples_family.ra.nitrifiers.melt <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
      Family == "Nitrosomonadaceae" ~ "AOB",
      Family == "Chromatiaceae" ~ "AOB",
      Family == "Nitrosocaldaceae" ~ "AOA",
      Family == "Nitrosopumilaceae" ~ "AOA",
      Family == "Nitrososphaeraceae" ~ "AOA",
      Family == "Candidatus Nitrosocaldaceae" ~ "AOA",
      Family == "Nitrospiraceae" ~ "NOB",
      Family == "Ectothiorhodospiraceae" ~ "NOB",
      Family == "Nitrobacteraceae" ~ "NOB",
      Family == "Gallionellaceae" ~ "NOB",
      Family == "Nitrospinaceae" ~ "NOB",
      TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Family) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  # bright pink/red
  AOB = "#1E88E5",  # strong blue
  NOB = "#FFC107"   # vivid amber/yellow
)
#Make hues based on families within each ammonia-nitrite oxidizing group
palette_nitrifiers_family_df <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt %>% 
  distinct(Family, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Family) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    #shade = seq(-0.1, 0.1, length.out = n()),
    shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_family_df

#Set up final palette
palette_nitrifiers_family <- setNames(
  palette_nitrifiers_family_df$color,
  palette_nitrifiers_family_df$Family)
palette_nitrifiers_family

##Apply the function to obtain top orders (n=15)
# top_nitrifying_families <- top_taxa_legend(phyloseq.bacteria.samples_family.ra.nitrifiers.melt, 
#                                            n = 5)
# top_nitrifying_families <- c("Nitrosopumilaceae", #AOA
#                              "Chromatiaceae",#AOB
#                              "Nitrosomonadaceae",#AOB
#                              "Nitrobacteraceae",#NOB
#                              "Nitrospiraceae",#NOB
#                              "Ectothiorhodospiraceae")#NOB
# top_nitrifying_families

#Plot
RA_family_enclosures_nit_plot_datenum <- ggplot(phyloseq.bacteria.samples_family.ra.nitrifiers.melt,
                                            aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_family,
                    #breaks = top_nitrifying_families
                    ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.5),  # x, y inside plot
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22, face = "bold"),
    legend.key.size = unit(0.7, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_enclosures_nit_plot_datenum

###Plot RA of just Nitrosopumilaceae between phases ########
phyloseq.bacteria.samples_family.ra.nitrifiers #9 nitrifying taxa, 216 samples

#Out of this overall communities object, select only Nitrosopumilaceae
phyloseq.bacteria.samples_family.ra.nitrosopumilaceae <- subset_taxa(phyloseq.bacteria.samples_family.ra, 
                                                                    Family == "Nitrosopumilaceae") 
phyloseq.bacteria.samples_family.ra.nitrosopumilaceae <- subset_samples(phyloseq.bacteria.samples_family.ra.nitrosopumilaceae, 
                                                                       sample_sums(phyloseq.bacteria.samples_family.ra.nitrosopumilaceae) > 0)
phyloseq.bacteria.samples_family.ra.nitrosopumilaceae #Nitrosopumilaceae (1 taxa) in 216 samples 

#Melt to long format
phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt <- psmelt(phyloseq.bacteria.samples_family.ra.nitrosopumilaceae)

#Going to make a single column for abbreviated phases, so I can color by this column 
phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt <-  phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt%>%
  mutate(System = ifelse(Enclosure == "H21", "Naive", "Established"), 
         System = factor(System, levels = c("Naive", "Established")))%>%
  mutate(Date_num_colors_systems = case_when(
  System == "Naive" & !is.na(Date_num_phase_naive_abbrv) ~
    paste0(Date_num_phase_naive_abbrv, "_", System),
  
  System == "Established" & !is.na(Date_num_phase_established_abbrv) ~
    paste0(Date_num_phase_established_abbrv, "_", System),
  
  TRUE ~ NA_character_
))%>%
  mutate(Date_num_colors_systems = factor(Date_num_colors_systems, 
                                          levels = c("L_Naive", 
                                                     "T1_Naive",
                                                     "E1_Naive", 
                                                     "T2_Naive", 
                                                     "T3_Naive", 
                                                     "E2_Naive",
                                                     "P_Naive", 
                                                     "L_Established", 
                                                     "T1_Established", 
                                                     "E_Established", 
                                                     "P_Established"
                                                     )))

labels_nitrosopumilaceae <- phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt %>%
  distinct(Date_num_colors_systems, Date_num_phase_abbrv) %>%
  mutate(Date_num_phase_abbrv = as.character(Date_num_phase_abbrv)) %>%
  tibble::deframe()

#Get wilcox stats
#Have to add another variable for each system to get their respective stats 
phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt_2 <- phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt%>%
  mutate(Date_num_phase_naive_abbrv_stats = paste0(Date_num_phase_naive_abbrv, "_Naive"), 
         Date_num_phase_naive_abbrv_stats = factor(Date_num_phase_naive_abbrv_stats, 
                                                   levels = c(
                                                     "L_Naive", 
                                                     "T1_Naive",
                                                     "E1_Naive", 
                                                     "T2_Naive", 
                                                     "T3_Naive", 
                                                     "E2_Naive",
                                                     "P_Naive"))
         )%>%
  mutate(Date_num_phase_established_abbrv_stats = paste0(Date_num_phase_established_abbrv, "_Established"), 
         Date_num_phase_established_abbrv_stats = factor(Date_num_phase_established_abbrv_stats, 
                                                   levels = c(
                                                     "L_Established", 
                                                     "T1_Established", 
                                                     "E_Established", 
                                                     "P_Established")))
         
#Get stats for naive system - only comparisons in order of time (L vs T1, T1 vs E1, E1 vs T2, T2 vs T3, T3 vs E2, E2 vs P)
stat_wilcox_test_nitrosopumilaceae_H21 <- phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt_2 %>%
  filter(Enclosure == "H21")%>%
  wilcox_test(
    Abundance ~ Date_num_phase_naive_abbrv_stats,
    comparisons = list(
      c("L_Naive", "T1_Naive"),
      c("T1_Naive", "E1_Naive"),
      c("E1_Naive", "T2_Naive"),
      c("T2_Naive", "T3_Naive"),
      c("T3_Naive", "E2_Naive"),
      c("E2_Naive", "P_Naive")
    ),
    p.adjust.method = "BH"
  ) %>%
  filter(p.adj <= 0.05) 

stat_wilcox_test_nitrosopumilaceae_df_H21 <- stat_wilcox_test_nitrosopumilaceae_H21 %>%
  add_xy_position(x = "Date_num_phase_naive_abbrv_stats", 
                  step.increase = 0.04, 
                  scales = "free"
                  )%>%
  mutate(System = "Naive")

#Get stats for established system - only comparisons in order of time (L vs T1, T1 vs E, E vs P)
stat_wilcox_test_nitrosopumilaceae_P1 <- 
  phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt_2 %>%
  filter(Enclosure == "P1") %>%
  wilcox_test(
    Abundance ~ Date_num_phase_established_abbrv_stats,
    comparisons = list(
      c("L_Established", "T1_Established"),
      c("T1_Established", "E_Established"),
      c("E_Established", "P_Established")
    ),
    p.adjust.method = "BH"
  ) %>%
  filter(p.adj <= 0.05)

stat_wilcox_test_nitrosopumilaceae_df_P1 <- stat_wilcox_test_nitrosopumilaceae_P1 %>%
  add_xy_position(x = "Date_num_phase_established_abbrv_stats",
                  step.increase = 0.04)%>%
  mutate(System = "Established")

stat_wilcox_test_nitrosopumilaceae_df_both_H21_P1 <- bind_rows(stat_wilcox_test_nitrosopumilaceae_df_P1, 
                                                               stat_wilcox_test_nitrosopumilaceae_df_H21)%>%
  mutate(System = factor(System, levels = c("Naive", "Established")))

#Plot
nitrosopumilaceae.RA.phases.boxplot <- ggplot(phyloseq.bacteria.samples_family.ra.nitrosopumilaceae.melt, 
                                  aes(x = Date_num_colors_systems, 
                                      y= Abundance, 
                                      color = Date_num_colors_systems, 
                                      fill = Date_num_colors_systems)) +
  theme_bw() +
  facet_wrap(~System, 
             scales = "free")+
  labs(x = "Phase", 
       y = expression(italic(Nitrosopumilaceae) ~ "RA (%)"),
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_x_discrete(labels = labels_nitrosopumilaceae)+
  scale_fill_manual(values = phases_naive_established_palette)+
  scale_color_manual(values = phases_naive_established_palette)+
  scale_y_continuous(expand= c(0.03,0,0.1,0)) +
  theme(
    legend.position = "none",
    # legend.text = element_text(size = 20),
    # legend.title = element_text(size = 22, face = "bold"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 32),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    axis.title = element_text(size = 24, colour = "black", face = "bold"),
    axis.text = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.5))+
  stat_pvalue_manual(stat_wilcox_test_nitrosopumilaceae_df_both_H21_P1,
                     label = "p.adj.signif",
                     hide.ns = TRUE,
                     tip.length = 0.01)
nitrosopumilaceae.RA.phases.boxplot

ggsave("nitrosopumilaceae_RA_between_phases.png", 
       nitrosopumilaceae.RA.phases.boxplot, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 16)

## SPECIES#######
phyloseq.bacteria.samples_species.ra #160,112  species

#Out of this overall communities object, select only nitrifiers 
phyloseq.bacteria.samples_species.ra.nitrifiers <- subset_taxa(phyloseq.bacteria.samples_species.ra, 
                                                               Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                                                                 Family == "Chromatiaceae" | # AOB (Nitrosococcus genus)
                                                                 Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales
                                                                 Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales
                                                                 Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales
                                                                 Genus == "JBMCYV01" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is JBMCYV01 sp04919959
                                                                 Order == "Nitrosomirales" | #AOA
                                                                 Order == "Candidatus Nitrosocaldales" | #AOA
                                                                 Family == "Nitrospiraceae" | # NOB/Commamox; some!
                                                                 Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                                                                 Family == "Nitrobacteraceae" | #NOB
                                                                 Family == "Gallionellaceae" | #NOB none
                                                                 Family == "Nitrospinaceae") # NOB; some, plus a new one!
phyloseq.bacteria.samples_species.ra.nitrifiers <- subset_samples(phyloseq.bacteria.samples_species.ra.nitrifiers, 
                                                                       sample_sums(phyloseq.bacteria.samples_species.ra.nitrifiers) > 0)
phyloseq.bacteria.samples_species.ra.nitrifiers #1437 nitrifying species in 216 samples 


#Melt to plot 
phyloseq.bacteria.samples_species.ra.nitrifiers.melt <- psmelt(phyloseq.bacteria.samples_species.ra.nitrifiers)


##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.bacteria.samples_species.ra.nitrifiers.melt <- phyloseq.bacteria.samples_species.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Chromatiaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "Candidatus Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "Ectothiorhodospiraceae" ~ "NOB",
    Family == "Nitrobacteraceae" ~ "NOB",
    Family == "Gallionellaceae" ~ "NOB",
    Family == "Nitrospinaceae" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Species) %>%
  mutate(Species = factor(Species, levels = unique(Species)))

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  # bright pink/red
  AOB = "#1E88E5",  # strong blue
  NOB = "#FFC107"   # vivid amber/yellow
)
#Make hues based on species within each ammonia-nitrite oxidizing group
palette_nitrifiers_species_df <- phyloseq.bacteria.samples_species.ra.nitrifiers.melt %>% 
  distinct(Species, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Species) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    #shade = seq(-0.1, 0.1, length.out = n()),
    shade = seq(-0.5, 0.5, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()

#Set up final palette
palette_nitrifiers_species <- setNames(
  palette_nitrifiers_species_df$color,
  palette_nitrifiers_species_df$Species)
palette_nitrifiers_species

#Which are the most abundant species?
top_nitrifying_species <- top_taxa_legend(phyloseq.bacteria.samples_species.ra.nitrifiers.melt, 
                                          taxlevel = "Species", 
                                          n = 10)
top_nitrifying_species

#Factor 
top_nitrifying_species <- factor(top_nitrifying_species, levels = c("Nitrosopumilus maritimus",
                                                                   "Candidatus Nitrosopumilus sp. SW",
                                                                   "Nitrosopumilus piranensis", 
                                                                   "unclassified Nitrosopumilus",
                                                                   "Candidatus Nitrosopumilus koreensis",
                                                                   "Nitrosopumilus sp.",
                                                                   "Nitrosopumilus cobalaminigenes" ,
                                                                   "Nitrosopumilus adriaticus",
                                                                   "unclassified Bradyrhizobium", 
                                                                   "Nitrospira sp."))
#Plot
RA_species_enclosures_nit_plot_datenum <- ggplot(phyloseq.bacteria.samples_species.ra.nitrifiers.melt,
                                          aes(x=factor(Date_num), y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
  # #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_species,
                    breaks = top_nitrifying_species) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.5),  # x, y inside plot
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(0.7, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_species_enclosures_nit_plot_datenum


#NITRIFIERS ONLY######
any(sample_sums(nitrifiers)== 0) ## no samples with 0 OTUs

nitrifiers.ra <- transform_sample_counts(nitrifiers, 
                                         function(x) x/sum(x)*100) ##Relative abundance from normalized data

nitrifiers.ra #1437 taxa and 216 samples

###FAMILY#######
nitrifiers.ra.family <- tax_glom(nitrifiers.ra, taxrank = "Family", NArm = F)
nitrifiers.ra.family #9 families 
##Melt 
nitrifiers.ra.family.melt <- psmelt(nitrifiers.ra.family)

####SUPPLEMENTARY TABLE 7 #####
stable7.1 <-  phyloseq.bacteria.samples_family.ra.nitrifiers.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Family, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  rename(Phase = Date_num_phase_abbrv, 
         `Mean Relative Abundance (%) Within Overall Microbial Community` = mean_abun,  
         `Standard Deviation (%) Within Overall Microbial Community` = sd_abun, 
         `Min Relative Abundance (%) Within Overall Microbial Community` = min_abun, 
         `Max Relative Abundance (%) Within Overall Microbial Community` = max_abun)
stable7.1


##Which are the top most abundant taxa by group? 
stable7.2 <-  nitrifiers.ra.family.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Family, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  rename(Phase = Date_num_phase_abbrv, 
         `Mean Relative Abundance (%) Within Only Nitrifying Community` = mean_abun,  
         `Standard Deviation (%) Within Only Nitrifying Community` = sd_abun, 
         `Min Relative Abundance (%) Within Only Nitrifying Community` = min_abun, 
         `Max Relative Abundance (%) Within Only Nitrifying Community` = max_abun)
stable7.2

stable7 <- merge(stable7.1, stable7.2, by = c("System", "Family", "Phase"))%>%
  arrange(System,  Phase, desc(`Mean Relative Abundance (%) Within Only Nitrifying Community`))
  
write_xlsx(stable7, 
           "SupplementaryTable7.xlsx")

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
nitrifiers.ra.family.melt <- nitrifiers.ra.family.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Chromatiaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "Candidatus Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "Ectothiorhodospiraceae" ~ "NOB",
    Family == "Nitrobacteraceae" ~ "NOB",
    Family == "Gallionellaceae" ~ "NOB",
    Family == "Nitrospinaceae" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Family) %>%
  mutate(Family = factor(Family, levels = unique(Family)))

##Create color palette
#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  # bright pink/red
  AOB = "#1E88E5",  # strong blue
  NOB = "#FFC107"   # vivid amber/yellow
)
#Make hues based on families within each ammonia-nitrite oxidizing group
palette_nitrifiers_only_family_df <- nitrifiers.ra.family.melt %>% 
  distinct(Family, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Family) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    #shade = seq(-0.1, 0.1, length.out = n()),
    shade = seq(0.01, 0.6, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_only_family_df

#Set up final palette
palette_nitrifiers_only_family <- setNames(
  palette_nitrifiers_only_family_df$color,
  palette_nitrifiers_only_family_df$Family)
palette_nitrifiers_only_family

#Plot
RA_enclosures_nitrifiers_only_family.plot <- ggplot(nitrifiers.ra.family.melt,
                                               aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  geom_vline(data = line_breaks_phases,
             aes(xintercept = Date_num),
             linetype = "dashed",
             color = "black",
             alpha = 0.8) +
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_only_family) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.6),  # x, y inside plot
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(0.7, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_enclosures_nitrifiers_only_family.plot


#####SPECIES#######
nitrifiers.ra.species <- tax_glom(nitrifiers.ra, taxrank = "Species", NArm = F)
nitrifiers.ra.species #1437 species and 216 samples

#Merge low abun species
nitrifiers.ra.species.filt <- merge_low_abundance_grouped_ra(nitrifiers.ra.species, 
                                                             "Enclosure", 
                                                             level = "Species", threshold = 0.5)
nitrifiers.ra.species.filt #25 Species over 0.5% mean RA
nitrifiers.ra.species.filt.melt <- psmelt(nitrifiers.ra.species.filt)%>%
  mutate(Species = factor(Species, 
                          levels = c(setdiff(Species, 
                                             unique(grep("Others", Species, value = TRUE))), 
                                     unique(grep("Others", Species, value = TRUE)))))##Factoring the Species column so that "Others.." is the last category
levels(nitrifiers.ra.species.filt.melt$Species) ##ok

##Which are the top most abundant taxa by group? 
top_nitrifier_species_list <- nitrifiers.ra.species.filt.melt %>%
  group_by(Enclosure, Species) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(Enclosure,  desc(mean_abun))

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
nitrifiers.ra.species.filt.melt <- nitrifiers.ra.species.filt.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Chromatiaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "Candidatus Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "Ectothiorhodospiraceae" ~ "NOB",
    Family == "Nitrobacteraceae" ~ "NOB",
    Family == "Gallionellaceae" ~ "NOB",
    Family == "Nitrospinaceae" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Species) %>%
  mutate(Species = factor(Species, levels = unique(Species)))

####SUPPLEMENTARY TABLE 8 #####
##Which are the top most abundant taxa by group? 
stable8 <-  nitrifiers.ra.species.filt.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Species, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  rename(Phase = Date_num_phase_abbrv, 
         `Mean Relative Abundance (%) Within Only Nitrifying Community` = mean_abun,  
         `Standard Deviation Within Only Nitrifying Community` = sd_abun, 
         `Min Relative Abundance (%) Within Only Nitrifying Community` = min_abun, 
         `Max Relative Abundance (%) Within Only Nitrifying Community` = max_abun)
stable8
write_xlsx(stable8, 
           "SupplementaryTable8.xlsx")


#Top species for the legend
top_nitrifying_only_species <- top_taxa_legend(nitrifiers.ra.species.filt.melt, 
                                               taxlevel = "Species",
                                               n = 12)
top_nitrifying_only_species <- factor(top_nitrifying_only_species,
                                      levels = c("unclassified Nitrosopumilus",
                                                 "Nitrosopumilus sp028278985",
                                                 "Nitrosopumilus maritimus",
                                                 "unclassified Nitrosopumilaceae",
                                                 "Nitrosopumilus sp000746765",
                                                 "Nitrosopumilus sp006740685",
                                                 "Nitrosopumilus sp016125975",
                                                 "Nitrosopumilus sp003702525",
                                                 "unclassified Nitrospira_D", 
                                                 "Nitrospira_D sp029865225", 
                                                 "Nitrosomonas sp021604405",
                                                 "Others <0.5% RA" ))
#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  # bright pink/red
  AOB = "#1E88E5",  # strong blue
  NOB = "#FFC107"   # vivid amber/yellow
)

#Make hues based on species within each ammonia-nitrite oxidizing group
palette_nitrifiers_only_species_df <-
  nitrifiers.ra.species.filt.melt %>%
  distinct(Species, Nitrifying_group) %>%
  filter(!is.na(Nitrifying_group)) %>%
  group_by(Nitrifying_group) %>%
  arrange(Species) %>%
  mutate(
    color = {
      n_species <- n()
      group <- first(Nitrifying_group)
      
      # if (group == "AOA") {
      #   colorRampPalette(c(
      #     "#FCE4EC",  # very light pink
      #     "#880E4F",  # wine
      #     "#4A148C",   # deep purple
      #     "#FF4081", # hot pink
      #     "#C2185B"  # deep magenta
      #   ))(n_species)
        
        if (group == "AOA") {
        colorRampPalette(c(
          "#FCE4EC",  # very light pink
          "#880E4F",  # wine
          "#4A148C",   # deep purple
          "#FF4081", # hot pink
          #"#D50000",  # red (clear shift toward red)
          "#B71C1C", 
          "#EF9A9A",
          #"#E53935",
          "#880E4F"   # deep wine/magenta (dark endpoint)
        ))(n_species)
        
      } else if (group == "AOB") {
        colorRampPalette(c(
          "#81D4FA",
          "#1E88E5",
          "#1565C0",
          "#0D47A1"
        ))(n_species)
        
      } else {
        colorRampPalette(c(
          "#FFD54F",  # saturated yellow
          "#6D4C41",   # brown
          "#F57C00",  # orange
          "#FFB300",  # amber
          "#BF360C" # burnt orange
        ))(n_species)
      }
    }[row_number()]
  ) %>%
  ungroup()

#Set up final palette
palette_nitrifiers_only_species <- setNames(
  palette_nitrifiers_only_species_df$color,
  palette_nitrifiers_only_species_df$Species)
palette_nitrifiers_only_species
palette_nitrifiers_only_species$'Others <0.5% RA' <- "grey95"


#Plot
RA_enclosures_nitrifiers_only_species.plot <- ggplot(nitrifiers.ra.species.filt.melt, 
                                                       aes(x=factor(Date_num), y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_col(color = "black")+
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
  # #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(
    values = palette_nitrifiers_only_species,
    breaks = top_nitrifying_only_species,
    #labels = function(x) str_wrap(x, width = 20), 
    drop = FALSE
  )+
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    # legend.position = "right",
    legend.position = c(1.085, 0.5),  # x, y inside plot
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_enclosures_nitrifiers_only_species.plot


######Together with alpha div of nitrifiers as well as family level in the overall community#######
#Have to edit legend and x axis on the family_level plot
alpha_div_nit_wq_date_num_factor_other_metadata_2 <- alpha_div_nit_wq_date_num_factor_other_metadata +
  theme(
    legend.position = c(0.93, 1.3),  
    axis.text.x = element_blank(),
        axis.title.x = element_blank())
RA_family_enclosures_nit_plot_datenum_2 <- RA_family_enclosures_nit_plot_datenum + 
  theme(
    legend.position = c(1.08, 0.5),  
    legend.text = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 20)) 
RA_enclosures_nitrifiers_only_species.plot_2 <- RA_enclosures_nitrifiers_only_species.plot +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 20), 
    axis.title.x = element_text(size = 30))

#Final plot
figure_alpha_div_nit_species_ra_time_copper <-
  alpha_div_nit_wq_date_num_factor_other_metadata_2 /
  RA_family_enclosures_nit_plot_datenum_2  /
  RA_enclosures_nitrifiers_only_species.plot_2 +
  plot_layout(heights = c(1.4, 0.6, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_div_nit_species_ra_time_copper

#Saving figure
ggsave("figure_alpha_div_nit_species_ra_time_copper.png", 
       figure_alpha_div_nit_species_ra_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 18, 
       width = 26)




#Correlation of Nitrosopumilaceae with Copper levels#########
####H21########
phyloseq.bacteria.samples_family.ra.AOA.melt.H21 <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
  filter(Family == "Nitrosopumilaceae")%>%
  filter(Enclosure == "H21")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month),
         Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%  # convert to factor 
  filter(!Collection_Month %in% c("2023-09", "2024-03"))%>%
  filter(!is.na(Copper_level_mg_L),
         !is.na(Collection_Date))%>%
  mutate(Attempt = case_when(is.na(Attempt) ~ "No CuSO4 Addition",
                             Attempt == "1" ~ "PHASE 1",
                             Attempt == "2" ~ "PHASE 2",
                             Attempt == "3" ~ "PHASE 3",
                             TRUE ~ as.character(Attempt)))%>%
  mutate(Attempt = factor(Attempt, 
                          levels = c("PHASE 1", "PHASE 2", "PHASE 3", 
                                     "No CuSO4 Addition")))%>%
  mutate(Collection_Month = case_when(Collection_Month == "2023-10"~ "Oct-2023", 
                                      Collection_Month == "2023-11"~ "Nov-2023", 
                                      Collection_Month == "2023-12"~ "Dec-2023", 
                                      Collection_Month == "2024-01"~ "Jan-2024", 
                                      Collection_Month == "2024-02"~ "Feb-2024", 
                                      Collection_Month == "2024-03"~ "Mar-2024",
                                      TRUE ~ as.character(Collection_Month)))%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("Oct-2023",
                                                                "Nov-2023", 
                                                                "Dec-2023", 
                                                                "Jan-2024", 
                                                                "Feb-2024", 
                                                                "Mar-2024")))

copper_AOA_relationship_plot_H21 <- ggplot(phyloseq.bacteria.samples_family.ra.AOA.melt.H21,
       aes(x = Copper_level_mg_L, 
           y = Abundance)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_color_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Nitrosopumilaceae RA (%)",
       x = "Copper levels (mg/L)", 
       title = "NAIVE", 
       color = "Collection Month") +
  theme_bw() +
  geom_smooth(method="loess", 
              se=TRUE,
              color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(method = "spearman",
    label.x.npc = "left",
    label.y.npc = "bottom",
    size = 8) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
copper_AOA_relationship_plot_H21
ggsave("copper_AOA_relationship_plot_H21.png", 
       copper_AOA_relationship_plot_H21, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)


#####GAM MODEL#######
gam_model_nit_AOA_H21 <- gam(Abundance ~ s(Copper_level_mg_L, by = Collection_Month) + Collection_Month,
                         data = phyloseq.bacteria.samples_family.ra.AOA.melt.H21)
summary(gam_model_nit_AOA_H21) ##No effect of enclosure, but copper effect did vary between enclosures 
plot(gam_model_nit_AOA_H21, pages = 1, shade = TRUE)

#####Spearman correlation#####
cor.test(x = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Abundance, 
         y = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Copper_level_mg_L, 
         method = 'spearman') #Significant for naive one

H21_pcor_AOA <- pcor.test(x = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Abundance,
                          y = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Copper_level_mg_L,
                          z = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Date_num,
                          method = "pearson")

####P1######
phyloseq.bacteria.samples_family.ra.AOA.melt.P1 <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
  filter(Family == "Nitrosopumilaceae")%>%
  filter(Enclosure == "P1")%>%
  filter(Collection_Date > "2023-06-01")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month))%>%  # convert to factor 
  mutate(Attempt = case_when(is.na(Attempt) ~ "No CuSO4 Addition",
                             Attempt == "1" ~ "PHASE 1",
                             Attempt == "2" ~ "PHASE 2",
                             Attempt == "3" ~ "PHASE 3",
                             TRUE ~ as.character(Attempt)))%>%
  mutate(Attempt = factor(Attempt, 
                          levels = c("PHASE 1", "PHASE 2", "PHASE 3", 
                                     "No CuSO4 Addition")))%>%
  mutate(Collection_Month = case_when(Collection_Month == "2023-11"~ "Nov-2023", 
                                      Collection_Month == "2023-12"~ "Dec-2023", 
                                      Collection_Month == "2024-01"~ "Jan-2024", 
                                      Collection_Month == "2024-02"~ "Feb-2024", 
                                      Collection_Month == "2024-03"~ "Mar-2024",
                                      Collection_Month == "2024-04"~ "Apr-2024",
                                      TRUE ~ as.character(Collection_Month)))%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("Nov-2023", 
                                                                "Dec-2023", 
                                                                "Jan-2024", 
                                                                "Feb-2024", 
                                                                "Mar-2024", 
                                                                "Apr-2024")))
 
#Plot 
copper_AOA_relationship_plot_P1 <- ggplot(phyloseq.bacteria.samples_family.ra.AOA.melt.P1,
      aes(x = Copper_level_mg_L, 
          y = Abundance)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_color_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Nitrosopumilaceae RA (%)",
       x = "Copper levels (mg/L)", 
       title = "ESTABLISHED", 
       color = "Collection Month") +
  theme_bw() +
  geom_smooth(method="loess", 
              se=TRUE,
              color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 8) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
copper_AOA_relationship_plot_P1
ggsave("copper_AOA_relationship_plot_P1.png", 
       copper_AOA_relationship_plot_P1, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)

#Modelling established Nitrosopumilaceae (AOA)
phyloseq.bacteria.samples_family.ra.AOA.melt.P1.clean <- phyloseq.bacteria.samples_family.ra.AOA.melt.P1 %>%
  filter(!is.na(Copper_level_mg_L),
         !is.na(Collection_Date))%>%
  mutate(Enclosure = "P1")

model_lm_nit_AOA_P1 <- lm(Abundance ~ Copper_level_mg_L + I(Copper_level_mg_L^2) +
                         Collection_Month +
                         Copper_level_mg_L:Collection_Month +
                         I(Copper_level_mg_L^2):Collection_Month,
                       data = phyloseq.bacteria.samples_family.ra.AOA.melt.P1.clean)
summary(model_lm_nit_AOA_P1)

#Confidence Intervals
confint(model_lm_nit_AOA_P1)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole.
Anova(model_lm_nit_AOA_P1, type = "III") 

# Add fitted values to your data
phyloseq.bacteria.samples_family.ra.AOA.melt.P1.clean$fitted <- fitted(model_lm_nit_AOA_P1)

# Plot actual vs fitted
ggplot(phyloseq.bacteria.samples_family.ra.AOA.melt.P1.clean, aes(x = fitted, y = Abundance)) +
  geom_point(color = "steelblue", size = 2) +     # points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # 1:1 line
  theme_minimal() +
  labs(
    x = "Fitted Shannon",
    y = "Observed Shannon",
    title = "Fitted vs Observed Shannon Diversity")

gam_model_nit_H21 <- gam(Shannon ~ s(Copper_level_mg_L, by = Collection_Month) + Collection_Month,
                         data = alpha_div_nit_meta_clean_H21)
summary(gam_model_nit_H21) ##No effect of enclosure, but copper effect did vary between enclosures 
plot(gam_model_nit_H21, pages = 1, shade = TRUE)

##Trying to plot? Need copper levels to predict at 
# copper_seq <- seq(
#   min(alpha_div_nit_meta_clean_H21_filt$Copper_level_mg_L, na.rm = TRUE),
#   max(alpha_div_nit_meta_clean_H21_filt$Copper_level_mg_L, na.rm = TRUE),
#   length.out = 20)

emmeans_h21 <- emmeans(model_lm_nit_H21, ~ Copper_level_mg_L + I(Copper_level_mg_L^2) | Collection_Month,
                       #at = list(Copper_level_mg_L = c(0, 0.05, 0.1, 0.15, 0.2)))%>%  # values to predict
                       at = list(Copper_level_mg_L = unique(alpha_div_nit_meta_clean_H21_filt$Copper_level_mg_L)))
#at = list(Copper_level_mg_L = copper_seq))%>%
emmeans_h21


emmeans_h21_filt <- emmeans_h21 %>%
  data.frame()%>%
  filter(
    !(Collection_Month == "2023-10" & Copper_level_mg_L > 0.02),
    !(Collection_Month == "2023-11" & Copper_level_mg_L > 0.17),
    !(Collection_Month == "2023-12" & Copper_level_mg_L > 0.10),
    !(Collection_Month == "2024-01" & Copper_level_mg_L > 0.22),
    !(Collection_Month == "2024-02" & (Copper_level_mg_L > 0.22 | Copper_level_mg_L < 0.06)))


emmeans_h21_plot <- emmeans_h21_filt %>%
  ggplot(aes(x = Copper_level_mg_L, y = emmean, color = Collection_Month, fill = Collection_Month)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), 
              alpha = 0.2, 
              color = NA) +
  geom_point(aes(x = Copper_level_mg_L,
                 y = Shannon),
             alpha = 0.35,
             size = 3,
             data = alpha_div_nit_meta_clean_H21_filt) +#raw data
  #geom_point(size = 4, shape = 20) + ##emmean
  # geom_errorbar(aes(ymin = `lower.CL`, ymax = `upper.CL`), 
  #               #position = position_dodge(width = 0.5), 
  #               width = 0.03,
  #               linewidth = 0.3) + #error bars for confidence intervals
  facet_grid(~ Collection_Month,
             scales = "free", labeller = as_labeller(c( "2023-10" = "Oct 2023",
                                                        "2023-11" = "Nov 2023", 
                                                        "2023-12" = "Dec 2023",
                                                        "2024-01" = "Jan 2024",
                                                        "2024-02" = "Feb 2024")))+
  theme_bw() +
  labs(title= "NITRIFYING TAXA",
       y = "Shannon's Diversity", 
       x = "Copper Levels (mg/L)") +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text = element_text(colour = "white", size = 28, face = "bold"),
    axis.title = element_text(colour = "black", size = 20),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.text.x = element_text(colour = "black", size = 20, 
                               angle = 45, 
                               vjust = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 30, face = "bold", hjust = 0.5))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) 
emmeans_h21_plot


###Correlation of Nitrobacteraceae with Copper levels#########
####H21########
phyloseq.bacteria.samples_family.ra.NOB.melt.H21 <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
  filter(Family == "Nitrobacteraceae")%>% #"Nitrobacteraceae"
  filter(Enclosure == "H21")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month),
         Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%  # convert to factor 
  filter(!Collection_Month %in% c("2023-09", "2024-03"))%>%
  filter(!is.na(Copper_level_mg_L),
         !is.na(Collection_Date))%>%
  mutate(Attempt = case_when(is.na(Attempt) ~ "No CuSO4 Addition",
                             Attempt == "1" ~ "PHASE 1",
                             Attempt == "2" ~ "PHASE 2",
                             Attempt == "3" ~ "PHASE 3",
                             TRUE ~ as.character(Attempt)))%>%
  mutate(Attempt = factor(Attempt, 
                          levels = c("PHASE 1", "PHASE 2", "PHASE 3", 
                                     "No CuSO4 Addition")))%>%
  mutate(Collection_Month = case_when(Collection_Month == "2023-10"~ "Oct-2023", 
                                      Collection_Month == "2023-11"~ "Nov-2023", 
                                      Collection_Month == "2023-12"~ "Dec-2023", 
                                      Collection_Month == "2024-01"~ "Jan-2024", 
                                      Collection_Month == "2024-02"~ "Feb-2024", 
                                      Collection_Month == "2024-03"~ "Mar-2024",
                                      TRUE ~ as.character(Collection_Month)))%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("Oct-2023",
                                                                "Nov-2023", 
                                                                "Dec-2023", 
                                                                "Jan-2024", 
                                                                "Feb-2024", 
                                                                "Mar-2024")))
#Plot
copper_NOB_relationship_plot_H21 <- ggplot(phyloseq.bacteria.samples_family.ra.NOB.melt.H21,
                                           aes(x = Copper_level_mg_L, 
                                               y = Abundance)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_color_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Nitrobacteraceae RA (%)",
       x = "Copper levels (mg/L)", 
       title = "NAIVE", 
       color = "Collection Month") +
  theme_bw() +
  geom_smooth(method="loess", 
              se=TRUE,
              color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(method = "spearman",
           label.x.npc = 0.8,
           label.y.npc = 0.8,
           size = 8) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
copper_NOB_relationship_plot_H21
ggsave("copper_NOB_relationship_plot_H21.png", 
       copper_NOB_relationship_plot_H21, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)

####P1######
phyloseq.bacteria.samples_family.ra.NOB.melt.P1 <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
  filter(Family == "Nitrobacteraceae")%>% #"Nitrobacteraceae"
  filter(Enclosure == "P1")%>%
  filter(Collection_Date > "2023-06-01")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month))%>%  # convert to factor 
  mutate(Attempt = case_when(is.na(Attempt) ~ "No CuSO4 Addition",
                             Attempt == "1" ~ "PHASE 1",
                             Attempt == "2" ~ "PHASE 2",
                             Attempt == "3" ~ "PHASE 3",
                             TRUE ~ as.character(Attempt)))%>%
  mutate(Attempt = factor(Attempt, 
                          levels = c("PHASE 1", "PHASE 2", "PHASE 3", 
                                     "No CuSO4 Addition")))%>%
  mutate(Collection_Month = case_when(Collection_Month == "2023-11"~ "Nov-2023", 
                                      Collection_Month == "2023-12"~ "Dec-2023", 
                                      Collection_Month == "2024-01"~ "Jan-2024", 
                                      Collection_Month == "2024-02"~ "Feb-2024", 
                                      Collection_Month == "2024-03"~ "Mar-2024",
                                      Collection_Month == "2024-04"~ "Apr-2024",
                                      TRUE ~ as.character(Collection_Month)))%>%
  mutate(Collection_Month = factor(Collection_Month, levels = c("Nov-2023", 
                                                                "Dec-2023", 
                                                                "Jan-2024", 
                                                                "Feb-2024", 
                                                                "Mar-2024", 
                                                                "Apr-2024")))

#Plot
copper_NOB_relationship_plot_P1 <- ggplot(phyloseq.bacteria.samples_family.ra.NOB.melt.P1,
                                          aes(x = Copper_level_mg_L, 
                                              y = Abundance)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_color_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Nitrobacteraceae RA (%)",
       x = "Copper levels (mg/L)", 
       title = "ESTABLISHED", 
       color = "Collection Month") +
  theme_bw() +
  geom_smooth(method="loess", 
              se=TRUE,
              color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(
    method = "spearman",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 8) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
copper_NOB_relationship_plot_P1
ggsave("copper_NOB_relationship_plot_P1.png", 
       copper_NOB_relationship_plot_P1, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)



#BRAY CURTIS####
##OVERALL COMMUNITIES#########
###DISTANCES IN NAIVE (H21) SYSTEM####
phyloseq.bacteria.samples_H21 #This phyloseq object only has samples from the naive system
#First, need to order samples according to Collection Date
str(phyloseq.bacteria.samples_H21@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for naive system (H21)
phyloseq.bacteria.samples_H21_metadata <- data.frame(sample_data(phyloseq.bacteria.samples_H21))

# Order by Date
phyloseq.bacteria.samples_H21_metadata <- phyloseq.bacteria.samples_H21_metadata[order
                                                                                             (phyloseq.bacteria.samples_H21_metadata$Collection_Date), 
                                                                                             ]
phyloseq.bacteria.samples_H21_metadata$SampleID #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0302 (March 02, 2024)

#Extract OTU table for naive system (H21), then reorder according to date as in metadata
phyloseq.bacteria.samples_H21_otu <- data.frame(phyloseq.bacteria.samples_H21@otu_table)
phyloseq.bacteria.samples_H21_otu <- phyloseq.bacteria.samples_H21_otu %>%
  select(phyloseq.bacteria.samples_H21_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
phyloseq.bacteria.samples_H21.ordered <- phyloseq.bacteria.samples_H21
colnames(phyloseq.bacteria.samples_H21.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
phyloseq.bacteria.samples_H21.ordered@otu_table <- otu_table(phyloseq.bacteria.samples_H21_otu, taxa_are_rows = T)
phyloseq.bacteria.samples_H21.ordered@sam_data <- sample_data(phyloseq.bacteria.samples_H21_metadata) 
colnames(phyloseq.bacteria.samples_H21.ordered@otu_table) #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0302 (March 02, 2024)

#Now, normalize counts (Relative abundance)
phyloseq.bacteria.samples_H21.ordered_RA <- transform_sample_counts(phyloseq.bacteria.samples_H21.ordered, 
                                                                  function(x) x/sum(x)*100)
sample_sums(phyloseq.bacteria.samples_H21.ordered_RA)

#Calculate BC distances
phyloseq.bacteria.samples_H21.ordered_RA.bray <- vegdist(t(phyloseq.bacteria.samples_H21.ordered_RA@otu_table), method = "bray")
phyloseq.bacteria.samples_H21.ordered_RA.bray

#Make into square matrix
phyloseq.bacteria.samples_H21.ordered_RA.bray_matrix <- as.matrix(phyloseq.bacteria.samples_H21.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_H21_BC_dist <- expand.grid(
  i = 1:nrow(phyloseq.bacteria.samples_H21.ordered_RA.bray_matrix),
  j = 1:nrow(phyloseq.bacteria.samples_H21.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_H21_BC_dist) #8464 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_H21_BC_dist <- all_pairs_comparisons_H21_BC_dist[all_pairs_comparisons_H21_BC_dist$i < all_pairs_comparisons_H21_BC_dist$j, ]
nrow(all_pairs_comparisons_H21_BC_dist) #4186 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_H21_BC_dist$BC_dist <- mapply(function(i, j) {
  phyloseq.bacteria.samples_H21.ordered_RA.bray_matrix[i, j]
  }, all_pairs_comparisons_H21_BC_dist$i, all_pairs_comparisons_H21_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_H21_BC_dist$time_diff <- as.numeric(
  phyloseq.bacteria.samples_H21_metadata$Collection_Date[all_pairs_comparisons_H21_BC_dist$j] - phyloseq.bacteria.samples_H21_metadata$Collection_Date[all_pairs_comparisons_H21_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_H21_BC_dist <- all_pairs_comparisons_H21_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_H21 <- rownames(phyloseq.bacteria.samples_H21.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_H21_BC_dist <- all_pairs_comparisons_H21_BC_dist%>%
  mutate(sample_i = sample_ids_H21[all_pairs_comparisons_H21_BC_dist$i], 
         sample_j = sample_ids_H21[all_pairs_comparisons_H21_BC_dist$j])
all_pairs_comparisons_H21_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_H21 <- all_pairs_comparisons_H21_BC_dist %>%
  filter(j == i + 1)
first_diff_df_H21$sample_j <- factor(first_diff_df_H21$sample_j, levels = first_diff_df_H21$sample_j)

#Add some more metadata for plotting 
first_diff_df_H21 <- first_diff_df_H21 %>%
  mutate(
    Collection_Date_i = phyloseq.bacteria.samples_H21_metadata$Collection_Date[i],
    Collection_Date_j = phyloseq.bacteria.samples_H21_metadata$Collection_Date[j],
    Date_num_i = phyloseq.bacteria.samples_H21_metadata$Date_num[i],
    Date_num_j = phyloseq.bacteria.samples_H21_metadata$Date_num[j], 
    Enclosure = "H21"
    )


###DISTANCES IN ESTABLISHED (P1) SYSTEM####
phyloseq.bacteria.samples_P1 #This phyloseq object only has samples from the established system
#First, need to order samples according to Collection Date
str(phyloseq.bacteria.samples_P1@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for established system (P1)
phyloseq.bacteria.samples_P1_metadata <- data.frame(sample_data(phyloseq.bacteria.samples_P1))

# Order by Date
phyloseq.bacteria.samples_P1_metadata <- phyloseq.bacteria.samples_P1_metadata[order
                                                                                             (phyloseq.bacteria.samples_P1_metadata$Collection_Date), 
]
phyloseq.bacteria.samples_P1_metadata$SampleID #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0430 (April 30, 2024)

#Extract OTU table for established system (P1), then reorder according to date as in metadata
phyloseq.bacteria.samples_P1_otu <- data.frame(phyloseq.bacteria.samples_P1@otu_table)
phyloseq.bacteria.samples_P1_otu <- phyloseq.bacteria.samples_P1_otu %>%
  select(phyloseq.bacteria.samples_P1_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
phyloseq.bacteria.samples_P1.ordered <- phyloseq.bacteria.samples_P1
colnames(phyloseq.bacteria.samples_P1.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
phyloseq.bacteria.samples_P1.ordered@otu_table <- otu_table(phyloseq.bacteria.samples_P1_otu, taxa_are_rows = T)
phyloseq.bacteria.samples_P1.ordered@sam_data <- sample_data(phyloseq.bacteria.samples_P1_metadata) 
colnames(phyloseq.bacteria.samples_P1.ordered@otu_table) #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0430 (April 30, 2024)

#Now, normalize counts (Relative abundance)
phyloseq.bacteria.samples_P1.ordered_RA <- transform_sample_counts(phyloseq.bacteria.samples_P1.ordered, 
                                                                          function(x) x/sum(x)*100)
sample_sums(phyloseq.bacteria.samples_P1.ordered_RA)

#Calculate BC distances
phyloseq.bacteria.samples_P1.ordered_RA.bray <- vegdist(t(phyloseq.bacteria.samples_P1.ordered_RA@otu_table), method = "bray")
phyloseq.bacteria.samples_P1.ordered_RA.bray

#Make into square matrix
phyloseq.bacteria.samples_P1.ordered_RA.bray_matrix <- as.matrix(phyloseq.bacteria.samples_P1.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_P1_BC_dist <- expand.grid(
  i = 1:nrow(phyloseq.bacteria.samples_P1.ordered_RA.bray_matrix),
  j = 1:nrow(phyloseq.bacteria.samples_P1.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_P1_BC_dist) #15376 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_P1_BC_dist <- all_pairs_comparisons_P1_BC_dist[all_pairs_comparisons_P1_BC_dist$i < all_pairs_comparisons_P1_BC_dist$j, ]
nrow(all_pairs_comparisons_P1_BC_dist) #7626 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_P1_BC_dist$BC_dist <- mapply(function(i, j) {
  phyloseq.bacteria.samples_P1.ordered_RA.bray_matrix[i, j]
}, all_pairs_comparisons_P1_BC_dist$i, all_pairs_comparisons_P1_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_P1_BC_dist$time_diff <- as.numeric(
  phyloseq.bacteria.samples_P1_metadata$Collection_Date[all_pairs_comparisons_P1_BC_dist$j] - phyloseq.bacteria.samples_P1_metadata$Collection_Date[all_pairs_comparisons_P1_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_P1_BC_dist <- all_pairs_comparisons_P1_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_P1 <- rownames(phyloseq.bacteria.samples_P1.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_P1_BC_dist <- all_pairs_comparisons_P1_BC_dist%>%
  mutate(sample_i = sample_ids_P1[all_pairs_comparisons_P1_BC_dist$i], 
         sample_j = sample_ids_P1[all_pairs_comparisons_P1_BC_dist$j])
all_pairs_comparisons_P1_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_P1 <- all_pairs_comparisons_P1_BC_dist %>%
  filter(j == i + 1)
first_diff_df_P1$sample_j <- factor(first_diff_df_P1$sample_j, levels = first_diff_df_P1$sample_j)

#Add some more metadata for plotting 
first_diff_df_P1 <- first_diff_df_P1 %>%
  mutate(
    Collection_Date_i = phyloseq.bacteria.samples_P1_metadata$Collection_Date[i],
    Collection_Date_j = phyloseq.bacteria.samples_P1_metadata$Collection_Date[j],
    Date_num_i = phyloseq.bacteria.samples_P1_metadata$Date_num[i],
    Date_num_j = phyloseq.bacteria.samples_P1_metadata$Date_num[j], 
    Enclosure = "P1"
  )


###JOINING DISTANCES IN ESTABLISHED (P1) AND NAIVE (H21) SYSTEMS####
first_diff_df_all <- bind_rows(first_diff_df_H21, 
                               first_diff_df_P1)
first_diff_df_all

##Add a row for date_num_j to be "1" just so it is plotted, but it will have no values
first_diff_df_all <- first_diff_df_all %>%
  add_row(
    Date_num_j = 1,
    Enclosure = "H21",
    i = 0,
    j = 1)%>%
  add_row(
    Date_num_j = 1,
    Enclosure = "P1",
    i = 0,
    j = 1)
first_diff_df_all

#### PLOT#######
first_diff_BC_H21_P1_plot <- ggplot(first_diff_df_all, aes(x = factor(Date_num_j))) +
  #BC distance lines
  geom_line(aes(y = BC_dist, color = "Bray-Curtis (BC) distance", group = Enclosure)) +
  geom_point(aes(y = BC_dist, color = "Bray-Curtis (BC) distance")) +
  #Rate of BC distance change (BC distance divided by days between dates)
  geom_point(aes(y = bc_rate_of_change, color = "Rate of change (BC/time)")) +
  geom_line(aes(y = bc_rate_of_change, color = "Rate of change (BC/time)", group = Enclosure),
            linetype = "dashed") +
  #Facet
  facet_grid(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  #Legend
  scale_color_manual(
    name = "Metric",
    values = c(
      "Bray-Curtis (BC) distance" = "black",
      "Rate of change (BC/time)" = "blue"
    ),
    labels = function(x) stringr::str_wrap(x, width = 15)
  )+
  labs(x = "Date",
       y = "BC Distances") +
  theme_minimal() +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),

      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme_bw()+
  #Legend
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,      # bigger dots
        linewidth = 1.5  # thicker lines 
      )))+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.5),  # x, y inside plot
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.title.y = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
first_diff_BC_H21_P1_plot


##NITRIFYING COMMUNITIES#########
nitrifiers_all #nitrifying taxa in the overall community

###DISTANCES IN NAIVE (H21) SYSTEM####
nitrifiers_all_H21 <- subset_samples(nitrifiers_all, Enclosure == "H21")
nitrifiers_all_H21 <- prune_taxa(taxa_sums(nitrifiers_all_H21) > 0, 
                                                  nitrifiers_all_H21)
nitrifiers_all_H21 #1436 taxa and 92 samples from the naive system

#First, need to order samples according to Collection Date
str(nitrifiers_all_H21@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for naive system (H21)
nitrifiers_all_H21_metadata <- data.frame(sample_data(nitrifiers_all_H21))

# Order by Date
nitrifiers_all_H21_metadata <- nitrifiers_all_H21_metadata[order(nitrifiers_all_H21_metadata$Collection_Date), ]
nitrifiers_all_H21_metadata$SampleID #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0302 (March 02, 2024)

#Extract OTU table for naive system (H21), then reorder according to date as in metadata
nitrifiers_all_H21_otu <- data.frame(nitrifiers_all_H21@otu_table)
nitrifiers_all_H21_otu <- nitrifiers_all_H21_otu %>%
  select(nitrifiers_all_H21_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
nitrifiers_all_H21.ordered <- nitrifiers_all_H21
colnames(nitrifiers_all_H21.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
nitrifiers_all_H21.ordered@otu_table <- otu_table(nitrifiers_all_H21_otu, taxa_are_rows = T)
nitrifiers_all_H21.ordered@sam_data <- sample_data(nitrifiers_all_H21_metadata) 
colnames(nitrifiers_all_H21.ordered@otu_table) #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0302 (March 02, 2024)

#Now, normalize counts (Relative abundance)
nitrifiers_all_H21.ordered_RA <- transform_sample_counts(nitrifiers_all_H21.ordered, 
                                                                          function(x) x/sum(x)*100)
sample_sums(nitrifiers_all_H21.ordered_RA)

#Calculate BC distances
nitrifiers_all_H21.ordered_RA.bray <- vegdist(t(nitrifiers_all_H21.ordered_RA@otu_table), method = "bray")
nitrifiers_all_H21.ordered_RA.bray

#Make into square matrix
nitrifiers_all_H21.ordered_RA.bray_matrix <- as.matrix(nitrifiers_all_H21.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_H21_nit_BC_dist <- expand.grid(
  i = 1:nrow(nitrifiers_all_H21.ordered_RA.bray_matrix),
  j = 1:nrow(nitrifiers_all_H21.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_H21_nit_BC_dist) #8464 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_H21_nit_BC_dist <- all_pairs_comparisons_H21_nit_BC_dist[all_pairs_comparisons_H21_nit_BC_dist$i < all_pairs_comparisons_H21_nit_BC_dist$j, ]
nrow(all_pairs_comparisons_H21_nit_BC_dist) #4186 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_H21_nit_BC_dist$BC_dist <- mapply(function(i, j) {
  nitrifiers_all_H21.ordered_RA.bray_matrix[i, j]
}, all_pairs_comparisons_H21_nit_BC_dist$i, all_pairs_comparisons_H21_nit_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_H21_nit_BC_dist$time_diff <- as.numeric(
  nitrifiers_all_H21_metadata$Collection_Date[all_pairs_comparisons_H21_nit_BC_dist$j] - nitrifiers_all_H21_metadata$Collection_Date[all_pairs_comparisons_H21_nit_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_H21_nit_BC_dist <- all_pairs_comparisons_H21_nit_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_H21_nit <- rownames(nitrifiers_all_H21.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_H21_nit_BC_dist <- all_pairs_comparisons_H21_nit_BC_dist%>%
  mutate(sample_i = sample_ids_H21_nit[all_pairs_comparisons_H21_nit_BC_dist$i], 
         sample_j = sample_ids_H21_nit[all_pairs_comparisons_H21_nit_BC_dist$j])
all_pairs_comparisons_H21_nit_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_H21_nit <- all_pairs_comparisons_H21_nit_BC_dist %>%
  filter(j == i + 1)
first_diff_df_H21_nit$sample_j <- factor(first_diff_df_H21_nit$sample_j, levels = first_diff_df_H21_nit$sample_j)

#Add some more metadata for plotting 
first_diff_df_H21_nit <- first_diff_df_H21_nit %>%
  mutate(
    Collection_Date_i = nitrifiers_all_H21_metadata$Collection_Date[i],
    Collection_Date_j = nitrifiers_all_H21_metadata$Collection_Date[j],
    Date_num_i = nitrifiers_all_H21_metadata$Date_num[i],
    Date_num_j = nitrifiers_all_H21_metadata$Date_num[j], 
    Enclosure = "H21"
  )


###DISTANCES IN ESTABLISHED (P1) SYSTEM####
nitrifiers_all_P1 <- subset_samples(nitrifiers_all, Enclosure == "P1")
nitrifiers_all_P1 <- prune_taxa(taxa_sums(nitrifiers_all_P1) > 0, 
                                 nitrifiers_all_P1)
nitrifiers_all_P1 #1434 taxa and 124 samples from the established system

#First, need to order samples according to Collection Date
str(nitrifiers_all_P1@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for established system (P1)
nitrifiers_all_P1_metadata <- data.frame(sample_data(nitrifiers_all_P1))

# Order by Date
nitrifiers_all_P1_metadata <- nitrifiers_all_P1_metadata[order(nitrifiers_all_P1_metadata$Collection_Date), ]
nitrifiers_all_P1_metadata$SampleID #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0430 (April 30, 2024)

#Extract OTU table for established system (P1), then reorder according to date as in metadata
nitrifiers_all_P1_otu <- data.frame(nitrifiers_all_P1@otu_table)
nitrifiers_all_P1_otu <- nitrifiers_all_P1_otu %>%
  select(nitrifiers_all_P1_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
nitrifiers_all_P1.ordered <- nitrifiers_all_P1
colnames(nitrifiers_all_P1.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
nitrifiers_all_P1.ordered@otu_table <- otu_table(nitrifiers_all_P1_otu, taxa_are_rows = T)
nitrifiers_all_P1.ordered@sam_data <- sample_data(nitrifiers_all_P1_metadata) 
colnames(nitrifiers_all_P1.ordered@otu_table) #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0430 (April 30, 2024)

#Now, normalize counts (Relative abundance)
nitrifiers_all_P1.ordered_RA <- transform_sample_counts(nitrifiers_all_P1.ordered, 
                                                                         function(x) x/sum(x)*100)
sample_sums(nitrifiers_all_P1.ordered_RA)

#Calculate BC distances
nitrifiers_all_P1.ordered_RA.bray <- vegdist(t(nitrifiers_all_P1.ordered_RA@otu_table), method = "bray")
nitrifiers_all_P1.ordered_RA.bray

#Make into square matrix
nitrifiers_all_P1.ordered_RA.bray_matrix <- as.matrix(nitrifiers_all_P1.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_P1_nit_BC_dist <- expand.grid(
  i = 1:nrow(nitrifiers_all_P1.ordered_RA.bray_matrix),
  j = 1:nrow(nitrifiers_all_P1.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_P1_nit_BC_dist) #15376 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_P1_nit_BC_dist <- all_pairs_comparisons_P1_nit_BC_dist[all_pairs_comparisons_P1_nit_BC_dist$i < all_pairs_comparisons_P1_nit_BC_dist$j, ]
nrow(all_pairs_comparisons_P1_nit_BC_dist) #7626 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_P1_nit_BC_dist$BC_dist <- mapply(function(i, j) {
  nitrifiers_all_P1.ordered_RA.bray_matrix[i, j]
}, all_pairs_comparisons_P1_nit_BC_dist$i, all_pairs_comparisons_P1_nit_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_P1_nit_BC_dist$time_diff <- as.numeric(
  nitrifiers_all_P1_metadata$Collection_Date[all_pairs_comparisons_P1_nit_BC_dist$j] - nitrifiers_all_P1_metadata$Collection_Date[all_pairs_comparisons_P1_nit_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_P1_nit_BC_dist <- all_pairs_comparisons_P1_nit_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_P1_nit <- rownames(nitrifiers_all_P1.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_P1_nit_BC_dist <- all_pairs_comparisons_P1_nit_BC_dist%>%
  mutate(sample_i = sample_ids_P1_nit[all_pairs_comparisons_P1_nit_BC_dist$i], 
         sample_j = sample_ids_P1_nit[all_pairs_comparisons_P1_nit_BC_dist$j])
all_pairs_comparisons_P1_nit_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_P1_nit <- all_pairs_comparisons_P1_nit_BC_dist %>%
  filter(j == i + 1)
first_diff_df_P1_nit$sample_j <- factor(first_diff_df_P1_nit$sample_j, 
                                        levels = first_diff_df_P1_nit$sample_j)

#Add some more metadata for plotting 
first_diff_df_P1_nit <- first_diff_df_P1_nit %>%
  mutate(
    Collection_Date_i = nitrifiers_all_P1_metadata$Collection_Date[i],
    Collection_Date_j = nitrifiers_all_P1_metadata$Collection_Date[j],
    Date_num_i = nitrifiers_all_P1_metadata$Date_num[i],
    Date_num_j = nitrifiers_all_P1_metadata$Date_num[j], 
    Enclosure = "P1"
  )


###JOINING DISTANCES IN ESTABLISHED (P1) AND NAIVE (H21) SYSTEMS####
first_diff_df_nit <- bind_rows(first_diff_df_H21_nit, 
                               first_diff_df_P1_nit)
first_diff_df_nit

##Add a row for date_num_j to be "1" just so it is plotted, but it will have no values
first_diff_df_nit <- first_diff_df_nit %>%
  add_row(
    Date_num_j = 1,
    Enclosure = "H21",
    i = 0,
    j = 1)%>%
  add_row(
    Date_num_j = 1,
    Enclosure = "P1",
    i = 0,
    j = 1)
first_diff_df_nit

#### PLOT#######
first_diff_nit_BC_H21_P1_plot <- ggplot(first_diff_df_nit, aes(x = factor(Date_num_j))) +
  #BC distance lines
  geom_line(aes(y = BC_dist, color = "Bray-Curtis (BC) distance", group = Enclosure)) +
  geom_point(aes(y = BC_dist, color = "Bray-Curtis (BC) distance")) +
  #Rate of BC distance change (BC distance divided by days between dates)
  geom_point(aes(y = bc_rate_of_change, color = "Rate of change (BC/time)")) +
  geom_line(aes(y = bc_rate_of_change, color = "Rate of change (BC/time)", group = Enclosure),
            linetype = "dashed") +
  #Facet
  facet_grid(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  #Legend
  scale_color_manual(
    name = "Metric",
    values = c(
      "Bray-Curtis (BC) distance" = "black",
      "Rate of change (BC/time)" = "blue"
    ),
    labels = function(x) stringr::str_wrap(x, width = 15)
  )+
  labs(x = "Date",
       y = "BC Distances") +
  theme_minimal() +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme_bw()+
  #Legend
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,      # bigger dots
        linewidth = 1.5  # thicker lines 
      )))+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.5),  # x, y inside plot
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.title.y = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
first_diff_nit_BC_H21_P1_plot


###CALCULATE CHANGE POINT ########
#At what points does the behavior of microbiome change (as measured by BC-diff) shift?
#This is basically segmenting the time series into different dynamical regimes


####OVERALL COMMUNITIES#####
#####NAIVE SYSTEM#####
#Pull the pairwise BC distances. Vector of day-to-day change n community composition
first_diff_df_H21_BC <- first_diff_df_H21 %>%
  pull(BC_dist)

#Calculate Changepoint 
#Find points where the mean and/or variance of this time series changes
#Looking for changes from high BC-diff to low BC-diff (mean change) or steady to noisy (variance change)
#PELT Alogrithm (best set of change points across the entire series)
#BIC (Bayesian Information Criterion): Don’t add a change point unless it improves the model 
changepoint_overall_communities_h21 <- cpt.meanvar(first_diff_df_H21_BC, 
                                                   method = "PELT", 
                                                   penalty = "MBIC", 
                                                   #minseglen = 2
                                                   )

changepoint_overall_communities_h21 <- cpt.meanvar(first_diff_df_H21_BC,
                       method = "PELT",
                       penalty = "Manual",
                       pen.value = 0.5,
                       minseglen = 7)
#Pull the indeces where changepoints happen
changepoint_overall_communities_h21_indeces <- cpts(changepoint_overall_communities_h21)

#Filter the indeces in the main dataframe of day-to-day variances.
changepoint_overall_communities_h21_days <- first_diff_df_H21 %>%
  filter(i %in% changepoint_overall_communities_h21_indeces) %>% #i is the index
  select(i, Date_num_i, Date_num_j) #Date_num_j is the day+1 change
changepoint_overall_communities_h21_days

gam_fit_H21_overrall <- gamm(BC_dist ~ s(Date_num_i),
                             correlation = corAR1(form = ~ Date_num_i),
                             data = first_diff_df_H21)
plot(gam_fit_H21_overrall$gam, shade = TRUE)

##ORDINATION####
###NAIVE SYSTEM######
#Object with bray curtis distances
phyloseq.bacteria.samples_H21.ordered_RA.bray
#Ordination
phyloseq.bacteria.samples_H21.ordered_RA.bray.ord <- metaMDS(phyloseq.bacteria.samples_H21.ordered_RA.bray, 
                                                 k=2, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
#Metadata 
phyloseq.bacteria.samples_H21_metadata <- data.frame(sample_data(phyloseq.bacteria.samples_H21.ordered_RA))

###Centroids by Phase#####
#Simple ordination plot
phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot <- ordiplot(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs <- scores(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs),
                                                Copper_level_mg_L = phyloseq.bacteria.samples_H21_metadata$Copper_level_mg_L, 
                                                SampleID = phyloseq.bacteria.samples_H21_metadata$SampleID, 
                                                Date_num = phyloseq.bacteria.samples_H21_metadata$Date_num,
                                                Collection_Date = phyloseq.bacteria.samples_H21_metadata$Collection_Date, 
                                                Ammonia_level = phyloseq.bacteria.samples_H21_metadata$Ammonia_mg_L, 
                                                Date_num_phase = phyloseq.bacteria.samples_H21_metadata$Date_num_phase, 
                                                Date_num_phase_abbrv = phyloseq.bacteria.samples_H21_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.bacteria.samples.H21.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Date_num_phase_abbrv, 
                                                     data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs, 
                                                     FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_segs <- merge(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs, 
                             setNames(phyloseq.bacteria.samples.H21.ra.dates.bray.cent.phase, 
                                      c("Date_num_phase_abbrv","cMDS1","cMDS2")),
                             by = 'Date_num_phase_abbrv', 
                             sort = F)

####PERMANOVA - Just phase####
set.seed(98)
phase_h21_BC_adonis  <- adonis2(phyloseq.bacteria.samples_H21.ordered_RA.bray ~ Date_num_phase_abbrv, 
                                phyloseq.bacteria.samples_H21_metadata, 
                                #strata = phyloseq.bacteria.samples_H21_metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_h21_BC_adonis #32.6% of the variation is due to Enclosure, p = 1e-04

##PERMDISPS
# Run the betadisper function, average distance to centroid
bray.h21.phase.disp <- betadisper(phyloseq.bacteria.samples_H21.ordered_RA.bray, 
                                  phyloseq.bacteria.samples_H21_metadata$Date_num_phase_abbrv)
bray.h21.phase.disp 
##Then test by permuting
set.seed(98)
bray.h21.phase.permdisp <- permutest(bray.h21.phase.disp, permutations = 9999)
bray.h21.phase.permdisp ##S, p = 1e-04


# Extract R2 and p-values
R2_adonis_h21_phase <- phase_h21_BC_adonis$R2[1] 
pvalue_adonis_h21_phase <- phase_h21_BC_adonis$`Pr(>F)`[1]
# R2_mirkat_h21_phase <- phase_h21_BC_mirkat$R2
# pvalue_mirkat_h21_phase<-  phase_h21_BC_mirkat$p_values


#PLOT
h21_phase_BC_beta_div <- ggplot(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_segs) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Naive System Phases", fill = "Naive System Phases") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Date_num_phase_abbrv, 
                                     colour = Date_num_phase_abbrv), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10, show.legend = F) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = naive_phase_palette)+
  scale_fill_manual(values = naive_phase_palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  annotate("text", x = -1, y = -0.9,
           label = "PERMANOVA\nPhase",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1, y = -0.85,
           label = paste("R² = ", round(R2_adonis_h21_phase* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_h21_phase, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
h21_phase_BC_beta_div

###ESTABLISHED SYSTEM######
#Object with bray curtis distances
phyloseq.bacteria.samples_P1.ordered_RA.bray
#Ordination
phyloseq.bacteria.samples_P1.ordered_RA.bray.ord <- metaMDS(phyloseq.bacteria.samples_P1.ordered_RA.bray, 
                                                             k=2, try = 50, 
                                                             trymax = 1000,
                                                             autotransform = F)
#Metadata 
phyloseq.bacteria.samples_P1_metadata <- data.frame(sample_data(phyloseq.bacteria.samples_P1.ordered_RA))

###Centroids by Phase#####
#Simple ordination plot
phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot <- ordiplot(phyloseq.bacteria.samples_P1.ordered_RA.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_scrs <- scores(phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_scrs),
                                                                     Copper_level_mg_L = phyloseq.bacteria.samples_P1_metadata$Copper_level_mg_L, 
                                                                     SampleID = phyloseq.bacteria.samples_P1_metadata$SampleID, 
                                                                     Date_num = phyloseq.bacteria.samples_P1_metadata$Date_num,
                                                                     Collection_Date = phyloseq.bacteria.samples_P1_metadata$Collection_Date, 
                                                                     Ammonia_level = phyloseq.bacteria.samples_P1_metadata$Ammonia_mg_L, 
                                                                     Date_num_phase = phyloseq.bacteria.samples_P1_metadata$Date_num_phase, 
                                                                     Date_num_phase_abbrv = phyloseq.bacteria.samples_P1_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.bacteria.samples.P1.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Date_num_phase_abbrv, 
                                                                    data = phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_scrs, 
                                                                    FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_segs <- merge(phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_scrs, 
                                                                     setNames(phyloseq.bacteria.samples.P1.ra.dates.bray.cent.phase, 
                                                                              c("Date_num_phase_abbrv","cMDS1","cMDS2")),
                                                                     by = 'Date_num_phase_abbrv', 
                                                                     sort = F)

####PERMANOVA - Just phase####
set.seed(98)
phase_P1_BC_adonis  <- adonis2(phyloseq.bacteria.samples_P1.ordered_RA.bray ~ Date_num_phase_abbrv, 
                                phyloseq.bacteria.samples_P1_metadata, 
                                #strata = phyloseq.bacteria.samples_P1_metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_P1_BC_adonis #15.2% of the variation is due to Phase, p = 1e-04

##PERMDISPS
# Run the betadisper function, average distance to centroid
bray.P1.phase.disp <- betadisper(phyloseq.bacteria.samples_P1.ordered_RA.bray, 
                                  phyloseq.bacteria.samples_P1_metadata$Date_num_phase_abbrv)
bray.P1.phase.disp 
##Then test by permuting
set.seed(98)
bray.P1.phase.permdisp <- permutest(bray.P1.phase.disp, permutations = 9999)
bray.P1.phase.permdisp ##S, p = 0.001


# Extract R2 and p-values
R2_adonis_P1_phase <- phase_P1_BC_adonis$R2[1] 
pvalue_adonis_P1_phase <- phase_P1_BC_adonis$`Pr(>F)`[1]
# R2_mirkat_P1_phase <- phase_P1_BC_mirkat$R2
# pvalue_mirkat_P1_phase<-  phase_P1_BC_mirkat$p_values


#PLOT
P1_phase_BC_beta_div <- ggplot(phyloseq.bacteria.samples_P1.ordered_RA.bray.ord_plot_segs) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Established System Phases", fill = "Established System Phases") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Date_num_phase_abbrv, 
                                     colour = Date_num_phase_abbrv), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = naive_phase_palette)+
  scale_fill_manual(values = naive_phase_palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  annotate("text", x = -1, y = -0.9,
           label = "PERMANOVA\nPhase",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1, y = -0.85,
           label = paste("R² = ", round(R2_adonis_P1_phase* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_P1_phase, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
P1_phase_BC_beta_div




#Naive system
#Days 0 -10 (2023-10-07 through 2023-10-17): "First" copper exposure, very low levels of copper (phase 1)
#Days 11 - 28 (2023-10-18 through 2023-11-04): Down time between "first" copper exposure and start of actual first exposure to high copper levels
#Day 29 - 56 (2023-11-05 through 2023-12-02): Actual first copper exposure (phase 2)
#Day 57-86 (2023-12-03 through 2024-01-01): Downtime between first (phase 2) and second (phase 3) copper exposure
#Day 87-131 (2024-01-02 through 2024-02-15): Second copper exposure (phase 3)
#Day 132 - 147 (2024-02-16 through 2024-03-02): Downtime after phase 3 

#Established system
#Days 0 - 50 (2023-11-14 through 2024-01-03): Phase 1, conservative copper dosing (phase 1)
#Days 51 - 65 (2024-01-04 through 2024-01-18): Phase 2, increasing copper dosing (phase 2)
#Days 66 - 98 (2024-01-19 through 2024-02-20): Phase 3, back to copper dosage as needed
#Day 99 - 168 (2024-02-21 through 2024-04-30) Final, after copper treatment finished

#Days since start variable
phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs_3 <- 
  phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs %>%
  mutate(Date_num = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07")),
    Enclosure == "P1"  ~ as.numeric(Collection_Date - as.Date("2023-11-14"))
  )) %>%
  #Separate if wanting to do separate day scales
  mutate(Date_num_naive = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07"))),
    Date_num_established = case_when(
      Enclosure == "P1" ~ as.numeric(Collection_Date - as.Date("2023-11-14"))))%>%
  mutate(Date_num_phase_naive = case_when(
    Enclosure == "H21" & Date_num >= 0   & Date_num <= 10  ~ "Phase 1 (Day 0 - 10)",
    Enclosure == "H21" & Date_num >= 11  & Date_num <= 28  ~ "Downtime Phase 1 - 2 (Day 11 - 28)",
    Enclosure == "H21" & Date_num >= 29  & Date_num <= 56  ~ "Phase 2 (Day 29 - 56)",
    Enclosure == "H21" & Date_num >= 57  & Date_num <= 86  ~ "Downtime Phase 2 - 3 (Day 57 - 86)",
    Enclosure == "H21" & Date_num >= 87  & Date_num <= 131 ~ "Phase 3 (Day 87 - 131)",
    Enclosure == "H21" & Date_num >= 132 & Date_num <= 147 ~ "Post-Treatment completion (Day 132 - 147)"
  )) %>%
  mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
                                       levels = c("Phase 1 (Day 0 - 10)", 
                                                  "Downtime Phase 1 - 2 (Day 11 - 28)", 
                                                  "Phase 2 (Day 29 - 56)", 
                                                  "Downtime Phase 2 - 3 (Day 57 - 86)", 
                                                  "Phase 3 (Day 87 - 131)", 
                                                  "Post-Treatment completion (Day 132 - 147)")))%>%
  mutate(Date_num_phase_established = case_when(
    Enclosure == "P1" & Date_num >= 0   & Date_num <= 50  ~ "Phase 1 (Day 0 - 50)",
    Enclosure == "P1" & Date_num >= 51  & Date_num <= 65  ~ "Phase 2 (Day 51 - 65)",
    Enclosure == "P1" & Date_num >= 66  & Date_num <= 98  ~ "Phase 3 (Day 66 - 98)",
    Enclosure == "P1" & Date_num >= 99  & Date_num <= 168 ~ "Post-Treatment completion (Day 99 - 168)"
  ))%>%
  mutate(Date_num_phase_established = factor(Date_num_phase_established, 
                                       levels = c("Phase 1 (Day 0 - 50)",
                                                  "Phase 2 (Day 51 - 65)",
                                                  "Phase 3 (Day 66 - 98)",
                                                  "Post-Treatment completion (Day 99 - 168)")))
#Shape by phase 
enclosure_BC_beta_div_3 <- ggplot(phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs_3) + 
  theme_bw() +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2,
  #                  color = Enclosure), show.legend = F)+
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
                                     colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure, shape = Date_num_phase_naive), 
             size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(15, 3, 16, 4, 17, 18), name = "Naive system\nPhase") +
  guides(
    shape = guide_legend(override.aes = list(size = 7), ncol = 1)) +
  new_scale("shape") +
  geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure, shape = Date_num_phase_established), 
             size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Established system\nPhase") +
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure), 
            colour= "white", size = 6, fontface = "bold") +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "System", 
       fill = "System") +
  scale_color_manual(values = enclosure.palette,
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  scale_fill_manual(values = enclosure.palette, 
                    labels = c("H21" = "Naive", 
                               "P1" = "Established"))+
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7), ncol = 1)) +
  annotate("text", x = -1, y = -0.9,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1, y = -0.85,
           label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_beta_div_3

##Colored by days since start
enclosure_BC_beta_div_4 <- ggplot(phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs_3) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       shape = "System", 
       color = "Day") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, group = Enclosure), 
               #color = "black",
               alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, 
                 colour = Date_num,
                 shape = Enclosure), 
             size = 5, 
             alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure), size = 10) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure),
            colour= "white", size = 5, fontface = "bold") +
  scale_shape_manual(values = c(15,16),
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  scale_color_viridis_c(option = "viridis")+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 20),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  annotate("text", x = -0.9, y = -0.9,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -0.9, y = -0.8,
           label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_beta_div_4
ggsave("enclosure_BC_beta_div_4.png", 
       enclosure_BC_beta_div_4, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 15)


#Colored by days since start, separate scales for each system
enclosure_BC_beta_div_5 <- ggplot(phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs_3) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", 
       title= "MICROBIOME", 
       shape = "System") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
  stat_ellipse(
    geom = "polygon",
    aes(x = MDS1, y = MDS2, fill = Enclosure, color = Enclosure, group = Enclosure),
    alpha = 0.2,
    lty = 2,
    linewidth = 1,
    level = 0.95,
    show.legend = FALSE
  )+
  #Problem with layering. If adding centroids here, dots go on top. 
  geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure,
                 color = Enclosure), size = 10,
             show.legend = FALSE) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure),
            colour= "white", size = 5, fontface = "bold") +
  scale_fill_manual(values = c(
    "#fc8d62",
    "#8da0cb"
  ))+
  scale_color_manual(values = c(
    "#fc8d62",
    "#8da0cb"
  ))+
  new_scale_color()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  #Naive system individual dots
  geom_point(data = phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs_3 %>%
               filter(Enclosure == "H21"), aes(x=MDS1, y=MDS2, 
                 colour = Date_num_naive,
                 shape = Enclosure), 
             size = 5) + # individuals
  scale_color_gradient2(low = "#fad2c2", 
                        mid = "#fc8d62", 
                        high = "#ad431a",
                        midpoint = 74,
                        limits = c(0,150),
                        breaks = c(0, 50, 100, 150))+
  # scale_color_viridis_c(option = "viridis",
  #                       limits = c(0,150),
  #                       breaks = c(0, 50, 100, 150))+
  labs(color = "Days (Naive system)")+
  #Established system individual dots. 
  new_scale_color()+
  geom_point(data = phyloseq.bacteria.samples.ra.dates.bray.enclosure.segs_3 %>%
               filter(Enclosure == "P1"), aes(x=MDS1, y=MDS2, 
                 colour = Date_num_established,
                 shape = Enclosure), 
             size = 5) + # individuals
  scale_color_gradient2(low = "#c6cfe6",
                        mid = "#8da0cb", 
                        high = "#012983", 
                        midpoint = 84,
                        limits = c(0,170),
                        breaks = c(0, 50, 100, 170))+
  labs(color = "Days (Established system)")+
  # scale_color_viridis_c(option = "rocket", 
  #                       limits = c(0,170),
  #                       breaks = c(0, 50, 100, 170))+
  scale_shape_manual(values = c(15,16),
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 18),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    color = guide_colorbar(order = 3),
    #color = guide_legend(order = 3),
    shape = guide_legend(override.aes = list(size = 7),
                         order = 1)) +
  # new_scale_color()+
  # geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure, 
  #                color = Enclosure), size = 10,
  #            show.legend = FALSE) + # centroids +
  # geom_text(aes (x= cMDS1, y = cMDS2,
  #                label= Enclosure),
  #           colour= "white", size = 5, fontface = "bold") +
  # scale_color_manual(values = c(
  #   "#fc8d62",
  #   "#8da0cb"
  # ))+
  annotate("text", x = -0.9, y = -0.9,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -0.9, y = -0.8,
           label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_beta_div_5
ggsave("enclosure_BC_beta_div_5.png", 
       enclosure_BC_beta_div_5, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 15)

#POSTER BIMS 2026
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Posters/enclosure_BC_beta_div.png", 
       enclosure_BC_beta_div+
         theme(plot.title = element_text(size = 50, face = "bold")), 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 19)

# Run envfit on your NMDS object 
# set.seed(98)
# fit <- envfit(phyloseq.bacteria.samples.ra.dates.bray.ord, 
#               phyloseq.bacteria.samples.df$Copper_level_mg_L, 
#               permutations = 999)
# 
# # Extract coordinates for the copper vector
# vec <- as.data.frame(fit$vectors$arrows * fit$vectors$r)  # scale by r
# colnames(vec) <- c("xend", "yend")
# vec$label <- rownames(vec)  # will be "Copper_level_mg_L"
# 
# enclosure_BC_beta_div_copper <- enclosure_BC_beta_div +
#   geom_segment(
#     data = vec,
#     aes(x = 0, y = 0, xend = xend, yend = yend),
#     arrow = arrow(length = unit(0.3, "cm")),  # adds arrowhead
#     colour = "red",
#     size = 1.2
#   )
# enclosure_BC_beta_div_copper




###BOTH NAIVE AND ESTABLISHED SYSTEM######
#All samples
phyloseq.bacteria.samples
#Relative abundances
phyloseq.bacteria.samples.ra

#Bray curtis distances
phyloseq.bacteria.samples_RA.bray <- vegdist(t(phyloseq.bacteria.samples.ra@otu_table), method = "bray")

#Ordination
phyloseq.bacteria.samples_RA.bray.ord <- metaMDS(phyloseq.bacteria.samples_RA.bray, 
                                                             k=2, try = 50, 
                                                             trymax = 1000,
                                                             autotransform = F)
#Metadata 
phyloseq.bacteria.samples_metadata <- data.frame(sample_data(phyloseq.bacteria.samples))

###Centroids by System and Phase#####
#Simple ordination plot
phyloseq.bacteria.samples_RA.bray.ord_plot <- ordiplot(phyloseq.bacteria.samples_RA.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples_RA.bray.ord_plot_scrs <- scores(phyloseq.bacteria.samples_RA.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.bacteria.samples_RA.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.bacteria.samples_RA.bray.ord_plot_scrs),
                                                         Enclosure = phyloseq.bacteria.samples_metadata$Enclosure,
                                                                     Copper_level_mg_L = phyloseq.bacteria.samples_metadata$Copper_level_mg_L, 
                                                                     SampleID = phyloseq.bacteria.samples_metadata$SampleID, 
                                                                     Date_num = phyloseq.bacteria.samples_metadata$Date_num,
                                                                     Collection_Date = phyloseq.bacteria.samples_metadata$Collection_Date, 
                                                                     Ammonia_level = phyloseq.bacteria.samples_metadata$Ammonia_mg_L, 
                                                                     Date_num_phase = phyloseq.bacteria.samples_metadata$Date_num_phase, 
                                                                     Date_num_phase_abbrv = phyloseq.bacteria.samples_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.bacteria.samples.H21.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Enclosure + Date_num_phase_abbrv, 
                                                                    data = phyloseq.bacteria.samples_RA.bray.ord_plot_scrs, 
                                                                    FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples_RA.bray.ord_plot_segs <- merge(phyloseq.bacteria.samples_RA.bray.ord_plot_scrs, 
                                                                     setNames(phyloseq.bacteria.samples.H21.ra.dates.bray.cent.phase, 
                                                                              c("Enclosure", "Date_num_phase_abbrv","cMDS1","cMDS2")),
                                                                     by = c('Enclosure', 'Date_num_phase_abbrv'), 
                                                                     sort = F)

####PERMANOVA - Enclosure + Phase####
set.seed(98)
phase_enclosure_BC_adonis  <- adonis2(phyloseq.bacteria.samples_RA.bray ~ Enclosure + Date_num_phase_abbrv, 
                                phyloseq.bacteria.samples_metadata, 
                                #strata = phyloseq.bacteria.samples_metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_enclosure_BC_adonis 

# ##PERMDISPS
# # Run the betadisper function, average distance to centroid
# bray.phase.enclosure.disp <- betadisper(phyloseq.bacteria.samples_RA.bray, 
#                                   phyloseq.bacteria.samples_metadata$Date_num_phase_abbrv)
# bray.phase.enclosure.disp 
# ##Then test by permuting
# set.seed(98)
# bray.h21.phase.permdisp <- permutest(bray.h21.phase.disp, permutations = 9999)
# bray.h21.phase.permdisp ##S, p = 1e-04


# # Extract R2 and p-values
# R2_adonis_phase_enclosure <- phase_BC_adonis$R2[1] 
# pvalue_adonis_phase_enclosure <- phase_BC_adonis$`Pr(>F)`[1]
# R2_mirkat_phase <- phase_BC_mirkat$R2
# pvalue_mirkat_phase<-  phase_BC_mirkat$p_values


#PLOT
phase_enclosure_BC_beta_div <- ggplot(phyloseq.bacteria.samples_RA.bray.ord_plot_segs) + 
  theme_bw() +
  facet_grid(~Enclosure)+
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Phase", fill = "Phase") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Date_num_phase_abbrv, 
                                     colour = Date_num_phase_abbrv), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  # scale_color_manual(values = naive_phase_palette)+
  # scale_fill_manual(values = naive_phase_palette)+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) 
  # annotate("text", x = -1, y = -0.9,
  #          label = "PERMANOVA\nPhase",
  #          hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  # annotate("text", x = -1, y = -0.85,
  #          label = paste("R² = ", round(R2_adonis_phase* 100, 1), "%",
  #                        "\np = ", round(pvalue_adonis_phase, 4)),
  #          hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
phase_enclosure_BC_beta_div


###Centroids by copper quartiles facetted by enclosure #####
## BC
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.ra.dates.bray.cent.copper_quart <- aggregate(
  cbind(MDS1, MDS2) ~ Enclosure + Copper_quartile,
  data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs,
  FUN = mean)

#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.ra.dates.bray.copper_quart.segs <- merge(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs, 
                                                          setNames(phyloseq.bacteria.samples.ra.dates.bray.cent.copper_quart, 
                                                                   c("Enclosure", "Copper_quartile",
                                                                     "cMDS1","cMDS2")),
                                                          by = c("Enclosure", "Copper_quartile"), 
                                                          sort = F)
##PERMANOVA - Just copper###
set.seed(98)
copper_quart_BC_adonis  <- adonis2(phyloseq.bacteria.samples.ra.dates.bray ~ Copper_quartile, 
                                phyloseq.bacteria.samples.df, 
                                strata = phyloseq.bacteria.samples.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
copper_quart_BC_adonis 


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_adonis  <- adonis2(phyloseq.bacteria.samples.ra.dates.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.df, 
#                                 strata = phyloseq.bacteria.samples.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_adonis  <- adonis2(phyloseq.bacteria.samples.ra.dates.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.df , 
#                                 strata = phyloseq.bacteria.samples.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS#
#COPPER quartile#
# Run the betadisper function, average distance to centroid
bray.copper.quart.disp <- betadisper(phyloseq.bacteria.samples.ra.dates.bray,
                                  phyloseq.bacteria.samples.df$Copper_quartile)
bray.copper.quart.disp
##Then test by permuting
set.seed(98)
bray.copper.quart.permdisp <- permutest(bray.copper.quart.disp, permutations = 9999)
bray.copper.quart.permdisp ##S, p = 1e-04

# Extract R2 and p-values
R2_adonis_copper_quart <- copper_quart_BC_adonis$R2[1] 
pvalue_adonis_copper_quart<-  copper_quart_BC_adonis$`Pr(>F)`[1]
# R2_adonis_copper <- enclosure_BC_adonis$R2[2] 
# pvalue_adonis_copper <- enclosure_BC_adonis$`Pr(>F)`[2]

#### PLOTS
## BC
copper_quart_BC_beta_div <- ggplot(phyloseq.bacteria.samples.ra.dates.bray.copper_quart.segs) + 
  theme_bw() +
  facet_grid(~Enclosure,
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Copper (mg/L) Quartile", 
       fill = "Copper (mg/L) Quartile") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, 
                                     fill= Copper_quartile.abvr, 
                                     colour = Copper_quartile.abvr), alpha = 0.2, 
               lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Copper_quartile.abvr), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Copper_quartile.abvr), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Copper_quartile.abvr), colour= "white", 
            size = 6, fontface = "bold") +
  scale_color_manual(values = plasma_quartiles)+
  scale_fill_manual(values = plasma_quartiles)+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold"))
  # annotate("text", x = 0.4, y = 0.9,
  #          label = "Enclosure",
  #          hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  # annotate("text", x =0.4, y = 0.9,
  #          label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
  #                        "\np = ", round(pvalue_adonis_enclosure, 4)),
  #          hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
copper_quart_BC_beta_div
ggsave("copper_quart_BC_beta_div.png", 
       copper_quart_BC_beta_div, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 14)


#make DF from metadata
phyloseq.bacteria.samples.df <- as(phyloseq.bacteria.samples.ra@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date),
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.df 


###Continuous copper facetted by enclosure #####
#No samples without copper levels:
phyloseq.bacteria.samples.ra.copper <- subset_samples(phyloseq.bacteria.samples.ra.dates, 
                                                                !is.na(Copper_level_mg_L))
phyloseq.bacteria.samples.ra.copper <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.ra.copper)> 0, 
                                                            phyloseq.bacteria.samples.ra.copper)
phyloseq.bacteria.samples.ra.copper #33481 taxa and 204 samples

####H21 #######
phyloseq.bacteria.samples.ra.copper.H21 <- subset_samples(phyloseq.bacteria.samples.ra.copper, 
                                                                    Enclosure == "H21")
phyloseq.bacteria.samples.ra.copper.H21 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.ra.copper.H21)> 0, 
                                                                phyloseq.bacteria.samples.ra.copper.H21)
phyloseq.bacteria.samples.ra.copper.H21 #33418 taxa and 82 samples 

##BC 
phyloseq.bacteria.samples.ra.copper.H21.bray <- vegdist(t(phyloseq.bacteria.samples.ra.copper.H21@otu_table), 
                                                                  method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.copper.H21.df <- as(phyloseq.bacteria.samples.ra.copper.H21@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date),
    #Date number 
    Date_num = as.numeric(Collection_Date - min(Collection_Date)),
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.copper.H21.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_H21 <- D2K(
  as.matrix(phyloseq.bacteria.samples.ra.copper.H21.bray))

# Covariates (date_num)
datenum_covariate_H21 <- phyloseq.bacteria.samples.copper.H21.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_H21 

# # Covariates (month)
# month_covariate_H21 <- model.matrix(~ Collection_Month, 
#                                         data = phyloseq.bacteria.samples.copper.H21.df)[, -1]
# 

#Run model
mirkat_H21_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.copper.H21.df$Copper_level_mg_L,
  X = datenum_covariate_H21,
  Ks = bray_kernel_H21,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_H21_copper_cont #R2 6.9% , p_values = 0.0001017485

#####PERMANOVA #######
set.seed(98)
adonis_H21_copper_cont  <- adonis2(phyloseq.bacteria.samples.ra.copper.H21.bray ~ Copper_level_mg_L, 
                                       phyloseq.bacteria.samples.copper.H21.df, 
                                       strata = phyloseq.bacteria.samples.copper.H21.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_H21_copper_cont # R2 13.5%, p value 2e-04 *

##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.copper.cont.disp.H21 <- betadisper(phyloseq.bacteria.samples.ra.copper.H21.bray,
                                            phyloseq.bacteria.samples.copper.H21.df$Copper_quartile)
bray.copper.cont.disp.H21
##Then test by permuting
set.seed(98)
bray.copper.cont.permdisp.H21 <- permutest(bray.copper.cont.disp.H21, permutations = 9999)
bray.copper.cont.permdisp.H21 ##S, p =0.0016

####P1 #######
phyloseq.bacteria.samples.ra.copper.P1 <- subset_samples(phyloseq.bacteria.samples.ra.copper, 
                                                                   Enclosure == "P1")
phyloseq.bacteria.samples.ra.copper.P1 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.ra.copper.P1)> 0, 
                                                               phyloseq.bacteria.samples.ra.copper.P1)
phyloseq.bacteria.samples.ra.copper.P1 #28427 taxa and 122 samples

##BC 
phyloseq.bacteria.samples.ra.copper.P1.bray <- vegdist(t(phyloseq.bacteria.samples.ra.copper.P1@otu_table), 
                                                                 method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.copper.P1.df <- as(phyloseq.bacteria.samples.ra.copper.P1@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date),
    #Date number 
    Date_num = as.numeric(Collection_Date - min(Collection_Date)),
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.copper.P1.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_P1 <- D2K(
  as.matrix(phyloseq.bacteria.samples.ra.copper.P1.bray))

# Covariates (date_num)
datenum_covariate_P1 <- phyloseq.bacteria.samples.copper.P1.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_P1 

# # Covariates (month)
# month_covariate_P1 <- model.matrix(~ Collection_Month, 
#                                        data = phyloseq.bacteria.samples.copper.P1.df)[, -1]
# 
# month_covariate_P1

#Run model
mirkat_P1_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.copper.P1.df$Copper_level_mg_L,
  X = datenum_covariate_P1,
  Ks = bray_kernel_P1,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_P1_copper_cont #p val 1.320553e-08, R2 2.5%

#####PERMANOVA #####
set.seed(98)
adonis_P1_copper_cont  <- adonis2(phyloseq.bacteria.samples.ra.copper.P1.bray ~ Copper_level_mg_L, 
                                      phyloseq.bacteria.samples.copper.P1.df, 
                                      strata = phyloseq.bacteria.samples.copper.P1.df$Collection_Month,
                                      by = "margin",
                                      permutations = 9999)
adonis_P1_copper_cont # R2 10.4%, p value 0.021 *


####MIRKAT GLMM FOR COPPER #####
#make DF from metadata
phyloseq.bacteria.samples.copper.all.df <- as(phyloseq.bacteria.samples.ra.copper@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date))%>%
  group_by(Enclosure)%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  mutate(
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.copper.all.df

#Bray-Curtis distances
phyloseq.bacteria.samples.ra.copper.bray <- vegdist(t(phyloseq.bacteria.samples.ra.copper@otu_table), 
                                                              method = "bray")
phyloseq.bacteria.samples.ra.copper.bray

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.ra.copper.bray))

#Jaccard distances
phyloseq.bacteria.samples.ra.copper.jac <- vegdist(t(phyloseq.bacteria.samples.ra.copper@otu_table), 
                                                             method = "jaccard")
phyloseq.bacteria.samples.ra.copper.jac

jac_kernel_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.ra.copper.jac))

# Covariates (date_num)
datenum_covariate_all <- phyloseq.bacteria.samples.copper.all.df%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_all

#Id for repeated measures 
Enclosure_d_mirkat_all <- phyloseq.bacteria.samples.copper.all.df %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_mirkat_all

#Run model - dont think it works. just two enclosures and enclosure id is the same as the outcome 
mirkat_enclosure_glm <- GLMMMiRKAT(
  y = phyloseq.bacteria.samples.copper.all.df$Copper_level_mg_L,
  X = datenum_covariate_all,
  Ks = list(bray_kernel_all, jac_kernel_all),
  id = Enclosure_d_mirkat_all,
  model  = "gaussian",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE,                    # random intercept only
  nperm  = 5000)                   # (not n.perm)
mirkat_enclosure_glm #0.00019996 omnibus_p 

###PLOT FOR BOTH ENCLOSURES#######
#### Extract R2 and p-values
#####H21
R2_mirkat_copper_cont_H21 <- mirkat_H21_copper_cont$R2
pvalue_mirkat_copper_cont_H21<- mirkat_H21_copper_cont$p_values

##PERMANOVA
R2_adonis_copper_cont_H21 <- adonis_H21_copper_cont$R2[1]
pvalue_adonis_copper_cont_H21<- adonis_H21_copper_cont$`Pr(>F)`[1]

#####P1
##Mirkat
R2_mirkat_copper_cont_P1 <- mirkat_P1_copper_cont$R2
pvalue_mirkat_copper_cont_P1<- mirkat_P1_copper_cont$p_values

##PERMANOVA
R2_adonis_copper_cont_P1 <- adonis_P1_copper_cont$R2[1]
pvalue_adonis_copper_cont_P1<- adonis_P1_copper_cont$`Pr(>F)`[1]
#### PLOTS
## BC
copper_cont_BC_beta_div <- ggplot(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs) + 
  theme_bw() +
  facet_wrap(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, 
  #                                    fill= Copper_quartile.abvr, 
  #                                    colour = Copper_quartile.abvr), alpha = 0.2, 
  #              lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, color = Copper_level_mg_L), size = 5, alpha = 0.8) + # individuals
  scale_color_viridis_c(option = "plasma")+
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Copper levels (mg/L)") +
  theme(legend.position = "bottom",
        legend.direction = "vertical", 
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold")) +
  #H21
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.5, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.7, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "P1"),
            aes(x = 0.2, y = -1.1, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "P1"),
            aes(x = 0.2, y = -1.4, 
                label = paste("R² = ", round(R2_adonis_copper_cont_P1* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_P1, 6))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")
copper_cont_BC_beta_div
ggsave("copper_cont_BC_beta_div.png", 
       copper_cont_BC_beta_div, 
       device = "png", 
       dpi = 600, 
       height = 9, 
       width = 16)


##Add collection month as shape
copper_cont_BC_beta_div_2 <- ggplot(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs) + 
  theme_bw() +
  facet_wrap(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, 
  #                                    fill= Copper_quartile.abvr, 
  #                                    colour = Copper_quartile.abvr), alpha = 0.2, 
  #              lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, color = Copper_level_mg_L,
                 shape = Collection_Month), size = 5, alpha = 0.8) + # individuals
  scale_color_viridis_c(option = "plasma")+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8))+
  # scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8), 
  #                    labels = c("2023-10" = "Oct 2023",
  #                               "2023-11" = "Nov 2023",
  #                               "2023-12" = "Dec 2023",
  #                               "2024-01" = "Jan 2024",
  #                               "2024-02" = "Feb 2024",
  #                               "2024-03" = "Mar 2024",
  #                               "2024-04" = "Apr 2024",
  #                               "2024-05" = "May 2024"))+
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Copper levels (mg/L)",
       shape = "Collection Month") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold")) +
  #H21
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.5, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.7, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.8, 
                label = paste("R² = ", round(R2_adonis_copper_cont_P1* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_P1, 6))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  guides(color = guide_colorbar(direction = "vertical"))
copper_cont_BC_beta_div_2
ggsave("copper_cont_BC_beta_div_2.png", 
       copper_cont_BC_beta_div_2, 
       device = "png", 
       dpi = 600, 
       height = 9, 
       width = 16)

#POSTER BIMS 2026
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Posters/copper_cont_BC_beta_div.png", 
       copper_cont_BC_beta_div + 
         theme(legend.position = "none",
               plot.title = element_text(colour = "black", size = 50, face = "bold"),
               strip.text = element_text(colour = "white", size = 45, face = "bold")), 
       device = "png", 
       dpi = 600, 
       height = 9.5, 
       width = 15.5)



#Adding date_num and phases
phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs_2 <- phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs %>%
  mutate(Date_num = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07")),
    Enclosure == "P1"  ~ as.numeric(Collection_Date - as.Date("2023-11-14"))
  )) %>%
  #Separate if wanting to do separate day scales
  mutate(Date_num_naive = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07"))),
    Date_num_established = case_when(
      Enclosure == "P1" ~ as.numeric(Collection_Date - as.Date("2023-11-14"))))%>%
  mutate(Date_num_phase_naive = case_when(
    Enclosure == "H21" & Date_num >= 0   & Date_num <= 10  ~ "Phase 1 (Day 0 - 10)",
    Enclosure == "H21" & Date_num >= 11  & Date_num <= 28  ~ "Downtime Phase 1 - 2 (Day 11 - 28)",
    Enclosure == "H21" & Date_num >= 29  & Date_num <= 56  ~ "Phase 2 (Day 29 - 56)",
    Enclosure == "H21" & Date_num >= 57  & Date_num <= 86  ~ "Downtime Phase 2 - 3 (Day 57 - 86)",
    Enclosure == "H21" & Date_num >= 87  & Date_num <= 131 ~ "Phase 3 (Day 87 - 131)",
    Enclosure == "H21" & Date_num >= 132 & Date_num <= 147 ~ "Post-treatment (Day 132 - 147)"
  )) %>%
  mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
                                       levels = c("Phase 1 (Day 0 - 10)", 
                                                  "Downtime Phase 1 - 2 (Day 11 - 28)", 
                                                  "Phase 2 (Day 29 - 56)", 
                                                  "Downtime Phase 2 - 3 (Day 57 - 86)", 
                                                  "Phase 3 (Day 87 - 131)", 
                                                  "Post-treatment (Day 132 - 147)")))%>%
  mutate(Date_num_phase_established = case_when(
    Enclosure == "P1" & Date_num >= 0   & Date_num <= 50  ~ "Phase 1 (Day 0 - 50)",
    Enclosure == "P1" & Date_num >= 51  & Date_num <= 65  ~ "Phase 2 (Day 51 - 65)",
    Enclosure == "P1" & Date_num >= 66  & Date_num <= 98  ~ "Phase 3 (Day 66 - 98)",
    Enclosure == "P1" & Date_num >= 99  & Date_num <= 168 ~ "Post-treatment (Day 99 - 168)"
  ))%>%
  mutate(Date_num_phase_established = factor(Date_num_phase_established, 
                                             levels = c("Phase 1 (Day 0 - 50)",
                                                        "Phase 2 (Day 51 - 65)",
                                                        "Phase 3 (Day 66 - 98)",
                                                        "Post-treatment (Day 99 - 168)")))


#Adding phases
copper_cont_BC_beta_div_3 <- ggplot(phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs_2) + 
  theme_bw() +
  facet_wrap(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
       color = "Copper levels (mg/L)") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, colour = Copper_level_mg_L, shape = Date_num_phase_naive), 
             size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(15, 3, 16, 4, 17, 18), name = "Naive system") +
  guides(
    shape = guide_legend(override.aes = list(size = 7), ncol = 1)) +
  new_scale("shape") +
  geom_point(aes(x=MDS1, y=MDS2, colour = Copper_level_mg_L, 
                 shape = Date_num_phase_established), 
             size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Established system") +
  scale_color_viridis_c(option = "plasma")+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 20),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  #H21
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.8, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples_H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.83, 
                label = paste("R² = ", round(R2_adonis_copper_cont_P1* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_P1, 6))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  guides(color = guide_colorbar(direction = "vertical",
                                order = 1))
copper_cont_BC_beta_div_3 
ggsave("copper_cont_BC_beta_div_3.png", 
       copper_cont_BC_beta_div_3, 
       device = "png", 
       dpi = 600, 
       height = 9, 
       width = 16)



##Nitirfiers within overall community#######
##Nitrifying taxa - filtered by dates too####
phyloseq.bacteria.samples.nitifiers.ra <- subset_taxa(phyloseq.bacteria.samples.ra.dates, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                            Family == "Chromatiaceae" | # no lineages
                            Family == "Nitrosopumilaceae" | # AOA; some!
                            Family == "Nitrososphaeraceae" | # no lineages
                            Order == "Candidatus Nitrosomirales" | # no lineages
                            Order == "Candidatus Nitrosocaldales" | # no lineages
                            Family == "Nitrospiraceae" | # NOB/Commamox; some!
                            Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                            Family == "Nitrobacteraceae" | # none
                            Family == "Gallionellaceae" | # none
                            Family == "Nitrospinaceae") # NOB; some, plus a new one!
phyloseq.bacteria.samples.nitifiers.ra #827 taxa and 218 samples
phyloseq.bacteria.samples.nitifiers.ra <- subset_samples(phyloseq.bacteria.samples.nitifiers.ra, sample_sums(phyloseq.bacteria.samples.nitifiers.ra) > 0)
phyloseq.bacteria.samples.nitifiers.ra #827 taxa and 218 samples 

##BC 
phyloseq.bacteria.samples.nitifiers.ra.bray <- vegdist(t(phyloseq.bacteria.samples.nitifiers.ra@otu_table), method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.nitifiers.df <- as(phyloseq.bacteria.samples.nitifiers.ra@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date),
    #Date number 
    Date_num = as.numeric(Collection_Date - min(Collection_Date)),
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.nitifiers.df 

###ORDINATION####
set.seed(98)
phyloseq.bacteria.samples.nitifiers.ra.bray.ord <- metaMDS(phyloseq.bacteria.samples.nitifiers.ra.bray, 
                                                 k=3, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
###Centroids by enclosure #####
## BC
#Simple ordination plot
phyloseq.bacteria.samples.nitifiers.ra.bray.plot <- ordiplot(phyloseq.bacteria.samples.nitifiers.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples.nitifiers.ra.bray.scrs <- scores(phyloseq.bacteria.samples.nitifiers.ra.bray.plot, display = "sites")
#Add metadata to coordinates
phyloseq.bacteria.samples.nitifiers.ra.bray.scrs <- cbind(as.data.frame(phyloseq.bacteria.samples.nitifiers.ra.bray.scrs),
                                                Copper_level_mg_L = phyloseq.bacteria.samples.nitifiers.df$Copper_level_mg_L, 
                                                SampleID = phyloseq.bacteria.samples.nitifiers.df$SampleID, 
                                                Collection_Date = phyloseq.bacteria.samples.nitifiers.df$Collection_Date, 
                                                Enclosure = phyloseq.bacteria.samples.nitifiers.df$Enclosure, 
                                                Attempt = phyloseq.bacteria.samples.nitifiers.df$Attempt, 
                                                Copper_quartile = phyloseq.bacteria.samples.nitifiers.df$Copper_quartile,
                                                Copper_quartile.abvr = phyloseq.bacteria.samples.nitifiers.df$Copper_quartile.abvr)%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.nitifiers.ra.bray.cent.enclosure <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
                                                              data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs, 
                                                              FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs <- merge(phyloseq.bacteria.samples.nitifiers.ra.bray.scrs, 
                                                          setNames(phyloseq.bacteria.samples.nitifiers.ra.bray.cent.enclosure, 
                                                                   c("Enclosure","cMDS1","cMDS2")),
                                                          by = 'Enclosure', 
                                                          sort = F)

####MIRKAT- Just enclosure#####
#Dataframe with metadata
df_metadata_nit_mirkat <- data.frame(sample_data(phyloseq.bacteria.samples.nitifiers.ra))
#BC distance matrix row order
labs_nit <- rownames(as.matrix(phyloseq.bacteria.samples.nitifiers.ra.bray))  
#Align metadata df to distance order
df_nit_mirkat <- df_metadata_nit_mirkat[labs_nit, , drop = FALSE] 

#Going to get kernels out of different distance matrices
#Bray curtis 
phyloseq.bacteria.samples.ra.dates.bray

#Jaccard distances
phyloseq.bacteria.samples.ra.dates.jac <- vegdist(t(phyloseq.bacteria.samples.ra.dates@otu_table), method = "jaccard")
phyloseq.bacteria.samples.ra.dates.jac

#Aitchison 
##Going to take out samples from P1 from april and may 2023 and september from H21 (raw counts object)
phyloseq.bacteria.samples.dates <- subset_samples(phyloseq.bacteria.samples, 
                                                     Collection_Date > "2023-10-01")
phyloseq.bacteria.samples.dates #33490 taxa and 218 samples
#Take out taxa sums = 0
phyloseq.bacteria.samples.dates  <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates) > 0, 
                                               phyloseq.bacteria.samples.dates)
phyloseq.bacteria.samples.dates  #33490 taxa and 218 samples
#Calculate distance
phyloseq.bacteria.samples.dates_clr <- microbiome::transform(phyloseq.bacteria.samples.dates, "clr") #convert raw counts to clr
phyloseq.bacteria.samples.dates_clr_dist_aitch <- dist(t(otu_table(phyloseq.bacteria.samples.dates_clr)), method = "euclidean") #calculate euclidean distances

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
#Bray curtis
bray_kernel_nit <- D2K(
  as.matrix(phyloseq.bacteria.samples.nitifiers.ra.bray))
#Jaccards
jac_kernel_nit <- D2K(
  as.matrix(phyloseq.bacteria.samples.ra.dates.jac))
# Aitchison
aitch_kernel_nit <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates_clr_dist_aitch))

#Modelling enclosure, have to make it dichotomus (used to be phyloseq.bacteria.samples.nitifiers.df)
Enclosure_d_nit_mirkat <- df_nit_mirkat %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_nit_mirkat

# # Covariates (month)
# month_covariate_nit <- model.matrix(~ Collection_Month, 
#                                 data = phyloseq.bacteria.samples.nitifiers.df)[, -1]
# 
# month_covariate_nit

# Covariates (date_num)
datenum_covariate_nit <- df_nit_mirkat%>%
  group_by(Enclosure)%>%
  #Date number 
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_nit 

#Run model (Are microbial communities associated with enclosure, after controlling for time?)
set.seed(98)
mirkat_enclosure_nit <- MiRKAT::MiRKAT(
  y = Enclosure_d_nit_mirkat,
  X = datenum_covariate_nit,
  Ks = list(bray_kernel_nit, jac_kernel_nit, aitch_kernel_nit),
  out_type = "D", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_enclosure_nit #R2 13.8% , p_values = 0



####PERMANOVA - Just enclosure#####
set.seed(98)
enclosure_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.bray ~ Enclosure, 
                                phyloseq.bacteria.samples.nitifiers.df, 
                                strata = phyloseq.bacteria.samples.nitifiers.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
enclosure_BC_nit_adonis #21% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.nitifiers.df, 
#                                 strata = phyloseq.bacteria.samples.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_nit_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.nitifiers.df , 
#                                 strata = phyloseq.bacteria.samples.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_nit_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS
#ENCLOSURE#
# Run the betadisper function, average distance to centroid
bray.enclosure.nit.disp <- betadisper(phyloseq.bacteria.samples.nitifiers.ra.bray, 
                                  phyloseq.bacteria.samples.nitifiers.df$Enclosure)
bray.enclosure.nit.disp
##Then test by permuting
set.seed(98)
bray.enclosure.nit.permdisp <- permutest(bray.enclosure.nit.disp, permutations = 9999)
bray.enclosure.nit.permdisp ##S, p = 6e-04

# #COPPER#
# # Run the betadisper function, average distance to centroid
# bray.copper.disp <- betadisper(phyloseq.bacteria.samples.nitifiers.ra.bray, 
#                                   phyloseq.bacteria.samples.nitifiers.df$Copper_level_mg_L)
# bray.copper.disp
# ##Then test by permuting
# set.seed(98)
# bray.copper.permdisp <- permutest(bray.copper.disp, permutations = 9999)
# bray.copper.permdisp ##S, p = 0.044

# Extract R2 and p-values
R2_adonis_enclosure_nit <- enclosure_BC_nit_adonis$R2[1] 
pvalue_adonis_enclosure_nit<-  enclosure_BC_nit_adonis$`Pr(>F)`[1]
R2_mirkat_enclosure_nit <- mirkat_enclosure_nit$R2 
pvalue_mirkat_enclosure_nit<-  mirkat_enclosure_nit$p_values

#### PLOTS#######
## BC
enclosure_BC_nit_beta_div <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       color = "System", fill = "System") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
                                     colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = enclosure.palette,
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  scale_fill_manual(values = enclosure.palette, 
                    labels = c("H21" = "Naive", 
                               "P1" = "Established"))+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  annotate("text", x = -1.2, y = -0.65,
           label = "MiRKAT\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x =-1.2, y = -0.6,
           label = paste("R² = ", round(R2_mirkat_enclosure_nit* 100, 1), "%",
                         "\np = ", round(pvalue_mirkat_enclosure_nit, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_nit_beta_div
ggsave("enclosure_BC_nit_beta_div.png", 
       enclosure_BC_nit_beta_div, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 14)

##Including collection month as shape 
enclosure_BC_nit_beta_div_2 <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       shape = "Collection Month", 
       color = "System", 
       fill = "System") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2,
  #                  color = Enclosure), show.legend = F)+
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
                                     colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, 
                 colour = Enclosure,
                 shape = Collection_Month), 
             size = 5, 
             alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure), 
            colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = enclosure.palette,
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  scale_fill_manual(values = enclosure.palette, 
                    labels = c("H21" = "Naive", 
                               "P1" = "Established"))+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8))+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  annotate("text", x = -1.2, y = -0.65,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x =-1.2, y = -0.6,
           label = paste("R² = ", round(R2_adonis_enclosure_nit * 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure_nit, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_nit_beta_div_2
ggsave("enclosure_BC_nit_beta_div_2.png", 
       enclosure_BC_nit_beta_div_2, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 15)

#POSTER BIMS 2026
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Posters/enclosure_BC_nit_beta_div.png", 
       enclosure_BC_nit_beta_div+
         theme(plot.title = element_text(size = 50, face = "bold")), 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 19)

#Add date num and phase
phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs_2 <- phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs%>%
mutate(Date_num = case_when(
  Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07")),
  Enclosure == "P1"  ~ as.numeric(Collection_Date - as.Date("2023-11-14"))
)) %>%
  #Separate if wanting to do separate day scales
  mutate(Date_num_naive = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07"))),
    Date_num_established = case_when(
      Enclosure == "P1" ~ as.numeric(Collection_Date - as.Date("2023-11-14"))))%>%
  mutate(Date_num_phase_naive = case_when(
    Enclosure == "H21" & Date_num >= 0   & Date_num <= 10  ~ "Phase 1 (Day 0 - 10)",
    Enclosure == "H21" & Date_num >= 11  & Date_num <= 28  ~ "Downtime Phase 1 - 2 (Day 11 - 28)",
    Enclosure == "H21" & Date_num >= 29  & Date_num <= 56  ~ "Phase 2 (Day 29 - 56)",
    Enclosure == "H21" & Date_num >= 57  & Date_num <= 86  ~ "Downtime Phase 2 - 3 (Day 57 - 86)",
    Enclosure == "H21" & Date_num >= 87  & Date_num <= 131 ~ "Phase 3 (Day 87 - 131)",
    Enclosure == "H21" & Date_num >= 132 & Date_num <= 147 ~ "Post-Treatment completion (Day 132 - 147)"
  )) %>%
  mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
                                       levels = c("Phase 1 (Day 0 - 10)", 
                                                  "Downtime Phase 1 - 2 (Day 11 - 28)", 
                                                  "Phase 2 (Day 29 - 56)", 
                                                  "Downtime Phase 2 - 3 (Day 57 - 86)", 
                                                  "Phase 3 (Day 87 - 131)", 
                                                  "Post-Treatment completion (Day 132 - 147)")))%>%
  mutate(Date_num_phase_established = case_when(
    Enclosure == "P1" & Date_num >= 0   & Date_num <= 50  ~ "Phase 1 (Day 0 - 50)",
    Enclosure == "P1" & Date_num >= 51  & Date_num <= 65  ~ "Phase 2 (Day 51 - 65)",
    Enclosure == "P1" & Date_num >= 66  & Date_num <= 98  ~ "Phase 3 (Day 66 - 98)",
    Enclosure == "P1" & Date_num >= 99  & Date_num <= 168 ~ "Post-Treatment completion (Day 99 - 168)"
  ))%>%
  mutate(Date_num_phase_established = factor(Date_num_phase_established, 
                                             levels = c("Phase 1 (Day 0 - 50)",
                                                        "Phase 2 (Day 51 - 65)",
                                                        "Phase 3 (Day 66 - 98)",
                                                        "Post-Treatment completion (Day 99 - 168)")))

#Color by day num
enclosure_BC_nit_beta_div_3 <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs_2) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       shape = "System", 
       color = "Day") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, group = Enclosure), 
               #color = "black",
               alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, 
                 colour = Date_num,
                 shape = Enclosure), 
             size = 5, 
             alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure), size = 10) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure),
            colour= "white", size = 5, fontface = "bold") +
  scale_shape_manual(values = c(15,16),
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  scale_color_viridis_c(option = "viridis")+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 20),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  annotate("text", x = -1.3, y = 0.51,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1.3, y = 0.58,
           label = paste("R² = ", round(R2_adonis_enclosure_nit * 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure_nit, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_nit_beta_div_3 

ggsave("enclosure_BC_nit_beta_div_3.png", 
       enclosure_BC_nit_beta_div_3, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 15)

#Color by day num, different scales per system
enclosure_BC_nit_beta_div_4 <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs_2) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", 
       title= "NITRIFYING TAXA", 
       shape = "System") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
  stat_ellipse(
    geom = "polygon",
    aes(x = MDS1, y = MDS2, fill = Enclosure, color = Enclosure, group = Enclosure),
    alpha = 0.2,
    lty = 2,
    linewidth = 1,
    level = 0.95,
    show.legend = FALSE
  )+
  #Problem with layering. If adding centroids here, dots go on top. 
  geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure,
                 color = Enclosure), size = 10,
             show.legend = FALSE) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Enclosure),
            colour= "white", size = 5, fontface = "bold") +
  scale_fill_manual(values = c(
    "#fc8d62",
    "#8da0cb"
  ))+
  scale_color_manual(values = c(
    "#fc8d62",
    "#8da0cb"
  ))+
  new_scale_color()+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  #Naive system individual dots
  geom_point(data = phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs_2%>%
               filter(Enclosure == "H21"), aes(x=MDS1, y=MDS2, 
                                               colour = Date_num_naive,
                                               shape = Enclosure), 
             size = 5) + # individuals
  scale_color_gradient2(low = "#fad2c2", 
                        mid = "#fc8d62", 
                        high = "#ad431a",
                        midpoint = 74,
                        limits = c(0,150),
                        breaks = c(0, 50, 100, 150))+
  # scale_color_viridis_c(option = "viridis",
  #                       limits = c(0,150),
  #                       breaks = c(0, 50, 100, 150))+
  labs(color = "Days (Naive system)")+
  #Established system individual dots. 
  new_scale_color()+
  geom_point(data = phyloseq.bacteria.samples.nitifiers.ra.bray.enclosure.segs_2 %>%
               filter(Enclosure == "P1"), aes(x=MDS1, y=MDS2, 
                                              colour = Date_num_established,
                                              shape = Enclosure), 
             size = 5) + # individuals
  scale_color_gradient2(low = "#c6cfe6",
                        mid = "#8da0cb", 
                        high = "#012983", 
                        midpoint = 84,
                        limits = c(0,170),
                        breaks = c(0, 50, 100, 170))+
  labs(color = "Days (Established system)")+
  # scale_color_viridis_c(option = "rocket", 
  #                       limits = c(0,170),
  #                       breaks = c(0, 50, 100, 170))+
  scale_shape_manual(values = c(15,16),
                     labels = c("H21" = "Naive", 
                                "P1" = "Established"))+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 18),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    color = guide_colorbar(order = 3),
    #color = guide_legend(order = 3),
    shape = guide_legend(override.aes = list(size = 7),
                         order = 1)) +
  # new_scale_color()+
  # geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure,
  #                color = Enclosure), size = 10,
  #            show.legend = FALSE) + # centroids +
  # geom_text(aes (x= cMDS1, y = cMDS2,
  #                label= Enclosure),
  #           colour= "white", size = 5, fontface = "bold") +
  # scale_color_manual(values = c(
  #   "#fc8d62",
  #   "#8da0cb"
  # ))+
  annotate("text", x = -1.3, y = 0.51,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1.3, y = 0.58,
           label = paste("R² = ", round(R2_adonis_enclosure_nit * 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure_nit, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_nit_beta_div_4
ggsave("enclosure_BC_nit_beta_div_4.png", 
       enclosure_BC_nit_beta_div_4, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 15)

# Run envfit on your NMDS object 
# set.seed(98)
# fit <- envfit(phyloseq.bacteria.samples.nitifiers.ra.bray.ord, 
#               phyloseq.bacteria.samples.nitifiers.df$Copper_level_mg_L, 
#               permutations = 999)
# 
# # Extract coordinates for the copper vector
# vec <- as.data.frame(fit$vectors$arrows * fit$vectors$r)  # scale by r
# colnames(vec) <- c("xend", "yend")
# vec$label <- rownames(vec)  # will be "Copper_level_mg_L"
# 
# enclosure_BC_nit_beta_div_copper <- enclosure_BC_nit_beta_div +
#   geom_segment(
#     data = vec,
#     aes(x = 0, y = 0, xend = xend, yend = yend),
#     arrow = arrow(length = unit(0.3, "cm")),  # adds arrowhead
#     colour = "red",
#     size = 1.2
#   )
# enclosure_BC_nit_beta_div_copper

###Centroids by copper quartiles facetted by enclosure #####
## BC
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.nitifiers.ra.bray.cent.copper_quart <- aggregate(
  cbind(MDS1, MDS2) ~ Enclosure + Copper_quartile,
  data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs,
  FUN = mean)

#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.nitifiers.ra.bray.copper_quart.segs <- merge(phyloseq.bacteria.samples.nitifiers.ra.bray.scrs, 
                                                             setNames(phyloseq.bacteria.samples.nitifiers.ra.bray.cent.copper_quart, 
                                                                      c("Enclosure", "Copper_quartile",
                                                                        "cMDS1","cMDS2")),
                                                             by = c("Enclosure", "Copper_quartile"), 
                                                             sort = F)
##PERMANOVA - Just copper###
# set.seed(98)
# copper_quart_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.bray ~ Copper_quartile, 
#                                    phyloseq.bacteria.samples.nitifiers.df, 
#                                    strata = phyloseq.bacteria.samples.nitifiers.df$Collection_Month,
#                                    by = "margin",
#                                    permutations = 9999)
# copper_quart_BC_nit_adonis #30.1% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.nitifiers.df, 
#                                 strata = phyloseq.bacteria.samples.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_nit_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.nitifiers.df , 
#                                 strata = phyloseq.bacteria.samples.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_nit_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS#
#COPPER quartile#
# Run the betadisper function, average distance to centroid
bray.nit.copper.quart.disp <- betadisper(phyloseq.bacteria.samples.nitifiers.ra.bray,
                                     phyloseq.bacteria.samples.nitifiers.df$Copper_quartile)
bray.nit.copper.quart.disp
##Then test by permuting
set.seed(98)
bray.nit.copper.quart.permdisp <- permutest(bray.nit.copper.quart.disp, permutations = 9999)
bray.nit.copper.quart.permdisp ##S, p = 0.0017

# Extract R2 and p-values
# R2_adonis_copper_quart_nit <- copper_quart_BC_nit_adonis$R2[1] 
# pvalue_adonis_copper_quart_nit<-  copper_quart_BC_nit_adonis$`Pr(>F)`[1]
# R2_adonis_copper <- enclosure_BC_nit_adonis$R2[2] 
# pvalue_adonis_copper <- enclosure_BC_nit_adonis$`Pr(>F)`[2]

#### PLOTS
## BC
copper_quart_BC_nit_beta_div <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.copper_quart.segs) + 
  theme_bw() +
  facet_grid(~Enclosure,
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       color = "Copper (mg/L) Quartile", 
       fill = "Copper (mg/L) Quartile") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, 
                                     fill= Copper_quartile.abvr, 
                                     colour = Copper_quartile.abvr), alpha = 0.2, 
               lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Copper_quartile.abvr), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Copper_quartile.abvr), size = 10) + # centroids
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Copper_quartile.abvr), colour= "white", 
            size = 6, fontface = "bold") +
  scale_color_manual(values = plasma_quartiles)+
  scale_fill_manual(values = plasma_quartiles)+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold"))
# annotate("text", x = 0.4, y = 0.9,
#          label = "Enclosure",
#          hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
# annotate("text", x =0.4, y = 0.9,
#          label = paste("R² = ", round(R2_adonis_enclosure_nit* 100, 1), "%",
#                        "\np = ", round(pvalue_adonis_enclosure_nit, 4)),
#          hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
copper_quart_BC_nit_beta_div
ggsave("copper_quart_BC_nit_beta_div.png", 
       copper_quart_BC_nit_beta_div, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 14)


###Continuous copper facetted by enclosure #####
#No samples without copper levels:
phyloseq.bacteria.samples.nitifiers.ra.copper <- subset_samples(phyloseq.bacteria.samples.nitifiers.ra, 
                                                                !is.na(Copper_level_mg_L))
phyloseq.bacteria.samples.nitifiers.ra.copper <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.nitifiers.ra.copper)> 0, 
                                                            phyloseq.bacteria.samples.nitifiers.ra.copper)
phyloseq.bacteria.samples.nitifiers.ra.copper #827 taxa and 204 samples

####H21 #######
phyloseq.bacteria.samples.nitifiers.ra.copper.H21 <- subset_samples(phyloseq.bacteria.samples.nitifiers.ra.copper, 
                                                                Enclosure == "H21")
phyloseq.bacteria.samples.nitifiers.ra.copper.H21 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.nitifiers.ra.copper.H21)> 0, 
                                                                phyloseq.bacteria.samples.nitifiers.ra.copper.H21)
phyloseq.bacteria.samples.nitifiers.ra.copper.H21 # 824 taxa and 82 samples 

##BC 
phyloseq.bacteria.samples.nitifiers.ra.copper.H21.bray <- vegdist(t(phyloseq.bacteria.samples.nitifiers.ra.copper.H21@otu_table), 
                                                              method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.nitifiers.copper.H21.df <- as(phyloseq.bacteria.samples.nitifiers.ra.copper.H21@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date),
    #Date number 
    Date_num = as.numeric(Collection_Date - min(Collection_Date)),
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.nitifiers.copper.H21.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_H21 <- D2K(
  as.matrix(phyloseq.bacteria.samples.nitifiers.ra.copper.H21.bray))

#Jaccard distances
phyloseq.bacteria.samples.nitifiers.ra.copper.H21.jac <- vegdist(t(phyloseq.bacteria.samples.nitifiers.ra.copper.H21@otu_table), 
                                                                 method = "jaccard")
phyloseq.bacteria.samples.nitifiers.ra.copper.H21.jac

jac_kernel_nit_H21 <- D2K(
  as.matrix(phyloseq.bacteria.samples.nitifiers.ra.copper.H21.jac))

# # Covariates (month)
# month_covariate_nit_H21 <- model.matrix(~ Collection_Month, 
#                                     data = phyloseq.bacteria.samples.nitifiers.copper.H21.df)[, -1]
# 
# month_covariate_nit_H21

# Covariates (date_num)
datenum_covariate_H21_nit <- phyloseq.bacteria.samples.nitifiers.copper.H21.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_H21_nit

#Collection month for random intercepts
collection_month_H21_nit <- phyloseq.bacteria.samples.nitifiers.copper.H21.df%>%
  pull(Collection_Month)%>%
  as.matrix()
collection_month_H21_nit

#Run model
mirkat_nit_H21_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.nitifiers.copper.H21.df$Copper_level_mg_L,
  X = datenum_covariate_H21_nit,
  Ks = bray_kernel_nit_H21,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_nit_H21_copper_cont #R2 3.1% , p_values = 0.004

#####PERMANOVA #######
set.seed(98)
adonis_nit_H21_copper_cont  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.copper.H21.bray ~ Copper_level_mg_L, 
                                       phyloseq.bacteria.samples.nitifiers.copper.H21.df, 
                                       strata = phyloseq.bacteria.samples.nitifiers.copper.H21.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_nit_H21_copper_cont # R2 9.7%, p value 0.027*

##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.nit.copper.cont.disp.H21 <- betadisper(phyloseq.bacteria.samples.nitifiers.ra.copper.H21.bray,
                                           phyloseq.bacteria.samples.nitifiers.copper.H21.df$Copper_quartile)
bray.nit.copper.cont.disp.H21
##Then test by permuting
set.seed(98)
bray.nit.copper.cont.permdisp.H21 <- permutest(bray.nit.copper.cont.disp.H21, permutations = 9999)
bray.nit.copper.cont.permdisp.H21 ##NS, p = 0.65

####P1 #######
phyloseq.bacteria.samples.nitifiers.ra.copper.P1 <- subset_samples(phyloseq.bacteria.samples.nitifiers.ra.copper, 
                                                                    Enclosure == "P1")
phyloseq.bacteria.samples.nitifiers.ra.copper.P1 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.nitifiers.ra.copper.P1)> 0, 
                                                                phyloseq.bacteria.samples.nitifiers.ra.copper.P1)
phyloseq.bacteria.samples.nitifiers.ra.copper.P1 #624 taxa and 122 samples

##BC 
phyloseq.bacteria.samples.nitifiers.ra.copper.P1.bray <- vegdist(t(phyloseq.bacteria.samples.nitifiers.ra.copper.P1@otu_table), 
                                                                  method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.nitifiers.copper.P1.df <- as(phyloseq.bacteria.samples.nitifiers.ra.copper.P1@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date),
    #Date number 
    Date_num = as.numeric(Collection_Date - min(Collection_Date)),
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.nitifiers.copper.P1.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_P1 <- D2K(
  as.matrix(phyloseq.bacteria.samples.nitifiers.ra.copper.P1.bray))

# # Covariates (month)
# month_covariate_nit_P1 <- model.matrix(~ Collection_Month, 
#                                         data = phyloseq.bacteria.samples.nitifiers.copper.P1.df)[, -1]
# 
# month_covariate_nit_P1

# Covariates (date_num)
datenum_covariate_P1_nit <- phyloseq.bacteria.samples.nitifiers.copper.P1.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_P1_nit

#Run model
mirkat_nit_P1_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.nitifiers.copper.P1.df$Copper_level_mg_L,
  X = datenum_covariate_P1_nit,
  Ks = bray_kernel_nit_P1,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_nit_P1_copper_cont #p val 0 , R2 20.6%

#####PERMANOVA #####
set.seed(98)
adonis_nit_P1_copper_cont  <- adonis2(phyloseq.bacteria.samples.nitifiers.ra.copper.P1.bray ~ Copper_level_mg_L, 
                                       phyloseq.bacteria.samples.nitifiers.copper.P1.df, 
                                       strata = phyloseq.bacteria.samples.nitifiers.copper.P1.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_nit_P1_copper_cont # R2 26.9%, p value 1e-04*
##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.nit.copper.cont.disp.P1 <- betadisper(phyloseq.bacteria.samples.nitifiers.ra.copper.P1.bray,
                                           phyloseq.bacteria.samples.nitifiers.copper.P1.df$Copper_quartile)
bray.nit.copper.cont.disp.P1
##Then test by permuting
set.seed(98)
bray.nit.copper.cont.permdisp.P1 <- permutest(bray.nit.copper.cont.disp.P1, permutations = 9999)
bray.nit.copper.cont.permdisp.P1 ##S, p = 1e-04 

###MIRKAT GLMM FOR COPPER #####
#make DF from metadata
phyloseq.bacteria.samples.nitifiers.copper.all.df <- as(phyloseq.bacteria.samples.nitifiers.ra.copper@sam_data, "data.frame") %>%
  mutate(
    # Collection date as Date object
    Collection_Date = as.Date(Collection_Date))%>%
  group_by(Enclosure)%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  mutate(
    # Collection date as factor
    Collection_Date_fact = factor(Collection_Date), 
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month),  # convert to factor for PERMANOVA
    # Copper quartiles
    Copper_quartile = factor(
      ntile(Copper_level_mg_L, 4),
      labels = c("Q1_low", "Q2_midlow", "Q3_midhigh", "Q4_high")),
    Copper_quartile.abvr = case_when(
      Copper_quartile == "Q1_low" ~ "Q1", 
      Copper_quartile == "Q2_midlow" ~ "Q2", 
      Copper_quartile == "Q3_midhigh" ~ "Q3", 
      Copper_quartile == "Q4_high" ~ "Q4"))
phyloseq.bacteria.samples.nitifiers.copper.all.df

#Bray-Curtis distances
phyloseq.bacteria.samples.nitifiers.ra.copper.bray <- vegdist(t(phyloseq.bacteria.samples.nitifiers.ra.copper@otu_table), 
                                                              method = "bray")
phyloseq.bacteria.samples.nitifiers.ra.copper.bray

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.nitifiers.ra.copper.bray))

#Jaccard distances
phyloseq.bacteria.samples.nitifiers.ra.copper.jac <- vegdist(t(phyloseq.bacteria.samples.nitifiers.ra.copper@otu_table), 
                                                             method = "jaccard")
phyloseq.bacteria.samples.nitifiers.ra.copper.jac

jac_kernel_nit_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.nitifiers.ra.copper.jac))

# Covariates (date_num)
datenum_covariate_nit_all <- phyloseq.bacteria.samples.nitifiers.copper.all.df%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_nit_all

#Id for repeated measures 
Enclosure_d_nit_mirkat_all <- phyloseq.bacteria.samples.nitifiers.copper.all.df %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_nit_mirkat_all

#Run model - dont think it works. just two enclosures and enclosure id is the same as the outcome 
mirkat_enclosure_nit_glm <- GLMMMiRKAT(
  y = phyloseq.bacteria.samples.nitifiers.copper.all.df$Copper_level_mg_L,
  X = datenum_covariate_nit_all,
  Ks = list(bray_kernel_nit_all, jac_kernel_nit_all),
  id = Enclosure_d_nit_mirkat_all,
  model  = "gaussian",               # "gaussian" if y is continuous; "poisson" for counts
  slope  = FALSE,                    # random intercept only
  nperm  = 5000)                   # (not n.perm)
mirkat_enclosure_nit_glm #0.00019996 omnibus_p 

###PLOT FOR BOTH ENCLOSURES#######
#### Extract R2 and p-values
#####H21
##Mirkat
R2_mirkat_copper_cont_nit_H21 <- mirkat_nit_H21_copper_cont$R2
pvalue_mirkat_copper_cont_nit_H21<- mirkat_nit_H21_copper_cont$p_values

##PERMANOVA
R2_adonis_copper_cont_nit_H21 <- adonis_nit_H21_copper_cont$R2[1]
pvalue_adonis_copper_cont_nit_H21<- adonis_nit_H21_copper_cont$`Pr(>F)`[1]

#####P1
##Mirkat
R2_mirkat_copper_cont_nit_P1 <- mirkat_nit_P1_copper_cont$R2
pvalue_mirkat_copper_cont_nit_P1<- mirkat_nit_P1_copper_cont$p_values

##PERMANOVA
R2_adonis_copper_cont_nit_P1 <- adonis_nit_P1_copper_cont$R2[1]
pvalue_adonis_copper_cont_nit_P1<- adonis_nit_P1_copper_cont$`Pr(>F)`[1]

#### PLOTS
## BC
copper_cont_BC_nit_beta_div <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.scrs) + 
  theme_bw() +
  facet_wrap(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, 
  #                                    fill= Copper_quartile.abvr, 
  #                                    colour = Copper_quartile.abvr), alpha = 0.2, 
  #              lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, color = Copper_level_mg_L), size = 5, alpha = 0.8) + # individuals
  scale_color_viridis_c(option = "plasma")+
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       color = "Copper levels (mg/L)") +
  theme(legend.position = "bottom",
        legend.direction = "vertical", 
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold")) +
  #H21
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
    aes(x = -1, y = 0.6, label = "PERMANOVA\nCopper levels (mg/L)"),
    hjust = 0.5, vjust = -0.5,
    size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.4, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0.13, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_P1* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_P1, 6))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")
copper_cont_BC_nit_beta_div
ggsave("copper_cont_BC_nit_beta_div.png", 
       copper_cont_BC_nit_beta_div, 
       device = "png", 
       dpi = 600, 
       height = 9, 
       width = 16)

##Addinc collection month as shape 
copper_cont_BC_nit_beta_div_2 <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.scrs) + 
  theme_bw() +
  facet_wrap(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  # stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, 
  #                                    fill= Copper_quartile.abvr, 
  #                                    colour = Copper_quartile.abvr), alpha = 0.2, 
  #              lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, color = Copper_level_mg_L,
                 shape = Collection_Month), size = 5, alpha = 0.8) + # individuals
  scale_color_viridis_c(option = "plasma")+
  scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8))+
  # scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8), 
  #                    labels = c("2023-10" = "Oct 2023",
  #                               "2023-11" = "Nov 2023",
  #                               "2023-12" = "Dec 2023",
  #                               "2024-01" = "Jan 2024",
  #                               "2024-02" = "Feb 2024",
  #                               "2024-03" = "Mar 2024",
  #                               "2024-04" = "Apr 2024",
  #                               "2024-05" = "May 2024"))+
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       color = "Copper levels (mg/L)",
       shape = "Collection Month") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold")) +
  #H21
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = -0.65, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = -0.88, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = -0.2, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = -0.38, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_P1* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_P1, 6))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  guides(color = guide_colorbar(direction = "vertical"))
copper_cont_BC_nit_beta_div_2
ggsave("copper_cont_BC_nit_beta_div_2.png", 
       copper_cont_BC_nit_beta_div_2, 
       device = "png", 
       dpi = 600, 
       height = 9, 
       width = 16)

#POSTER BIMS 2026
ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Posters/copper_cont_BC_nit_beta_div.png", 
       copper_cont_BC_nit_beta_div + 
         theme(legend.position = "none",
               plot.title = element_text(colour = "black", size = 50, face = "bold"),
               strip.text = element_text(colour = "white", size = 45, face = "bold")), 
       device = "png", 
       dpi = 600, 
       height = 9.5, 
       width = 15.5)



#Adding date_num and phases
phyloseq.bacteria.samples.nitifiers.ra.bray.scrs_2 <- phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
  mutate(Date_num = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07")),
    Enclosure == "P1"  ~ as.numeric(Collection_Date - as.Date("2023-11-14"))
  )) %>%
  #Separate if wanting to do separate day scales
  mutate(Date_num_naive = case_when(
    Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07"))),
    Date_num_established = case_when(
      Enclosure == "P1" ~ as.numeric(Collection_Date - as.Date("2023-11-14"))))%>%
  mutate(Date_num_phase_naive = case_when(
    Enclosure == "H21" & Date_num >= 0   & Date_num <= 10  ~ "Phase 1 (Day 0 - 10)",
    Enclosure == "H21" & Date_num >= 11  & Date_num <= 28  ~ "Downtime Phase 1 - 2 (Day 11 - 28)",
    Enclosure == "H21" & Date_num >= 29  & Date_num <= 56  ~ "Phase 2 (Day 29 - 56)",
    Enclosure == "H21" & Date_num >= 57  & Date_num <= 86  ~ "Downtime Phase 2 - 3 (Day 57 - 86)",
    Enclosure == "H21" & Date_num >= 87  & Date_num <= 131 ~ "Phase 3 (Day 87 - 131)",
    Enclosure == "H21" & Date_num >= 132 & Date_num <= 147 ~ "Post-treatment (Day 132 - 147)"
  )) %>%
  mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
                                       levels = c("Phase 1 (Day 0 - 10)", 
                                                  "Downtime Phase 1 - 2 (Day 11 - 28)", 
                                                  "Phase 2 (Day 29 - 56)", 
                                                  "Downtime Phase 2 - 3 (Day 57 - 86)", 
                                                  "Phase 3 (Day 87 - 131)", 
                                                  "Post-treatment (Day 132 - 147)")))%>%
  mutate(Date_num_phase_established = case_when(
    Enclosure == "P1" & Date_num >= 0   & Date_num <= 50  ~ "Phase 1 (Day 0 - 50)",
    Enclosure == "P1" & Date_num >= 51  & Date_num <= 65  ~ "Phase 2 (Day 51 - 65)",
    Enclosure == "P1" & Date_num >= 66  & Date_num <= 98  ~ "Phase 3 (Day 66 - 98)",
    Enclosure == "P1" & Date_num >= 99  & Date_num <= 168 ~ "Post-treatment (Day 99 - 168)"
  ))%>%
  mutate(Date_num_phase_established = factor(Date_num_phase_established, 
                                             levels = c("Phase 1 (Day 0 - 50)",
                                                        "Phase 2 (Day 51 - 65)",
                                                        "Phase 3 (Day 66 - 98)",
                                                        "Post-treatment (Day 99 - 168)")))
##Adding phases of copper exposure 
copper_cont_BC_nit_beta_div_3 <- ggplot(phyloseq.bacteria.samples.nitifiers.ra.bray.scrs_2) + 
  theme_bw() +
  facet_wrap(~Enclosure,
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING TAXA", 
       color = "Copper levels (mg/L)") +
  # geom_segment(aes(x=cMDS1, y=cMDS2,
  #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_point(aes(x=MDS1, y=MDS2, colour = Copper_level_mg_L, shape = Date_num_phase_naive), 
             size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(15, 3, 16, 4, 17, 18), name = "Naive system") +
  guides(
    shape = guide_legend(override.aes = list(size = 7), ncol = 1)) +
  new_scale("shape") +
  geom_point(aes(x=MDS1, y=MDS2, colour = Copper_level_mg_L, shape = Date_num_phase_established), 
             size = 5, alpha = 0.8) +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Established system") +
  scale_color_viridis_c(option = "plasma")+
  theme(legend.position = "bottom",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 20),
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 38, face = "bold"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))) +
  #H21
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.3, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0.1, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = -0.05, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_P1* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_P1, 6))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  guides(color = guide_colorbar(direction = "vertical",
                                order= 1))
copper_cont_BC_nit_beta_div_3
ggsave("copper_cont_BC_nit_beta_div_3.png", 
       copper_cont_BC_nit_beta_div_3, 
       device = "png", 
       dpi = 600, 
       height = 9, 
       width = 16)


#NITRIFIERS######
#No samples without copper levels:
nitrifiers.copper <- subset_samples(nitrifiers, !is.na(Copper_level_mg_L))
nitrifiers.copper <- prune_taxa(taxa_sums(nitrifiers.copper)> 0, 
                                               nitrifiers.copper)

nitrifiers.copper.ra <- transform_sample_counts(nitrifiers.copper, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data

#make DF from metadata
nitrifiers.copper.df <- as(nitrifiers.copper@sam_data, "data.frame") %>%
  mutate(
    # Convert Collection_Date to Date if it isn’t already
    Collection_Date = as.Date(Collection_Date),
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month)  # convert to factor for PERMANOVA
  )

nitrifiers.copper.df 

###BRAY CURTIS#####
nitrifiers.ra.bray <- vegdist(t(nitrifiers.ra@otu_table), method = "bray")
nitrifiers.ra.bray

# ###ORDINATION####
# set.seed(98)
# nitrifiers.ra.bray.ord <- metaMDS(nitrifiers.ra.bray, k=3, try = 50, 
#                                                  trymax = 1000,
#                                                  autotransform = F)
# #### ADDING CENTROIDS FOR PLOTTING
# ## BC
# #Simple ordination plot
# nitrifiers.ra.bray.plot <- ordiplot(nitrifiers.ra.bray.ord$points)
# 
# #Now, extract coordinates
# nitrifiers.ra.bray.scrs <- scores(nitrifiers.ra.bray.plot, display = "sites")
# #Add metadata to coordinates
# nitrifiers.ra.bray.scrs <- cbind(as.data.frame(nitrifiers.ra.bray.scrs),
#                                                 Copper_level_mg_L = nitrifiers.copper.df$Copper_level_mg_L, 
#                                                 SampleID = nitrifiers.copper.df$SampleID, 
#                                                 Collection_Date = nitrifiers.copper.df$Collection_Date,
#                                                 Collection_Month = nitrifiers.copper.df$Collection_Month, 
#                                                 Enclosure = nitrifiers.copper.df$Enclosure)
# ##Calculate centroids according to Enclosure
# nitrifiers.ra.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
#                                                     data = nitrifiers.ra.bray.scrs, 
#                                                     FUN = mean) 
# #Merge centroids with coordinates and metadata
# nitrifiers.ra.bray.segs <- merge(nitrifiers.ra.bray.scrs, 
#                                                 setNames(nitrifiers.ra.bray.cent, c("Enclosure","cMDS1","cMDS2")),
#                                                 by = 'Enclosure', 
#                                                 sort = F)
# 
# 
# ##PERMANOVA###
# set.seed(98)
# nit_enclosure_copper_BC_adonis  <- adonis2(nitrifiers.ra.bray ~ Enclosure +
#                                          Copper_level_mg_L, 
#                                        nitrifiers.copper.df, 
#                                        strata = nitrifiers.copper.df$Collection_Month,
#                                        by = "margin",
#                                        permutations = 9999)
# nit_enclosure_copper_BC_adonis #4.8% of the variation is due to Enclosure, p = 1e-04
# #3.5% of the variation is due to copper levels, p = 8e-04
# 
# 
# #With interaction
# set.seed(98)
# nit_enclosure_copper_interac_BC_adonis  <- adonis2(nitrifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                nitrifiers.copper.df , 
#                                                strata = nitrifiers.copper.df$Collection_Month,
#                                                by = "margin",
#                                                permutations = 9999)
# nit_enclosure_copper_interac_BC_adonis
# #Enclosure:Copper_level_mg_L interaction significant (p = 1e-04). 7.3% of variation
# 
# ##PERMDISPS
# #ENCLOSURE#
# # Run the betadisper function, average distance to centroid
# bray.enclosure_nit.disp <- betadisper(nitrifiers.ra.bray, 
#                                   nitrifiers.copper.df$Enclosure)
# bray.enclosure_nit.disp
# ##Then test by permuting
# set.seed(98)
# bray.enclosure_nit.permdisp <- permutest(bray.enclosure_nit.disp, permutations = 9999)
# bray.enclosure_nit.permdisp ##NS, p = 0.65
# 
# #COPPER#
# # Run the betadisper function, average distance to centroid
# bray.copper_nit.disp <- betadisper(nitrifiers.ra.bray, 
#                                nitrifiers.copper.df$Copper_level_mg_L)
# bray.copper_nit.disp
# ##Then test by permuting
# set.seed(98)
# bray.copper_nit.permdisp <- permutest(bray.copper_nit.disp, permutations = 9999)
# bray.copper_nit.permdisp ##NS, p = 0.49
# 
# # Extract R2 and p-values
# R2_adonis_enclosure_nit <- nit_enclosure_copper_BC_adonis$R2[1] 
# pvalue_adonis_enclosure_nit <-  nit_enclosure_copper_BC_adonis$`Pr(>F)`[1]
# R2_adonis_copper_nit <- nit_enclosure_copper_BC_adonis$R2[2] 
# pvalue_adonis_copper_nit <- nit_enclosure_copper_BC_adonis$`Pr(>F)`[2]
# 
# R2_adonis_enclosurebynit <- nit_enclosure_copper_interac_BC_adonis$R2[1] 
# pvalue_adonis_enclosurebynit <-  nit_enclosure_copper_interac_BC_adonis$`Pr(>F)`[1]
# 
# #### PLOTS
# ## BC
# enclosure_nit_BC_beta_div <- ggplot(nitrifiers.ra.bray.segs) + theme_bw() +
#   labs(x="NMDS1", y="NMDS2", title= "NITRIFYING COMMUNITIES", shape = "Enclosure") +
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
#                                      colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
#   geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure), size = 5, alpha = 0.8) + # individuals
#   geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) + # centroids
#   #geom_text(aes (x= MDS1, y = MDS2,label= Collection_Date), colour= "black", size = 2.8, fontface = "bold") +
#   geom_text(aes (x= cMDS1, y = cMDS2,label= Enclosure), colour= "white", size = 2.8, fontface = "bold") +
#   scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
#   scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
#   theme(legend.position = "bottom",
#         legend.title = element_blank(),
#         legend.text = element_text(colour = "black", size = 22),
#         plot.title = element_text(size = 36),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth = 1),
#         axis.ticks = element_line(colour = "black", linewidth = 0.75),
#         axis.title = element_text(size = 28),
#         axis.text = element_text(size = 20, colour = "black"))+
#   guides(
#     shape = guide_legend(override.aes = list(size = 7))
#   )+
#   annotate("text", x = 0.8, y = 1.7, 
#            label = "Enclosure", 
#            hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
#   annotate("text", x = 0.8, y = 1.7, 
#            label = paste("R² = ", round(R2_adonis_enclosure_nit* 100, 1), "%",
#                          "\np = ", round(pvalue_adonis_enclosure_nit, 4)), 
#            hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+# Annotate R² and p-values
#   annotate("text", x = 0.8, y = 1.3, 
#            label = "Copper levels in water", 
#            hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Copper)
#   annotate("text", x = 0.8, y = 1.3, 
#            label = paste("R² = ", round(R2_adonis_copper_nit* 100, 1), "%",
#                          "\np = ", round(pvalue_adonis_copper_nit, 4)), 
#            hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values 
#   annotate("text", x = 0.8, y = 0.9, 
#            label = "Enclosure:Copper levels in water", 
#            hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Copper:Enclosure)
#   annotate("text", x = 0.8, y = 0.9, 
#            label = paste("R² = ", round(R2_adonis_enclosurebynit* 100, 1), "%",
#                          "\np = ", round(pvalue_adonis_enclosurebynit, 4)), 
#            hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values 
# 
# enclosure_nit_BC_beta_div
# 
# #Inclusing colection month as shape 
# enclosure_nit_BC_beta_div_2 <- ggplot(nitrifiers.ra.bray.segs) + 
#   theme_bw() +
#   labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
#        shape = "Collection Month", 
#        color = "System", 
#        fill = "System") +
#   # geom_segment(aes(x=cMDS1, y=cMDS2,
#   #                  xend= MDS1, yend = MDS2,
#   #                  color = Enclosure), show.legend = F)+
#   stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
#                                      colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   geom_point(aes(x=MDS1, y=MDS2, 
#                  colour = Enclosure,
#                  shape = Collection_Month), 
#              size = 5, 
#              alpha = 0.8) + # individuals
#   geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) + # centroids +
#   geom_text(aes (x= cMDS1, y = cMDS2,
#                  label= Enclosure), 
#             colour= "white", size = 6, fontface = "bold") +
#   scale_color_manual(values = enclosure.palette,
#                      labels = c("H21" = "Naive", 
#                                 "P1" = "Established"))+
#   scale_fill_manual(values = enclosure.palette, 
#                     labels = c("H21" = "Naive", 
#                                "P1" = "Established"))+
#   scale_shape_manual(values = c(15, 16, 17, 18, 3, 4, 8))+
#   theme(legend.position = "bottom",
#         legend.title = element_text(colour = "black", 
#                                     size = 22,
#                                     face = "bold"),
#         legend.text = element_text(colour = "black", size = 22),
#         plot.title = element_text(size = 45),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth = 1),
#         axis.ticks = element_line(colour = "black", linewidth = 0.75),
#         axis.title = element_text(size = 30),
#         axis.text = element_text(size = 26, colour = "black"))+
#   guides(
#     shape = guide_legend(override.aes = list(size = 7))) +
#   annotate("text", x = -1, y = -0.9,
#            label = "MiRKAT\nSystem",
#            hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
#   annotate("text", x = -1, y = -0.85,
#            label = paste("R² = ", round(R2_mirkat_enclosure* 100, 1), "%",
#                          "\np = ", round(pvalue_mirkat_enclosure, 4)),
#            hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
# enclosure_nit_BC_beta_div_2

# Run envfit on your NMDS object 
# set.seed(98)
# fit <- envfit(nitrifiers.ra.bray.ord, 
#               nitrifiers.copper.df$Copper_level_mg_L, 
#               permutations = 999)
# 
# # Extract coordinates for the copper vector
# vec <- as.data.frame(fit$vectors$arrows * fit$vectors$r)  # scale by r
# colnames(vec) <- c("xend", "yend")
# vec$label <- rownames(vec)  # will be "Copper_level_mg_L"
# 
# enclosure_BC_beta_div_copper <- enclosure_BC_beta_div +
#   geom_segment(
#     data = vec,
#     aes(x = 0, y = 0, xend = xend, yend = yend),
#     arrow = arrow(length = unit(0.3, "cm")),  # adds arrowhead
#     colour = "red",
#     size = 1.2
#   )


#RESISTOME#######
##Import count tables#######
ARGcounts <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplus_dev_branch/SNPconfirmed_AMR_analytic_matrix.csv')
# ARGcounts <- readr::read_csv(
#   '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplus_dev_branch_80gf/SNPconfirmed_AMR_analytic_matrix.csv')
colnames(ARGcounts)
rownames(ARGcounts)

##Annotations - downloaded from AMRplusplus#######
tax.table.ARG <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplus_dev_branch/Annotations/megares_annotations_v3.00.csv')
#There is a multi-compount category, fixing it to "Multi-compound"

#Edit to make compatible with phyloseq
tax.table.ARG <- tax.table.ARG%>%
  column_to_rownames(var = "header")%>% #make "header' rownames so it matches with the OTU table rownames
  rename_with(~ str_to_title(.))%>% #want the annotation cols with the first letter capitalized
  as.matrix() ##make into matrix so it is compatible with tax_table function from phyloseq



##OTU table ########
otu_table_ARG <- ARGcounts%>%
  column_to_rownames(var = "gene_accession")
#Just with samples part of the study
otu_table_ARG <- otu_table_ARG%>%
  select(-all_of(dropping_samples)) 
ncol(otu_table_ARG) #220 = 217 samples + 3 mocks 


##PHYLOSEQ####
OTU_ARG <-phyloseq::otu_table(otu_table_ARG, 
                              taxa_are_rows = TRUE)
TAX_ARG <-phyloseq::tax_table(tax.table.ARG)
phyloseq_ARG <- phyloseq(OTU_ARG, TAX_ARG, sampledata_phyloseq) 
phyloseq_ARG ##8732 taxa and 220 samples

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU_ARG), metadata$SampleID) #No



##PREPROCESSING ####
phyloseq_ARG #8733 taxa and 220 samples
min(sample_sums(phyloseq_ARG)) # 6797 (H21_0120)
max(sample_sums(phyloseq_ARG)) # 1,987,181  (ZymoMock1a_S142) 
mean(sample_sums(phyloseq_ARG)) #278,138.7
median(sample_sums(phyloseq_ARG)) #227,647.5
sort(sample_sums(phyloseq_ARG))

###Zymo#####
### pulling out samples from ZYMOs and EB, NTC and those samples with low OTUs
phyloseq_ARG.controls <- subset_samples(phyloseq_ARG, 
                                        grepl("NTC|EB|Zymo", sample_names(phyloseq_ARG)))
phyloseq_ARG.controls <- prune_taxa(taxa_sums(phyloseq_ARG.controls) > 0, phyloseq_ARG.controls) 
phyloseq_ARG.controls #4057 taxa and 3 samples (only zymos)

###Samples#####
##New phyloseq of just samples
phyloseq_ARG.samples <- subset_samples(phyloseq_ARG, 
                                       !grepl("NTC|EB|Zymo", sample_names(phyloseq_ARG)))
phyloseq_ARG.samples #8732 taxa and 217 samples 
#Taking out those with low counts:
sort(sample_sums(phyloseq_ARG.samples)) #H21_0120 has low counts, dropping it here as i did with kraken taxonomy
phyloseq_ARG.samples <- prune_samples(sample_sums(phyloseq_ARG.samples) > 10000, phyloseq_ARG.samples) 
phyloseq_ARG.samples <- prune_taxa(taxa_sums(phyloseq_ARG.samples) > 0, phyloseq_ARG.samples)
phyloseq_ARG.samples #8333 taxa and 216 samples (dropped H21_0120)
sort(sample_sums(phyloseq_ARG.samples)) #OK

#####H21####
phyloseq_ARG.samples_H21 <- subset_samples(phyloseq_ARG.samples, Enclosure == "H21")
phyloseq_ARG.samples_H21 <- prune_taxa(taxa_sums(phyloseq_ARG.samples_H21) > 0, 
                                             phyloseq_ARG.samples_H21)
phyloseq_ARG.samples_H21 #8304 taxa and 92 samples
range(phyloseq_ARG.samples_H21@sam_data$Collection_Date)#OK, "2023-10-09" through "2024-03-02"

#####P1####
phyloseq_ARG.samples_P1 <- subset_samples(phyloseq_ARG.samples, Enclosure == "P1")
phyloseq_ARG.samples_P1 <- prune_taxa(taxa_sums(phyloseq_ARG.samples_P1) > 0, 
                                            phyloseq_ARG.samples_P1)
phyloseq_ARG.samples_P1 #8269 taxa and 124 samples 
range(phyloseq_ARG.samples_P1@sam_data$Collection_Date)#OK, "2023-11-14" through "2024-04-30"

#Am I missing samples that have a kraken taxonomic classification?
setdiff(sample_names(phyloseq_ARG.samples), sample_names(phyloseq.bacteria.samples)) #Nope


###PERCENTAGE OF READS ALIGNED TO MEGARES GENES########
#Trimmed reads
cuso4_trimmed_read_counts_samples_metadata$Num_Reads_Forward_Trimmed_Paired

#Total aligned reads
total_aligned_reads <- sample_sums(phyloseq_ARG.samples)
total_aligned_reads_df <- data.frame(
  SampleID = names(total_aligned_reads),
  total_aligned_reads = as.numeric(total_aligned_reads))

#Merge with trimmed reads and metadata
Percentage_reads_aligned_MEGARes <- left_join(
  cuso4_trimmed_read_counts_samples_metadata,
  total_aligned_reads_df, 
  by = "SampleID")%>%
  #Calculate percentage reads aligned 
  mutate(Percentage_reads_aligned = (total_aligned_reads / Num_Reads_Forward_Trimmed_Paired) * 100)
nrow(Percentage_reads_aligned_MEGARes) #216 samples

#Descriptive stats per group
Percentage_reads_aligned_MEGARes %>%
  group_by(Enclosure)%>%
  summarise(mean_percentage_aligned_reads = mean(Percentage_reads_aligned),
            sd_percentage_aligned_reads = sd(Percentage_reads_aligned), 
            min_percentage_aligned_reads = min(Percentage_reads_aligned), 
            max_percentage_aligned_reads = max(Percentage_reads_aligned))
# Enclosure mean_percentage_aligned_reads sd_percentage_aligned_reads min_percentage_aligned_reads max_percentage_aligned_reads
# H21                               0.517                       0.156                        0.293                        0.982
# P1                                0.480                       0.173                        0.231                        1.13 

####Plot - Naive vs Established####
MEGARes_aligned_read_percentages_P1vsH21<- ggplot(Percentage_reads_aligned_MEGARes, 
                                                     aes(x = Enclosure, 
                                                         y= Percentage_reads_aligned, 
                                                         color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Percentage (%) Aligned Reads", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 15, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) +
  geom_pwc (method = "wilcox_test",
            label = "Wilcoxon, p = {p}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
MEGARes_aligned_read_percentages_P1vsH21

#####SUPPLEMENTARY FIGURE 1######
sfigure1 <- cowplot::plot_grid(sequencing_depth_P1vsH21+
                                 theme(axis.title.y = element_text(size = 16)),  
                               kraken2_classified_read_percentages_P1vsH21+
                                 theme(axis.title.y = element_text(size = 16)),
                               bacteria_archaea_samplesums_P1vsH21+
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.title.y = element_text(size = 16),
                                       legend.position = "none"),
                               MEGARes_aligned_read_percentages_P1vsH21,
                               rel_heights = c(0.55, 0.55, 0.55, 0.6),
                               align = "v", 
                               ncol = 1, 
                               labels = "AUTO", label_size = 22)
sfigure1
ggsave("SupplementaryFigure1.png", 
       sfigure1, 
       device = "png", 
       dpi = 500, 
       width = 8, height =16)

####Plot - Naive and Established over time####
MEGARes_aligned_read_percentages_P1andH21_overtime<- ggplot(Percentage_reads_aligned_MEGARes, 
                                    aes(x = factor(Date_num), 
                                        y= Percentage_reads_aligned, 
                                        color = Enclosure)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Percentage (%) Aligned Reads", 
       x = "Day") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.01,0,0.1,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_color_manual(values=enclosure.palette)+
  theme(
        plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "none",
        # legend.text = element_text(size = 20),
        # legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title = element_text(size = 28, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text= element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5)) 
MEGARes_aligned_read_percentages_P1andH21_overtime

##Read counts over time
read_counts_overtime <- cowplot::plot_grid(sequencing_depth_P1andH21_overtime +
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.title.y = element_text(size = 16)),
                               kraken2_classified_read_percentages_P1andH21_overtime +
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       strip.text = element_blank(),
                                       axis.title.y = element_text(size = 14)),
                               bacteria_archaea_samplesums_P1andH21_overtime+
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       strip.text = element_blank(),
                                       axis.title.y = element_text(size = 16)),
                               MEGARes_aligned_read_percentages_P1andH21_overtime+
                                 theme(strip.text = element_blank(),
                                       axis.title.y = element_text(size = 14)),
                               align = "v", 
                               ncol = 1, 
                               labels = "AUTO", label_size = 22, 
                               rel_heights = c(0.35, 0.3, 0.3, 0.35)) # Adjust heights as needed)
read_counts_overtime

##Group level#### 
phyloseq_ARG.samples.group <- tax_glom(phyloseq_ARG.samples, taxrank = "Group", NArm = F)
phyloseq_ARG.samples.group #1373 groups and 216 samples

##TAX GLOMMING - RAW COUNTS##### 
phyloseq_ARG.samples.type <- tax_glom(phyloseq_ARG.samples.group, taxrank = "Type", NArm = F) 
phyloseq_ARG.samples.type #4 types (216 samples)

phyloseq_ARG.samples.class <- tax_glom(phyloseq_ARG.samples.group, taxrank = "Class", NArm = F) # classes
phyloseq_ARG.samples.class #57 classes (216 samples)

phyloseq_ARG.samples.mechanism <- tax_glom(phyloseq_ARG.samples.group, taxrank = "Mechanism", NArm = F) 
phyloseq_ARG.samples.mechanism #214 mechanisms (216 samples)

##RESISTOME ALPHA DIVERSITY ######
alpha_div_ARG1 <- phyloseq::estimate_richness(phyloseq_ARG.samples.group, 
                                          measures = c("Observed", "Shannon")) # richness, diversity
alpha_div_ARG2 <- microbiome::evenness(phyloseq_ARG.samples.group, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness
# combine metrics with metadata
alpha_div_ARG <- cbind(alpha_div_ARG1, alpha_div_ARG2)%>%
  rownames_to_column(var = "SampleID")
alpha_div_ARG
metadata_alpha_div_ARG_df <- data.frame(phyloseq_ARG.samples.group@sam_data)
alpha_div_ARG_meta <- left_join(metadata_alpha_div_ARG_df, 
                            alpha_div_ARG, 
                            by = "SampleID") 
alpha_div_ARG_meta # metadata and div metrics
nrow(alpha_div_ARG_meta) #216 samples 


#Pivot to long format 
alpha_div_ARG_meta_long <- 
  alpha_div_ARG_meta %>%
  pivot_longer(cols = c(Observed, Shannon, pielou),  
               names_to = "alpha_div_ARG_metric", 
               values_to = "alpha_div_ARG_value") 
alpha_div_ARG_meta_long

##Factoring alpha div metrics
alpha_div_ARG_meta_long$alpha_div_ARG_metric <- factor(alpha_div_ARG_meta_long$alpha_div_ARG_metric,
                                               levels = c("Observed","pielou", "Shannon"))

####Alpha diversity indexes and Water quality over time#####
#Editing alpha diversity and metadata dataframe
alpha_div_ARG_wq_time_long <- alpha_div_ARG_meta %>%
  mutate(Copper_keep = Copper_level_mg_L) %>%   #keep a column as copper to be able to color the plot 
  pivot_longer(cols = c("Copper_level_mg_L",
                        "Temperature_F",
                        "Chlorine_mg_L", 
                        "pH_spu", 
                        "Ammonia_mg_L",
                        "Nitrite_mg_L",
                        "Nitrate_UV_mg_L", 
                        "Salinity_ppt",
                        "Alkalinity_mg_L", 
                        "Shannon",
                        "Observed",
                        "pielou"),
               names_to = "Index",
               values_to = "Index_value")%>%
  mutate(Index = factor(Index, levels = c(
    "Observed",
    "Shannon",
    "pielou",
    "Copper_level_mg_L",
    "Ammonia_mg_L",
    "Nitrite_mg_L",
    "Nitrate_UV_mg_L",
    "pH_spu", 
    "Salinity_ppt", 
    "Temperature_F",
    "Chlorine_mg_L", 
    "Alkalinity_mg_L"
  )))

#####Adding vertical dashed line for when "phase" changes#######
#Vertical lines to be placed
line_breaks_phases <- data.frame(
  Enclosure = c("H21", "H21", "H21", "H21", "H21", "H21","H21",
                "H21", "H21", "H21", "H21", "H21", "H21","H21",
                "P1", "P1", "P1", "P1", "P1", "P1", "P1", "P1"), 
  Date_num = c("1", "27", "29", "38", "39", "51", "52", "81", "86", "108", "109", "135", "136", "146",
               "1", "53", "54", "65", "66", "104", "105", "169"))
line_breaks_phases

#####Plot - Just  Shannon, Copper and Ammonia levels - Date number since start of sampling as as.factor####
alpha_div_ARG_wq_date_num_factor <- ggplot(alpha_div_ARG_wq_time_long%>%
                                         filter(Index %in% c("Copper_level_mg_L",
                                                             "Shannon",
                                                             "Ammonia_mg_L")),
                                       aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  scale_x_discrete(
    drop = TRUE,
    expand = expansion(mult = c(0.03, 0.03)),
    breaks = function(x) {
      x_num <- sort(unique(as.numeric(x)))
      
      # targets up to 120 only
      targets <- c(1, seq(30, 120, by = 30))
      
      closest <- unique(sapply(targets, function(t) {
        x_num[which.min(abs(x_num - t))]
      }))
      
      # add max separately
      final_vals <- unique(c(closest, max(x_num)))
      
      as.character(final_vals)
    },
    labels = function(x) {
      x
    }
  )+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "RESISTOME\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)")))+
  theme(
    legend.position = c(0.85, 1.7),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.y  = element_text(colour = "white", size = 30, 
                                 face = "bold", angle = 0),
    strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_ARG_wq_date_num_factor

#####Plot - Other Water quality measures- Date number since start of sampling as as.factor#####
alpha_div_ARG_wq_date_num_factor_other_metadata <- ggplot(alpha_div_ARG_wq_time_long%>%
                                                        filter(Index %in% c("Copper_level_mg_L",
                                                                            "Shannon",
                                                                            "Ammonia_mg_L",
                                                                            #"Chlorine_mg_L", 
                                                                            #"Alkalinity_mg_L",
                                                                            "Temperature_F",
                                                                            "pH_spu",
                                                                            "Salinity_ppt"
                                                        )),
                                                      aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black", 
  #            alpha = 0.8) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  scale_x_discrete(
    drop = TRUE,
    expand = expansion(mult = c(0.03, 0.03)),
    breaks = function(x) {
      x_num <- sort(unique(as.numeric(x)))
      
      # targets up to 120 only
      targets <- c(1, seq(30, 120, by = 30))
      
      closest <- unique(sapply(targets, function(t) {
        x_num[which.min(abs(x_num - t))]
      }))
      
      # add max separately
      final_vals <- unique(c(closest, max(x_num)))
      
      as.character(final_vals)
    },
    labels = function(x) {
      x
    }
  )+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title = "RESISTOME\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)",
                                      "Temperature_F"= "Temperature (F)",
                                      "Salinity_ppt" = "Salinity (ppt)", 
                                      "pH_spu" = "pH (spu)"
                                      #"Chlorine_mg_L" = "Chlorine (mg/L)",
                                      #"Alkalinity_mg_L" = "Alkalinity (mg/L)"
             )))+
  theme(
    legend.position = c(0.95, 1.3),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.y  = element_text(colour = "white", size = 30, 
                                 face = "bold", angle = 0),
    strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(colour = "black", size = 18),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_ARG_wq_date_num_factor_other_metadata


##RELATIVE ABUNDANCE - RESISTOME####
any(sample_sums(phyloseq_ARG.samples.group)== 0) ## no samples with 0 OTUs

###TAX GLOMMING - RA NORMALIZED COUNTS##### 
#Group level
phyloseq_ARG.samples.ra.group <- transform_sample_counts(phyloseq_ARG.samples.group, function(x) x/sum(x)*100) 
phyloseq_ARG.samples.ra.group #1372 groups (216 samples)

#Type
phyloseq_ARG.samples.ra.type <- tax_glom(phyloseq_ARG.samples.ra.group, taxrank = "Type", NArm = F) 
phyloseq_ARG.samples.ra.type #4 types (216 samples)

#Class
phyloseq_ARG.samples.ra.class <- tax_glom(phyloseq_ARG.samples.ra.group, taxrank = "Class", NArm = F) # classes
phyloseq_ARG.samples.ra.class #57 classes (216 samples)

#Mechanism
phyloseq_ARG.samples.ra.mechanism <- tax_glom(phyloseq_ARG.samples.ra.group, taxrank = "Mechanism", NArm = F) 
phyloseq_ARG.samples.ra.mechanism #214 mechanisms (216 samples)


###RA PLOT GROUP OVERALL RESISTOME ####
phyloseq_ARG.samples.ra.group #1372 taxa and 216 samples

# #Merge low abundance groups into one
phyloseq_ARG.samples.ra.group_filt <- merge_low_abundance_ARG_ra(phyloseq_ARG.samples.ra.group,
                                                                              "Enclosure",
                                                                              level = "Group",
                                                                              threshold = 0.1)
phyloseq_ARG.samples.ra.group_filt #321 gene groups with mean RA above 0.1%

#Melt to plot 
phyloseq_ARG.samples.ra.group_filt.melt <- psmelt(phyloseq_ARG.samples.ra.group_filt)
#For some reason the function is leaving the Others group as NA, fixing that...
phyloseq_ARG.samples.ra.group_filt.melt <-
  phyloseq_ARG.samples.ra.group_filt.melt %>%
  mutate(Group = ifelse(is.na(Group), "Others <0.1% RA", Group))

#Factoring so "Others" group is last
phyloseq_ARG.samples.ra.group_filt.melt <- phyloseq_ARG.samples.ra.group_filt.melt%>%
  mutate(Group = factor(Group,
                        levels = c(setdiff(Group,
                                           unique(grep("Others", Group, value = TRUE))),
                                   unique(grep("Others", Group, value = TRUE)))))##Factoring the Group column so that "Others.." is the last category

##Create color palette
palette_gene_groups <- distinctColorPalette(length(unique(phyloseq_ARG.samples.ra.group_filt.melt$Group)))
copper_filt_names <- unique(phyloseq_ARG.samples.ra.group_filt.melt$Group)# Create a named vector for the palette, where the names correspond to phlyum names
palette_gene_groups <- setNames((palette_gene_groups)[1:length(palette_gene_groups)], copper_filt_names)
#order.filt.palette <- unname(alphabet2())
palette_gene_groups$'Others <0.1% RA' <- "grey95"
#Plot
RA_enclosures_ARG_genegroup.plot <- ggplot(phyloseq_ARG.samples.ra.group_filt.melt,
                                                  aes(x=factor(Date_num), y= Abundance, fill = Group)) +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_fill_manual(values = palette_gene_groups)+
  # #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  guides(fill=guide_legend(title.position="top", ncol = 2))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.08, 0.5),  # x, y inside plot
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 22),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_enclosures_ARG_genegroup.plot 

###RA PLOT CLASS LEVEL COPPER ####
phyloseq_ARG.samples.ra.class #58 classes and 216 samples

#Out of this overall communities object, select only copper 
phyloseq_ARG.samples.ra.class.copper <- subset_taxa(phyloseq_ARG.samples.ra.class, 
                                                          Class == "Copper resistance") 

phyloseq_ARG.samples.ra.class.copper # 1 class and 216 samples 

#Melt to plot 
phyloseq_ARG.samples.ra.class.copper.melt <- psmelt(phyloseq_ARG.samples.ra.class.copper)

#Plot  
RA_enclosures_ARG_copper.plot <- ggplot(phyloseq_ARG.samples.ra.class.copper.melt,
                                        aes(x=factor(Date_num), y= Abundance, fill = Class)) +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  # #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme_bw()+
  theme(legend.text = element_text(size = 20),
        legend.position = "right",
        #legend.position = "none",
        legend.title = element_text(size = 24, face = "bold"),
        legend.key.size = unit(0.7, "cm"),
        strip.background = element_rect(fill = "black"),
        #strip.text.x  = element_text(colour = "white", size = 45, face = "bold"),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        #axis.title = element_blank(),
        axis.title.x = element_text(colour = "black", size = 22),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_text(colour = "black", size = 20,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 18),
        axis.ticks.x = element_line(colour = "black", linewidth = 1),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5))
RA_enclosures_ARG_copper.plot

###RA PLOT GROUP LEVEL COPPER ####
phyloseq_ARG.samples.ra.group #1372 taxa and 216 samples

#Out of this overall communities object, select only copper 
phyloseq_ARG.samples.ra.group.copper <- subset_taxa(phyloseq_ARG.samples.ra.group, 
                                                          Class == "Copper resistance") 

phyloseq_ARG.samples.ra.group.copper # 58 gene groups and 216 samples 

# #Merge low abundance groups into one
phyloseq_ARG.samples.ra.group.copper_filt <- merge_low_abundance_ARG_ra(phyloseq_ARG.samples.ra.group.copper,
                                                                                  "Enclosure",
                                                                                  level = "Group",
                                                                                  threshold = 0.1)
phyloseq_ARG.samples.ra.group.copper_filt #19 gene groups with mean RA above 0.1%

#Melt to plot 
phyloseq_ARG.samples.ra.group.copper_filt.melt <- psmelt(phyloseq_ARG.samples.ra.group.copper_filt)
#For some reason the function is leaving the Others group as NA, fixing that...
phyloseq_ARG.samples.ra.group.copper_filt.melt <-
  phyloseq_ARG.samples.ra.group.copper_filt.melt %>%
  mutate(Group = ifelse(is.na(Group), "Others <0.1% RA", Group))

#Factoring so "Others" group is last
phyloseq_ARG.samples.ra.group.copper_filt.melt <- phyloseq_ARG.samples.ra.group.copper_filt.melt%>%
  mutate(Group = factor(Group,
                        levels = c(setdiff(Group,
                                           unique(grep("Others", Group, value = TRUE))),
                                   unique(grep("Others", Group, value = TRUE)))))##Factoring the Group column so that "Others.." is the last category


####SUPPLEMENTARY TABLE 6 #####
stable6 <-  phyloseq_ARG.samples.ra.group.copper_filt.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Group, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  rename(Phase = Date_num_phase_abbrv, 
         `Mean Relative Abundance (%)` = mean_abun,  
         `Standard Deviation (%)` = sd_abun, 
         `Min Relative Abundance (%)` = min_abun, 
         `Max Relative Abundance (%)` = max_abun)
stable6
write_xlsx(stable6, 
           "SupplementaryTable6.xlsx")

##Create color palette
palette_copper_gene_groups <- distinctColorPalette(length(unique(phyloseq_ARG.samples.ra.group.copper_filt.melt$Group)))
copper_filt_names <- unique(phyloseq_ARG.samples.ra.group.copper_filt.melt$Group)# Create a named vector for the palette, where the names correspond to phlyum names
palette_copper_gene_groups <- setNames((palette_copper_gene_groups)[1:length(palette_copper_gene_groups)], copper_filt_names)
#order.filt.palette <- unname(alphabet2())
palette_copper_gene_groups$'Others <0.1% RA' <- "grey95"

# #Just top most abundant groups to include in legend
# top_taxa_copper_groups <- top_taxa_legend(phyloseq_ARG.samples.ra.group.copper.melt, taxlevel = "Group", n = 23)

#Copper gene groups palette
# palette_copper_gene_groups <- c(
#   "#1B9E77", "#D95F02", "#7570B3",  "#66A61E","#FDBF6F",
#   "#E6AB02", "#666666", "#1F78B4", "#A6761D", "#B2DF8A",
#   "#FB9A99", "#E7298A", "#6A3D9A", "#CAB2D6", "#FFFF99"
# )
#Plot
RA_enclosures_ARG_copper_genegroup.plot <- ggplot(phyloseq_ARG.samples.ra.group.copper_filt.melt,
                                                    aes(x=factor(Date_num), y= Abundance, fill = Group)) +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_fill_manual(values = palette_copper_gene_groups,
                    labels = function(x) str_wrap(x, width = 15)) +
  # #Scale x, want to keep 1 and the closest 30-multiple, plus max date
  # scale_x_discrete(
  #   drop = TRUE,
  #   expand = expansion(mult = c(0.03, 0.03)),
  #   breaks = function(x) {
  #     x_num <- sort(unique(as.numeric(x)))
  #     
  #     # targets up to 120 only
  #     targets <- c(1, seq(30, 120, by = 30))
  #     
  #     closest <- unique(sapply(targets, function(t) {
  #       x_num[which.min(abs(x_num - t))]
  #     }))
  #     
  #     # add max separately
  #     final_vals <- unique(c(closest, max(x_num)))
  #     
  #     as.character(final_vals)
  #   },
  #   labels = function(x) {
  #     x
  #   }
  # )+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135", "146"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","169"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  guides(fill=guide_legend(title.position="top", ncol = 2))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.08, 0.5),  # x, y inside plot
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    strip.background = element_rect(fill = "black"),
    panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    strip.text.x  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title.x = element_text(colour = "black", size = 22),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_enclosures_ARG_copper_genegroup.plot 

######PLOT - Together with alpha div of nitrifying community, as well as nitrifying community RA at the family level#######
#Have to edit figures that will go into the final plot
alpha_div_nit_wq_date_num_factor_other_metadata_3 <- alpha_div_nit_wq_date_num_factor_other_metadata +
  theme(
    legend.position = c(0.7, 1.6),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 20))
RA_family_enclosures_nit_plot_datenum_3 <- RA_family_enclosures_nit_plot_datenum + 
  theme(
    legend.position = c(1.09, 0.5), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 20)) 
RA_enclosures_nitrifiers_only_species.plot_3 <- RA_enclosures_nitrifiers_only_species.plot +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 20), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18, face = "bold")
    )
RA_enclosures_ARG_copper_genegroup.plot

#Final plot
figure_alpha_div_nit_species_ra_time_copperARG <-
  alpha_div_nit_wq_date_num_factor_other_metadata_3 /
  RA_family_enclosures_nit_plot_datenum_3  /
  RA_enclosures_nitrifiers_only_species.plot_3 +
  RA_enclosures_ARG_copper_genegroup.plot +
  plot_layout(heights = c(1.1, 0.4, 0.4, 0.4))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_div_nit_species_ra_time_copperARG

#Saving figure
ggsave("figure_alpha_div_nit_species_ra_time_copperARG.png", 
       figure_alpha_div_nit_species_ra_time_copperARG, 
       device = "png", 
       dpi = 600, 
       height = 20, 
       width = 26)

######PLOT - Together with alpha div of nitrifying community#######
#Have to edit figures that will go into the final plot
alpha_div_nit_wq_date_num_factor_other_metadata_4 <- alpha_div_nit_wq_date_num_factor_other_metadata +
  theme(
    legend.position = c(0.7, 1.6),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 18)
    )
RA_family_enclosures_nit_plot_datenum_4 <- RA_family_enclosures_nit_plot_datenum + 
  theme(
    legend.position = c(1.09, 0.5), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 17),
    axis.text.y = element_text(size = 17)
    ) 

RA_enclosures_ARG_copper_genegroup.plot_4 <- RA_enclosures_ARG_copper_genegroup.plot +
  guides(fill=guide_legend(title.position="top", ncol = 2))+
  theme(
    legend.position = c(1.08, 0.5), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.key.size =  unit(0.7, "cm"),
    axis.title.y = element_text(size = 17)
    )

#Final plot
figure_alpha_div_nit_species_ra_time_copper_ARG_nospec <-
  alpha_div_nit_wq_date_num_factor_other_metadata_4 /
  RA_family_enclosures_nit_plot_datenum_4  /
  RA_enclosures_ARG_copper_genegroup.plot_4 +
  plot_layout(ncol =1, heights = c(1.4, 0.6, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_div_nit_species_ra_time_copper_ARG_nospec

#Saving figure
ggsave("figure_alpha_div_nit_species_ra_time_copper_ARG_nospec.png", 
       figure_alpha_div_nit_species_ra_time_copper_ARG_nospec, 
       device = "png", 
       dpi = 600, 
       height = 15, 
       width = 26)

######PLOT - Together with alpha div of Overall community#######
#Have to edit figures that will go into the final plot
alpha_div_wq_date_num_factor_other_metadata_6 <- alpha_div_wq_date_num_factor_other_metadata +
  theme(
    legend.position = c(0.93, 1.5),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15),
    strip.text.y = element_text(size = 25) 
    )
RA_family_enclosures_overall_plot_datenum_6 <- RA_family_enclosures_overall_plot_datenum  +
  theme(legend.position = c(1.06, 0.5), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 17),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 16, face = "bold"))
RA_enclosures_ARG_copper_genegroup.plot_6 <- RA_enclosures_ARG_copper_genegroup.plot +
  theme(
    legend.position = c(1.08, 0.5), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 17),
    axis.text.y = element_text(size = 17))
  
#Final plot
figure_alpha_overall_div_family_time_copper_ARG <- 
  alpha_div_wq_date_num_factor_other_metadata_6 /
  RA_family_enclosures_overall_plot_datenum_6 /
  RA_enclosures_ARG_copper_genegroup.plot_6 +
  plot_layout(ncol = 1, heights = c(1.3, 0.8, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_overall_div_family_time_copper_ARG
ggsave("figure_alpha_overall_div_family_time_copper_ARG.png", 
       figure_alpha_overall_div_family_time_copper_ARG, 
       device = "png", 
       dpi = 600, 
       height = 15, 
       width = 26)


####JOIN BETA DIV WITH RESISTOME AND OTHER PLOTS#######
alpha_div_wq_date_num_factor_other_metadata_5 <- alpha_div_wq_date_num_factor_other_metadata +
  theme(
    legend.position = c(0.93, 1.3),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 17),
    strip.text.y = element_text(size = 30) 
  )
RA_family_enclosures_overall_plot_datenum_5 <- RA_family_enclosures_overall_plot_datenum  +
  theme(legend.position = c(1.08, 0.45), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"))

first_diff_BC_H21_P1_plot_1 <- first_diff_BC_H21_P1_plot +
  theme(axis.text.x = element_blank(),
        legend.position = c(1.07, 0.5),  # x, y inside plot
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))

RA_enclosures_ARG_copper_genegroup.plot_5 <- RA_enclosures_ARG_copper_genegroup.plot +
  theme(
    legend.position = c(1.08, 0.5), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20))

#Final plot
figure_alpha_div_overall_species_ra_time_BC <-
  alpha_div_wq_date_num_factor_other_metadata_5 /
  RA_family_enclosures_overall_plot_datenum_5  /
  first_diff_BC_H21_P1_plot_1 / 
  RA_enclosures_ARG_copper_genegroup.plot_5 +
  plot_layout(heights = c(0.9, 0.6, 0.6, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_div_overall_species_ra_time_BC
ggsave("figure_alpha_div_overall_species_ra_time_BC.png", 
       figure_alpha_div_overall_species_ra_time_BC, 
       device = "png", 
       dpi = 600, 
       height = 23, 
       width = 26)


####FIGURE 1 #####
#Metadata and alpha diversity measures of overall community
alpha_div_wq_date_num_factor_other_metadata_7 <- alpha_div_wq_date_num_factor_other_metadata +
  theme(
    plot.title = element_blank(),
    legend.position = c(0.98, 1.22),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15),
    strip.text.y = element_text(size = 30))

#RA plot of overall microbial community
RA_family_enclosures_overall_plot_datenum_7 <- RA_family_enclosures_overall_plot_datenum  +
  scale_fill_manual(values = family_named_palette,
                    breaks = top_families,
                    #labels = function(x) str_wrap(x, width = 20)
  ) +
  theme(legend.position = c(1.084, 0.47), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 21, face = "bold"))

#RA plot of ARGs
RA_enclosures_ARG_copper_genegroup.plot_7 <- RA_enclosures_ARG_copper_genegroup.plot +
  theme(
    legend.position = c(1.08, 0.5), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 21, face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

#RA plot of nitrifiers within the overall community at the family level
RA_family_enclosures_nit_plot_datenum_3 <- RA_family_enclosures_nit_plot_datenum + 
  theme(
    legend.position = c(1.075, 0.5),  
    legend.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 18),
    legend.key.size = unit(0.5, "cm"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20)) 

#RA plot of nitrifiers at the species level only within nitrifiers
RA_enclosures_nitrifiers_only_species.plot_3 <- RA_enclosures_nitrifiers_only_species.plot +
  theme(
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 21, face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20), 
    axis.title.x = element_text(size = 30))


#Final plot
figure1 <-
  alpha_div_wq_date_num_factor_other_metadata_7 /
  RA_family_enclosures_overall_plot_datenum_7 / 
  RA_enclosures_ARG_copper_genegroup.plot_7 /
  RA_family_enclosures_nit_plot_datenum_3 /
  RA_enclosures_nitrifiers_only_species.plot_3 +
  plot_layout(heights = c(0.7, 0.4, 0.4, 0.4, 0.4))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure1
ggsave("Figure1.png", 
       figure1, 
       device = "png", 
       dpi = 600, 
       height = 23, 
       width = 26)


#JOIN BETA DIV OF NITRIFIERS WITH RESISTOME AND OTHER PLOTS 
#Join with other plots 
alpha_div_nit_wq_date_num_factor_other_metadata_5 <- alpha_div_nit_wq_date_num_factor_other_metadata +
  theme(
    legend.position = c(0.93, 1.3),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 17),
    strip.text.y = element_text(size = 30) 
  )
RA_family_enclosures_nit_plot_datenum_5 <- RA_family_enclosures_nit_plot_datenum +
  theme(legend.position = c(1.08, 0.6), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))

RA_enclosures_nitrifiers_only_species.plot_5 <- RA_enclosures_nitrifiers_only_species.plot+
  theme(legend.position = c(1.07, 0.6), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.key.size = unit(0.6, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = "bold"))

first_diff_nit_BC_H21_P1_plot_1 <- first_diff_nit_BC_H21_P1_plot +
  theme(axis.text.x = element_blank(),
        legend.position = c(1.07, 0.5),  # x, y inside plot
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))

RA_enclosures_ARG_copper_genegroup.plot_5 <- RA_enclosures_ARG_copper_genegroup.plot +
  theme(
    legend.position = c(1.08, 0.5), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 20))

#Final plot
figure_alpha_div_nit_species_ARG_ra_time_BC <-
  alpha_div_nit_wq_date_num_factor_other_metadata_5 /
  RA_enclosures_nitrifiers_only_species.plot_5  /
  first_diff_nit_BC_H21_P1_plot_1 / 
  RA_enclosures_ARG_copper_genegroup.plot_5 +
  plot_layout(heights = c(0.9, 0.6, 0.6, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_div_nit_species_ARG_ra_time_BC
ggsave("figure_alpha_div_nit_species_ARG_ra_time_BC.png", 
       figure_alpha_div_nit_species_ARG_ra_time_BC, 
       device = "png", 
       dpi = 600, 
       height = 23, 
       width = 26)


#Final plot - with family level
figure_alpha_div_nit_family_ra_time_BC <-
  alpha_div_nit_wq_date_num_factor_other_metadata_5 /
  RA_family_enclosures_nit_plot_datenum_5 /
  RA_enclosures_nitrifiers_only_species.plot_5  /
  first_diff_nit_BC_H21_P1_plot_1 +
  plot_layout(heights = c(0.9, 0.6, 0.6, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_div_nit_family_ra_time_BC
ggsave("figure_alpha_div_nit_family_ra_time_BC.png", 
       figure_alpha_div_nit_family_ra_time_BC, 
       device = "png", 
       dpi = 600, 
       height = 23, 
       width = 26)

####Correlation of Rhodobacteraceae with RA of Copper resistance class RA#########
phyloseq.bacteria.samples_family.ra #6755 families 

#Out of this overall communities object, select only rhodobacteraceae
phyloseq.bacteria.samples_family.ra.rhodobacteraceae <- subset_taxa(phyloseq.bacteria.samples_family.ra, 
                                                                    Family == "Rhodobacteraceae") 
phyloseq.bacteria.samples_family.ra.rhodobacteraceae <- subset_samples(phyloseq.bacteria.samples_family.ra.rhodobacteraceae, 
                                                                       sample_sums(phyloseq.bacteria.samples_family.ra.rhodobacteraceae) > 0)
phyloseq.bacteria.samples_family.ra.rhodobacteraceae #Rhodobacteraceae (1 taxa) in 216 samples 

#Melt to long format
phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt <- psmelt(phyloseq.bacteria.samples_family.ra.rhodobacteraceae)

#Join with copper resistance class RA
#Copper resistance class melted object (long format)
phyloseq_ARG.samples.ra.class.copper.melt_join <- phyloseq_ARG.samples.ra.class.copper.melt %>%
  rename(Abundance_copper = Abundance)
phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt <- phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt %>%
  rename(Abundance_rhodobacteraceae = Abundance)

#Join 
rhodobacteraceae_copper.ra.melt <- merge(phyloseq.bacteria.samples_family.ra.rhodobacteraceae.melt%>%
                                           select(SampleID, Abundance_rhodobacteraceae, Enclosure), 
                                         phyloseq_ARG.samples.ra.class.copper.melt_join%>%
                                           select(SampleID, Abundance_copper, Enclosure), 
                                         by = c("SampleID", "Enclosure"))


#Plot
copperARG_rhodobacteraceae_correlation_plot<- ggplot(rhodobacteraceae_copper.ra.melt,
                                                       aes(y = Abundance_copper, 
                                                           x = Abundance_rhodobacteraceae,
                                                          fill = Enclosure, color = Enclosure)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Enclosure)) +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  scale_color_manual(values = enclosure.palette)+
  scale_fill_manual(values = enclosure.palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(x = expression(italic(Rhodobacteraceae) ~ "RA (%)"), 
       y = "Copper ARG Class RA (%)") +
  theme_bw() +
  geom_smooth(method="lm",
              se=TRUE,
              #color = "black", 
              linewidth = 0.6,
              alpha = 0.3) +
  stat_cor(method = "spearman",
           label.x.npc = "center",
           label.y.npc = "bottom",
           color = "black",
           size = 8) +
  theme(legend.position = "none",
        legend.text = element_text(size = 28, vjust = 0.5),
        legend.title = element_text(size = 28, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 38, face = "bold"),
        axis.title = element_text(colour = "black", size = 35),
        axis.text.x = element_text(colour = "black", size = 30),
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))
copperARG_rhodobacteraceae_correlation_plot
ggsave("copperARG_rhodobacteraceae_correlation_plot.png", 
       copperARG_rhodobacteraceae_correlation_plot, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 16)
