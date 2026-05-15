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
library(phyloseq); library (tidyverse); library(ggplot2);  library(stringr); 
library(dplyr);library(metagMisc); library(metagenomeSeq); library(vegan); library(cowplot);
library(ggdendro); library(pairwiseAdonis); library(randomcoloR); library(ggpubr); library(ppcor)
library(ggsignif); library (ANCOMBC);library(maaslin3); library (UpSetR); library(MicrobiotaProcess); library(microbiome)
library(ggtext); library(ggnewscale); library(rstatix); library(ggrepel); library(ggh4x); library(svglite);
library(lmerTest); library(mgcv); library(rmcorr); library("emmeans"); library(patchwork); library(colorspace)
library(MiRKAT);library(writexl)
library(pals)


##Source functions
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbundanceOthersPercentage.R')
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/top_taxa_legend_updated.R")
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/fill_taxonomy_updated.R')
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbun_group_microbiome.R')


#Importing data from kraken output nt_core - counts will be classified reads#### 
counts <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/Paired_end_mode_updated_20260513/Conf_01/kraken_analytic_matrix.conf_0.1.csv')
##Separating into taxonomy levels
counts_separated_tax <- counts %>%
  separate(taxa, 
           into = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|",#splits the strings by the "|" symbol.
           fill = "right") #fill = "right", missing components are added as "NA" to the right (last columns) instead of to the left 

##Extracting just taxonomy  (columns 1:8 are taxonomy, the rest are counts)
tax.table<- counts_separated_tax %>%
  dplyr::select(1:8)
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
setdiff(sample_names(OTU), metadata$SampleID) #Yes, "H21_1021a" and "H21_1021b". These will not be included in the phyloseq object.

#Are there samples in metadata that don't have sequencing data?
setdiff(metadata$SampleID, sample_names(OTU)) #Yes, "P1_1126", "P1_1203", 
#"P1_1216", "P1_1225", "P1_1228", "P1_0104", "P1_0112", "P1_0115", 
#"P1_0205", "P1_0212", "P1_0218", "H21_1021", "H21_1122b"

#COLOR PALETTES#####
enclosure.palette <- c("H21" = "#fc8d62",  
                       "P1"  = "#8da0cb" )

reads.palette <- c("Raw" ="#E69F00", 
                   "Trimmed" = "#0072B2")

#PREPROCESSING ####
phyloseq #69,067 taxa and 240 samples 
      
##Selecting only Bacteria/Archaea####
phyloseq.bacteria <- subset_taxa(phyloseq, Domain=="Archaea" | Domain=="Bacteria")
phyloseq.bacteria #23725 taxa and 240 samples

##Selecting only viruses######
phyloseq.viruses <- subset_taxa(phyloseq, Domain=="Viruses")
taxanames_viruses <- c("Kingdom", "Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species") ##they have a different classification system, updating it here
colnames(phyloseq.viruses@tax_table) <- taxanames_viruses #replacing col names of the tax_table for new ones
colnames(phyloseq.viruses@tax_table) #OK taxonomy ranks
phyloseq.viruses #14570 taxa and 240 samples

##Selecting only eukaryota #####
phyloseq.eukaryota <- subset_taxa(phyloseq, Domain=="Eukaryota")
colnames(phyloseq.eukaryota@tax_table) ##These are OK taxonomy ranks
phyloseq.eukaryota #30,772 taxa and 240 samples

#WORKING ON BACTERIA/ARCHAEA ONLY####
# some QC checks of the "classified" reads per samples
min(sample_sums(phyloseq.bacteria)) # 0 (P1_0308)
max(sample_sums(phyloseq.bacteria)) # 69,452,267  (H21_0119) 
mean(sample_sums(phyloseq.bacteria)) #9,058,722
median(sample_sums(phyloseq.bacteria)) # 7,137,108
sort(sample_sums(phyloseq.bacteria))

##Zymo and controls#####
### Getting samples from ZYMOs and EB, NTC
phyloseq.bacteria.controls <- subset_samples(phyloseq.bacteria, 
  grepl("NTC|EB|Zymo", sample_names(phyloseq.bacteria)))
phyloseq.bacteria.controls <- prune_taxa(taxa_sums(phyloseq.bacteria.controls) > 0, phyloseq.bacteria.controls) 
phyloseq.bacteria.controls #1315 taxa, 15 samples(NTC, EB and Zymos)

##Samples#####
##New phyloseq of just samples
phyloseq.bacteria.samples <- subset_samples(phyloseq.bacteria, 
                                             !grepl("NTC|EB|Zymo", sample_names(phyloseq.bacteria)))
phyloseq.bacteria.samples #23725 taxa and 225 samples
sort(sample_sums(phyloseq.bacteria.samples))
#Taking out those with low counts
phyloseq.bacteria.samples <- prune_samples(sample_sums(phyloseq.bacteria.samples) > 100000, phyloseq.bacteria.samples) 
phyloseq.bacteria.samples <- prune_taxa(taxa_sums(phyloseq.bacteria.samples) > 0, phyloseq.bacteria.samples) 
phyloseq.bacteria.samples #23595 taxa and 222 samples  (dropped P1_0308, H21_0109, and H21_0120)
sort(sample_sums(phyloseq.bacteria.samples)) #OK


###Dropping samples before copper dosage started####
#What's the range of dates 
range(phyloseq.bacteria.samples@sam_data$Collection_Date)#"2023-04-20" "2024-04-30"

#Actual Copper dosing starts from 10/09/2023 (sampling ends on 03/02/2024) for naive system
#Actual copper dosing starts from 11/14/2023 (sampling ends on 04/30/2024) for established system
phyloseq.bacteria.samples.dates <- subset_samples(phyloseq.bacteria.samples, Collection_Date > "2023-10-05")
phyloseq.bacteria.samples.dates
phyloseq.bacteria.samples.dates <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates) > 0, phyloseq.bacteria.samples.dates) 
phyloseq.bacteria.samples.dates #23,476 taxa and 217 samples
setdiff(sample_names(phyloseq.bacteria.samples), sample_names(phyloseq.bacteria.samples.dates)) #Dropped "H21_0912" "H21_1005" "P1_0420"  "P1_0427"  "P1_0504" 

#Also, have H21_1202a and 1202b. These were taken on Dec 2nd, 2023. A is before backwash, B is afterbackwash. 
#Have more reliable metadata for H21_1202a (before backwash)
phyloseq.bacteria.samples.dates <- subset_samples(phyloseq.bacteria.samples.dates, 
                                                  SampleID != "H21_1202b")
phyloseq.bacteria.samples.dates #23,476 taxa and 216 samples (Dopped H21_1202b)
phyloseq.bacteria.samples.dates <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates) > 0, phyloseq.bacteria.samples.dates) 
phyloseq.bacteria.samples.dates #23,461 taxa and 216 samples
setdiff(sample_names(phyloseq.bacteria.samples), sample_names(phyloseq.bacteria.samples.dates)) 

###H21####
phyloseq.bacteria.samples.dates_H21 <- subset_samples(phyloseq.bacteria.samples.dates, Enclosure == "H21")
phyloseq.bacteria.samples.dates_H21 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates_H21) > 0, 
                                                  phyloseq.bacteria.samples.dates_H21)
phyloseq.bacteria.samples.dates_H21 #21,812 taxa and 92 samples
range(phyloseq.bacteria.samples.dates_H21@sam_data$Collection_Date)#OK, now "2023-10-09" through "2024-03-02"

###P1####
phyloseq.bacteria.samples.dates_P1 <- subset_samples(phyloseq.bacteria.samples.dates, Enclosure == "P1")
phyloseq.bacteria.samples.dates_P1 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates_P1) > 0, 
                                                 phyloseq.bacteria.samples.dates_P1)
phyloseq.bacteria.samples.dates_P1 #19724 taxa and 124 samples
range(phyloseq.bacteria.samples.dates_P1@sam_data$Collection_Date)#OK, now "2023-11-14" through "2024-04-30"


##NITRIFYING TAXA####
nitrifiers_all <- subset_taxa(phyloseq.bacteria.samples.dates, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
nitrifiers_all #449 taxa and 217 samples
nitrifiers <- subset_samples(nitrifiers_all, sample_sums(nitrifiers_all) > 0)
nitrifiers #449 taxa and 217 samples 

##QC checks again
min(sample_sums(phyloseq.bacteria.samples.dates)) #2,434,204 (H21_1101)
max(sample_sums(phyloseq.bacteria.samples.dates)) #35,431,988 (H21_0119) 
mean(sample_sums(phyloseq.bacteria.samples.dates)) #9,095,353
median(sample_sums(phyloseq.bacteria.samples.dates)) #7,622,284
sort(sample_sums(phyloseq.bacteria.samples.dates)) 

#COMPARING SEQUENCING DEPTHS#######
cuso4_raw_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Raw/CuSO4_clean_raw_read_counts.csv')
str(cuso4_raw_read_counts)

#Keeping just those samples that I'm analyzing in phyloseq.bacteria.samples.dates
sampleIDs_phyloseq.bacteria.samples.dates <- phyloseq.bacteria.samples.dates@sam_data$SampleID 
length(sampleIDs_phyloseq.bacteria.samples.dates)#216 samples

#Filtering just phyloseq.bacteria.samples.dates SampleIDs, adding metadata
cuso4_raw_read_counts_samples_metadata <- cuso4_raw_read_counts %>%
  filter(SampleID %in% sampleIDs_phyloseq.bacteria.samples.dates)%>%
  dplyr::left_join(metadata, by = "SampleID")
nrow(cuso4_raw_read_counts_samples_metadata) #Ok, 216 samples

###Established vs Naive####
sequencing_depth_P1vsH21<- ggplot(cuso4_raw_read_counts_samples_metadata, 
                                  aes(x = Enclosure, y= Num_Reads_Forward_Raw, 
                                      color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Paired-end Reads", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
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

#COMPARING TRIMMED READS#######
cuso4_trimmed_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Trimmed/CuSO4_clean_trimmed_read_counts.csv')
str(cuso4_trimmed_read_counts)

#Filtering just phyloseq.bacteria.samples.dates SampleIDs, adding metadata
cuso4_trimmed_read_counts_samples_metadata <- cuso4_trimmed_read_counts %>%
  filter(SampleID %in% sampleIDs_phyloseq.bacteria.samples.dates)%>%
  dplyr::left_join(cuso4_raw_read_counts_samples_metadata, by = "SampleID")
nrow(cuso4_trimmed_read_counts_samples_metadata) #Ok, 217 samples

###Established vs Naive####
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

###Established and Naive, Raw and Trimmed over time####
cuso4_trimmed_read_counts_samples_metadata_long <- 
  cuso4_trimmed_read_counts_samples_metadata %>%
  pivot_longer(cols = c(Num_Reads_Forward_Trimmed_Paired,
             Num_Reads_Forward_Raw),
    names_to = "Read_Status",
    values_to = "Num_Paired_Reads") %>%
  mutate(
    Read_Status = dplyr::recode(Read_Status,
                                "Num_Reads_Forward_Trimmed_Paired" = "Trimmed",
                                "Num_Reads_Forward_Raw" = "Raw"))
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
  scale_color_manual(values=reads.palette)+
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "bottom",
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


###Established and Naive, Raw vs Trimmed####
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

#COMPARING CLASSIFIED READS BY KRAKEN#######
kraken_unclassified_reads <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/Paired_end_mode_updated_20260513/Conf_01/unclassifieds_kraken_analytic_matrix.conf_0.1.csv')

#Filtering just samples included in phyloseq.bacteria.samples.dates, adding metadata, calculating percentage classified
kraken_unclassified_reads_samples_metadata <- kraken_unclassified_reads %>%
  filter(SampleID %in% sampleIDs_phyloseq.bacteria.samples.dates)%>%
  dplyr::left_join(cuso4_trimmed_read_counts_samples_metadata, by = "SampleID")%>%
  rename(Kraken2_Input_PairedEnd_Reads = Total, 
         Kraken2_Unclassified_PairedEnd_Reads = NumberUnclassified, 
         Kraken2_Unclassified_Percentage_Reads = PercentUnclassified)%>%
  mutate(Kraken2_Classified_Percentage_Reads = (100 - Kraken2_Unclassified_Percentage_Reads))
nrow(kraken_unclassified_reads_samples_metadata) #Ok, 216 samples

###Kraken2 Classified Percentages Established vs Naive####
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

#COMPARING SAMPLE SUMS (CLASSIFIED READS FROM KRAKEN)#######
##ALL TAXA#####
sample.sums <- sample_sums(phyloseq.bacteria.samples.dates) #making a sample sums object
phyloseq.bacteria.samples.dates.samplessums.df <- cbind(phyloseq.bacteria.samples.dates@sam_data, 
                                      sample.sums) #combining sample sums with metaphyloseq
phyloseq.bacteria.samples.dates.samplessums.df
phyloseq.bacteria.samples.dates.samplessums.df$sampleID <- rownames(phyloseq.bacteria.samples.dates.samplessums.df) ##making a sampleID column


###Established vs Naive####
bacteria_archaea_samplesums_P1vsH21<- ggplot(phyloseq.bacteria.samples.dates.samplessums.df, 
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
wilcox_test(phyloseq.bacteria.samples.dates.samplessums.df, sample.sums~Enclosure) #S. p = 1.98e-13

###SUPPLEMENTARY FIGURE 1######
sfigure1 <- cowplot::plot_grid(sequencing_depth_P1vsH21, 
                               kraken2_classified_read_percentages_P1vsH21,
                               bacteria_archaea_samplesums_P1vsH21,
                               align = "v", 
                               ncol = 1, 
                               labels = "AUTO", label_size = 22)
sfigure1
ggsave("SupplementaryFigure1.svg", 
       sfigure1, 
       device = "svg", width = 8, height =16)

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
wilcox_test(nitrifiers.samplesums.df, sample.sums.nit~Enclosure) #S. p = 0.00228

#RELATIVE ABUNDANCE####
any(sample_sums(phyloseq.bacteria.samples.dates)== 0) ## no samples with 0 OTUs
phyloseq.bacteria.samples.dates.ra <- transform_sample_counts(phyloseq.bacteria.samples.dates, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data
##CLASSIFICATION PERCENTAGES AT DIFFERENT TAXONOMIC LEVELS####
###PHYLUM######
phyloseq.bacteria.samples.dates_phylum.ra <- tax_glom(phyloseq.bacteria.samples.dates.ra, taxrank = "Phylum", NArm = F) 
phyloseq.bacteria.samples.dates_phylum.ra #1815 phyla and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_phylum.ra)[, "Phylum"])) #1815 phyla (so No duplicates)

Unknown_phylum_abundance <- phyloseq.bacteria.samples.dates_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Phylum, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_phylum_abundance ##0.199% abundance by Unknown Phyla

Unclassified_phylum_abundance <- phyloseq.bacteria.samples.dates_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##12.6% abundance by Unclassified Phyla

Classified_phylum_abundance <- phyloseq.bacteria.samples.dates_phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##87.2% abundance by Classified Phyla

##Checking on excel
write.csv(phyloseq.bacteria.samples.dates_phylum.ra@otu_table, "phylum_otus.csv")
write.csv(phyloseq.bacteria.samples.dates_phylum.ra@tax_table, "phylum_taxa.csv")  

#How many unclassified?
phyloseq.bacteria.samples.dates_phylum.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples.dates_phylum.ra)
phyloseq.bacteria.samples.dates_phylum.unclassified.ra #9 unclassified Phyla

#How many unknown?
phyloseq.bacteria.samples.dates_phylum.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples.dates_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples.dates_phylum.ra)
phyloseq.bacteria.samples.dates_phylum.unknown.ra #1690 "unknown" Phyla

#Keep just classified Phyla
phyloseq.bacteria.samples.dates_phylum.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples.dates_phylum.ra)
phyloseq.bacteria.samples.dates_phylum.classified.ra ##116 classified (not unknown or unclassified) Phyla

###CLASS#####
phyloseq.bacteria.samples.dates_class.ra <- tax_glom(phyloseq.bacteria.samples.dates.ra, taxrank = "Class", NArm = F) 
phyloseq.bacteria.samples.dates_class.ra #2865 taxa and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_class.ra)[, "Class"])) #2865 classes (so No duplicates)

Unknown_class_abundance <- phyloseq.bacteria.samples.dates_class.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Class, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_class_abundance #0.287% Abundance by Unknown classes

Unclassified_class_abundance <- phyloseq.bacteria.samples.dates_class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##21.5% Abundance by Unclassified Classes

Classified_class_abundance <- phyloseq.bacteria.samples.dates_class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##78.2% Abundance by Classified classes

##Checking on excel
write.csv(phyloseq.bacteria.samples.dates_class.ra@otu_table, "class_otus.csv")
write.csv(phyloseq.bacteria.samples.dates_class.ra@tax_table, "class_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples.dates_class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_class.ra)[, "Class"]),
  phyloseq.bacteria.samples.dates_class.ra)
phyloseq.bacteria.samples.dates_class.unclassified.ra #70 unclassified classes

#How many unknown?
phyloseq.bacteria.samples.dates_class.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples.dates_class.ra)[, "Class"]),
  phyloseq.bacteria.samples.dates_class.ra)
phyloseq.bacteria.samples.dates_class.unknown.ra #2650 "unknown" classes

#Keep just classified Classes
phyloseq.bacteria.samples.dates_class.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_class.ra)[, "Class"]),
  phyloseq.bacteria.samples.dates_class.ra)
phyloseq.bacteria.samples.dates_class.classified.ra #145 classified classes

###ORDER######
phyloseq.bacteria.samples.dates_order.ra <- tax_glom(phyloseq.bacteria.samples.dates.ra, taxrank = "Order", NArm = F) 
phyloseq.bacteria.samples.dates_order.ra #3400  orders

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_order.ra)[, "Order"])) #3396 orders (3 duplicates)
order_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples.dates_order.ra)[, "Order"])
unique(order_taxa_vec[duplicated(order_taxa_vec)]) 
# "Mycoplasmoidales" "Candidatus Cenarchaeales"  "Candidatus Fermentimicrarchaeales" are the duplicaes

Unknown_order_abundance <- phyloseq.bacteria.samples.dates_order.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Order, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_order_abundance ##0.375% abundance by Unknown Orders

Unclassified_order_abundance <- phyloseq.bacteria.samples.dates_order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##36.7% abundance by Unclassified Orders

Classified_order_abundance <- phyloseq.bacteria.samples.dates_order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##62.9% abundance by Classified orders

#Checking on excel
write.csv(phyloseq.bacteria.samples.dates_order.ra@otu_table, "order_otus.csv")
write.csv(phyloseq.bacteria.samples.dates_order.ra@tax_table, "order_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples.dates_order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_order.ra)[, "Order"]),
  phyloseq.bacteria.samples.dates_order.ra)
phyloseq.bacteria.samples.dates_order.unclassified.ra #131 unclassified orders

#How many unknown?
phyloseq.bacteria.samples.dates_order.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples.dates_order.ra)[, "Order"]),
  phyloseq.bacteria.samples.dates_order.ra)
phyloseq.bacteria.samples.dates_order.unknown.ra #2943 "unknown" orders

#Keep just classified Orders
phyloseq.bacteria.samples.dates_order.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_order.ra)[, "Order"]),
  phyloseq.bacteria.samples.dates_order.ra)
phyloseq.bacteria.samples.dates_order.classified.ra #326 classified orders
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_order.classified.ra)[, "Order"])) ##322 classified orders (unique - without duplicates)

###FAMILY######
phyloseq.bacteria.samples.dates_family.ra <- tax_glom(phyloseq.bacteria.samples.dates.ra, taxrank = "Family", NArm = F) 
phyloseq.bacteria.samples.dates_family.ra #4313 families
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.ra)[, "Family"])) #4307 taxa (5 duplicates)
family_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.ra)[, "Family"])
unique(family_taxa_vec[duplicated(family_taxa_vec)]) 
#"Mycoplasmoidaceae", "unclassified Mycoplasmoidales", "Metamycoplasmataceae", "Candidatus Cenarchaeaceae"         
#"Candidatus Fermentimicrarchaeaceae"

Unknown_family_abundance <- phyloseq.bacteria.samples.dates_family.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_family_abundance #0.506% abundance by Unknown Families

Unclassified_family_abundance <- phyloseq.bacteria.samples.dates_family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##61.9% abundance by Unclassified Families

Classified_family_abundance <- phyloseq.bacteria.samples.dates_family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##37.6% abundance by Classified Families

#Checking on excel
write.csv(phyloseq.bacteria.samples.dates_family.ra@otu_table, "family_otus.csv")
write.csv(phyloseq.bacteria.samples.dates_family.ra@tax_table, "family_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples.dates_family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.ra)[, "Family"]),
  phyloseq.bacteria.samples.dates_family.ra)
phyloseq.bacteria.samples.dates_family.unclassified.ra #270 unclassified families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.unclassified.ra)[, "Family"])) ##269 classified families (unique - without duplicates)

"unclassified Mycoplasmoidales"  
#How many unknown?
phyloseq.bacteria.samples.dates_family.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.ra)[, "Family"]),
  phyloseq.bacteria.samples.dates_family.ra)
phyloseq.bacteria.samples.dates_family.unknown.ra #3298 "unknown" families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.unknown.ra)[, "Family"]))#3298 "unknown" taxa (unique - without duplicates)

#Keep just classified Families
phyloseq.bacteria.samples.dates_family.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.ra)[, "Family"]),
  phyloseq.bacteria.samples.dates_family.ra)
phyloseq.bacteria.samples.dates_family.classified.ra #745 classified families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_family.classified.ra)[, "Family"]))#740 classified families (unique - without duplicates)

###GENUS ######
phyloseq.bacteria.samples.dates_genus.ra <- tax_glom(phyloseq.bacteria.samples.dates.ra, taxrank = "Genus", NArm = F) 
phyloseq.bacteria.samples.dates_genus.ra #7365 genera
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.ra)[, "Genus"])) #7288 taxa (77 duplicates)
genus_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.ra)[, "Genus"])
unique(genus_taxa_vec[duplicated(genus_taxa_vec)]) #41 duplicated unique ones

Unknown_genus_abundance <- phyloseq.bacteria.samples.dates_genus.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_genus_abundance ##0.742%  abundance by unknown genera

Unclassified_genus_abundance <- phyloseq.bacteria.samples.dates_genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance)) #Sum across OTUs
Unclassified_genus_abundance ##80.1% abundance by unclassified genera

Classified_genus_abundance <- phyloseq.bacteria.samples.dates_genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##19.2% abundance by Classified Genera

#Checking on excel
write.csv(phyloseq.bacteria.samples.dates_genus.ra@otu_table, "genus_otus.csv")
write.csv(phyloseq.bacteria.samples.dates_genus.ra@tax_table, "genus_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples.dates_genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples.dates_genus.ra)
phyloseq.bacteria.samples.dates_genus.unclassified.ra #648 unclassified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.unclassified.ra)[, "Genus"])) ##645 unclassified genera (unique - without duplicates)

#How many unknown?
phyloseq.bacteria.samples.dates_genus.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples.dates_genus.ra)
phyloseq.bacteria.samples.dates_genus.unknown.ra #3755 "unknown" genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.unknown.ra)[, "Genus"])) ##3755 unknown genera (unique - without duplicates)


#Keep just classified Genera
phyloseq.bacteria.samples.dates_genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples.dates_genus.ra)
phyloseq.bacteria.samples.dates_genus.classified.ra #2962 classified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_genus.classified.ra)[, "Genus"])) ##2888 classified genera (unique - without duplicates)


###SPECIES######
phyloseq.bacteria.samples.dates.ra ##23461 Species- OTUs
phyloseq.bacteria.samples.dates_species.ra <- phyloseq.bacteria.samples.dates.ra
phyloseq.bacteria.samples.dates_species.ra #23461 Species

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.ra)[, "Species"])) #33365 species (125 duplicates)
species_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.ra)[, "Species"])
unique(species_taxa_vec[duplicated(species_taxa_vec)]) #167 duplicated unique ones

Unknown_species_abundance <- phyloseq.bacteria.samples.dates_species.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Species, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_species_abundance ##0.00000616%  abundance by unknown species

Unclassified_species_abundance <- phyloseq.bacteria.samples.dates_species.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Species, ignore.case = TRUE)) %>%  # Filter unclassified <tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_species_abundance ##86.6% abundance by unclassified species

Classified_species_abundance <- phyloseq.bacteria.samples.dates_species.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Species, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_species_abundance ##13.4% abundance by Classified Species

#Checking on excel
write.csv(phyloseq.bacteria.samples.dates_species.ra@otu_table, "species_otus.csv")
write.csv(phyloseq.bacteria.samples.dates_species.ra@tax_table, "species_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples.dates_species.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.ra)[, "Species"]),
  phyloseq.bacteria.samples.dates_species.ra)
phyloseq.bacteria.samples.dates_species.unclassified.ra #2001 unclassified species
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.unclassified.ra)[, "Species"])) ##1981 unclassified species (unique - without duplicates)

#How many unknown?
phyloseq.bacteria.samples.dates_species.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.ra)[, "Species"]),
  phyloseq.bacteria.samples.dates_species.ra)
phyloseq.bacteria.samples.dates_species.unknown.ra #4 "unknown" species
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.unknown.ra)[, "Species"])) ##4 unknown species (unique - without duplicates)


#Keep just classified Genera
phyloseq.bacteria.samples.dates_species.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.ra)[, "Species"]),
  phyloseq.bacteria.samples.dates_species.ra)
phyloseq.bacteria.samples.dates_species.classified.ra #21,456 classified species
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples.dates_species.classified.ra)[, "Species"])) ##21,159 classified species (unique - without duplicates)


###SUPPLEMENTARY TABLE 1 ####
stable1 <- data.frame(
  "Taxonomic level" = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  # Total taxa at each level
  "Number of Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples.dates_phylum.ra)),
    length(taxa_names(phyloseq.bacteria.samples.dates_class.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_order.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_family.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_genus.ra)),
    length(taxa_names(phyloseq.bacteria.samples.dates_species.ra))),

  
  # Unclassified taxa
  "Number of Unclassified Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples.dates_phylum.unclassified.ra)),
    length(taxa_names(phyloseq.bacteria.samples.dates_class.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_order.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_family.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_genus.unclassified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_species.unclassified.ra))),
  
  # Unknown taxa
  "Number of Unknown Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples.dates_phylum.unknown.ra)),
    length(taxa_names(phyloseq.bacteria.samples.dates_class.unknown.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_order.unknown.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_family.unknown.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_genus.unknown.ra)), 
    NA),
  
  # Classified taxa
  "Number of Classified Taxa" = c(
    length(taxa_names(phyloseq.bacteria.samples.dates_phylum.classified.ra)),
    length(taxa_names(phyloseq.bacteria.samples.dates_class.classified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_order.classified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_family.classified.ra)), 
    length(taxa_names(phyloseq.bacteria.samples.dates_genus.classified.ra)),
    length(taxa_names(phyloseq.bacteria.samples.dates_species.classified.ra))),
  
  # Relative abundances
  "Mean Relative Abundance (%) of Unclassified Taxa Across Samples" = c(
    Unclassified_phylum_abundance$Unclassified_sum, 
    Unclassified_class_abundance$Unclassified_sum,
    Unclassified_order_abundance$Unclassified_sum,
    Unclassified_family_abundance$Unclassified_sum,
    Unclassified_genus_abundance$Unclassified_sum,
    Unclassified_species_abundance$Unclassified_sum
  ),
  "Mean Relative Abundance (%) of Unknown Taxa Across Samples" = c(
    Unknown_phylum_abundance$Unknown_sum, 
    Unknown_class_abundance$Unknown_sum, 
    Unknown_order_abundance$Unknown_sum,
    Unknown_family_abundance$Unknown_sum,
    Unknown_genus_abundance$Unknown_sum,
    NA),
  "Mean Relative Abundance (%) of Classified Taxa Across Samples" = c(
    Classified_phylum_abundance$Classified_sum,
    Classified_class_abundance$Classified_sum, 
    Classified_order_abundance$Classified_sum,
    Classified_family_abundance$Classified_sum,
    Classified_genus_abundance$Classified_sum,
    Classified_species_abundance$Classified_sum),
  check.names = FALSE
) %>%
  mutate(`Percentage of Classified Taxa` = 
           (`Number of Classified Taxa` / `Number of Taxa`) * 100)
stable1
#Make into excel file
write_xlsx(stable1, 
          "SupplementaryTable1.xlsx")

###SUPPLEMENTARY TABLE 2######
#CSVs from classifications at each taxonomic level

#ALPHA DIVERSITY ######
## ALL COMMUNITIES#####
alpha_div1 <- phyloseq::estimate_richness(phyloseq.bacteria.samples.dates, 
                                          measures = c("Observed", "Shannon")) # richness, diversity
alpha_div2 <- microbiome::evenness(phyloseq.bacteria.samples.dates, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness

#Combine alpha div metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div

#Metadata and div metrics
alpha_div_meta <- cbind(phyloseq.bacteria.samples.dates@sam_data, 
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
metadata_phyloseq.bacteria.samples.dates <- data.frame(phyloseq.bacteria.samples.dates@sam_data)
#When do phases switch?
#Naive
metadata_phyloseq.bacteria.samples.dates %>%
  filter(Enclosure == "H21") %>%
  group_by(Date_num_phase_naive) %>%
  summarise(
    min_naive = min(Date_num_naive, na.rm = TRUE),
    max_naive = max(Date_num_naive, na.rm = TRUE)
  )

#Established
metadata_phyloseq.bacteria.samples.dates %>%
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
  labs(title = "MICROBIOME\n  ",
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
nitrifiers #449 taxa and 216 samples
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

####Adding vertical dashed line for when "phase" changes#####
metadata_phyloseq.bacteria.samples.dates <- data.frame(phyloseq.bacteria.samples.dates@sam_data)
#When do phases switch?
#Naive
metadata_phyloseq.bacteria.samples.dates %>%
  filter(Enclosure == "H21") %>%
  group_by(Date_num_phase_naive) %>%
  summarise(
    min_naive = min(Date_num_naive, na.rm = TRUE),
    max_naive = max(Date_num_naive, na.rm = TRUE)
  )

#Established
metadata_phyloseq.bacteria.samples.dates %>%
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
  labs(title = "NITRIFYING TAXA\n  ",
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)")))+
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
alpha_div_nit_wq_date_num_factor_other_metadata

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
phyloseq.bacteria.samples.dates_order.ra #3400 taxa and 216 samples 

#Grouping the low abundance orders into one category
phyloseq.bacteria.samples.dates.order.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.samples.dates_order.ra, 
                                                                        "Enclosure", 
                                                                        level = "Order", 
                                                                        threshold = 0.5)
phyloseq.bacteria.samples.dates.order.filt #20 orders over 0.5% mean RA
phyloseq.bacteria.samples.dates.order.filt.melt <- psmelt(phyloseq.bacteria.samples.dates.order.filt)%>%
  mutate(Order = factor(Order, 
                         levels = c(setdiff(Order, 
                                            unique(grep("Others", Order, value = TRUE))), 
                                    unique(grep("Others", Order, value = TRUE)))))##Factoring the Order column so that "Others.." is the last category
levels(phyloseq.bacteria.samples.dates.order.filt.melt$Order) ##ok

##Create color palette
order.filt.palette <- distinctColorPalette(length(unique(phyloseq.bacteria.samples.dates.order.filt.melt$Order)))
#order.filt.palette <- unname(alphabet2())
order_filt_names <- unique(phyloseq.bacteria.samples.dates.order.filt.melt$Order)# Create a named vector for the palette, where the names correspond to phlyum names
order_named_palette <- setNames((order.filt.palette)[1:length(order_filt_names)], order_filt_names)
order_named_palette$'Others <0.5% RA' <- "grey95"
order_named_palette$'Flavobacteriales' <-  "#63A184"
order_named_palette$'Rhodobacterales' <- "#E3B199"
order_named_palette$'unclassified Alphaproteobacteria' <- "dodgerblue"
##Apply the function to obtain top orders (n=15)
top_orders <- top_taxa_legend(phyloseq.bacteria.samples.dates.order.filt.melt, 
                              taxlevel = "Order", n = 15)
top_orders

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_order_enclosures_overall_plot_datenum <- ggplot(phyloseq.bacteria.samples.dates.order.filt.melt,
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
  scale_fill_manual(values = order_named_palette,
                    breaks = top_orders,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
        #legend.position = "right",
        legend.position = c(1.09, 0.5),  # x, y inside plot
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

####  Together with the alpha div measures #######
figure_alpha_overall_div_order_time_copper <- alpha_div_wq_date_num_factor_other_metadata + 
  theme(axis.text.x = element_blank())+
  RA_order_enclosures_overall_plot_datenum  +
  plot_layout(ncol = 1, heights = c(1.1, 0.9))
figure_alpha_overall_div_order_time_copper
ggsave("figure_alpha_overall_div_order_time_copper.png", 
       figure_alpha_overall_div_order_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 15, 
       width = 26)

####  Together with the alpha div measures - making smaller to fit together with antimicrobial resistance#######
figure_alpha_overall_div_order_time_copper_smaller <- alpha_div_wq_date_num_factor + 
  theme(axis.text.x = element_blank())+
  RA_order_enclosures_overall_plot_datenum  + 
  theme(legend.position = "none")+
  theme(legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"))+
  plot_layout(ncol = 1, heights = c(0.8, 1.2))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_overall_div_order_time_copper_smaller
ggsave("figure_alpha_overall_div_order_time_copper_smaller.png", 
       figure_alpha_overall_div_order_time_copper_smaller, 
       device = "png", 
       dpi = 600, 
       height = 12, 
       width = 26)

## FAMILY #####
phyloseq.bacteria.samples.dates_family.ra #4313 families and 216 samples 

phyloseq.bacteria.samples.dates.family.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.samples.dates_family.ra, 
                                                                             "Enclosure", 
                                                                             level = "Family", 
                                                                             threshold = 0.5)
phyloseq.bacteria.samples.dates.family.filt #27 families over 0.5% mean RA
phyloseq.bacteria.samples.dates.family.filt.melt <- psmelt(phyloseq.bacteria.samples.dates.family.filt)%>%
  mutate(Family = factor(Family, 
                        levels = c(setdiff(Family, 
                                           unique(grep("Others", Family, value = TRUE))), 
                                   unique(grep("Others", Family, value = TRUE)))))##Factoring the Family column so that "Others.." is the last category
levels(phyloseq.bacteria.samples.dates.family.filt.melt$Family) ##ok

##Create color palette - based on families within the same order
palette_family_level_df <- phyloseq.bacteria.samples.dates.family.filt.melt %>% 
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
family_named_palette$'Others <0.5% RA' <- "grey95"

##Apply the function to obtain top familys (n=15)
top_families <- top_taxa_legend(phyloseq.bacteria.samples.dates.family.filt.melt, 
                              taxlevel = "Family", n = 15)
top_families

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_family_enclosures_overall_plot_datenum <- ggplot(phyloseq.bacteria.samples.dates.family.filt.melt,
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
  scale_fill_manual(values = family_named_palette,
                    breaks = top_families,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.08, 0.5),  # x, y inside plot
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

####  Together with the alpha div measures #######
figure_alpha_overall_div_family_time_copper <- alpha_div_wq_date_num_factor_other_metadata + 
  theme(axis.text.x = element_blank())+
  RA_family_enclosures_overall_plot_datenum  + 
  plot_layout(ncol = 1, heights = c(1.1, 0.9))
figure_alpha_overall_div_family_time_copper
ggsave("figure_alpha_overall_div_family_time_copper.png", 
       figure_alpha_overall_div_family_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 15, 
       width = 26)


##RHOBACTERALES only##########
phyloseq.bacteria.samples.dates_order.ra #5883 orders 

#Out of this overall communities object, select only rhodobacterales 
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales <- subset_taxa(phyloseq.bacteria.samples.dates_order.ra, 
                                                              Order == "Rhodobacterales") # NOB; some, plus a new one!
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales <- subset_samples(phyloseq.bacteria.samples.dates_order.ra.rhodobacterales, 
                                                                 sample_sums(phyloseq.bacteria.samples.dates_order.ra.rhodobacterales) > 0)
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales #Rhodobacterales (1 taxa) in 218 samples 


#Melt to plot 
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt <- psmelt(phyloseq.bacteria.samples.dates_order.ra.rhodobacterales)


###Correlation of rhodobacterales with Copper levels#########
####H21########
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt.H21 <- phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt%>%
  filter(Order == "Rhodobacterales")%>%
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
copper_rhodobacterales_relationship_plot_H21 <- ggplot(phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt.H21,
                                           aes(x = Copper_level_mg_L, 
                                               y = Abundance)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Collection_Month)) +
  facet_wrap(~Attempt,
             scales = "free_y",
             ncol = 1)+
  scale_fill_manual(values = month_palette)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(y = "Rhodobacterales RA (%)",
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
copper_rhodobacterales_relationship_plot_H21

ggsave("copper_rhodobacterales_relationship_plot_H21.png", 
       copper_rhodobacterales_relationship_plot_H21, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)


####P1########
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt.P1 <- phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt%>%
  filter(Enclosure == "P1")%>%
  filter(Order == "Rhodobacterales")%>%
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
phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt.P1

#Plot
copper_rhodobacterales_relationship_plot_P1 <- ggplot(phyloseq.bacteria.samples.dates_order.ra.rhodobacterales.melt.P1,
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
  labs(y = "Rhodobacterales RA (%)",
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
        axis.text.y = element_text(colour = "black", size = 25),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 50, face = "bold"))+
  guides(color = guide_legend(override.aes = list(size = 7),
                              nrow = 1))
copper_rhodobacterales_relationship_plot_P1
ggsave("copper_rhodobacterales_relationship_plot_P1.png", 
       copper_rhodobacterales_relationship_plot_P1, 
       device = "png", 
       dpi = 600, 
       height = 11, 
       width = 23)


#NITRIFIERS WITHIN THE OVERALL COMMUNITY#####
## FAMILY #######
phyloseq.bacteria.samples.dates_family.ra #4313 families

##Which families are nitrifiers? 
nitrifiers.melt <- psmelt(nitrifiers)
unique(nitrifiers.melt$Family) #"Nitrosopumilaceae", "Nitrobacteraceae","Nitrospinaceae"             
#"Nitrospiraceae", "Chromatiaceae", "Ectothiorhodospiraceae"     
#"Nitrosomonadaceae", "Gallionellaceae", "Nitrososphaeraceae"         
#"Candidatus Nitrosocaldaceae"

#Out of this overall communities object, select only nitrifiers 
phyloseq.bacteria.samples.dates_family.ra.nitrifiers <- subset_taxa(phyloseq.bacteria.samples.dates_family.ra, 
                                                              Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                                                                Family == "Chromatiaceae" | # no lineages
                                                                Family == "Nitrosopumilaceae" | # AOA; some!
                                                                Family == "Nitrososphaeraceae" | # no lineages
                                                                Family == "Candidatus Nitrosocaldaceae" | # no lineages
                                                                Family == "Nitrospiraceae" | # NOB/Commamox; some!
                                                                Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                                                                Family == "Nitrobacteraceae" | # none
                                                                Family == "Gallionellaceae" | # none
                                                                Family == "Nitrospinaceae") # NOB; some, plus a new one!
phyloseq.bacteria.samples.dates_family.ra.nitrifiers <- subset_samples(phyloseq.bacteria.samples.dates_family.ra.nitrifiers, 
                                                                 sample_sums(phyloseq.bacteria.samples.dates_family.ra.nitrifiers) > 0)
phyloseq.bacteria.samples.dates_family.ra.nitrifiers #10 nitrifying families in 216 samples 


#Melt to plot 
phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt <- psmelt(phyloseq.bacteria.samples.dates_family.ra.nitrifiers)


##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt <- phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
      Family == "Nitrosomonadaceae" ~ "AOB",
      Family == "Chromatiaceae" ~ "AOB",
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
palette_nitrifiers_family_df <- phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt %>% 
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
top_nitrifying_families <- top_taxa_legend(phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt, 
                                           n = 5)
# top_nitrifying_families <- c("Nitrosopumilaceae", #AOA
#                              "Chromatiaceae",#AOB
#                              "Nitrosomonadaceae",#AOB
#                              "Nitrobacteraceae",#NOB
#                              "Nitrospiraceae",#NOB
#                              "Ectothiorhodospiraceae")#NOB
# top_nitrifying_families

#Plot
RA_family_enclosures_nit_plot_datenum <- ggplot(phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt,
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
  scale_fill_manual(values = palette_nitrifiers_family,
                    breaks = top_nitrifying_families
                    ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.5),  # x, y inside plot
    legend.text = element_text(size = 21),
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

####  Together with the alpha div measures #######
figure_alpha_nit_div_family_time_copper <- alpha_div_nit_wq_date_num_factor_other_metadata + 
  theme(axis.text.x = element_blank())+
  RA_family_enclosures_nit_plot_datenum  + 
  plot_layout(ncol = 1, heights = c(1.1, 0.9))
figure_alpha_nit_div_family_time_copper
ggsave("figure_alpha_nit_div_family_time_copper.png", 
       figure_alpha_nit_div_family_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 15, 
       width = 26)

## SPECIES#######
phyloseq.bacteria.samples.dates_species.ra #23461 species


#Out of this overall communities object, select only nitrifiers 
phyloseq.bacteria.samples.dates_species.ra.nitrifiers <- subset_taxa(phyloseq.bacteria.samples.dates_species.ra, 
                                                                    Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                                                                      Family == "Chromatiaceae" | # no lineages
                                                                      Family == "Nitrosopumilaceae" | # AOA; some!
                                                                      Family == "Nitrososphaeraceae" | # no lineages
                                                                      Family == "Candidatus Nitrosocaldaceae" | # no lineages
                                                                      Family == "Nitrospiraceae" | # NOB/Commamox; some!
                                                                      Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                                                                      Family == "Nitrobacteraceae" | # none
                                                                      Family == "Gallionellaceae" | # none
                                                                      Family == "Nitrospinaceae") # NOB; some, plus a new one!
phyloseq.bacteria.samples.dates_species.ra.nitrifiers <- subset_samples(phyloseq.bacteria.samples.dates_species.ra.nitrifiers, 
                                                                       sample_sums(phyloseq.bacteria.samples.dates_species.ra.nitrifiers) > 0)
phyloseq.bacteria.samples.dates_species.ra.nitrifiers #449 nitrifying species in 216 samples 


#Melt to plot 
phyloseq.bacteria.samples.dates_species.ra.nitrifiers.melt <- psmelt(phyloseq.bacteria.samples.dates_species.ra.nitrifiers)


##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.bacteria.samples.dates_species.ra.nitrifiers.melt <- phyloseq.bacteria.samples.dates_species.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Chromatiaceae" ~ "AOB",
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
palette_nitrifiers_species_df <- phyloseq.bacteria.samples.dates_species.ra.nitrifiers.melt %>% 
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
top_nitrifying_species <- top_taxa_legend(phyloseq.bacteria.samples.dates_species.ra.nitrifiers.melt, 
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
RA_species_enclosures_nit_plot_datenum <- ggplot(phyloseq.bacteria.samples.dates_species.ra.nitrifiers.melt,
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

####  Together with the alpha div measures #######
figure_alpha_nit_div_species_time_copper <- alpha_div_nit_wq_date_num_factor_other_metadata + 
  theme(axis.text.x = element_blank())+
  RA_species_enclosures_nit_plot_datenum  + 
  plot_layout(ncol = 1, heights = c(1.1, 0.9))
figure_alpha_nit_div_family_time_copper
ggsave("figure_alpha_nit_div_species_time_copper.png", 
       figure_alpha_nit_div_species_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 15, 
       width = 26)


#NITRIFIERS ONLY######
any(sample_sums(nitrifiers)== 0) ## no samples with 0 OTUs

nitrifiers.ra <- transform_sample_counts(nitrifiers, 
                                         function(x) x/sum(x)*100) ##Relative abundance from normalized data

nitrifiers.ra #449 taxa and 216 samples

###FAMILY#######
nitrifiers.ra.family <- tax_glom(nitrifiers.ra, taxrank = "Family", NArm = F)
nitrifiers.ra.family #10 families 
##Melt 
nitrifiers.ra.family.melt <- psmelt(nitrifiers.ra.family)

##Which are the top most abundant taxa by group? 
nitrifiers.ra.family.melt %>%
  group_by(Enclosure, Family) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(Enclosure,  desc(mean_abun))%>%
  print(n=40)


##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
nitrifiers.ra.family.melt <- nitrifiers.ra.family.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Chromatiaceae" ~ "AOB",
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
  scale_fill_manual(values = palette_nitrifiers_only_family) +
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
RA_enclosures_nitrifiers_only_family.plot


#####SPECIES#######
nitrifiers.ra.species <- tax_glom(nitrifiers.ra, taxrank = "Species", NArm = F)
nitrifiers.ra.species #449 species and 216 samples

#Merge low abun species
nitrifiers.ra.species.filt <- merge_low_abundance_grouped_ra(nitrifiers.ra.species, 
                                                             "Enclosure", 
                                                             level = "Species", threshold = 0.5)
nitrifiers.ra.species.filt #16 Species over 0.5% mean RA
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


#Top species for the legend
top_nitrifying_only_species <- top_taxa_legend(nitrifiers.ra.species.filt.melt, 
                                               taxlevel = "Species",
                                               n = 8)
top_nitrifying_only_species <- factor(top_nitrifying_only_species, 
                                      levels = c("unclassified Nitrosopumilus", 
                                                 "Nitrosopumilus maritimus", 
                                                 "Candidatus Nitrosopumilus sp. SW", 
                                                 "unclassified Nitrosopumilaceae",
                                                 "Nitrosopumilus piranensis", 
                                                 "unclassified Bradyrhizobium", 
                                                 "Candidatus Nitronauta litoralis", 
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
      
      if (group == "AOA") {
        colorRampPalette(c(
          "#FCE4EC",  # very light pink
          "#880E4F",  # wine
          "#4A148C",   # deep purple
          "#FF4081", # hot pink
          "#C2185B"  # deep magenta
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
  palette_nitrifiers_only_species_df $color,
  palette_nitrifiers_only_species_df $Species)
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
  geom_vline(data = line_breaks_phases,
             aes(xintercept = Date_num),
             linetype = "dashed",
             color = "black",
             alpha = 0.8) +
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
  scale_fill_manual(
    values = palette_nitrifiers_only_species,
    breaks = top_nitrifying_only_species,
    labels = function(x) str_wrap(x, width = 20), 
    drop = FALSE
  )+
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.07, 0.5),  # x, y inside plot
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
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
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
RA_family_enclosures_nit_plot_datenum_2 <- RA_family_enclosures_nit_plot_datenum + 
  theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17),
        axis.text.y = element_text(size = 12)) 
RA_enclosures_nitrifiers_only_species.plot_2 <- RA_enclosures_nitrifiers_only_species.plot +
  theme(
    axis.title.y = element_text(size = 17),
    axis.text.y = element_text(size = 12))

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
       height = 15, 
       width = 26)


#Correlation of Nitrosopumilaceae with Copper levels#########
####H21########
phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21 <- phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt%>%
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

copper_AOA_relationship_plot_H21 <- ggplot(phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21,
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
                         data = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21)
summary(gam_model_nit_AOA_H21) ##No effect of enclosure, but copper effect did vary between enclosures 
plot(gam_model_nit_AOA_H21, pages = 1, shade = TRUE)

#####Spearman correlation#####
cor.test(x = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21$Abundance, 
         y = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21$Copper_level_mg_L, 
         method = 'spearman') #Significant for naive one

H21_pcor_AOA <- pcor.test(x = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21$Abundance,
                          y = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21$Copper_level_mg_L,
                          z = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.H21$Date_num,
                          method = "pearson")

####P1######
phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1 <- phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt%>%
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
copper_AOA_relationship_plot_P1 <- ggplot(phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1,
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
phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1.clean <- phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1 %>%
  filter(!is.na(Copper_level_mg_L),
         !is.na(Collection_Date))%>%
  mutate(Enclosure = "P1")

model_lm_nit_AOA_P1 <- lm(Abundance ~ Copper_level_mg_L + I(Copper_level_mg_L^2) +
                         Collection_Month +
                         Copper_level_mg_L:Collection_Month +
                         I(Copper_level_mg_L^2):Collection_Month,
                       data = phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1.clean)
summary(model_lm_nit_AOA_P1)

#Confidence Intervals
confint(model_lm_nit_AOA_P1)

#Anova type3 - instead of testing each coefficient individually, Type III ANOVA tests the factor as a whole.
Anova(model_lm_nit_AOA_P1, type = "III") 

# Add fitted values to your data
phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1.clean$fitted <- fitted(model_lm_nit_AOA_P1)

# Plot actual vs fitted
ggplot(phyloseq.bacteria.samples.dates_family.ra.AOA.melt.P1.clean, aes(x = fitted, y = Abundance)) +
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
phyloseq.bacteria.samples.dates_family.ra.NOB.melt.H21 <- phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt%>%
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
copper_NOB_relationship_plot_H21 <- ggplot(phyloseq.bacteria.samples.dates_family.ra.NOB.melt.H21,
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
phyloseq.bacteria.samples.dates_family.ra.NOB.melt.P1 <- phyloseq.bacteria.samples.dates_family.ra.nitrifiers.melt%>%
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
copper_NOB_relationship_plot_P1 <- ggplot(phyloseq.bacteria.samples.dates_family.ra.NOB.melt.P1,
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



#BRAY CURTIS#####
##Going to take out samples from P1 from april and may 2023 and september from H21
phyloseq.bacteria.samples.dates.ra.dates <- subset_samples(phyloseq.bacteria.samples.dates.ra, 
                                                     Collection_Date > "2023-10-01")
phyloseq.bacteria.samples.dates.ra.dates #33490 taxa and 218 samples 
#Take out taxa sums = 0
phyloseq.bacteria.samples.dates.ra.dates <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.ra.dates) > 0, 
                                                 phyloseq.bacteria.samples.dates.ra.dates)
phyloseq.bacteria.samples.dates.ra.dates # 33490 taxa and 218 samples 

#BC distances
phyloseq.bacteria.samples.dates.ra.dates.bray <- vegdist(t(phyloseq.bacteria.samples.dates.ra.dates@otu_table), method = "bray")
phyloseq.bacteria.samples.dates.ra.dates.bray

#make DF from metadata
phyloseq.bacteria.samples.dates.df <- as(phyloseq.bacteria.samples.dates.ra.dates@sam_data, "data.frame") %>%
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
      Copper_quartile == "Q4_high" ~ "Q4")
    )%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>% #Date as number since start 
  mutate()
phyloseq.bacteria.samples.dates.df 

##ORDINATION####
set.seed(98)
phyloseq.bacteria.samples.dates.ra.dates.bray.ord <- metaMDS(phyloseq.bacteria.samples.dates.ra.dates.bray, 
                                                 k=3, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
###Centroids by enclosure #####
## BC
#Simple ordination plot
phyloseq.bacteria.samples.dates.ra.dates.bray.plot <- ordiplot(phyloseq.bacteria.samples.dates.ra.dates.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples.dates.ra.dates.bray.scrs <- scores(phyloseq.bacteria.samples.dates.ra.dates.bray.plot, display = "sites")
#Add metadata to coordinates
phyloseq.bacteria.samples.dates.ra.dates.bray.scrs <- cbind(as.data.frame(phyloseq.bacteria.samples.dates.ra.dates.bray.scrs),
                                                Copper_level_mg_L = phyloseq.bacteria.samples.dates.df$Copper_level_mg_L, 
                                                SampleID = phyloseq.bacteria.samples.dates.df$SampleID, 
                                                Collection_Date = phyloseq.bacteria.samples.dates.df$Collection_Date, 
                                                Enclosure = phyloseq.bacteria.samples.dates.df$Enclosure, 
                                                Attempt = phyloseq.bacteria.samples.dates.df$Attempt, 
                                                Copper_quartile = phyloseq.bacteria.samples.dates.df$Copper_quartile,
                                                Copper_quartile.abvr = phyloseq.bacteria.samples.dates.df$Copper_quartile.abvr)%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.dates.ra.dates.bray.cent.enclosure <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
                                                     data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs, 
                                                     FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs <- merge(phyloseq.bacteria.samples.dates.ra.dates.bray.scrs, 
                             setNames(phyloseq.bacteria.samples.dates.ra.dates.bray.cent.enclosure, 
                                      c("Enclosure","cMDS1","cMDS2")),
                             by = 'Enclosure', 
                             sort = F)
####MIRKAT- Just enclosure#####
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.ra.dates.bray))

#Modelling enclosure, have to make it dichotomus 
Enclosure_d_mirkat <- phyloseq.bacteria.samples.dates.df%>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_mirkat

# Covariates (date_num)
datenum_covariate <- phyloseq.bacteria.samples.dates.df%>%
  group_by(Enclosure)%>%
  #Date number 
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate 


#Run model
set.seed(98)
mirkat_enclosure <- MiRKAT(
  y = Enclosure_d_mirkat,
  X = datenum_covariate,
  Ks = bray_kernel,
  out_type = "D", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_enclosure #R2 29% , p_values = 0

####PERMANOVA - Just enclosure####
set.seed(98)
enclosure_BC_adonis  <- adonis2(phyloseq.bacteria.samples.dates.ra.dates.bray ~ Enclosure, 
                                       phyloseq.bacteria.samples.dates.df, 
                                       strata = phyloseq.bacteria.samples.dates.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
enclosure_BC_adonis #32.6% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_adonis  <- adonis2(phyloseq.bacteria.samples.dates.ra.dates.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.dates.df, 
#                                 strata = phyloseq.bacteria.samples.dates.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_adonis  <- adonis2(phyloseq.bacteria.samples.dates.ra.dates.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.dates.df , 
#                                 strata = phyloseq.bacteria.samples.dates.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS
#ENCLOSURE#
# Run the betadisper function, average distance to centroid
bray.enclosure.disp <- betadisper(phyloseq.bacteria.samples.dates.ra.dates.bray, 
                                  phyloseq.bacteria.samples.dates.df$Enclosure)
bray.enclosure.disp
##Then test by permuting
set.seed(98)
bray.enclosure.permdisp <- permutest(bray.enclosure.disp, permutations = 9999)
bray.enclosure.permdisp ##S, p = 1e-04

# #COPPER#
# # Run the betadisper function, average distance to centroid
# bray.copper.disp <- betadisper(phyloseq.bacteria.samples.dates.ra.dates.bray, 
#                                   phyloseq.bacteria.samples.dates.df$Copper_level_mg_L)
# bray.copper.disp
# ##Then test by permuting
# set.seed(98)
# bray.copper.permdisp <- permutest(bray.copper.disp, permutations = 9999)
# bray.copper.permdisp ##S, p = 0.044

# Extract R2 and p-values (mirkat)
R2_adonis_enclosure <- enclosure_BC_adonis$R2[1] 
pvalue_adonis_enclosure <- enclosure_BC_adonis$`Pr(>F)`[1]
R2_mirkat_enclosure <- mirkat_enclosure$R2
pvalue_mirkat_enclosure<-  mirkat_enclosure$p_values


#### PLOTS
## BC
enclosure_BC_beta_div <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
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
  annotate("text", x = -1, y = -0.9,
           label = "MiRKAT\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1, y = -0.85,
           label = paste("R² = ", round(R2_mirkat_enclosure* 100, 1), "%",
                         "\np = ", round(pvalue_mirkat_enclosure, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_beta_div
ggsave("enclosure_BC_beta_div.png", 
       enclosure_BC_beta_div, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 14)


##Colored by collection month and shape by ecnlosure 
enclosure_BC_beta_div_2 <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
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
  annotate("text", x = -1, y = -0.9,
           label = "PERMANOVA\nSystem",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -1, y = -0.85,
           label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
enclosure_BC_beta_div_2
ggsave("enclosure_BC_beta_div_2.png", 
       enclosure_BC_beta_div_2, 
       device = "png", 
       dpi = 600, 
       height = 8, 
       width = 15)

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
phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs_3 <- 
  phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs %>%
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
enclosure_BC_beta_div_3 <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs_3) + 
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
enclosure_BC_beta_div_4 <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs_3) + 
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
enclosure_BC_beta_div_5 <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs_3) + 
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
  geom_point(data = phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs_3 %>%
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
  geom_point(data = phyloseq.bacteria.samples.dates.ra.dates.bray.enclosure.segs_3 %>%
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
# fit <- envfit(phyloseq.bacteria.samples.dates.ra.dates.bray.ord, 
#               phyloseq.bacteria.samples.dates.df$Copper_level_mg_L, 
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

###Centroids by copper quartiles facetted by enclosure #####
## BC
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.dates.ra.dates.bray.cent.copper_quart <- aggregate(
  cbind(MDS1, MDS2) ~ Enclosure + Copper_quartile,
  data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs,
  FUN = mean)

#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.dates.ra.dates.bray.copper_quart.segs <- merge(phyloseq.bacteria.samples.dates.ra.dates.bray.scrs, 
                                                          setNames(phyloseq.bacteria.samples.dates.ra.dates.bray.cent.copper_quart, 
                                                                   c("Enclosure", "Copper_quartile",
                                                                     "cMDS1","cMDS2")),
                                                          by = c("Enclosure", "Copper_quartile"), 
                                                          sort = F)
##PERMANOVA - Just copper###
set.seed(98)
copper_quart_BC_adonis  <- adonis2(phyloseq.bacteria.samples.dates.ra.dates.bray ~ Copper_quartile, 
                                phyloseq.bacteria.samples.dates.df, 
                                strata = phyloseq.bacteria.samples.dates.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
copper_quart_BC_adonis 


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_adonis  <- adonis2(phyloseq.bacteria.samples.dates.ra.dates.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.dates.df, 
#                                 strata = phyloseq.bacteria.samples.dates.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_adonis  <- adonis2(phyloseq.bacteria.samples.dates.ra.dates.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.dates.df , 
#                                 strata = phyloseq.bacteria.samples.dates.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS#
#COPPER quartile#
# Run the betadisper function, average distance to centroid
bray.copper.quart.disp <- betadisper(phyloseq.bacteria.samples.dates.ra.dates.bray,
                                  phyloseq.bacteria.samples.dates.df$Copper_quartile)
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
copper_quart_BC_beta_div <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.copper_quart.segs) + 
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
phyloseq.bacteria.samples.dates.df <- as(phyloseq.bacteria.samples.dates.ra@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.df 


###Continuous copper facetted by enclosure #####
#No samples without copper levels:
phyloseq.bacteria.samples.dates.ra.copper <- subset_samples(phyloseq.bacteria.samples.dates.ra.dates, 
                                                                !is.na(Copper_level_mg_L))
phyloseq.bacteria.samples.dates.ra.copper <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.ra.copper)> 0, 
                                                            phyloseq.bacteria.samples.dates.ra.copper)
phyloseq.bacteria.samples.dates.ra.copper #33481 taxa and 204 samples

####H21 #######
phyloseq.bacteria.samples.dates.ra.copper.H21 <- subset_samples(phyloseq.bacteria.samples.dates.ra.copper, 
                                                                    Enclosure == "H21")
phyloseq.bacteria.samples.dates.ra.copper.H21 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.ra.copper.H21)> 0, 
                                                                phyloseq.bacteria.samples.dates.ra.copper.H21)
phyloseq.bacteria.samples.dates.ra.copper.H21 #33418 taxa and 82 samples 

##BC 
phyloseq.bacteria.samples.dates.ra.copper.H21.bray <- vegdist(t(phyloseq.bacteria.samples.dates.ra.copper.H21@otu_table), 
                                                                  method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.dates.copper.H21.df <- as(phyloseq.bacteria.samples.dates.ra.copper.H21@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.copper.H21.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_H21 <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.ra.copper.H21.bray))

# Covariates (date_num)
datenum_covariate_H21 <- phyloseq.bacteria.samples.dates.copper.H21.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_H21 

# # Covariates (month)
# month_covariate_H21 <- model.matrix(~ Collection_Month, 
#                                         data = phyloseq.bacteria.samples.dates.copper.H21.df)[, -1]
# 

#Run model
mirkat_H21_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.dates.copper.H21.df$Copper_level_mg_L,
  X = datenum_covariate_H21,
  Ks = bray_kernel_H21,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_H21_copper_cont #R2 6.9% , p_values = 0.0001017485

#####PERMANOVA #######
set.seed(98)
adonis_H21_copper_cont  <- adonis2(phyloseq.bacteria.samples.dates.ra.copper.H21.bray ~ Copper_level_mg_L, 
                                       phyloseq.bacteria.samples.dates.copper.H21.df, 
                                       strata = phyloseq.bacteria.samples.dates.copper.H21.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_H21_copper_cont # R2 13.5%, p value 2e-04 *

##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.copper.cont.disp.H21 <- betadisper(phyloseq.bacteria.samples.dates.ra.copper.H21.bray,
                                            phyloseq.bacteria.samples.dates.copper.H21.df$Copper_quartile)
bray.copper.cont.disp.H21
##Then test by permuting
set.seed(98)
bray.copper.cont.permdisp.H21 <- permutest(bray.copper.cont.disp.H21, permutations = 9999)
bray.copper.cont.permdisp.H21 ##S, p =0.0016

####P1 #######
phyloseq.bacteria.samples.dates.ra.copper.P1 <- subset_samples(phyloseq.bacteria.samples.dates.ra.copper, 
                                                                   Enclosure == "P1")
phyloseq.bacteria.samples.dates.ra.copper.P1 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.ra.copper.P1)> 0, 
                                                               phyloseq.bacteria.samples.dates.ra.copper.P1)
phyloseq.bacteria.samples.dates.ra.copper.P1 #28427 taxa and 122 samples

##BC 
phyloseq.bacteria.samples.dates.ra.copper.P1.bray <- vegdist(t(phyloseq.bacteria.samples.dates.ra.copper.P1@otu_table), 
                                                                 method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.dates.copper.P1.df <- as(phyloseq.bacteria.samples.dates.ra.copper.P1@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.copper.P1.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_P1 <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.ra.copper.P1.bray))

# Covariates (date_num)
datenum_covariate_P1 <- phyloseq.bacteria.samples.dates.copper.P1.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_P1 

# # Covariates (month)
# month_covariate_P1 <- model.matrix(~ Collection_Month, 
#                                        data = phyloseq.bacteria.samples.dates.copper.P1.df)[, -1]
# 
# month_covariate_P1

#Run model
mirkat_P1_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.dates.copper.P1.df$Copper_level_mg_L,
  X = datenum_covariate_P1,
  Ks = bray_kernel_P1,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_P1_copper_cont #p val 1.320553e-08, R2 2.5%

#####PERMANOVA #####
set.seed(98)
adonis_P1_copper_cont  <- adonis2(phyloseq.bacteria.samples.dates.ra.copper.P1.bray ~ Copper_level_mg_L, 
                                      phyloseq.bacteria.samples.dates.copper.P1.df, 
                                      strata = phyloseq.bacteria.samples.dates.copper.P1.df$Collection_Month,
                                      by = "margin",
                                      permutations = 9999)
adonis_P1_copper_cont # R2 10.4%, p value 0.021 *


####MIRKAT GLMM FOR COPPER #####
#make DF from metadata
phyloseq.bacteria.samples.dates.copper.all.df <- as(phyloseq.bacteria.samples.dates.ra.copper@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.copper.all.df

#Bray-Curtis distances
phyloseq.bacteria.samples.dates.ra.copper.bray <- vegdist(t(phyloseq.bacteria.samples.dates.ra.copper@otu_table), 
                                                              method = "bray")
phyloseq.bacteria.samples.dates.ra.copper.bray

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.ra.copper.bray))

#Jaccard distances
phyloseq.bacteria.samples.dates.ra.copper.jac <- vegdist(t(phyloseq.bacteria.samples.dates.ra.copper@otu_table), 
                                                             method = "jaccard")
phyloseq.bacteria.samples.dates.ra.copper.jac

jac_kernel_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.ra.copper.jac))

# Covariates (date_num)
datenum_covariate_all <- phyloseq.bacteria.samples.dates.copper.all.df%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_all

#Id for repeated measures 
Enclosure_d_mirkat_all <- phyloseq.bacteria.samples.dates.copper.all.df %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_mirkat_all

#Run model - dont think it works. just two enclosures and enclosure id is the same as the outcome 
mirkat_enclosure_glm <- GLMMMiRKAT(
  y = phyloseq.bacteria.samples.dates.copper.all.df$Copper_level_mg_L,
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
copper_cont_BC_beta_div <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.scrs) + 
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
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.5, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.7, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = 0.2, y = -1.1, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
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
copper_cont_BC_beta_div_2 <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.scrs) + 
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
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.5, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.7, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
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
phyloseq.bacteria.samples.dates.ra.dates.bray.scrs_2 <- phyloseq.bacteria.samples.dates.ra.dates.bray.scrs %>%
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
copper_cont_BC_beta_div_3 <- ggplot(phyloseq.bacteria.samples.dates.ra.dates.bray.scrs_2) + 
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
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs_2 %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs_2 %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.8, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs_2 %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.ra.dates.bray.scrs_2 %>%
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
phyloseq.bacteria.samples.dates.nitifiers.ra <- subset_taxa(phyloseq.bacteria.samples.dates.ra.dates, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
phyloseq.bacteria.samples.dates.nitifiers.ra #827 taxa and 218 samples
phyloseq.bacteria.samples.dates.nitifiers.ra <- subset_samples(phyloseq.bacteria.samples.dates.nitifiers.ra, sample_sums(phyloseq.bacteria.samples.dates.nitifiers.ra) > 0)
phyloseq.bacteria.samples.dates.nitifiers.ra #827 taxa and 218 samples 

##BC 
phyloseq.bacteria.samples.dates.nitifiers.ra.bray <- vegdist(t(phyloseq.bacteria.samples.dates.nitifiers.ra@otu_table), method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.dates.nitifiers.df <- as(phyloseq.bacteria.samples.dates.nitifiers.ra@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.nitifiers.df 

###ORDINATION####
set.seed(98)
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.ord <- metaMDS(phyloseq.bacteria.samples.dates.nitifiers.ra.bray, 
                                                 k=3, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
###Centroids by enclosure #####
## BC
#Simple ordination plot
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.plot <- ordiplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs <- scores(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.plot, display = "sites")
#Add metadata to coordinates
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs <- cbind(as.data.frame(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs),
                                                Copper_level_mg_L = phyloseq.bacteria.samples.dates.nitifiers.df$Copper_level_mg_L, 
                                                SampleID = phyloseq.bacteria.samples.dates.nitifiers.df$SampleID, 
                                                Collection_Date = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Date, 
                                                Enclosure = phyloseq.bacteria.samples.dates.nitifiers.df$Enclosure, 
                                                Attempt = phyloseq.bacteria.samples.dates.nitifiers.df$Attempt, 
                                                Copper_quartile = phyloseq.bacteria.samples.dates.nitifiers.df$Copper_quartile,
                                                Copper_quartile.abvr = phyloseq.bacteria.samples.dates.nitifiers.df$Copper_quartile.abvr)%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.cent.enclosure <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
                                                              data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs, 
                                                              FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs <- merge(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs, 
                                                          setNames(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.cent.enclosure, 
                                                                   c("Enclosure","cMDS1","cMDS2")),
                                                          by = 'Enclosure', 
                                                          sort = F)

####MIRKAT- Just enclosure#####
#Dataframe with metadata
df_metadata_nit_mirkat <- data.frame(sample_data(phyloseq.bacteria.samples.dates.nitifiers.ra))
#BC distance matrix row order
labs_nit <- rownames(as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.bray))  
#Align metadata df to distance order
df_nit_mirkat <- df_metadata_nit_mirkat[labs_nit, , drop = FALSE] 

#Going to get kernels out of different distance matrices
#Bray curtis 
phyloseq.bacteria.samples.dates.ra.dates.bray

#Jaccard distances
phyloseq.bacteria.samples.dates.ra.dates.jac <- vegdist(t(phyloseq.bacteria.samples.dates.ra.dates@otu_table), method = "jaccard")
phyloseq.bacteria.samples.dates.ra.dates.jac

#Aitchison 
##Going to take out samples from P1 from april and may 2023 and september from H21 (raw counts object)
phyloseq.bacteria.samples.dates.dates <- subset_samples(phyloseq.bacteria.samples.dates, 
                                                     Collection_Date > "2023-10-01")
phyloseq.bacteria.samples.dates.dates #33490 taxa and 218 samples
#Take out taxa sums = 0
phyloseq.bacteria.samples.dates.dates  <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.dates) > 0, 
                                               phyloseq.bacteria.samples.dates.dates)
phyloseq.bacteria.samples.dates.dates  #33490 taxa and 218 samples
#Calculate distance
phyloseq.bacteria.samples.dates.dates_clr <- microbiome::transform(phyloseq.bacteria.samples.dates.dates, "clr") #convert raw counts to clr
phyloseq.bacteria.samples.dates.dates_clr_dist_aitch <- dist(t(otu_table(phyloseq.bacteria.samples.dates.dates_clr)), method = "euclidean") #calculate euclidean distances

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
#Bray curtis
bray_kernel_nit <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.bray))
#Jaccards
jac_kernel_nit <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.ra.dates.jac))
# Aitchison
aitch_kernel_nit <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.dates_clr_dist_aitch))

#Modelling enclosure, have to make it dichotomus (used to be phyloseq.bacteria.samples.dates.nitifiers.df)
Enclosure_d_nit_mirkat <- df_nit_mirkat %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_nit_mirkat

# # Covariates (month)
# month_covariate_nit <- model.matrix(~ Collection_Month, 
#                                 data = phyloseq.bacteria.samples.dates.nitifiers.df)[, -1]
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
enclosure_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.bray ~ Enclosure, 
                                phyloseq.bacteria.samples.dates.nitifiers.df, 
                                strata = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
enclosure_BC_nit_adonis #21% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.dates.nitifiers.df, 
#                                 strata = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_nit_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.dates.nitifiers.df , 
#                                 strata = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_nit_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS
#ENCLOSURE#
# Run the betadisper function, average distance to centroid
bray.enclosure.nit.disp <- betadisper(phyloseq.bacteria.samples.dates.nitifiers.ra.bray, 
                                  phyloseq.bacteria.samples.dates.nitifiers.df$Enclosure)
bray.enclosure.nit.disp
##Then test by permuting
set.seed(98)
bray.enclosure.nit.permdisp <- permutest(bray.enclosure.nit.disp, permutations = 9999)
bray.enclosure.nit.permdisp ##S, p = 6e-04

# #COPPER#
# # Run the betadisper function, average distance to centroid
# bray.copper.disp <- betadisper(phyloseq.bacteria.samples.dates.nitifiers.ra.bray, 
#                                   phyloseq.bacteria.samples.dates.nitifiers.df$Copper_level_mg_L)
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
enclosure_BC_nit_beta_div <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs) + theme_bw() +
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
enclosure_BC_nit_beta_div_2 <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs) + 
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
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs_2 <- phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs%>%
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
enclosure_BC_nit_beta_div_3 <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs_2) + 
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
enclosure_BC_nit_beta_div_4 <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs_2) + 
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
  geom_point(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs_2%>%
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
  geom_point(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.enclosure.segs_2 %>%
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
# fit <- envfit(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.ord, 
#               phyloseq.bacteria.samples.dates.nitifiers.df$Copper_level_mg_L, 
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
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.cent.copper_quart <- aggregate(
  cbind(MDS1, MDS2) ~ Enclosure + Copper_quartile,
  data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs,
  FUN = mean)

#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.copper_quart.segs <- merge(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs, 
                                                             setNames(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.cent.copper_quart, 
                                                                      c("Enclosure", "Copper_quartile",
                                                                        "cMDS1","cMDS2")),
                                                             by = c("Enclosure", "Copper_quartile"), 
                                                             sort = F)
##PERMANOVA - Just copper###
# set.seed(98)
# copper_quart_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.bray ~ Copper_quartile, 
#                                    phyloseq.bacteria.samples.dates.nitifiers.df, 
#                                    strata = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Month,
#                                    by = "margin",
#                                    permutations = 9999)
# copper_quart_BC_nit_adonis #30.1% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.bacteria.samples.dates.nitifiers.df, 
#                                 strata = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_nit_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_nit_adonis  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.bacteria.samples.dates.nitifiers.df , 
#                                 strata = phyloseq.bacteria.samples.dates.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_nit_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS#
#COPPER quartile#
# Run the betadisper function, average distance to centroid
bray.nit.copper.quart.disp <- betadisper(phyloseq.bacteria.samples.dates.nitifiers.ra.bray,
                                     phyloseq.bacteria.samples.dates.nitifiers.df$Copper_quartile)
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
copper_quart_BC_nit_beta_div <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.copper_quart.segs) + 
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
phyloseq.bacteria.samples.dates.nitifiers.ra.copper <- subset_samples(phyloseq.bacteria.samples.dates.nitifiers.ra, 
                                                                !is.na(Copper_level_mg_L))
phyloseq.bacteria.samples.dates.nitifiers.ra.copper <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.nitifiers.ra.copper)> 0, 
                                                            phyloseq.bacteria.samples.dates.nitifiers.ra.copper)
phyloseq.bacteria.samples.dates.nitifiers.ra.copper #827 taxa and 204 samples

####H21 #######
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21 <- subset_samples(phyloseq.bacteria.samples.dates.nitifiers.ra.copper, 
                                                                Enclosure == "H21")
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21)> 0, 
                                                                phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21)
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21 # 824 taxa and 82 samples 

##BC 
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.bray <- vegdist(t(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21@otu_table), 
                                                              method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df <- as(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_H21 <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.bray))

#Jaccard distances
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.jac <- vegdist(t(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21@otu_table), 
                                                                 method = "jaccard")
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.jac

jac_kernel_nit_H21 <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.jac))

# # Covariates (month)
# month_covariate_nit_H21 <- model.matrix(~ Collection_Month, 
#                                     data = phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df)[, -1]
# 
# month_covariate_nit_H21

# Covariates (date_num)
datenum_covariate_H21_nit <- phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_H21_nit

#Collection month for random intercepts
collection_month_H21_nit <- phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df%>%
  pull(Collection_Month)%>%
  as.matrix()
collection_month_H21_nit

#Run model
mirkat_nit_H21_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df$Copper_level_mg_L,
  X = datenum_covariate_H21_nit,
  Ks = bray_kernel_nit_H21,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_nit_H21_copper_cont #R2 3.1% , p_values = 0.004

#####PERMANOVA #######
set.seed(98)
adonis_nit_H21_copper_cont  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.bray ~ Copper_level_mg_L, 
                                       phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df, 
                                       strata = phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_nit_H21_copper_cont # R2 9.7%, p value 0.027*

##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.nit.copper.cont.disp.H21 <- betadisper(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.H21.bray,
                                           phyloseq.bacteria.samples.dates.nitifiers.copper.H21.df$Copper_quartile)
bray.nit.copper.cont.disp.H21
##Then test by permuting
set.seed(98)
bray.nit.copper.cont.permdisp.H21 <- permutest(bray.nit.copper.cont.disp.H21, permutations = 9999)
bray.nit.copper.cont.permdisp.H21 ##NS, p = 0.65

####P1 #######
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1 <- subset_samples(phyloseq.bacteria.samples.dates.nitifiers.ra.copper, 
                                                                    Enclosure == "P1")
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1 <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1)> 0, 
                                                                phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1)
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1 #624 taxa and 122 samples

##BC 
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1.bray <- vegdist(t(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1@otu_table), 
                                                                  method = "bray")

#make DF from metadata
phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df <- as(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_P1 <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1.bray))

# # Covariates (month)
# month_covariate_nit_P1 <- model.matrix(~ Collection_Month, 
#                                         data = phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df)[, -1]
# 
# month_covariate_nit_P1

# Covariates (date_num)
datenum_covariate_P1_nit <- phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_P1_nit

#Run model
mirkat_nit_P1_copper_cont <- MiRKAT(
  y = phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df$Copper_level_mg_L,
  X = datenum_covariate_P1_nit,
  Ks = bray_kernel_nit_P1,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_nit_P1_copper_cont #p val 0 , R2 20.6%

#####PERMANOVA #####
set.seed(98)
adonis_nit_P1_copper_cont  <- adonis2(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1.bray ~ Copper_level_mg_L, 
                                       phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df, 
                                       strata = phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_nit_P1_copper_cont # R2 26.9%, p value 1e-04*
##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.nit.copper.cont.disp.P1 <- betadisper(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.P1.bray,
                                           phyloseq.bacteria.samples.dates.nitifiers.copper.P1.df$Copper_quartile)
bray.nit.copper.cont.disp.P1
##Then test by permuting
set.seed(98)
bray.nit.copper.cont.permdisp.P1 <- permutest(bray.nit.copper.cont.disp.P1, permutations = 9999)
bray.nit.copper.cont.permdisp.P1 ##S, p = 1e-04 

###MIRKAT GLMM FOR COPPER #####
#make DF from metadata
phyloseq.bacteria.samples.dates.nitifiers.copper.all.df <- as(phyloseq.bacteria.samples.dates.nitifiers.ra.copper@sam_data, "data.frame") %>%
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
phyloseq.bacteria.samples.dates.nitifiers.copper.all.df

#Bray-Curtis distances
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.bray <- vegdist(t(phyloseq.bacteria.samples.dates.nitifiers.ra.copper@otu_table), 
                                                              method = "bray")
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.bray

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.bray))

#Jaccard distances
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.jac <- vegdist(t(phyloseq.bacteria.samples.dates.nitifiers.ra.copper@otu_table), 
                                                             method = "jaccard")
phyloseq.bacteria.samples.dates.nitifiers.ra.copper.jac

jac_kernel_nit_all <- D2K(
  as.matrix(phyloseq.bacteria.samples.dates.nitifiers.ra.copper.jac))

# Covariates (date_num)
datenum_covariate_nit_all <- phyloseq.bacteria.samples.dates.nitifiers.copper.all.df%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_nit_all

#Id for repeated measures 
Enclosure_d_nit_mirkat_all <- phyloseq.bacteria.samples.dates.nitifiers.copper.all.df %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_nit_mirkat_all

#Run model - dont think it works. just two enclosures and enclosure id is the same as the outcome 
mirkat_enclosure_nit_glm <- GLMMMiRKAT(
  y = phyloseq.bacteria.samples.dates.nitifiers.copper.all.df$Copper_level_mg_L,
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
copper_cont_BC_nit_beta_div <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs) + 
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
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
    aes(x = -1, y = 0.6, label = "PERMANOVA\nCopper levels (mg/L)"),
    hjust = 0.5, vjust = -0.5,
    size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.4, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0.13, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
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
copper_cont_BC_nit_beta_div_2 <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs) + 
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
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = -0.65, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = -0.88, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = -0.2, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
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
phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs_2 <- phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
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
copper_cont_BC_nit_beta_div_3 <- ggplot(phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs_2) + 
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
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.3, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0.1, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.bacteria.samples.dates.nitifiers.ra.bray.scrs %>%
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
#Import count tables#######
ARGcounts <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplus_dev_branch/SNPconfirmed_AMR_analytic_matrix.csv')
colnames(ARGcounts)
rownames(ARGcounts)

#Annotations - downloaded from AMRplusplus -REMEMBER TO UPDATE THIS ONCE YOU GET ACTUAL FILES#######
tax.table.ARG <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplus_dev_branch/Annotations/megares_annotations_v3.00.csv')

#Edit to make compatible with phyloseq
tax.table.ARG <- tax.table.ARG%>%
  column_to_rownames(var = "header")%>% #make "header' rownames so it matches with the OTU table rownames
  rename_with(~ str_to_title(.))%>% #want the annotation cols with the first letter capitalized
  as.matrix() ##make into matrix so it is compatible with tax_table function from phyloseq

#OTU table ########
otu_table_ARG <- ARGcounts%>%
  column_to_rownames(var = "gene_accession")

#PHYLOSEQ####
OTU_ARG <-phyloseq::otu_table(otu_table_ARG, 
                              taxa_are_rows = TRUE)
TAX_ARG <-phyloseq::tax_table(tax.table.ARG)
phyloseq_ARG <- phyloseq(OTU_ARG, TAX_ARG, sampledata_phyloseq) 
phyloseq_ARG ##8733 taxa and 226 samples

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU_ARG), metadata$SampleID) #Yes, "H21_1021a" and "H21_1021b"

#PREPROCESSING ####
phyloseq_ARG #8733 taxa and 226 samples
min(sample_sums(phyloseq_ARG)) # 6799 (H21_0120)
max(sample_sums(phyloseq_ARG)) # 2,001,381  (ZymoMock1a_S142) 
mean(sample_sums(phyloseq_ARG)) #277,444.7
median(sample_sums(phyloseq_ARG)) #227,141
sort(sample_sums(phyloseq_ARG))

##Zymo#####
### pulling out samples from ZYMOs and EB, NTC and those samples with low OTUs
phyloseq_ARG.controls <- subset_samples(phyloseq_ARG, 
                                        grepl("NTC|EB|Zymo", sample_names(phyloseq_ARG)))
phyloseq_ARG.controls <- prune_taxa(taxa_sums(phyloseq_ARG.controls) > 0, phyloseq_ARG.controls) 
phyloseq_ARG.controls #4058 taxa and 3 samples (only zymos)

##Samples#####
##New phyloseq of just samples
phyloseq_ARG.samples <- subset_samples(phyloseq_ARG, 
                                       !grepl("NTC|EB|Zymo", sample_names(phyloseq_ARG)))
phyloseq_ARG.samples #8733 taxa and 223 samples 
#Taking out those with low counts:
sort(sample_sums(phyloseq_ARG.samples)) #H21_0120 has low counts, droppin git here as i did with kraken taxonomy
phyloseq_ARG.samples <- prune_samples(sample_sums(phyloseq_ARG.samples) > 10000, phyloseq_ARG.samples) 
phyloseq_ARG.samples <- prune_taxa(taxa_sums(phyloseq_ARG.samples) > 0, phyloseq_ARG.samples)
phyloseq_ARG.samples #8335 taxa and 222 samples (dropped H21_0120)
sort(sample_sums(phyloseq_ARG.samples)) #OK

###Dropping samples before copper dosage started####
#What's the range of dates 
range(phyloseq_ARG.samples@sam_data$Collection_Date)#"2023-04-20" "2024-04-30"

#Actual Copper dosing starts from 10/09/2023 (sampling ends on 03/02/2024) for naive system
#Actual copper dosing starts from 11/14/2023 (sampling ends on 04/30/2024) for established system
phyloseq_ARG.samples.dates <- subset_samples(phyloseq_ARG.samples, Collection_Date > "2023-10-05")
phyloseq_ARG.samples.dates
phyloseq_ARG.samples.dates <- prune_taxa(taxa_sums(phyloseq_ARG.samples.dates) > 0, phyloseq_ARG.samples.dates) 
phyloseq_ARG.samples.dates #8334 taxa and 217 samples
setdiff(sample_names(phyloseq_ARG.samples), sample_names(phyloseq_ARG.samples.dates)) #Dropped "H21_0912" "H21_1005" "P1_0420"  "P1_0427"  "P1_0504" 

#Also, have H21_1202a and 1202b. These were taken on Dec 2nd, 2023. A is before backwash, B is afterbackwash. 
#Have more reliable metadata for H21_1202a (before backwash)
phyloseq_ARG.samples.dates <- subset_samples(phyloseq_ARG.samples.dates, 
                                             SampleID != "H21_1202b")
phyloseq_ARG.samples.dates #8334 taxa and 216 samples (Dopped H21_1202b)
phyloseq_ARG.samples.dates <- prune_taxa(taxa_sums(phyloseq_ARG.samples.dates) > 0, phyloseq_ARG.samples.dates) 
phyloseq_ARG.samples.dates #8334 taxa and 216 samples  

###H21####
phyloseq_ARG.samples.dates_H21 <- subset_samples(phyloseq_ARG.samples.dates, Enclosure == "H21")
phyloseq_ARG.samples.dates_H21 <- prune_taxa(taxa_sums(phyloseq_ARG.samples.dates_H21) > 0, 
                                             phyloseq_ARG.samples.dates_H21)
phyloseq_ARG.samples.dates_H21 #8305 taxa and 92 samples
range(phyloseq_ARG.samples.dates_H21@sam_data$Collection_Date)#OK, now "2023-10-09" through "2024-03-02"

###P1####
phyloseq_ARG.samples.dates_P1 <- subset_samples(phyloseq_ARG.samples.dates, Enclosure == "P1")
phyloseq_ARG.samples.dates_P1 <- prune_taxa(taxa_sums(phyloseq_ARG.samples.dates_P1) > 0, 
                                            phyloseq_ARG.samples.dates_P1)
phyloseq_ARG.samples.dates_P1 #8270 taxa and 124 samples 
range(phyloseq_ARG.samples.dates_P1@sam_data$Collection_Date)#OK, now "2023-11-14" through "2024-04-30"

#Am I missing samples that have a kraken taxonomic classification?
setdiff(sample_names(phyloseq_ARG.samples.dates), sample_names(phyloseq.bacteria.samples.dates)) #Nope

#Group level#### 
phyloseq_ARG.samples.dates.group <- tax_glom(phyloseq_ARG.samples.dates, taxrank = "Group", NArm = F)
phyloseq_ARG.samples.dates.group #1373 groups and 216 samples

###TAX GLOMMING - RAW COUNTS##### 
phyloseq_ARG.samples.dates.type <- tax_glom(phyloseq_ARG.samples.dates.group, taxrank = "Type", NArm = F) 
phyloseq_ARG.samples.dates.type #5 types (216 samples)

phyloseq_ARG.samples.dates.class <- tax_glom(phyloseq_ARG.samples.dates.group, taxrank = "Class", NArm = F) # classes
phyloseq_ARG.samples.dates.class #58 classes (216 samples)

phyloseq_ARG.samples.dates.mechanism <- tax_glom(phyloseq_ARG.samples.dates.group, taxrank = "Mechanism", NArm = F) 
phyloseq_ARG.samples.dates.mechanism #97 mechanisms (217 samples)




# enclosure_BC_beta_div_copper