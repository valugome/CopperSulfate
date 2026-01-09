# setwd,load libraries, source functions ####
setwd('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Analysis/Bacteria_archaea/Figures')

# install.packages("devtools")
# devtools::install_github("vmikk/metagMisc")
# # BiocManager::install("phyloseq")
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install.packages("svglite")
library(phyloseq); library (tidyverse); library(ggplot2);  library(stringr); 
library(dplyr);library(metagMisc); library(metagenomeSeq); library(vegan); library(cowplot);
library(ggdendro); library(pairwiseAdonis); library(randomcoloR); library(ggpubr); library(ppcor)
library(ggsignif); library (ANCOMBC);library(maaslin3); library (UpSetR); library(MicrobiotaProcess); library(microbiome)
library(ggtext); library(ggnewscale); library(rstatix); library(ggrepel); library(ggh4x); library(svglite);
library(lmerTest); library(mgcv); library(rmcorr); library("emmeans"); library(patchwork)

##Source functions
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbundanceOthersPercentage.R')
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/top_taxa_legend_updated.R")
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/fill_taxonomy_updated.R')
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbun_group_microbiome.R')


#Importing data from kraken output nt_core - counts will be classified reads#### 
counts <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/Paired_end_mode_updated_20251005/kraken_analytic_matrix.conf_0.0_cuso4.csv')

##Unclassified counts###
unclassified_counts <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/Paired_end_mode_updated_20251005/unclassifieds_kraken_analytic_matrix.conf_0.0_cuso4.csv')

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

###OTUs
##OTU table 
otu_table <- counts[, -1]%>% #Excludes the first column (taxonomy)
  mutate(OTU = paste0("OTU", 1:nrow(counts))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##make into matrix so it is compatible with otu_table function from phyloseq
otu_table

#IMPORT METADATA####
##Metadata####
#P1 enclosure (established)
metadata_P1 <- readr::read_csv("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/P1_Schooling_Copper_Treatment_Water_Sample_Metadata_clean.csv")
#H21 enclosure (naive)
metadata_H21 <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/H21_Copper_Treatment_Water_Sample_Metadata_clean.csv')

#P1 enclosure water quality metrics
water_quality_P1 <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/P1_copper_WQ_clean.csv')

#H21 enclosure water quality metrics 
water_quality_H21 <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/H21_WQ_090123_030824_clean.csv')

#H21 enclosure treatment metrics
treatment_dates_H21 <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/H21_Copper_Tx_12.2023_clean.csv')

##Clean up####
###P1 metadata####
#Data types 
str(metadata_P1) #OK. Numerics are ok. 

##Replace spaces, commas, () with "_" ([] characters inside these brackets will be replaced)
colnames(metadata_P1) <- stringr::str_replace_all(
  colnames(metadata_P1),
  "[ ,()/]", "_")%>%
  str_replace_all(., "__", "_")%>% ##replace double "__" with "_"
  str_replace_all(., "_$", "") ##Discard "_" at the end 
colnames(metadata_P1) #Ok

#Date check
str(metadata_P1) #OK. Only want to change Collection_Date from "chr' to date 
metadata_P1$Collection_Date 
metadata_P1$Collection_Date <- as.Date(metadata_P1$Collection_Date, format = "%m/%d/%y")
str(metadata_P1$Collection_Date) #Ok now 

##Add enclosure as column 
metadata_P1$Enclosure <- "P1"
#Attempt as factor 
metadata_P1$Attempt <- factor(metadata_P1$Attempt)

####P1 Water quality#####
str(water_quality_P1) #All characters. Need to fix those that are int or floats

##Fixing colnames 
colnames(water_quality_P1)
#New col names will have units
newcolnames_P1_WQ <- c("Request_Date", "Enclosure", 
                       "Temperature_F", "Chlorine_mg_L", 
                       "pH_spu", "Ammonia_mg_L", "Nitrite_mg_L", 
                       "Nitrate_UV_mg_L", "Salinity_ppt",
                       "Alkalinity_mg_L", "Calcium_mg_L",
                       "Phosphate_mg_L", "Copper_mg_L", "Duplicate_date_metadata", 
                       "Attempt")
colnames(water_quality_P1) <- newcolnames_P1_WQ
colnames(water_quality_P1) #OK

##Discard units from data
water_quality_P1 <- water_quality_P1 %>%
  mutate(across(!c(Request_Date, Enclosure, Duplicate_date_metadata, Attempt),
                ~ parse_number(.)))
str(water_quality_P1) #OK now

#Fixing date
str(water_quality_P1) #OK. Only want to change Collection_Date from "chr' to date 
water_quality_P1$Request_Date 
water_quality_P1$Request_Date <- as.Date(water_quality_P1$Request_Date, format = "%m/%d/%y")
str(water_quality_P1$Request_Date) #Ok now 

#Duplicate_date and Attempt as factors
water_quality_P1$Duplicate_date_metadata <- factor(water_quality_P1$Duplicate_date_metadata, levels = c("0", "1"))
water_quality_P1$Attempt <- factor(water_quality_P1$Attempt)
str(water_quality_P1$Duplicate_date_metadata) #Ok now 

##Add enclosure as column
water_quality_P1$Enclosure <- "P1"

##Fix duplicates (I checked metadata and then added "1" to those dates that matched the data on the metadata spreadsheet)
water_quality_P1_collapsed <- water_quality_P1 %>%
  filter(Duplicate_date_metadata == "1")
water_quality_P1_collapsed
any(duplicated(water_quality_P1_collapsed$Request_Date)) #Ok, no duplicated dates now

###H21 metadata#####
colnames(metadata_H21)
##Replace spaces, commas, (), / with "_" ([] characters inside these brackets will be replaced)
colnames(metadata_H21) <- stringr::str_replace_all(
  colnames(metadata_H21),
  "[ ,()/]", "_")%>%
  str_replace_all(., "__", "_")%>% ##replace double "__" with "_"
  str_replace_all(., "_$", "") ##Discard "_" at the end 
colnames(metadata_H21) #Ok

#Data types
str(metadata_H21) ##Copper_addition_mL and Copper_level_mg_Lhave to be num, Attempt, Backwash, Post_backwash, Water_change, Post_waterchange, Other_antiparasitic have to be chr
#Checking Copper addition
unique(metadata_H21$Copper_addition_mL) #has a weird "?". Will replace with NA
unique(metadata_H21$Copper_level_mg_L) #There is a weird "na", Will have to replace with actual NA

#Making changes: Copper_addition_mL and Copper_level_mg_Lhave to be num
metadata_H21 <- metadata_H21 %>%
  mutate(
    Copper_addition_mL = if_else(
      Copper_addition_mL == "?",
      NA_real_,
      parse_number(Copper_addition_mL)),
    Copper_level_mg_L = if_else(
      Copper_level_mg_L == "na",
      NA_real_,
      parse_number(Copper_level_mg_L)))
str(metadata_H21) #Ok, only need to fix date now

#Making changes: Attempt, Backwash, Post_backwash, Water_change, Post_waterchange, Other_antiparasitic have to be chr
metadata_H21 <- metadata_H21 %>%
  mutate(across(c("Attempt", "Backwash", "Post_backwash", "Water_change", "Post_waterchange", "Other_antiparasitic"),
                ~ factor(.)))
str(metadata_H21) #Ok, only need to fix date now

#Date check
str(metadata_H21) #OK. Only want to change Collection_Date from "chr' to date 
metadata_H21$Collection_Date 
metadata_H21$Collection_Date <- as.Date(metadata_H21$Collection_Date, format = "%m/%d/%y")
str(metadata_H21$Collection_Date) #Ok now 

##Add enclosure as column 
metadata_H21$Enclosure <- "H21"

####H21 Water quality#####
str(water_quality_H21) #All characters. Need to fix those that are int or floats
##Fixing colnames 
colnames(water_quality_H21)
#New col names will have units
newcolnames_H21_WQ <- c("Request_Date", "Enclosure", 
                       "Temperature_F", "pH_spu", 
                       "Ammonia_mg_L", "Nitrite_mg_L", 
                       "Nitrate_UV_mg_L", "Salinity_ppt",
                       "Phosphate_mg_L", "Copper_mg_L",
                       "Nitrite_IC_ppm", "Duplicate_date_metadata",
                       "WQ_comments")
colnames(water_quality_H21) <- newcolnames_H21_WQ
colnames(water_quality_H21) #OK

#Discard units from data
water_quality_H21 <- water_quality_H21 %>%
  mutate(across(!c(Request_Date, Enclosure,Duplicate_date_metadata, WQ_comments),
                ~ parse_number(.)))
str(water_quality_H21) #OK now

#Fixing date
str(water_quality_H21) #OK. Only want to change Collection_Date from "chr' to date 
water_quality_H21$Request_Date 
water_quality_H21$Request_Date <- as.Date(water_quality_H21$Request_Date, format = "%m/%d/%y")
str(water_quality_H21$Request_Date) #Ok now

##Add enclosure as column 
water_quality_H21$Enclosure <- "H21"

##Fix duplicates (I checked metadata and then added "1" to those dates that matched the data on the metadata spreadsheet)
water_quality_H21_collapsed <- water_quality_H21 %>%
  filter(Duplicate_date_metadata == "1")
water_quality_H21_collapsed
water_quality_H21_collapsed$Duplicate_date_metadata <- factor(water_quality_H21_collapsed$Duplicate_date_metadata)
any(base::duplicated(water_quality_H21_collapsed$Request_Date)) #Ok, no duplicated dates now

##Treatment metadata ####
str(treatment_dates_H21) #Ok
colnames(treatment_dates_H21)
#New col names will have units
newcolnames_H21_treatment <- c("Treatment_Date", "Treatment_Enclosure", 
                        "Treatment_Copper_reading_mg_L", "Treatment_Copper_target_mg_L", 
                        "Treatment_Copper_addition_mL", "Treatment_ammonia_reading_mg_L",
                        "Treatment_Attempt", "Backwash",
                        "Post_backwash", "Water_change", 
                        "Post_waterchange","Water_change_percent",
                        "Other_antiparasitic",
                        "Treatment_Comments")
colnames(treatment_dates_H21) <- newcolnames_H21_treatment
colnames(treatment_dates_H21) #OK

##Discard units from data
str(treatment_dates_H21) #ok, only Treatment_Copper_addition_mL and Treatment_ammonia_reading_mg_L needs to be fixed to be a num, and Attempt to a factor 

#Treatment_Copper_addition_mL and needs to be fixed, remove the unit "mL", and make Attempt and others a factor
treatment_dates_H21 <- treatment_dates_H21 %>%
  mutate(Treatment_Copper_addition_mL = parse_number(Treatment_Copper_addition_mL),
         Treatment_ammonia_reading_mg_L = as.numeric(Treatment_ammonia_reading_mg_L),
         across(c("Treatment_Attempt", "Backwash", "Post_backwash", "Water_change", "Post_waterchange", "Other_antiparasitic"),
                ~ factor(.)))
str(treatment_dates_H21) #OK now, only need to fix date now 

#Fixing date
str(treatment_dates_H21) #OK. Only want to change Collection_Date from "chr' to date 
treatment_dates_H21$Treatment_Date  
treatment_dates_H21$Treatment_Date  <- as.Date(treatment_dates_H21$Treatment_Date , format = "%m/%d/%y")
str(treatment_dates_H21$Treatment_Date)  #Ok now 

##Add enclosure as column 
treatment_dates_H21$Enclosure <- "H21"

##Checking for duplicates 
treatment_dates_H21 %>%
  count(Treatment_Date) %>%
  filter(n > 1) ##There are 8. These are usually a pre and post backwash and/or water change. That's why I added columns to identify them when merging metadata. 


##Joining metadata info#####
##Metadata and water quality
metadata_join_P1 <- left_join(metadata_P1, water_quality_P1_collapsed, 
                              by = c("Collection_Date" = "Request_Date", "Enclosure", "Attempt"))

##For h21 (naive) a bit more complicated treatment, have to add also treatment metadata
metadata_join_H21 <- left_join(metadata_H21, water_quality_H21_collapsed, 
                               by = c("Collection_Date" = "Request_Date", 
                                      "Enclosure"))%>%
  left_join(treatment_dates_H21, by = c("Collection_Date" = "Treatment_Date", 
                                        "Backwash","Post_backwash","Water_change",
                                        "Post_waterchange", 
                                        "Enclosure"))
str(metadata_join_H21)

#Handle duplicate columns
metadata_join_H21 <- metadata_join_H21 %>%
  mutate(Water_change_percent = coalesce(Water_change_percent.x,
                               Water_change_percent.y), 
         Other_antiparasitic.y = coalesce(Other_antiparasitic.x, 
                                          Other_antiparasitic.y), 
         Ammonia_mg_L = coalesce (Ammonia_mg_L, NH3_mg_L, Treatment_ammonia_reading_mg_L), 
         Copper_level_mg_L = coalesce(Copper_level_mg_L, Copper_mg_L, Treatment_Copper_reading_mg_L), 
         Copper_addition_mL = coalesce(Copper_addition_mL,Treatment_Copper_addition_mL),
         Attempt = coalesce(Attempt, Treatment_Attempt)) %>%
  dplyr::select(!c("Water_change_percent.x", 
            "Water_change_percent.y", 
            "Other_antiparasitic.y", 
            "Other_antiparasitic.x",
            "Treatment_ammonia_reading_mg_L", 
            "Copper_mg_L",
            "Treatment_Copper_reading_mg_L",
            "Treatment_Attempt", 
            ))
str(metadata_join_H21) #ok, good now. 

##Finally, add Zymos and negative controls
controls_and_zymo <- data.frame(SampleID = c("EB2_S186",
                                             "NTC10_S231",
                                             "NTC11_S155",
                                             "NTC1_S163",
                                             "NTC2_S7",
                                             "NTC3_S180",
                                             "NTC4_S22",
                                             "NTC5_S112",
                                             "NTC6_S42",
                                             "NTC7_S129",
                                             "NTC8_S60",
                                             "NTC9_S141",
                                             "ZymoMock1a_S142",
                                             "ZymoMock1_S238",
                                             "ZymoMock2_S78"))

#Final metadata 
metadata <- bind_rows(metadata_join_P1, 
                      metadata_join_H21, 
                      controls_and_zymo)

##Host free reads####
#Df obtained from running ('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Host_removal_stats/read_counts/Settingup_stats_HOSTREM.R')
#hostrem <- read.csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/Host_removal_stats/read_counts/HostRem_stats_reads.csv')

#Merge with metadata
# hostrem <- hostrem %>%
#   select(SampleID, Hostrem_output_total_num_seqs)%>%
#   left_join(metadata, by = "SampleID")%>%
#   rename(HostFree_Reads=Hostrem_output_total_num_seqs)#Hostrem_output_total_num_seqs is the total number of host free reads

##Making into phyloseq-compatible object
sampledata_phyloseq <- metadata %>%
  mutate(rows = SampleID)%>%
  column_to_rownames(var= "rows") %>%##Make sampleID column into row names, so they match sample_names() with OTU and TAX
  sample_data(metadata) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

#PHYLOSEQ####
##Sample H21_1102 I redid (H21_1102_re). Keeping whichever one has the highest counts
colSums(otu_table)
##H21_1102 has more counts. Keeping that one 
otu_table_filt <- otu_table%>%
  data.frame()%>%
  dplyr::select(-H21_1102_re)%>%
  as.matrix()

#Make phyloseq object
OTU <-phyloseq::otu_table(otu_table_filt, taxa_are_rows = TRUE)
TAX <-phyloseq::tax_table(filled_taxonomy_2)
phyloseq <- phyloseq(OTU, TAX, sampledata_phyloseq)

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU), metadata$SampleID) #Yes, "H21_1021a" and "H21_1021b"

#Are there samples in metadata that don't have sequencing data?
setdiff(metadata$SampleID, sample_names(OTU)) #Yes, "P1_1126", "P1_1203", 
#"P1_1216", "P1_1225", "P1_1228", "P1_0104", "P1_0112", "P1_0115", 
#"P1_0205", "P1_0212", "P1_0218", "H21_1021", "H21_1122b"

##H21####
phyloseq_H21 <- subset_samples(phyloseq, Enclosure == "H21")
phyloseq_H21 <- prune_taxa(taxa_sums(phyloseq_H21) > 0, phyloseq_H21)
phyloseq_H21 #187717 taxa and 97 samples

##P1####
phyloseq_P1 <- subset_samples(phyloseq, Enclosure == "P1")
phyloseq_P1 <- prune_taxa(taxa_sums(phyloseq_P1) > 0, phyloseq_P1)
phyloseq_P1 #126208 taxa and 128 samples  

#Color Palettes#####
enclosure.palette <- c("H21" = "#fc8d62",  
                       "P1"  = "#8da0cb" )
attempt.palette <- c("1" = "#0072B2", 
                     "2" = "#E69F00",
                     "3" = "#009E73")

#PREPROCESSING ####
phyloseq #187,717 taxa and 240 samples 
      
### Selecting only Bacteria/Archaea
phyloseq.bacteria <- subset_taxa(phyloseq, Domain=="Archaea" | Domain=="Bacteria")
phyloseq.bacteria #33613 taxa, 240 samples

##Selecting only viruses
phyloseq.viruses <- subset_taxa(phyloseq, Domain=="Viruses")
taxanames_viruses <- c("Kingdom", "Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species") ##they have a different classification system, updating it here
colnames(phyloseq.viruses@tax_table) <- taxanames_viruses #replacing col names of the tax_table for new ones
colnames(phyloseq.viruses@tax_table) #OK taxonomy ranks
phyloseq.viruses #48729 taxa and 240 samples

##Selecting only eukaryota 
phyloseq.eukaryota <- subset_taxa(phyloseq, Domain=="Eukaryota")
colnames(phyloseq.eukaryota@tax_table) ##These are OK taxonomy ranks
phyloseq.eukaryota #105,375 taxa and 238 samples

##WORKING ON BACTERIA/ARCHAEA ONLY
# some QC checks of the "classified" reads per samples
min(sample_sums(phyloseq.bacteria)) # 0 (P1_0308)
max(sample_sums(phyloseq.bacteria)) # 80,778,091  (H21_0119) 
mean(sample_sums(phyloseq.bacteria)) #19,485,561
median(sample_sums(phyloseq.bacteria)) # 16,466,566
sort(sample_sums(phyloseq.bacteria))

##Zymo and controls#####
### pulling out samples from ZYMOs and EB, NTC and those samples with low OTUs
phyloseq.bacteria.controls <- subset_samples(phyloseq.bacteria, 
  grepl("NTC|EB|Zymo", sample_names(phyloseq.bacteria)))
phyloseq.bacteria.controls <- prune_taxa(taxa_sums(phyloseq.bacteria.controls) > 0, phyloseq.bacteria.controls) 
phyloseq.bacteria.controls #11962 taxa, 15 samples (NTC, EB and Zymos)

##Samples#####
##New phyloseq of just samples
phyloseq.bacteria.samples <- subset_samples(phyloseq.bacteria, 
                                             !grepl("NTC|EB|Zymo", sample_names(phyloseq.bacteria)))
phyloseq.bacteria.samples #33,613 taxa and 225 samples
#Taking out those with low counts
phyloseq.bacteria.samples <- prune_samples(sample_sums(phyloseq.bacteria.samples) > 100000, phyloseq.bacteria.samples) 
phyloseq.bacteria.samples <- prune_taxa(taxa_sums(phyloseq.bacteria.samples) > 0, phyloseq.bacteria.samples) 
phyloseq.bacteria.samples #33,613 taxa and 223 samples (dropped P1_0308 and H21_0109)
sort(sample_sums(phyloseq.bacteria.samples)) #OK

##Nitrifying taxa####
nitrifiers <- subset_taxa(phyloseq.bacteria.samples, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
nitrifiers #827 taxa and 223 samples
nitrifiers <- subset_samples(nitrifiers, sample_sums(nitrifiers) > 0)
nitrifiers #827 taxa and 223 samples 

##QC checks again
min(sample_sums(phyloseq.bacteria.samples)) #172,269 (H21_0120)
max(sample_sums(phyloseq.bacteria.samples)) #80,778,091 (H21_0119) 
mean(sample_sums(phyloseq.bacteria.samples)) #20,219,572
median(sample_sums(phyloseq.bacteria.samples)) #17,030,801
sort(sample_sums(phyloseq.bacteria.samples)) 

#COMPARING SEQUENCING DEPTHS #######
unclassified_counts_metadata <- unclassified_counts%>%
  filter(!c(grepl("EB|NTC|Zymo", SampleID)))%>%
  mutate(Enclosure = ifelse(grepl("H21", SampleID), "H21", "P1"))

###Established vs Naive####
sequencing_depth_P1vsH21<- ggplot(unclassified_counts_metadata, 
                                             aes(x = Enclosure, y= Total, 
                                                 color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Reads per Sample", color = "Enclosure", fill = "Enclosure", title = "Sequencing Depth") +
  #facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "WATER"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
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
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
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

##Stats 
wilcox_test(unclassified_counts_metadata, Total~Enclosure) #S. p = 0.00155

ggsave("sequencing_depth_P1vsH21.svg", 
       plot = sequencing_depth_P1vsH21, 
       device = "svg", width = 14, height =8)


#COMPARING SAMPLE SUMS#######
##ALL#####
sample.sums <- sample_sums(phyloseq.bacteria.samples) #making a sample sums object
phyloseq.bacteria.samples.samplessums.df <- cbind(phyloseq.bacteria.samples@sam_data, 
                                      sample.sums) #combining sample sums with metaphyloseq
phyloseq.bacteria.samples.samplessums.df
phyloseq.bacteria.samples.samplessums.df$sampleID <- rownames(phyloseq.bacteria.samples.samplessums.df) ##making a sampleID column


###Established vs Naive####
bacteria_archaea_samplesums_P1vsH21<- ggplot(phyloseq.bacteria.samples.samplessums.df, 
                                  aes(x = Enclosure, y= sample.sums, 
                                      color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "OTUs per Sample", color = "Enclosure", fill = "Enclosure", title = "Bacterial - Archaeal Read Counts") +
  #facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "WATER"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
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
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
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
wilcox_test(phyloseq.bacteria.samples.samplessums.df, sample.sums~Enclosure) #S. p = 2.97e-13

ggsave("bacteria_archaea_samplesums_P1vsH21.svg", 
       plot = bacteria_archaea_samplesums_P1vsH21, 
       device = "svg", width = 14, height =8)

##NITRIFIERS#####
sample.sums.nit <- sample_sums(nitrifiers) #making a sample sums object
nitrifiers.samplesums.df <- cbind(nitrifiers@sam_data, 
                                      sample.sums.nit) #combining sample sums with metadata
nitrifiers.samplesums.df
nitrifiers.samplesums.df$SampleID <- rownames(nitrifiers.samplesums.df) ##making a sampleID column

###Established vs Naive####
bacteria_archaea_samplesums_P1vsH21_nit<- ggplot(nitrifiers.samplesums.df, 
                                  aes(x = Enclosure, y= sample.sums.nit, 
                                      color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "OTUs per Sample", color = "Enclosure", fill = "Enclosure", title = "Nitrifying Taxa - Read Counts") +
  #facet_grid(~sample_type, scales = "free",  labeller = as_labeller(c("Feces" = "FECES", "Water" = "WATER"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
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
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
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
bacteria_archaea_samplesums_P1vsH21_nit
##Stats 
wilcox_test(nitrifiers.samplesums.df, sample.sums.nit~Enclosure) #S. p = 0.00335

ggsave("bacteria_archaea_samplesums_P1vsH21_nit.svg", 
       plot = bacteria_archaea_samplesums_P1vsH21_nit, 
       device = "svg", width = 14, height =8)

#TSS (RA) ####
any(sample_sums(phyloseq.bacteria.samples)== 0) ## no samples with 0 OTUs
phyloseq.bacteria.samples.ra <- transform_sample_counts(phyloseq.bacteria.samples, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data
##CLASSIFICATION PERCENTAGES AT DIFFERENT LEVELS ####
##PHYLUM
phyloseq.bacteria.samples_phylum.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Phylum", NArm = F) 
phyloseq.bacteria.samples_phylum.ra #3742 taxa and 223 samples 

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"])) #3742 taxa (so No duplicates)

Unknown_phylum_abundance <- phyloseq.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Phylum, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_phylum_abundance ##0.48% abundance by Unknown Phyla

Unclassified_phylum_abundance <- phyloseq.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##3.77% abundance by Unclassified Phyla

Classified_phylum_abundance <- phyloseq.bacteria.samples_phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##95.7% abundance by Classified Phyla

##Checking on excel
write.csv(phyloseq.bacteria.samples_phylum.ra@otu_table, "phylum_otus.csv")
write.csv(phyloseq.bacteria.samples_phylum.ra@tax_table, "phylum_taxa.csv")  

#How many unclassified?
phyloseq.bacteria.samples_phylum.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples_phylum.ra)
phyloseq.bacteria.samples_phylum.unclassified.ra #9 unclassified Phyla

#How many unknown?
phyloseq.bacteria.samples_phylum.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples_phylum.ra)
phyloseq.bacteria.samples_phylum.unknown.ra #3612 "unknown" Phyla

#Keep just classified Phyla
phyloseq.bacteria.samples_phylum.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria.samples_phylum.ra)
phyloseq.bacteria.samples_phylum.classified.ra ##121 classified (not unknown or unclassified) Phyla

##CLASS
phyloseq.bacteria.samples_class.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Class", NArm = F) 
phyloseq.bacteria.samples_class.ra #5149 classes and 223 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"])) #5149 classes (so No duplicates)

Unknown_class_abundance <- phyloseq.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Class, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_class_abundance #1.33% Abundance by Unknown classes

Unclassified_class_abundance <- phyloseq.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##5.78% Abundance by Unclassified Classes

Classified_class_abundance <- phyloseq.bacteria.samples_class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##92.9% Abundance by Classified classes

##Checking on excel
write.csv(phyloseq.bacteria.samples_class.ra@otu_table, "class_otus.csv")
write.csv(phyloseq.bacteria.samples_class.ra@tax_table, "class_taxa.csv") 


#How many unclassified?
phyloseq.bacteria.samples_class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"]),
  phyloseq.bacteria.samples_class.ra)
phyloseq.bacteria.samples_class.unclassified.ra #70 unclassified classes

#How many unknown?
phyloseq.bacteria.samples_class.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"]),
  phyloseq.bacteria.samples_class.ra)
phyloseq.bacteria.samples_class.unknown.ra #4930 "unknown" classes

#Keep just classified Classes
phyloseq.bacteria.samples_class.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_class.ra)[, "Class"]),
  phyloseq.bacteria.samples_class.ra)
phyloseq.bacteria.samples_class.classified.ra #149 classified classes

##ORDER
phyloseq.bacteria.samples_order.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Order", NArm = F) 
phyloseq.bacteria.samples_order.ra #5909 orders
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"])) #5907 orders (2 duplicates)
order_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"])
unique(order_taxa_vec[duplicated(order_taxa_vec)]) 
#"Candidatus Fermentimicrarchaeales", "Candidatus Cenarchaeales"

Unknown_order_abundance <- phyloseq.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Order, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_order_abundance ##2.46% abundance by Unknown Orders

Unclassified_order_abundance <- phyloseq.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##9.19% abundance by Unclassified Orders

Classified_order_abundance <- phyloseq.bacteria.samples_order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##88.4% abundance by Classified orders

#Checking on excel
write.csv(phyloseq.bacteria.samples_order.ra@otu_table, "order_otus.csv")
write.csv(phyloseq.bacteria.samples_order.ra@tax_table, "order_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"]),
  phyloseq.bacteria.samples_order.ra)
phyloseq.bacteria.samples_order.unclassified.ra #129 unclassified orders

#How many unknown?
phyloseq.bacteria.samples_order.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"]),
  phyloseq.bacteria.samples_order.ra)
phyloseq.bacteria.samples_order.unknown.ra #5448 "unknown" orders

#Keep just classified Orders
phyloseq.bacteria.samples_order.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_order.ra)[, "Order"]),
  phyloseq.bacteria.samples_order.ra)
phyloseq.bacteria.samples_order.classified.ra #332 classified orders
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_order.classified.ra)[, "Order"])) ##330 classified orders (unique - without duplicates)

##FAMILY
phyloseq.bacteria.samples_family.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Family", NArm = F) 
phyloseq.bacteria.samples_family.ra #7028 families
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"])) #7026 taxa (2 duplicates)
family_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"])
unique(family_taxa_vec[duplicated(family_taxa_vec)]) 
#"Candidatus Fermentimicrarchaeales", "Candidatus Cenarchaeales"

Unknown_family_abundance <- phyloseq.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_family_abundance #3.85% abundance by Unknown Families

Unclassified_family_abundance <- phyloseq.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##13.2% abundance by Unclassified Families

Classified_family_abundance <- phyloseq.bacteria.samples_family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##82.9% abundance by Classified Families

#Checking on excel
write.csv(phyloseq.bacteria.samples_family.ra@otu_table, "family_otus.csv")
write.csv(phyloseq.bacteria.samples_family.ra@tax_table, "family_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"]),
  phyloseq.bacteria.samples_family.ra)
phyloseq.bacteria.samples_family.unclassified.ra #271 unclassified families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.unclassified.ra)[, "Family"])) ##271 classified families (unique - without duplicates)

#How many unknown?
phyloseq.bacteria.samples_family.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"]),
  phyloseq.bacteria.samples_family.ra)
phyloseq.bacteria.samples_family.unknown.ra #5968 "unknown" families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.unknown.ra)[, "Family"]))#5968 "unknown" taxa (unique - without duplicates)

#Keep just classified Families
phyloseq.bacteria.samples_family.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_family.ra)[, "Family"]),
  phyloseq.bacteria.samples_family.ra)
phyloseq.bacteria.samples_family.classified.ra #789 classified families
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_family.classified.ra)[, "Family"]))#787 classified families (unique - without duplicates)

##GENUS 
phyloseq.bacteria.samples_genus.ra <- tax_glom(phyloseq.bacteria.samples.ra, taxrank = "Genus", NArm = F) 
phyloseq.bacteria.samples_genus.ra #10690 genera
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"])) #10646 taxa (44 duplicates)
genus_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"])
unique(genus_taxa_vec[duplicated(genus_taxa_vec)]) #24 duplicated unique ones

Unknown_genus_abundance <- phyloseq.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_genus_abundance ##6.21%  abundance by unknown genera

Unclassified_genus_abundance <- phyloseq.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_genus_abundance ##20.0% abundance by unclassified genera

Classified_genus_abundance <- phyloseq.bacteria.samples_genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##73.8% abundance by Classified Genera

#Checking on excel
write.csv(phyloseq.bacteria.samples_genus.ra@otu_table, "genus_otus.csv")
write.csv(phyloseq.bacteria.samples_genus.ra@tax_table, "genus_taxa.csv") 

#How many unclassified?
phyloseq.bacteria.samples_genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples_genus.ra)
phyloseq.bacteria.samples_genus.unclassified.ra #657 unclassified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.unclassified.ra)[, "Genus"])) ##657 unclassified genera (unique - without duplicates)

#How many unknown?
phyloseq.bacteria.samples_genus.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples_genus.ra)
phyloseq.bacteria.samples_genus.unknown.ra #6711 "unknown" genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.unknown.ra)[, "Genus"])) ##6711 unknown genera (unique - without duplicates)


#Keep just classified Genera
phyloseq.bacteria.samples_genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria.samples_genus.ra)[, "Genus"]),
  phyloseq.bacteria.samples_genus.ra)
phyloseq.bacteria.samples_genus.classified.ra #3322 classified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria.samples_genus.classified.ra)[, "Genus"])) ## 3278 classified genera (unique - without duplicates)


#ALPHA DIVERSITY ######
## ALL COMMUNITIES#####
alpha_div1 <- phyloseq::estimate_richness(phyloseq.bacteria.samples, measures = c("Observed", "Shannon")) # richness, diversity
alpha_div2 <- microbiome::evenness(phyloseq.bacteria.samples, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness

# combine metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div

alpha_div_meta <- cbind(phyloseq.bacteria.samples@sam_data, 
                        alpha_div) %>%
  #rownames_to_column(var = "SampleID")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))%>% #group into collection months
  mutate(Collection_Month = factor(Collection_Month))%>% # convert to factor for stat tests
  group_by(Enclosure)%>%
  filter(Collection_Date > "2023-10-01")%>% #Filtering out those samples in P1 from april and may 2023 and september from H21
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()
alpha_div_meta # metadata and div metrics

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
  labs(title= "ALPHA DIVERSITY", color = "Enclosure", fill = "Enclosure") +
  facet_wrap(~alpha_div_metric, 
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  # geom_pwc (method = "wilcox_test",
  #           label = "Wilcoxon, p = {p}",
  #           step.increase = 0.1,
  #           size = 0.5,
  #           label.size = 5,
  #           tip.length = 0.02,
  #           hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_P1vsH21


###Time series#####
###Alpha diversity indeces#####
alpha_div_time <- ggplot(alpha_div_meta_long, 
                            aes(x = Date_num, y= alpha_div_value)) +
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
alpha_div_time

##Water quality levels over time 
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


####Water quality levels over time######
##Time series
alpha_div_wq_time <- ggplot(alpha_div_wq_time_long%>%
                              filter(Index %in% c("Copper_level_mg_L",
                                                  "Ammonia_mg_L",
                                                  "Nitrite_mg_L",
                                                  "Nitrate_UV_mg_L", 
                                                  "Shannon",
                                                  "Observed",
                                                  "pielou")),
                            aes(x = Collection_Date, y = Index_value, color = Copper_keep)) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = expansion(mult = c(0.03, 0.03)))+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title= "MICROBIOME", 
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Ammonia_mg_L" = "NH3\n(mg/L)",
                                      "Nitrite_mg_L" = "Nitrite\n(mg/L)",
                                      "Nitrate_UV_mg_L" = "Nitrate\n(mg/L)",
                                      "Salinity_ppt" = "Salinity\n(ppt)", 
                                      "pH_spu" = "pH\n(spu)",
                                      "Shannon" = "Shannon",
                                      "Observed" = "Richness\n(Observed)",
                                      "pielou" = "Evenness\n(Pielou's)")))+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 15, hjust = 0.5),
        legend.title = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text.y  = element_text(colour = "white", size = 10, face = "bold"),
        strip.text.x  = element_text(colour = "white", size = 26, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 14,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))
alpha_div_wq_time
ggsave("alpha_div_wq_time_overall.png",
       alpha_div_wq_time, 
       device = "png", 
       dpi = 600, 
       height = 7, 
       width = 15)

alpha_div_wq_time_2 <- ggplot(alpha_div_wq_time_long%>%
                                filter(Index %in% c("Copper_level_mg_L",
                                                    "Shannon"
                                )),
                              aes(x = Collection_Date, y = Index_value, color = Copper_keep)) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = expansion(mult = c(0.03, 0.03)))+
  scale_color_viridis_c(option = "plasma")+
  #scale_color_manual(values = attempt.palette)+
  # guides(color = guide_legend(override.aes = list(size = 7)))+
  theme_bw() +
  labs(title = "MICROBIOME",
    color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon")))+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 15, angle = 0, vjust = 0.5),
        legend.title = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text.y  = element_text(colour = "white", size = 22, face = "bold"),
        strip.text.x  = element_text(colour = "white", size = 30, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 16,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 18),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))
alpha_div_wq_time_2

####Copper levels over time######
##Time series
copper_time <- ggplot(alpha_div_meta,
                      aes(x = Collection_Date, y = Copper_level_mg_L)) +
  theme_bw() +
  #labs(title= "ALPHA DIVERSITY") +
  facet_grid(~Enclosure,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  #geom_line(size = 1, color = "black", aes(group = 1)) +
  geom_point(aes(color = Copper_level_mg_L), size = 3, shape = 18) +
  scale_color_viridis_c(option = "plasma")+
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
copper_time


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
alpha_div1_nit <- phyloseq::estimate_richness(nitrifiers, 
                                              measures = c("Observed", "Shannon")) # richness, diversity
alpha_div2_nit <- microbiome::evenness(nitrifiers, index = "pielou", 
                                   zeroes = TRUE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness

# combine metrics with metadata
alpha_div_nit <- cbind(alpha_div1_nit, alpha_div2_nit)
alpha_div_nit

##Add metadata
alpha_div_nit_meta <- cbind(nitrifiers@sam_data, 
                        alpha_div_nit) %>%
  #rownames_to_column(var = "SampleID")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))%>% #group into collection months
  # mutate(Collection_Month = factor(Collection_Month)) %>% # convert to factor for stat tests
  group_by(Enclosure)%>%
  filter(Collection_Date > "2023-10-01")%>% #Filtering out those samples in P1 from april and may 2023 and september from H21
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()  
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
  labs(title= "ALPHA DIVERSITY - NITRIFIERS", color = "Enclosure", fill = "Enclosure") +
  facet_wrap(~alpha_div_metric, 
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                      "Shannon" = "DIVERSITY\n(SHANNON)",
                                      "pielou" = "EVENNESS"))) +
  geom_boxplot(alpha = 0.1) +
  geom_point(size = 3, shape = 18) +
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
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
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

##Water quality levels over time 
alpha_div_nit_wq_time_long <- alpha_div_nit_meta%>%
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
alpha_div_nit_wq_time_long 

####Water quality levels over time######
##Time series
alpha_div_nit_wq_time <- ggplot(alpha_div_nit_wq_time_long%>%
                                  filter(Index %in% c("Copper_level_mg_L",
                                                      "Ammonia_mg_L",
                                                      "Nitrite_mg_L",
                                                      "Nitrate_UV_mg_L", 
                                                      "Shannon",
                                                      "Observed",
                                                      "pielou")),
                                aes(x = Collection_Date, y = Index_value, color = Copper_keep)) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = expansion(mult = c(0.03, 0.03)))+
  scale_color_viridis_c(option = "plasma")+
  theme_bw() +
  labs(title= "NITRIFYING TAXA", 
       color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Ammonia_mg_L" = "NH3\n(mg/L)",
                                      "Nitrite_mg_L" = "Nitrite\n(mg/L)",
                                      "Nitrate_UV_mg_L" = "Nitrate\n(mg/L)",
                                      "Salinity_ppt" = "Salinity\n(ppt)", 
                                      "pH_spu" = "pH\n(spu)",
                                      "Shannon" = "Shannon",
                                      "Observed" = "Richness\n(Observed)",
                                      "pielou" = "Evenness\n(Pielou's)")))+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 15, hjust = 0.5),
        legend.title = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text.y  = element_text(colour = "white", size = 10, face = "bold"),
        strip.text.x  = element_text(colour = "white", size = 26, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 14,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))
alpha_div_nit_wq_time
ggsave("alpha_div_wq_time_nitrifiers.png",
       alpha_div_nit_wq_time, 
       device = "png", 
       dpi = 600, 
       height = 7, 
       width = 15)

alpha_div_nit_wq_time_2 <- ggplot(alpha_div_nit_wq_time_long%>%
                                    filter(Index %in% c("Copper_level_mg_L",
                                                        "Shannon")),
                                  aes(x = Collection_Date, y = Index_value, color = Copper_keep)) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = expansion(mult = c(0.03, 0.03)))+
  scale_color_viridis_c(option = "plasma")+
  #scale_color_manual(values = attempt.palette)+
  # guides(color = guide_legend(override.aes = list(size = 7)))+
  theme_bw() +
  labs(title = "NITRIFYING TAXA", color = "Copper level (mg/L)") +
  facet_grid(Index~ Enclosure,
             scales = "free", 
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon")))+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 15, angle = 0, vjust = 0.5),
        legend.title = element_text(size = 15, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text.y  = element_text(colour = "white", size = 22, face = "bold"),
        strip.text.x  = element_text(colour = "white", size = 30, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 16,
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 18),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))
alpha_div_nit_wq_time_2

####Copper vs Shannon levels #####
##Time series
copper_shannon_nit <- ggplot(alpha_div_nit_meta,
                         aes(x = Copper_level_mg_L, 
                             y = Shannon, 
                             color = Collection_Month)) +
  geom_point(size = 3, shape = 18) +
  theme_bw() +
  labs(title= "Shannon's Diversity vs Copper levels (mg/L) \nNitrifiers", 
       y = "SHANNON'S INDEX",
       x = "Copper levels (mg/L)",
       color = "Month") +
  facet_grid(~Enclosure,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  #geom_smooth(method="loess", se=TRUE) +
  #scale_color_viridis_c(option = "plasma")+
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
copper_shannon_nit

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


#BETA DIV#####
##ALL TAXA######
###FAMILY 
phyloseq.bacteria.samples_order.ra #5909 orders

phyloseq.bacteria.samples.order.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria.samples_order.ra, 
                                                                        "Enclosure", 
                                                                        level = "Order", threshold = 0.5)
phyloseq.bacteria.samples.order.filt #30 orders over 0.5% mean RA
phyloseq.bacteria.samples.order.filt.melt <- psmelt(phyloseq.bacteria.samples.order.filt)%>%
  mutate(Order = factor(Order, 
                         levels = c(setdiff(Order, 
                                            unique(grep("Others", Order, value = TRUE))), 
                                    unique(grep("Others", Order, value = TRUE)))))##Factoring the Order column so that "Others.." is the last category
levels(phyloseq.bacteria.samples.order.filt.melt$Order) ##ok

##Create color palette
order.filt.palette <- distinctColorPalette(length(unique(phyloseq.bacteria.samples.order.filt.melt$Order)))
order_filt_names <- unique(phyloseq.bacteria.samples.order.filt.melt$Order)# Create a named vector for the palette, where the names correspond to phlyum names
order_named_palette <- setNames((order.filt.palette)[1:length(order_filt_names)], order_filt_names)
order_named_palette$'Others <0.5% RA' <- "grey95"

##Apply the function to obtain top orders (n=15)
top_orders <- top_taxa_legend(phyloseq.bacteria.samples.order.filt.melt, taxlevel = "Order", n = 15)
top_orders

#Plot
RA_enclosures_overall_plot <- ggplot(phyloseq.bacteria.samples.order.filt.melt%>%
                                                  filter(!Collection_Date < "2023-10-01"), 
                                                aes(x=Collection_Date, y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)",
       title = "MICROBIOME") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = expansion(mult = c(0.03, 0.03)))+
  scale_fill_manual(values = order_named_palette,
                    breaks = top_orders) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.7, "cm"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
        strip.text.x  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 18,
                                   vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(colour = "black", size = 22),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),
        axis.title.x = element_blank())
RA_enclosures_overall_plot 

figure_alpha_overall_div_time_copper <- alpha_div_wq_time_2 + RA_enclosures_overall_plot  + 
  plot_layout(ncol = 1, heights = c(0.85, 1.25))
figure_alpha_overall_div_time_copper
ggsave("figure_alpha_overall_div_time_copper.png", 
       figure_alpha_overall_div_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 23)

##NITRIFIERS#####
###FAMILY 
phyloseq.bacteria.samples_family.ra #7028 families

#Out of this overall communities object, select only nitrifiers 
phyloseq.bacteria.samples_family.ra.nitrifiers <- subset_taxa(phyloseq.bacteria.samples_family.ra, 
                                                              Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
                                                                Family == "Chromatiaceae" | # no lineages
                                                                Family == "Nitrosopumilaceae" | # AOA; some!
                                                                Family == "Nitrosophaeraceae" | # no lineages
                                                                Order == "Nitrosomirales" | # no lineages
                                                                Order == "Nitrosocaldales" | # no lineages
                                                                Family == "Nitrospiraceae" | # NOB/Commamox; some!
                                                                Family == "Ectothiorhodospiraceae" | #NOB; some, plus a new one!
                                                                Family == "Nitrobacteraceae" | # none
                                                                Family == "Gallionellaceae" | # none
                                                                Family == "Nitrospinaceae") # NOB; some, plus a new one!
phyloseq.bacteria.samples_family.ra.nitrifiers <- subset_samples(phyloseq.bacteria.samples_family.ra.nitrifiers, 
                                                                 sample_sums(phyloseq.bacteria.samples_family.ra.nitrifiers) > 0)
phyloseq.bacteria.samples_family.ra.nitrifiers #8 nitrifying families in 223 samples 


#Melt to plot 
phyloseq.bacteria.samples_family.ra.nitrifiers.melt <- psmelt(phyloseq.bacteria.samples_family.ra.nitrifiers)

##Create color palette
family.palette <- distinctColorPalette(length(unique(phyloseq.bacteria.samples_family.ra.nitrifiers.melt$Family)))
family_names <- unique(phyloseq.bacteria.samples_family.ra.nitrifiers.melt$Family)# Create a named vector for the palette, where the names correspond to family names
family_named_palette <- setNames((family.palette)[1:length(family_names)], family_names)
#phylum_named_palette$'Others <0.5% RA' <- "grey95"

#Plot
RA_enclosures_nitrifiers.plot <- ggplot(phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
                                                  filter(!Collection_Date < "2023-10-01"), 
                                                aes(x=Collection_Date, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    expand = expansion(mult = c(0.03, 0.03)))+
  scale_fill_manual(values = family_named_palette) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.8, "cm"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
        strip.text.x  = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 18,
                                   vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(colour = "black", size = 22),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.8),
        axis.title.x = element_blank())
RA_enclosures_nitrifiers.plot 

figure_alpha_nit_div_time_copper <- alpha_div_nit_wq_time_2 + RA_enclosures_nitrifiers.plot  + 
  plot_layout(ncol = 1, heights = c(1.4, 0.6))
figure_alpha_nit_div_time_copper

ggsave("figure_alpha_nit_div_time_copper.png", 
       figure_alpha_nit_div_time_copper, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 23)



#Correlation of Nitrosopumilaceae with Copper levels
#H21
phyloseq.bacteria.samples_family.ra.AOA.melt.H21 <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
  filter(OTU == "OTU33084")%>%
  filter(Enclosure == "H21")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month),
         Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%  # convert to factor 
  filter(!Collection_Month %in% c("2023-09", "2024-03"))%>%
  filter(!is.na(Copper_level_mg_L),
         !is.na(Collection_Date))

gam_model_nit_AOA_H21 <- gam(Abundance ~ s(Copper_level_mg_L, by = Collection_Month) + Collection_Month,
                         data = phyloseq.bacteria.samples_family.ra.AOA.melt.H21)
summary(gam_model_nit_AOA_H21) ##No effect of enclosure, but copper effect did vary between enclosures 
plot(gam_model_nit_AOA_H21, pages = 1, shade = TRUE)

###Spearman correlation#####
cor.test(x = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Abundance, 
         y = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Copper_level_mg_L, 
         method = 'spearman') #Significant for naive one

H21_pcor_AOA <- pcor.test(x = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Abundance,
                          y = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Copper_level_mg_L,
                          z = phyloseq.bacteria.samples_family.ra.AOA.melt.H21$Date_num,
                          method = "pearson")

#P1
phyloseq.bacteria.samples_family.ra.AOA.melt.P1 <- phyloseq.bacteria.samples_family.ra.nitrifiers.melt%>%
  filter(OTU == "OTU33084")%>%
  filter(Enclosure == "P1")%>%
  filter(Collection_Date > "2023-06-01")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"),
         Collection_Month = factor(Collection_Month))  # convert to factor 

ggplot(phyloseq.bacteria.samples_family.ra.AOA.melt.P1,
      aes(x = Copper_level_mg_L, 
          y = Abundance)) +
  geom_point(size = 3, shape = 18, aes(color = Collection_Month)) +
  facet_grid(~Collection_Month,
             scales = "free")+
  labs(y = "NITRIFIERS RA (%)",
       x = "Copper levels (mg/L)") +
  theme_bw() +
  geom_smooth(method="loess", se=TRUE) +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top") +
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

#Modelling naive Nitrosopumilaceae (AOA)
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
phyloseq.bacteria.samples_family.ra.AOA.melt.P1.clean$fitted <- fitted(model_lm_nit_AOA_H21)

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


#No samples without copper levels:
phyloseq.bacteria.samples.copper <- subset_samples(phyloseq.bacteria.samples, !is.na(Copper_level_mg_L))
phyloseq.bacteria.samples.copper <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.copper)> 0, 
                                               phyloseq.bacteria.samples.copper)

#make DF from metadata
phyloseq.bacteria.samples.df <- as(phyloseq.bacteria.samples.copper@sam_data, "data.frame") %>%
  mutate(
    # Convert Collection_Date to Date if it isn’t already
    Collection_Date = as.Date(Collection_Date),
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month)  # convert to factor for PERMANOVA
  )

phyloseq.bacteria.samples.df 

###BRAY CURTIS#####
phyloseq.bacteria.samples.ra.bray <- vegdist(t(phyloseq.bacteria.samples.ra@otu_table), method = "bray")
phyloseq.bacteria.samples.ra.bray

###ORDINATION####
set.seed(98)
phyloseq.bacteria.samples.ra.bray.ord <- metaMDS(phyloseq.bacteria.samples.ra.bray, k=3, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
#### ADDING CENTROIDS FOR PLOTTING
## BC
#Simple ordination plot
phyloseq.bacteria.samples.ra.bray.plot <- ordiplot(phyloseq.bacteria.samples.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.bacteria.samples.ra.bray.scrs <- scores(phyloseq.bacteria.samples.ra.bray.plot, display = "sites")
#Add metadata to coordinates
phyloseq.bacteria.samples.ra.bray.scrs <- cbind(as.data.frame(phyloseq.bacteria.samples.ra.bray.scrs),
                                                Copper_level_mg_L = phyloseq.bacteria.samples.df$Copper_level_mg_L, 
                                                SampleID = phyloseq.bacteria.samples.df$SampleID, 
                                                Collection_Date = phyloseq.bacteria.samples.df$Collection_Date, 
                                                Enclosure = phyloseq.bacteria.samples.df$Enclosure)
##Calculate centroids according to Enclosure
phyloseq.bacteria.samples.ra.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
                                                     data = phyloseq.bacteria.samples.ra.bray.scrs, 
                                                     FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.bacteria.samples.ra.bray.segs <- merge(phyloseq.bacteria.samples.ra.bray.scrs, 
                             setNames(phyloseq.bacteria.samples.ra.bray.cent, c("Enclosure","cMDS1","cMDS2")),
                             by = 'Enclosure', 
                             sort = F)


##PERMANOVA###
set.seed(98)
enclosure_copper_BC_adonis  <- adonis2(phyloseq.bacteria.samples.ra.bray ~ Enclosure + 
                                  Copper_level_mg_L, 
                                phyloseq.bacteria.samples.df, 
                                strata = phyloseq.bacteria.samples.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
enclosure_copper_BC_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
#4.1% of the variation is due to copper levels, p = 1e-04


#With interaction
set.seed(98)
enclosure_copper_interac_BC_adonis  <- adonis2(phyloseq.bacteria.samples.ra.bray ~ Enclosure*Copper_level_mg_L, 
                                               phyloseq.bacteria.samples.df , 
                                strata = phyloseq.bacteria.samples.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
enclosure_copper_interac_BC_adonis
#Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS
#ENCLOSURE#
# Run the betadisper function, average distance to centroid
bray.enclosure.disp <- betadisper(phyloseq.bacteria.samples.ra.bray, 
                                  phyloseq.bacteria.samples.df$Enclosure)
bray.enclosure.disp
##Then test by permuting
set.seed(98)
bray.enclosure.permdisp <- permutest(bray.enclosure.disp, permutations = 9999)
bray.enclosure.permdisp ##NS, p = 0.214

#COPPER#
# Run the betadisper function, average distance to centroid
bray.copper.disp <- betadisper(phyloseq.bacteria.samples.ra.bray, 
                                  phyloseq.bacteria.samples.df$Copper_level_mg_L)
bray.copper.disp
##Then test by permuting
set.seed(98)
bray.copper.permdisp <- permutest(bray.copper.disp, permutations = 9999)
bray.copper.permdisp ##S, p = 0.044

# Extract R2 and p-values
R2_adonis_enclosure <- enclosure_BC_adonis$R2[1] 
pvalue_adonis_enclosure<-  enclosure_BC_adonis$`Pr(>F)`[1]
R2_adonis_copper <- enclosure_BC_adonis$R2[2] 
pvalue_adonis_copper <- enclosure_BC_adonis$`Pr(>F)`[2]

#### PLOTS
## BC
enclosure_BC_beta_div <- ggplot(phyloseq.bacteria.samples.ra.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", shape = "Enclosure") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
                                     colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) + # centroids
  #geom_text(aes (x= MDS1, y = MDS2,label= Collection_Date), colour= "black", size = 2.8, fontface = "bold") +
  geom_text(aes (x= cMDS1, y = cMDS2,label= Enclosure), colour= "white", size = 2.8, fontface = "bold") +
  scale_color_manual(values = enclosure.palette)+
  scale_fill_manual(values = enclosure.palette)+
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7)),
    fill = "none",
    color = "none"
  )+
  annotate("text", x = 2, y = 1.2, 
           label = "Enclosure", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = 2, y = 1.2, 
           label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+# Annotate R² and p-values
  annotate("text", x = 2, y = 0.7, 
           label = "Copper levels in water", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (sample type)
  annotate("text", x = 2, y = 0.7, 
           label = paste("R² = ", round(R2_adonis_copper* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_copper, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")# Annotate R² and p-values
 
enclosure_BC_beta_div


# Run envfit on your NMDS object 
# set.seed(98)
# fit <- envfit(phyloseq.bacteria.samples.ra.bray.ord, 
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



##NITRIFIERS######
###TSS (RA) ####
any(sample_sums(nitrifiers)== 0) ## no samples with 0 OTUs

nitrifiers.ra <- transform_sample_counts(nitrifiers, 
                                         function(x) x/sum(x)*100) ##Relative abundance from normalized data

nitrifiers.ra

####RA PLOT #######
#####FAMILY#######
nitrifiers.ra.family <- tax_glom(nitrifiers.ra, taxrank = "Family", NArm = F)
nitrifiers.ra.family #8 families 
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

nitrifiers.ra.family.melt$Collection_Date <- factor(nitrifiers.ra.family.melt$Collection_Date)

##Create color palette
family.palette <- distinctColorPalette(length(unique(nitrifiers.ra.family.melt$Family)))
family_names <- unique(nitrifiers.ra.family.melt$Family)# Create a named vector for the palette, where the names correspond to family names
family_named_palette <- setNames((family.palette)[1:length(family_names)], family_names)
#phylum_named_palette$'Others <0.5% RA' <- "grey95"

#Plot
RA_enclosures_nitrifiers.plot <- ggplot(nitrifiers.ra.family.melt, aes(x=Collection_Date, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)") +
  facet_grid(~Enclosure, scales = "free")+
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = family_named_palette) +
  guides(fill=guide_legend(title.position="top", nrow = 3))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 22),
        legend.text = element_text(size = 22),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.75),
        axis.title.y = element_text(size = 25),
        axis.text.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 0.5, colour = "black"),
        axis.title.x = element_blank()) 
RA_enclosures_nitrifiers.plot
dendroRA.phylum.plot.2 <- plot_grid(dendro.bray.plot, dendroRA.phylum.plot, 
                                    align = "v", 
                                    ncol = 1, 
                                    rel_heights = c(0.3, 0.7))

#No samples without copper levels:
nitrifiers.copper <- subset_samples(nitrifiers, !is.na(Copper_level_mg_L))
nitrifiers.copper <- prune_taxa(taxa_sums(nitrifiers.copper)> 0, 
                                               nitrifiers.copper)

nitrifiers.copper.ra <- transform_sample_counts(nitrifiers.copper, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data

#make DF from metadata
nitrifiers.df <- as(nitrifiers.copper@sam_data, "data.frame") %>%
  mutate(
    # Convert Collection_Date to Date if it isn’t already
    Collection_Date = as.Date(Collection_Date),
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month)  # convert to factor for PERMANOVA
  )

nitrifiers.df 

###BRAY CURTIS#####
nitrifiers.ra.bray <- vegdist(t(nitrifiers.ra@otu_table), method = "bray")
nitrifiers.ra.bray

###ORDINATION####
set.seed(98)
nitrifiers.ra.bray.ord <- metaMDS(nitrifiers.ra.bray, k=3, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
#### ADDING CENTROIDS FOR PLOTTING
## BC
#Simple ordination plot
nitrifiers.ra.bray.plot <- ordiplot(nitrifiers.ra.bray.ord$points)

#Now, extract coordinates
nitrifiers.ra.bray.scrs <- scores(nitrifiers.ra.bray.plot, display = "sites")
#Add metadata to coordinates
nitrifiers.ra.bray.scrs <- cbind(as.data.frame(nitrifiers.ra.bray.scrs),
                                                Copper_level_mg_L = nitrifiers.df$Copper_level_mg_L, 
                                                SampleID = nitrifiers.df$SampleID, 
                                                Collection_Date = nitrifiers.df$Collection_Date, 
                                                Enclosure = nitrifiers.df$Enclosure)
##Calculate centroids according to Enclosure
nitrifiers.ra.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
                                                    data = nitrifiers.ra.bray.scrs, 
                                                    FUN = mean) 
#Merge centroids with coordinates and metadata
nitrifiers.ra.bray.segs <- merge(nitrifiers.ra.bray.scrs, 
                                                setNames(nitrifiers.ra.bray.cent, c("Enclosure","cMDS1","cMDS2")),
                                                by = 'Enclosure', 
                                                sort = F)


##PERMANOVA###
set.seed(98)
nit_enclosure_copper_BC_adonis  <- adonis2(nitrifiers.ra.bray ~ Enclosure +
                                         Copper_level_mg_L, 
                                       nitrifiers.df, 
                                       strata = nitrifiers.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
nit_enclosure_copper_BC_adonis #4.8% of the variation is due to Enclosure, p = 1e-04
#3.5% of the variation is due to copper levels, p = 8e-04


#With interaction
set.seed(98)
nit_enclosure_copper_interac_BC_adonis  <- adonis2(nitrifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
                                               nitrifiers.df , 
                                               strata = nitrifiers.df$Collection_Month,
                                               by = "margin",
                                               permutations = 9999)
nit_enclosure_copper_interac_BC_adonis
#Enclosure:Copper_level_mg_L interaction significant (p = 1e-04). 7.3% of variation

##PERMDISPS
#ENCLOSURE#
# Run the betadisper function, average distance to centroid
bray.enclosure_nit.disp <- betadisper(nitrifiers.ra.bray, 
                                  nitrifiers.df$Enclosure)
bray.enclosure_nit.disp
##Then test by permuting
set.seed(98)
bray.enclosure_nit.permdisp <- permutest(bray.enclosure_nit.disp, permutations = 9999)
bray.enclosure_nit.permdisp ##NS, p = 0.65

#COPPER#
# Run the betadisper function, average distance to centroid
bray.copper_nit.disp <- betadisper(nitrifiers.ra.bray, 
                               nitrifiers.df$Copper_level_mg_L)
bray.copper_nit.disp
##Then test by permuting
set.seed(98)
bray.copper_nit.permdisp <- permutest(bray.copper_nit.disp, permutations = 9999)
bray.copper_nit.permdisp ##NS, p = 0.49

# Extract R2 and p-values
R2_adonis_enclosure_nit <- nit_enclosure_copper_BC_adonis$R2[1] 
pvalue_adonis_enclosure_nit <-  nit_enclosure_copper_BC_adonis$`Pr(>F)`[1]
R2_adonis_copper_nit <- nit_enclosure_copper_BC_adonis$R2[2] 
pvalue_adonis_copper_nit <- nit_enclosure_copper_BC_adonis$`Pr(>F)`[2]

R2_adonis_enclosurebynit <- nit_enclosure_copper_interac_BC_adonis$R2[1] 
pvalue_adonis_enclosurebynit <-  nit_enclosure_copper_interac_BC_adonis$`Pr(>F)`[1]

#### PLOTS
## BC
enclosure_nit_BC_beta_div <- ggplot(nitrifiers.ra.bray.segs) + theme_bw() +
  labs(x="NMDS1", y="NMDS2", title= "NITRIFYING COMMUNITIES", shape = "Enclosure") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
                                     colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) + # centroids
  #geom_text(aes (x= MDS1, y = MDS2,label= Collection_Date), colour= "black", size = 2.8, fontface = "bold") +
  geom_text(aes (x= cMDS1, y = cMDS2,label= Enclosure), colour= "white", size = 2.8, fontface = "bold") +
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 22),
        plot.title = element_text(size = 36),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 20, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7))
  )+
  annotate("text", x = 0.8, y = 1.7, 
           label = "Enclosure", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = 0.8, y = 1.7, 
           label = paste("R² = ", round(R2_adonis_enclosure_nit* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosure_nit, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black")+# Annotate R² and p-values
  annotate("text", x = 0.8, y = 1.3, 
           label = "Copper levels in water", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Copper)
  annotate("text", x = 0.8, y = 1.3, 
           label = paste("R² = ", round(R2_adonis_copper_nit* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_copper_nit, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") +# Annotate R² and p-values 
  annotate("text", x = 0.8, y = 0.9, 
           label = "Enclosure:Copper levels in water", 
           hjust = 0.5, vjust = -0.5, size = 6, colour = "black", fontface = "bold") + ##annotate variable (Copper:Enclosure)
  annotate("text", x = 0.8, y = 0.9, 
           label = paste("R² = ", round(R2_adonis_enclosurebynit* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_enclosurebynit, 4)), 
           hjust = 0.5, vjust = 1.1, size = 6, colour = "black") # Annotate R² and p-values 

enclosure_nit_BC_beta_div


# Run envfit on your NMDS object 
# set.seed(98)
# fit <- envfit(nitrifiers.ra.bray.ord, 
#               nitrifiers.df$Copper_level_mg_L, 
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