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
library(lmerTest); library(mgcv); library(rmcorr)

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
                       "Phosphate_mg_L", "Copper_mg_L", "Duplicate_date_metadata")
colnames(water_quality_P1) <- newcolnames_P1_WQ
colnames(water_quality_P1) #OK

##Discard units from data
water_quality_P1 <- water_quality_P1 %>%
  mutate(across(!c(Request_Date, Enclosure, Duplicate_date_metadata),
                ~ parse_number(.)))
str(water_quality_P1) #OK now

#Fixing date
str(water_quality_P1) #OK. Only want to change Collection_Date from "chr' to date 
water_quality_P1$Request_Date 
water_quality_P1$Request_Date <- as.Date(water_quality_P1$Request_Date, format = "%m/%d/%y")
str(water_quality_P1$Request_Date) #Ok now 

#Duplicate_date as factor 
water_quality_P1$Duplicate_date_metadata <- factor(water_quality_P1$Duplicate_date_metadata, levels = c("0", "1"))
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
                              by = c("Collection_Date" = "Request_Date", "Enclosure"))

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

## Let's subset out our nitrifying taxa####
nitrifiers <- subset_taxa(phyloseq.bacteria.samples, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
nitrifiers #813 taxa and 223 samples
nitrifiers <- subset_samples(nitrifiers, sample_sums(nitrifiers) > 0)
nitrifiers #813 taxa and 223 samples 

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
alpha_div_meta_long$alpha_div_metric<- factor(alpha_div_meta_long$alpha_div_metric, levels = c("Observed","pielou", "Shannon"))

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
alpha_div_wq_time_long <- alpha_div_meta%>%
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
                                     "pH_spu",
                                     "Ammonia_mg_L",
                                     "Nitrite_mg_L",
                                     "Nitrate_UV_mg_L", 
                                     "Salinity_ppt",
                                     "Shannon",
                                     "Observed",
                                     "pielou"
                                     )),
                      aes(x = Date_num, y = Index_value)) +
  geom_point(size = 3, shape = 18)+
  theme_bw() +
  labs(title= "BACTERIAL - ARCHAEAL COMMUNITIES") +
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
  #geom_line(size = 1, color = "black", aes(group = 1)) +
  #geom_point(aes(color = Copper_level_mg_L), size = 3, shape = 18) +
  scale_color_viridis_c(option = "plasma")+
    theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 16, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 8, angle = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_wq_time

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
                      aes(x = Copper_level_mg_L, y = Shannon)) +
  theme_bw() +
  labs(title= "Shannon's Diversity vs Copper levels (mg/L)", y = "SHANNON'S INDEX") +
  facet_grid(~Enclosure,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  geom_smooth(method="loess", se=TRUE) +
  #geom_line(size = 1, color = "black", aes(group = 1)) +
  geom_point(aes(color = Copper_level_mg_L), size = 3, shape = 18) +
  scale_color_viridis_c(option = "plasma")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12, angle = 45, vjust = 0.5),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 20),
        axis.text.x = element_text(colour = "black", size = 8, angle = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
copper_shannon


##Linear models ####
alpha_div_meta_clean <- alpha_div_meta %>%
    filter(!is.na(Shannon),
           !is.na(Copper_level_mg_L),
           !is.na(Collection_Date))%>%
  arrange(Enclosure, Collection_Date)%>%
  mutate(Enclosure = factor(Enclosure, levels = c("P1", "H21")))

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
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness

# combine metrics with metadata
alpha_div_nit <- cbind(alpha_div1_nit, alpha_div2_nit)
alpha_div_nit

##Add metadata
alpha_div_nit_meta <- cbind(nitrifiers@sam_data, 
                        alpha_div_nit) %>%
  #rownames_to_column(var = "SampleID")%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))%>% #group into collection months
  mutate(Collection_Month = factor(Collection_Month)) %>% # convert to factor for stat tests
  group_by(Enclosure)%>%
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
alpha_div_nit_wq_time <- ggplot(alpha_div_nit_wq_time_long%>%
                              filter(Index %in% c("Copper_level_mg_L",
                                                  "pH_spu",
                                                  "Ammonia_mg_L",
                                                  "Nitrite_mg_L",
                                                  "Nitrate_UV_mg_L", 
                                                  "Salinity_ppt",
                                                  "Shannon",
                                                  "Observed",
                                                  "pielou"
                              )),
                            aes(x = Date_num, y = Index_value)) +
  geom_point(size = 3, shape = 18)+
  theme_bw() +
  labs(title= "NITRIFIERS") +
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
  #geom_line(size = 1, color = "black", aes(group = 1)) +
  #geom_point(aes(color = Copper_level_mg_L), size = 3, shape = 18) +
  scale_color_viridis_c(option = "plasma")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 16, face = "bold"),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = "black", size = 8, angle = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        plot.title = element_text(colour = "black", size = 30, face = "bold"))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_nit_wq_time


###Copper vs Shannon levels #####
##Time series
copper_shannon_nit <- ggplot(alpha_div_nit_meta,
                         aes(x = Copper_level_mg_L, y = Shannon)) +
  theme_bw() +
  labs(title= "Shannon's Diversity vs Copper levels (mg/L) \nNitrifiers", y = "SHANNON'S INDEX") +
  facet_grid(~Enclosure,
             scales = "free",
             #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  geom_smooth(method="loess", se=TRUE) +
  #geom_line(size = 1, color = "black", aes(group = 1)) +
  geom_text(aes(label = SampleID), vjust = -0.5, size = 3, angle = 90)+
  geom_point(aes(color = Copper_level_mg_L), size = 3, shape = 18) +
  scale_color_viridis_c(option = "plasma")+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12, angle = 45, vjust = 0.5),
        legend.title = element_text(size = 22, face = "bold"),
        strip.background = element_rect(fill = "black"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", size = 20),
        axis.text.x = element_text(colour = "black", size = 8, angle = 45),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
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
  mutate(Enclosure = factor(Enclosure, levels = c("P1", "H21")),
         Collection_Date = factor(Collection_Date))

#Per enclosure
alpha_div_nit_meta_clean_H21 <- alpha_div_nit_meta_clean%>%
  filter(Enclosure =="H21")
alpha_div_nit_meta_clean_P1 <- alpha_div_nit_meta_clean%>%
  filter(Enclosure =="P1")

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


####LME#####
##Collection_month: trying to accounts for repeated measurements or clustering by month: each month may have its own baseline Shannon diversity.
#Tells me the average effect of copper, enclosure, and their interaction on Shannon diversity across all months.
#This assumes a relationship that is aprox linear!
model_lme_nit <- lmer(Shannon ~ Copper_level_mg_L * Enclosure + (1 | Collection_Month),
                  data = alpha_div_nit_meta_clean)
summary(model_lme_nit) ##Effect of copper, enclosure, and their interaction 

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

# alpha_div_nit_meta_clean_P1$Date_num <- as.numeric(as.Date(alpha_div_nit_meta_clean_P1$Collection_Date))
# alpha_div_nit_meta_clean_H21$Date_num <- as.numeric(as.Date(alpha_div_nit_meta_clean_H21$Collection_Date))

#install.packages("ppcor")
library(ppcor)
#H21
pcor.test(x = alpha_div_nit_meta_clean_H21$Shannon,
          y = alpha_div_nit_meta_clean_H21$Copper_level_mg_L,
          z = alpha_div_nit_meta_clean_H21$Date_num,
          method = "spearman")
#P1
pcor.test(alpha_div_nit_meta_clean_P1$Shannon,
          alpha_div_nit_meta_clean_P1$Copper_level_mg_L,
          alpha_div_nit_meta_clean_P1["Date_num"],
          method = "spearman")


#BETA DIV#####
##ALL TAXA######
###TSS (RA) ####
any(sample_sums(phyloseq.bacteria.samples)== 0) ## no samples with 0 OTUs

#No samples without copper levels:
phyloseq.bacteria.samples.copper <- subset_samples(phyloseq.bacteria.samples, !is.na(Copper_level_mg_L))
phyloseq.bacteria.samples.copper <- prune_taxa(taxa_sums(phyloseq.bacteria.samples.copper)> 0, 
                                               phyloseq.bacteria.samples.copper)

phyloseq.bacteria.samples.ra <- transform_sample_counts(phyloseq.bacteria.samples.copper, 
                                    function(x) x/sum(x)*100) ##Relative abundance from normalized data

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

#No samples without copper levels:
nitrifiers.copper <- subset_samples(nitrifiers, !is.na(Copper_level_mg_L))
nitrifiers.copper <- prune_taxa(taxa_sums(nitrifiers.copper)> 0, 
                                               nitrifiers.copper)

nitrifiers.ra <- transform_sample_counts(nitrifiers.copper, 
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