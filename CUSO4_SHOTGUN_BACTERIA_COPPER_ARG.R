##THIS SCRIPT IS FOR ANALYZING COMMUNITIES CLASSIFIED BY KRAKEN. 
##READS CAME FROM THOSE ALIGNED TO COPPER RESISTANCE GENES IN AMRPLUSPLUS

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
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/R_functions/MergeLowAbun_group_ARG.R")


#Importing data from kraken output nt_core - counts will be classified reads#### 
counts <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/Copper_ARG_reads/Conf_01/kraken_analytic_matrix.conf_0.1.csv')
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
          "kraken_taxonomy_copper_ARG_reads.csv",
          row.names = F)

###OTU table #####
otu_table <- counts[, -1]%>% #Excludes the first column (taxonomy)
  mutate(OTU = paste0("OTU", 1:nrow(counts))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##make into matrix so it is compatible with otu_table function from phyloseq
otu_table


#IMPORT METADATA####
#This comes from an already clean metadata file with data for both systems, as well as positive and negative controls from the "Metadata_cleaning.R" script
metadata <- read_tsv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/Copper_ARG_reads/Conf_01/annotations_copper_ARG_groups.tsv')

##Making into phyloseq-compatible object
sampledata_phyloseq <- metadata %>%
  select(Type, Class, Mechanism, Group)%>%
  mutate(rows = Group)%>%
  column_to_rownames(var= "rows") %>%##Make sampleID column into row names, so they match sample_names() with OTU and TAX
  sample_data(metadata) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

#PHYLOSEQ####
#Make phyloseq object
OTU <-phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
TAX <-phyloseq::tax_table(filled_taxonomy_2)
phyloseq <- phyloseq(OTU, TAX, sampledata_phyloseq)


#PREPROCESSING ####
phyloseq #543 taxa and 58 samples (58 copper resistance gene groups)
      
##Selecting only Bacteria/Archaea####
phyloseq.bacteria <- subset_taxa(phyloseq, Domain=="Archaea" | Domain=="Bacteria")
phyloseq.bacteria #525 taxa and 58 samples

##Selecting only viruses######
phyloseq.viruses <- subset_taxa(phyloseq, Domain=="Viruses")
taxanames_viruses <- c("Kingdom", "Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species") ##they have a different classification system, updating it here
colnames(phyloseq.viruses@tax_table) <- taxanames_viruses #replacing col names of the tax_table for new ones
colnames(phyloseq.viruses@tax_table) #OK taxonomy ranks
phyloseq.viruses #1 taxa and 58 samples

##Selecting only eukaryota #####
phyloseq.eukaryota <- subset_taxa(phyloseq, Domain=="Eukaryota")
colnames(phyloseq.eukaryota@tax_table) ##These are OK taxonomy ranks
phyloseq.eukaryota #17 taxa and 58 samples

#WORKING ON BACTERIA/ARCHAEA ONLY####
# some QC checks of the "classified" reads per samples
min(sample_sums(phyloseq.bacteria)) # 14 (COPP)
max(sample_sums(phyloseq.bacteria)) # 2316  (COPR) 
mean(sample_sums(phyloseq.bacteria)) # 222
median(sample_sums(phyloseq.bacteria)) # 87.5
sort(sample_sums(phyloseq.bacteria))

##NITRIFYING TAXA####
nitrifiers_all <- subset_taxa(phyloseq.bacteria, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
nitrifiers_all #14 taxa and 58 samples
nitrifiers <- subset_samples(nitrifiers_all, 
                             sample_sums(nitrifiers_all) > 0)
nitrifiers #14 taxa and 34 samples (34 copper resistance groups in nitrifying taxa)

#COMPARING CLASSIFIED READS BY KRAKEN#######
kraken_unclassified_reads <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/Copper_ARG_reads/Conf_01/unclassifieds_kraken_analytic_matrix.conf_0.1.csv')

#Filtering just samples included in phyloseq.bacteria, adding metadata, calculating percentage classified
kraken_unclassified_reads_samples_metadata <- kraken_unclassified_reads %>%
  dplyr::left_join(metadata, by = c("SampleID" = "Group"))%>%
  rename(Kraken2_Input_PairedEnd_Reads = Total, 
         Kraken2_Unclassified_PairedEnd_Reads = NumberUnclassified, 
         Kraken2_Unclassified_Percentage_Reads = PercentUnclassified)%>%
  mutate(Kraken2_Classified_Percentage_Reads = (100 - Kraken2_Unclassified_Percentage_Reads))
nrow(kraken_unclassified_reads_samples_metadata) #Ok, 58 samples (Copper resistaance groups)

###Kraken2 Classified Percentages Established vs Naive####
kraken2_classified_read_percentages_cu_groups <- ggplot(kraken_unclassified_reads_samples_metadata, 
                               aes(x = SampleID, 
                                   y= Kraken2_Classified_Percentage_Reads, 
                                   color = Mechanism)) +
  theme_bw() +
  labs(y= "Percentage (%) Classified Reads", color = "Mechanism") +
  geom_point(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.05,0,0.1,0)) +
  theme(plot.title = element_text(colour = "black", size = 32, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        strip.background = element_rect(fill = "black"),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.title.y = element_text(size = 22, colour = "black"),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(colour = "black", size =16, angle = 90, 
                                  vjust = 0.5, hjust = 0.5)
  )
kraken2_classified_read_percentages_cu_groups


#RELATIVE ABUNDANCE####
any(sample_sums(phyloseq.bacteria)== 0) ## no samples with 0 OTUs
phyloseq.bacteria.ra <- transform_sample_counts(phyloseq.bacteria, 
                                                              function(x) x/sum(x)*100) ##Relative abundance from normalized data
##CLASSIFICATION PERCENTAGES AT DIFFERENT TAXONOMIC LEVELS####
###PHYLUM######
phyloseq.bacteria_phylum.ra <- tax_glom(phyloseq.bacteria.ra, taxrank = "Phylum", NArm = F) 
phyloseq.bacteria_phylum.ra #24 phyla and 58 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria_phylum.ra)[, "Phylum"])) #24 phyla (so No duplicates)

Unknown_phylum_abundance <- phyloseq.bacteria_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Phylum, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_phylum_abundance ##0.723% abundance by Unknown Phyla

Unclassified_phylum_abundance <- phyloseq.bacteria_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##21.8% abundance by Unclassified Phyla

Classified_phylum_abundance <- phyloseq.bacteria_phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##77.5% abundance by Classified Phyla

#How many unclassified?
phyloseq.bacteria_phylum.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria_phylum.ra)
phyloseq.bacteria_phylum.unclassified.ra #3 unclassified Phyla

#How many unknown?
phyloseq.bacteria_phylum.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria_phylum.ra)
phyloseq.bacteria_phylum.unknown.ra #3 "unknown" Phyla

#Keep just classified Phyla
phyloseq.bacteria_phylum.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria_phylum.ra)[, "Phylum"]),
  phyloseq.bacteria_phylum.ra)
phyloseq.bacteria_phylum.classified.ra ##18 classified (not unknown or unclassified) Phyla

###CLASS#####
phyloseq.bacteria_class.ra <- tax_glom(phyloseq.bacteria.ra, taxrank = "Class", NArm = F) 
phyloseq.bacteria_class.ra #45 taxa and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria_class.ra)[, "Class"])) #45 classes (so No duplicates)

Unknown_class_abundance <- phyloseq.bacteria_class.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Class, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_class_abundance #0.759% Abundance by Unknown classes

Unclassified_class_abundance <- phyloseq.bacteria_class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##30% Abundance by Unclassified Classes

Classified_class_abundance <- phyloseq.bacteria_class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##69.2% Abundance by Classified classes


#How many unclassified?
phyloseq.bacteria_class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria_class.ra)[, "Class"]),
  phyloseq.bacteria_class.ra)
phyloseq.bacteria_class.unclassified.ra #11 unclassified classes

#How many unknown?
phyloseq.bacteria_class.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria_class.ra)[, "Class"]),
  phyloseq.bacteria_class.ra)
phyloseq.bacteria_class.unknown.ra #8 "unknown" classes

#Keep just classified Classes
phyloseq.bacteria_class.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria_class.ra)[, "Class"]),
  phyloseq.bacteria_class.ra)
phyloseq.bacteria_class.classified.ra #26 classified classes

###ORDER######
phyloseq.bacteria_order.ra <- tax_glom(phyloseq.bacteria.ra, taxrank = "Order", NArm = F) 
phyloseq.bacteria_order.ra #91 orders

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria_order.ra)[, "Order"])) #91 orders (no duplicates)

Unknown_order_abundance <- phyloseq.bacteria_order.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Order, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_order_abundance ##0.759% abundance by Unknown Orders

Unclassified_order_abundance <- phyloseq.bacteria_order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##41% abundance by Unclassified Orders

Classified_order_abundance <- phyloseq.bacteria_order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##58.3% abundance by Classified orders

#How many unclassified?
phyloseq.bacteria_order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria_order.ra)[, "Order"]),
  phyloseq.bacteria_order.ra)
phyloseq.bacteria_order.unclassified.ra #21 unclassified orders

#How many unknown?
phyloseq.bacteria_order.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria_order.ra)[, "Order"]),
  phyloseq.bacteria_order.ra)
phyloseq.bacteria_order.unknown.ra #8 "unknown" orders

#Keep just classified Orders
phyloseq.bacteria_order.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria_order.ra)[, "Order"]),
  phyloseq.bacteria_order.ra)
phyloseq.bacteria_order.classified.ra #62 classified orders
length(unique(phyloseq::tax_table(phyloseq.bacteria_order.classified.ra)[, "Order"])) ##62 classified orders (unique - without duplicates)

###FAMILY######
phyloseq.bacteria_family.ra <- tax_glom(phyloseq.bacteria.ra, taxrank = "Family", NArm = F) 
phyloseq.bacteria_family.ra #162 families
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria_family.ra)[, "Family"])) #162 taxa (no duplicates)

Unknown_family_abundance <- phyloseq.bacteria_family.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_family_abundance #0.848% abundance by Unknown Families

Unclassified_family_abundance <- phyloseq.bacteria_family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##47.3% abundance by Unclassified Families

Classified_family_abundance <- phyloseq.bacteria_family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##51.9% abundance by Classified Families

#How many unclassified?
phyloseq.bacteria_family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria_family.ra)[, "Family"]),
  phyloseq.bacteria_family.ra)
phyloseq.bacteria_family.unclassified.ra #41 unclassified families
length(unique(phyloseq::tax_table(phyloseq.bacteria_family.unclassified.ra)[, "Family"])) ##41 classified families (unique - without duplicates)

#How many unknown?
phyloseq.bacteria_family.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria_family.ra)[, "Family"]),
  phyloseq.bacteria_family.ra)
phyloseq.bacteria_family.unknown.ra #12 "unknown" families
length(unique(phyloseq::tax_table(phyloseq.bacteria_family.unknown.ra)[, "Family"]))#12 "unknown" taxa (unique - without duplicates)

#Keep just classified Families
phyloseq.bacteria_family.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria_family.ra)[, "Family"]),
  phyloseq.bacteria_family.ra)
phyloseq.bacteria_family.classified.ra #109 classified families
length(unique(phyloseq::tax_table(phyloseq.bacteria_family.classified.ra)[, "Family"]))#109 classified families (unique - without duplicates)

###GENUS ######
phyloseq.bacteria_genus.ra <- tax_glom(phyloseq.bacteria.ra, taxrank = "Genus", NArm = F) 
phyloseq.bacteria_genus.ra #304 genera
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria_genus.ra)[, "Genus"])) #304 taxa (no duplicates)

Unknown_genus_abundance <- phyloseq.bacteria_genus.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_genus_abundance ##0.986%  abundance by unknown genera

Unclassified_genus_abundance <- phyloseq.bacteria_genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance)) #Sum across OTUs
Unclassified_genus_abundance ##56.5% abundance by unclassified genera

Classified_genus_abundance <- phyloseq.bacteria_genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##42.5% abundance by Classified Genera

#How many unclassified?
phyloseq.bacteria_genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria_genus.ra)[, "Genus"]),
  phyloseq.bacteria_genus.ra)
phyloseq.bacteria_genus.unclassified.ra #81 unclassified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria_genus.unclassified.ra)[, "Genus"])) ##81 unclassified genera (unique - without duplicates)

#How many unknown?
phyloseq.bacteria_genus.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq.bacteria_genus.ra)[, "Genus"]),
  phyloseq.bacteria_genus.ra)
phyloseq.bacteria_genus.unknown.ra #16 "unknown" genera
length(unique(phyloseq::tax_table(phyloseq.bacteria_genus.unknown.ra)[, "Genus"])) ##3755 unknown genera (unique - without duplicates)


#Keep just classified Genera
phyloseq.bacteria_genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria_genus.ra)[, "Genus"]),
  phyloseq.bacteria_genus.ra)
phyloseq.bacteria_genus.classified.ra #207 classified genera
length(unique(phyloseq::tax_table(phyloseq.bacteria_genus.classified.ra)[, "Genus"])) ##207 classified genera (unique - without duplicates)


###SPECIES######
phyloseq.bacteria.ra ##525 Species- OTUs
phyloseq.bacteria_species.ra <- phyloseq.bacteria.ra
phyloseq.bacteria_species.ra #525 Species

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.bacteria_species.ra)[, "Species"])) #525 species (no duplicates)

Unclassified_species_abundance <- phyloseq.bacteria_species.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Species, ignore.case = TRUE)) %>%  # Filter unclassified <tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_species_abundance ##72.2% abundance by unclassified species

Classified_species_abundance <- phyloseq.bacteria_species.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Species, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_species_abundance ##27.8% abundance by Classified Species

#How many unclassified?
phyloseq.bacteria_species.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.bacteria_species.ra)[, "Species"]),
  phyloseq.bacteria_species.ra)
phyloseq.bacteria_species.unclassified.ra #178 unclassified species
length(unique(phyloseq::tax_table(phyloseq.bacteria_species.unclassified.ra)[, "Species"])) ##178 unclassified species (unique - without duplicates)

#Keep just classified Genera
phyloseq.bacteria_species.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.bacteria_species.ra)[, "Species"]),
  phyloseq.bacteria_species.ra)
phyloseq.bacteria_species.classified.ra #347 classified species
length(unique(phyloseq::tax_table(phyloseq.bacteria_species.classified.ra)[, "Species"])) ##347 classified species (unique - without duplicates)

#RELATIVE ABUNDANCE######
#ALL TAXA######
## ORDER #####
phyloseq.bacteria_order.ra #91 taxa and 58 samples (ARG copper resistance groups)

#Grouping the low abundance orders into one category
phyloseq.bacteria.order.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria_order.ra, 
                                                               "Mechanism",
                                                                        level = "Order", 
                                                                        threshold = 0.7)
phyloseq.bacteria.order.filt #23 orders over 0.5% mean RA
phyloseq.bacteria.order.filt.melt <- psmelt(phyloseq.bacteria.order.filt)%>%
  mutate(Order = factor(Order, 
                         levels = c(setdiff(Order, 
                                            unique(grep("Others", Order, value = TRUE))), 
                                    unique(grep("Others", Order, value = TRUE)))))##Factoring the Order column so that "Others.." is the last category
levels(phyloseq.bacteria.order.filt.melt$Order) ##ok

##Create color palette
#order.filt.palette <- distinctColorPalette(length(unique(phyloseq.bacteria.order.filt.melt$Order)))
order.filt.palette <- unname(alphabet2())
order_filt_names <- unique(phyloseq.bacteria.order.filt.melt$Order)# Create a named vector for the palette, where the names correspond to phlyum names
order_named_palette <- setNames((order.filt.palette)[1:length(order_filt_names)], order_filt_names)
order_named_palette$'Others <0.7% RA' <- "grey95"
order_named_palette$'Flavobacteriales' <-  "#63A184"
order_named_palette$'Rhodobacterales' <- "#E3B199"
order_named_palette$'unclassified Alphaproteobacteria' <- "dodgerblue"
##Apply the function to obtain top orders (n=15)
top_orders <- top_taxa_legend(phyloseq.bacteria.order.filt.melt, 
                              taxlevel = "Order", n = 18)
top_orders

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_order_overall_copper_ARG_plot <- ggplot(phyloseq.bacteria.order.filt.melt,
                                     aes(x=Sample, 
                                         y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Copper ARG Group") +
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_fill_manual(values = order_named_palette,
                    breaks = top_orders,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
        legend.position = "right",
        # legend.position = c(1.09, 0.5),  # x, y inside plot
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18, face = "bold"),
        legend.key.size = unit(0.6, "cm"),
        panel.border = element_rect(colour = "black", linewidth= 1),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 20,
                                   angle = 90,
                                   vjust = 0.5, hjust = 0.5),
        axis.title = element_text(colour = "black", size = 22),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_order_overall_copper_ARG_plot

## FAMILY #####
phyloseq.bacteria_family.ra #162 families and 216 samples 

phyloseq.bacteria.family.filt <- merge_low_abundance_grouped_ra(phyloseq.bacteria_family.ra, 
                                                                             "Mechanism", 
                                                                             level = "Family", 
                                                                             threshold = 0.7)
phyloseq.bacteria.family.filt #33 families over 0.5% mean RA
phyloseq.bacteria.family.filt.melt <- psmelt(phyloseq.bacteria.family.filt)%>%
  mutate(Family = factor(Family, 
                        levels = c(setdiff(Family, 
                                           unique(grep("Others", Family, value = TRUE))), 
                                   unique(grep("Others", Family, value = TRUE)))))##Factoring the Family column so that "Others.." is the last category
levels(phyloseq.bacteria.family.filt.melt$Family) ##ok

##Create color palette - based on families within the same order
palette_family_level_df <- phyloseq.bacteria.family.filt.melt %>% 
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
family_named_palette$'Others <0.7% RA' <- "grey95"

##Apply the function to obtain top familys (n=15)
top_families <- top_taxa_legend(phyloseq.bacteria.family.filt.melt, 
                              taxlevel = "Family", n = 15)
top_families

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_family_overall_copper_ARG_plot <- ggplot(phyloseq.bacteria.family.filt.melt,
                                             aes(x=Sample, y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Copper ARG Groups") +
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_fill_manual(values = family_named_palette,
                    breaks = top_families,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    #legend.position = c(1.08, 0.5),  # x, y inside plot
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
                               angle = 90,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_overall_copper_ARG_plot

#NITRIFIERS WITHIN THE OVERALL COMMUNITY#####
## FAMILY #######
phyloseq.bacteria_family.ra #162 families

##Which families are nitrifiers? 
nitrifiers.melt <- psmelt(nitrifiers)
unique(nitrifiers.melt$Family) #""Nitrobacteraceae", "Nitrospinaceae", "Nitrosopumilaceae", "Chromatiaceae", "Nitrospiraceae"  

#Out of this overall communities object, select only nitrifiers 
phyloseq.bacteria_family.ra.nitrifiers <- subset_taxa(phyloseq.bacteria_family.ra, 
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
phyloseq.bacteria_family.ra.nitrifiers <- subset_samples(phyloseq.bacteria_family.ra.nitrifiers, 
                                                                 sample_sums(phyloseq.bacteria_family.ra.nitrifiers) > 0)
phyloseq.bacteria_family.ra.nitrifiers #5 nitrifying families in 34 samples (groups) 


#Melt to plot 
phyloseq.bacteria_family.ra.nitrifiers.melt <- psmelt(phyloseq.bacteria_family.ra.nitrifiers)

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.bacteria_family.ra.nitrifiers.melt <- phyloseq.bacteria_family.ra.nitrifiers.melt %>%
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
palette_nitrifiers_family_df <- phyloseq.bacteria_family.ra.nitrifiers.melt %>% 
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
top_nitrifying_families <- top_taxa_legend(phyloseq.bacteria_family.ra.nitrifiers.melt, 
                                           n = 5)
top_nitrifying_families

#Reorder by same linege
top_nitrifying_families <- c("Nitrobacteraceae",#NOB
                             "Nitrospinaceae", #NOB
                             "Nitrospiraceae", #NOB
                             "Nitrosopumilaceae", #AOA
                             "Chromatiaceae")#AOB
top_nitrifying_families

#Plot
RA_family_enclosures_nit_copper_ARG_plot <- ggplot(phyloseq.bacteria_family.ra.nitrifiers.melt,
                                            aes(x=Sample, 
                                                y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Copper ARG Group") +
  # facet_grid(~Enclosure, 
  #            scales = "free",
  #            labeller = as_labeller(c("P1" = "Established",
  #                                     "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_fill_manual(values = palette_nitrifiers_family,
                    breaks = top_nitrifying_families
                    ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    #legend.position = c(1.07, 0.5),  # x, y inside plot
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
                               angle = 90, 
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_enclosures_nit_copper_ARG_plot

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

