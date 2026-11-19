# Set working directory, load libraries####
setwd('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Analysis/Bacteria_archaea/Figures')

# Installs (if needed) and loads packages from CRAN, Bioconductor, GitHub
load_packages <- function(cran_pkgs = character(0),
                          bioc_pkgs = character(0),
                          github_pkgs = character(0)) {
  
  # Make sure BiocManager and remotes are available for installs
  if (length(bioc_pkgs) > 0 && !requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (length(github_pkgs) > 0 && !requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  #  CRAN packages 
  for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing CRAN package: ", pkg)
      install.packages(pkg)
    }
  }
  
  # Bioconductor packages 
  for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing Bioconductor package: ", pkg)
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }
  }
  
  # GitHub packages (named vector: pkg_name = "owner/repo") 
  for (pkg in names(github_pkgs)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing GitHub package: ", github_pkgs[[pkg]])
      remotes::install_github(github_pkgs[[pkg]])
    }
  }
  
  # Finally, load all
  all_pkgs <- c(cran_pkgs, bioc_pkgs, names(github_pkgs))
  invisible(lapply(all_pkgs, function(pkg) {
    library(pkg, character.only = TRUE)
  }))
  
  message("All packages loaded")
}

##Packages lists
#CRAN
cran_pkgs = c(
  "tidyverse", "ggplot2", "stringr", "dplyr", "vegan", "cowplot",
  "ggdendro", "randomcoloR", "ggpubr", "ggsignif", "UpSetR", 'zoo',
  "ggtext", "ggnewscale", "rstatix", "ggrepel", "ggh4x", "svglite",
  "writexl", "paletteer", "lme4", "lmerTest", "car", "emmeans", 'rmcorr','ppcor',
  "Polychrome", "colorspace", "devtools", "remotes", 'pals', 'gratia', 'purrr', 
  'mgcv', 'changepoint', 'RColorBrewer', 'patchwork'
)
#BioConductor
bioc_pkgs = c(
  "phyloseq", "metagenomeSeq", "ANCOMBC", "maaslin3",
  "MicrobiotaProcess", "microbiome"
)
#GitHub
github_pkgs = c(
  metagMisc      = "vmikk/metagMisc",
  pairwiseAdonis = "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",
  maaslin3       = "biobakery/maaslin3"  
)

#Load them with the function 
load_packages(cran_pkgs, bioc_pkgs, github_pkgs)

##Source functions
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbundanceOthersPercentage.R')
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/top_taxa_legend_updated.R")
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/fill_taxonomy_updated.R')
source('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/R_functions/MergeLowAbun_group_microbiome.R')
source("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/Feedlot_Lagoon_Project/R_functions/MergeLowAbun_group_ARG.R")

#KRAKEN2 TAXONOMY######
#Counts from core_nt kraken classifications###
counts_raw_ntcore <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/paired_end_ntcore_and_GTDB_20261027/Conf_05/Ntcore/kraken_analytic_matrix.conf_0.05.csv')

#Unclassified counts from core_nt kraken classifications##
unclassified_counts_ntcore <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/paired_end_ntcore_and_GTDB_20261027/Conf_05/Ntcore/unclassifieds_kraken_analytic_matrix.conf_0.05.csv')
  
#Counts of reads that were classified with ntcore as bacteria, archaea, just root, or just cellular organisms or left unclassified were afterwards classified using a GTDB kraken index##
counts_raw_GTDB <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/paired_end_ntcore_and_GTDB_20261027/Conf_05/GTDB/kraken_analytic_matrix_GTDB.conf_005.csv')

#Unclassified counts from GTDB kraken classifications##
unclassified_counts_GTDB <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/paired_end_ntcore_and_GTDB_20261027/Conf_05/GTDB/unclassifieds_kraken_analytic_matrix_GTDB.conf_005.csv')

##Separating into taxonomic ranks#####
###Nt core (has 8 classification ranks from Domain to Species)########
tax_ranks_separated_ntcore <- counts_raw_ntcore %>%
  separate(taxa, 
           into = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|",#splits the strings by the "|" symbol.
           fill = "right") #fill = "right", missing components are added as "NA" to the right (last columns) instead of to the left 
tax_ranks_separated_ntcore

###GTDB (has 7 classification ranks from Bacteria to Species)########
tax_ranks_separated_GTDB_tax <- counts_raw_GTDB %>%
  separate(taxa, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|",#splits the strings by the "|" symbol.
           fill = "right") #fill = "right", missing components are added as "NA" to the right (last columns) instead of to the left 
tax_ranks_separated_GTDB_tax

#Add rowname as "OTU#"
tax_ranks_separated_GTDB_tax <- tax_ranks_separated_GTDB_tax%>%
  mutate(OTU = paste0("GTDB_OTU", 1:nrow(.))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") ##Make OTU column into row names

##Kraken ntcore eukaryotic and viral counts###############
#From ntcore based classifications, keeping only eukaryotic and viral counts
ntcore_tax_eukaryota_viruses <- tax_ranks_separated_ntcore %>%
  dplyr::filter(Domain %in% c("Eukaryota", "Viruses"))%>%
  mutate(OTU = paste0("ntcore_OTU", 1:nrow(.))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") ##Make OTU column into row names
nrow(ntcore_tax_eukaryota_viruses) #96,695 eukaryotic and viral taxa out of 
nrow(tax_ranks_separated_ntcore) #125836 total taxa classified with ntcore

#Separate counts (eukarya, viral) from taxonomy
counts_ntcore_tax_eukaryota_viruses <- ntcore_tax_eukaryota_viruses %>%
  select(-c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#Tax table for (eukarya, viral)
tax_table_ntcore_eukaryota_viruses <- ntcore_tax_eukaryota_viruses%>%
  select(Domain, Kingdom, Phylum, Class, Order, Family, Genus, Species)

##Filling up taxonomy (actual NAs and string "NA"s in the taxonomy table)
filled_taxonomy_ntcore_eukaryota_viruses <- fill_taxonomy(tax_table_ntcore_eukaryota_viruses) ##apply the function to the taxonomy table
anyNA(filled_taxonomy_ntcore_eukaryota_viruses) ##OK, no NAs now 
grep("^NA$", filled_taxonomy_ntcore_eukaryota_viruses, value = T) ##OK, no "NA" strings now

#Now, keeping only GTDB-like ranks (Domain, Phylum, Class, Order, Family, Genus, Species)
filled_taxonomy_ntcore_eukaryota_viruses_GTDB_structure <- 
  filled_taxonomy_ntcore_eukaryota_viruses %>%
  dplyr::select(Domain, Phylum, Class, Order, Family, Genus, Species)
filled_taxonomy_ntcore_eukaryota_viruses_GTDB_structure
nrow(filled_taxonomy_ntcore_eukaryota_viruses_GTDB_structure) #96,695 eukaryotic and viral taxa, OK

##Kraken GTDB bacteria and archaeal counts #####
nrow(tax_ranks_separated_GTDB_tax) # 220,083 bacterial and archaeal taxa
#Separate counts (bacteria, archaea) from taxonomy
counts_GTDB_tax_bacteria_archaea <- tax_ranks_separated_GTDB_tax %>%
  dplyr::select(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
counts_GTDB_tax_bacteria_archaea

##Extracting just taxonomy (bacteria, archaea)
tax_table_GTDB_bacteria_archaea <- tax_ranks_separated_GTDB_tax %>%
  dplyr::select(Domain, Phylum, Class, Order, Family, Genus, Species)
tax_table_GTDB_bacteria_archaea

##Filling up taxonomy (actual NAs and string "NA"s in the taxonomy table)
filled_taxonomy_GTDB_bacteria_archaea <- fill_taxonomy(tax_table_GTDB_bacteria_archaea) ##apply the function to the taxonomy table
anyNA(filled_taxonomy_GTDB_bacteria_archaea) ##OK, no NAs now 
grep("^NA$", filled_taxonomy_GTDB_bacteria_archaea, value = T) ##OK, no "NA" strings now
nrow(filled_taxonomy_GTDB_bacteria_archaea) #220,083 bacterial and archaeal taxa, OK


#Making sure the total counts that went into being classified with GTDB (unclassified from ntcore, 
#classified as root, cellular organisms, bacteria, archaea) plus the ntcore counts
#of eukaryotic and viral organisms equal the total number of counts that went into 
#being classified (trimmomatic output)

#From the reports, calculated counts from each domain: 
domains_ntcore_counts <- read_delim('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/paired_end_ntcore_and_GTDB_20261027/Conf_05/Ntcore/table2_domains.tsv', delim = "\t")

#From the kraken reports, also calculated how many were root, unclassified, or just cellular organisms
unclassified_root_cellular_organisms <- read_delim('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Kraken2/paired_end_ntcore_and_GTDB_20261027/Conf_05/Ntcore/table1_unclassified_root_cellular_bacteria_archaea.tsv', 
                                                   delim = "\t")

#Total counts that went in from trimmomatic
total_trimmomatic_output_reads <- unclassified_counts_ntcore %>% #the 'total' here is how many reads from trimmomatic went into being classified initially with ntcore
  select(SampleID, Total)%>%
  rename(Trimmomatic_counts = Total)
total_trimmomatic_output_reads

#Now, total counts 'classified' by ntcore
total_ntcore_classified_reads <- tibble(
  SampleID = colnames(counts_raw_ntcore%>%select(-taxa)),
  Ntcore_classified_counts = colSums(counts_raw_ntcore%>%select(-taxa)),
)

#Now, the counts that went into being classified by GTDB
total_counts_into_GTDB <- unclassified_counts_GTDB %>%
  select(SampleID, Total) %>%
  rename(Counts_into_GTDB = Total)
total_counts_into_GTDB

#Now, the counts that are eukaryotic and viral from GTDB
total_eukaryotic_viral_counts <- tibble(
  SampleID = colnames(counts_ntcore_tax_eukaryota_viruses),
  Eukaryotic_Viral_counts = colSums(counts_ntcore_tax_eukaryota_viruses)
)

#Now merge them and see
checking_GTDB_ntcore_counts <- total_trimmomatic_output_reads %>%
  full_join(total_ntcore_classified_reads, by = "SampleID")%>%
  full_join(domains_ntcore_counts%>%select(domains_sum, sample), by = c("SampleID" = "sample"))%>%
  full_join(unclassified_root_cellular_organisms, by = c("SampleID" = "sample"))%>%
  full_join(total_counts_into_GTDB, by = "SampleID") %>%
  full_join(total_eukaryotic_viral_counts, by = "SampleID")%>%
  mutate(
    check_counts = Counts_into_GTDB + Eukaryotic_Viral_counts,
    counts_mismatch = check_counts == Trimmomatic_counts, 
    domain_plus_cellular_plus_root = domains_sum + root_only + cellular_only,
    bacteria_archaea_root_only_cellular_only_unclassified = bacteria + archaea + root_only + cellular_only + unclassified, 
    check_GTDB_counts = Counts_into_GTDB == bacteria_archaea_root_only_cellular_only_unclassified, 
    ntcore_more_tan_domains = Ntcore_classified_counts > domains_sum,
    gtdb_plus_ntcore_counts = Eukaryotic_Viral_counts + Counts_into_GTDB, 
    trimm_minus_counts = Trimmomatic_counts - gtdb_plus_ntcore_counts, 
    gtdb_plus_ntcore_counts_equals_trimmomatic_input = gtdb_plus_ntcore_counts == Trimmomatic_counts
    )
  
#Bind the two taxonomy tables (ntcore and GTDB) and make into a matrix######
filled_taxonomy_matrix<- bind_rows(filled_taxonomy_GTDB_bacteria_archaea,
                                   filled_taxonomy_ntcore_eukaryota_viruses_GTDB_structure) %>%
  as.matrix() ##convert into matrix for phyloseq
filled_taxonomy_matrix 
nrow(filled_taxonomy_matrix) #316,778 = 220,083 bacterial - archaeal taxa plus 96,695 eukaryotic and viral taxa

#Bind the two count tables (ntcore and GTDB)######
nrow(counts_GTDB_tax_bacteria_archaea) #220,083 bacterial and archaeal taxa
nrow(counts_ntcore_tax_eukaryota_viruses) #96,695 eukaryotic and viral taxa
counts_ntcore_GTDB <- bind_rows(counts_GTDB_tax_bacteria_archaea,
                                counts_ntcore_tax_eukaryota_viruses) 
counts_ntcore_GTDB
nrow(counts_ntcore_GTDB) #316,778 = 220,083 bacterial - archaeal taxa plus 96,695 eukaryotic and viral taxa


##There are some extra samples that will not be used for this project. Filtering those#####
dropping_samples <- c("H21_0912", "H21_1005", "P1_0420", "P1_0427", 
                      "P1_0504", "H21_1202b", "H21_1021a", "H21_1021b",
                      #Ozone in the established system was turned on 04/03/2023 at 10:00 AM. Going to drop those from 0404 on
                      "P1_0404", "P1_0407", "P1_0411", "P1_0414", "P1_0416", "P1_0423", "P1_0430",
                      #Protein Skimmer and UV were turned back on in the naive system on February 21st, 2024
                      "H21_0221", "H21_0223", "H21_0227", "H21_0228", "H21_0302"
)

#Dropping them from the count matrix
counts_ntcore_GTDB_filt <- counts_ntcore_GTDB %>%
  dplyr::select(-all_of(dropping_samples))
ncol(counts_ntcore_GTDB_filt) #222 = 207 samples + 3 mocks + 12 negative controls (11 NTC and one EB)

#OTU table (ntcore + GTDB)#####
otu_table <- counts_ntcore_GTDB_filt %>% 
  as.matrix() ##make into matrix so it is compatible with otu_table function from phyloseq
otu_table

#METADATA####
#This comes from an already clean metadata file with data for both systems, as well as positive and negative controls from the "Metadata_cleaning.R" script
metadata <- readr::read_csv("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Sample_metadata/metadata_all_systems_phyloseq_final.csv")

#Factor ordering
metadata <- metadata %>%
  # Convert Collection_Date to Date if it isn’t already
  mutate(
    Collection_Date  = as.Date(trimws(Collection_Date), format = "%Y-%m-%d"),
    Collection_Month = format(Collection_Date, "%Y-%m")
  )%>%
  mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
                                       levels = c("Low Copper Levels (Day 1-27)", 
                                                  "Transition Period 1 (Day 29-38)", 
                                                  "First Tx Copper Exposure (Day 39-51)", 
                                                  "Transition Period 2 (Day 52-81)", 
                                                  "Transition Period 3 (Day 86-108)", 
                                                  "Second Tx Copper Exposure (Day 109-135)"
                                                  #"Post-Treatment (Day 136-146)"
                                                  )))%>%
  mutate(Date_num_phase_naive_abbrv = factor(Date_num_phase_naive_abbrv, 
                                             levels = c("L",
                                                        "T1",
                                                        "E1",
                                                        "T2", 
                                                        "T3",
                                                        "E2" 
                                                        #"P"
                                                        )))%>%
  mutate(Date_num_phase_established = factor(Date_num_phase_established, 
                                             levels = c("Low Copper Levels (Day 1-53)",
                                                        "Transition Period 1 (Day 54-65)",
                                                        "Tx Copper Exposure (Day 66-104)",
                                                        "Post-Treatment (Day 105-141)")))%>%
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
                                   #"Post-Treatment (Day 136-146)", 
                                   #Established
                                   "Low Copper Levels (Day 1-53)",
                                   "Transition Period 1 (Day 54-65)",
                                   "Tx Copper Exposure (Day 66-104)",
                                   "Post-Treatment (Day 105-141)")))%>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = c("L", "T1", "E", 
                                                  "E1", "T2", "T3", "E2", 
                                                  "P")))

##Making into phyloseq-compatible object
sampledata_phyloseq <- metadata %>%
  mutate(rows = SampleID)%>%
  column_to_rownames(var= "rows") %>%##Make sampleID column into row names, so they match sample_names() with OTU and TAX
  sample_data(metadata) ##use phyloseq function sample_data() to make metadata into phyloseq sample data object

#PHYLOSEQ OBJECT####
#Make phyloseq object
OTU <-phyloseq::otu_table(otu_table, taxa_are_rows = TRUE)
TAX <-phyloseq::tax_table(filled_taxonomy_matrix)
phyloseq <- phyloseq(OTU, TAX, sampledata_phyloseq)

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU), metadata$SampleID) #No

#COLOR PALETTES#####
enclosure.palette <- c("H21" = "#fc8d62",  
                       "P1"  = "#8da0cb" )

reads.palette <- c("Raw (Sequencing Depth)" ="#E69F00", 
                   "QC - Trimmomatic" = "#0072B2")


naive_phase_palette <- c("L" = "#A3A33D",
                         "T1" = "#0CE3A3",
                         "E1" = "#7EDC8A",
                         "T2" = "#066E10",
                         "T3" = "#9E2307B4",
                         "E2" = "#DE1B1BE6",
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
phyloseq #316,778 taxa and 222 samples (207 actual samples, 3 controls, 12 NTC and EBs)

##CONTROLS####
### Getting samples from ZYMOs and EB, NTC
control_IDs <- grep("NTC|EB|Zymo", sample_names(phyloseq), value = T)
phyloseq.controls <- phyloseq::prune_samples(control_IDs, phyloseq)
phyloseq.controls <- prune_taxa(taxa_sums(phyloseq.controls) > 0, phyloseq.controls) 
phyloseq.controls #19,794 taxa, 15 samples (NTC, EB and Zymos)


###ZYMO MOCK COMMUNITIES#####
phyloseq.controls.zymo <- subset_samples(phyloseq.controls, 
                                                  grepl("Zymo", sample_names(phyloseq.controls)))
phyloseq.controls.zymo <- prune_taxa(taxa_sums(phyloseq.controls.zymo) > 0, 
                                              phyloseq.controls.zymo) 
phyloseq.controls.zymo #19,626 taxa, 3 samples

#Relative abundance
phyloseq.controls.zymo.ra <- transform_sample_counts(phyloseq.controls.zymo, 
                                                              function(x) x/sum(x)*100) ##Relative abundance 
####GENUS LEVEL ############
phyloseq.controls.zymo.ra.genus <- tax_glom(phyloseq.controls.zymo.ra, taxrank = "Genus", NArm = F)
phyloseq.controls.zymo.ra.genus #8989 genera and 3 samples

#Melt to long format at genus level
phyloseq.controls.zymo.ra.genus.melt <- psmelt(phyloseq.controls.zymo.ra.genus)

#What are the top genera
top_genera_zymo <- phyloseq.controls.zymo.ra.genus.melt %>%
  group_by(Family, Genus) %>%
  summarise(
    mean_ra = mean(Abundance, na.rm = TRUE),
    sd_ra   = sd(Abundance,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  #Rank by abundance, not alphabetically
  arrange(desc(mean_ra)) %>%
  mutate(
    `Mean Relative Abundance (%) ± SD` = paste0(
      round(mean_ra, 3),
      " ± ",
      round(sd_ra, 3)))%>%
  head(n = 30)%>%
  dplyr::select(Genus,`Mean Relative Abundance (%) ± SD`)
top_genera_zymo #Does not include cryptococcus or sacharomyces, since those are not bacteria or archaea in GTDB

#Plot - GENUS 
phyloseq.controls.zymo.ra.genus.filt <- merge_low_abundance_grouped_ra(phyloseq.controls.zymo.ra.genus, 
                                                                                variable = "Enclosure",
                                                                                level = "Genus", 
                                                                                threshold = 0.5)
phyloseq.controls.zymo.ra.genus.filt #15 genera over 0.5% mean RA

#Melt the object
phyloseq.controls.zymo.ra.genus.filt.melt <- psmelt(phyloseq.controls.zymo.ra.genus.filt)

#Order by family:
genus_order_zymo_controls <- phyloseq.controls.zymo.ra.genus.filt.melt %>%
  distinct(Family, Genus) %>%
  mutate(is_other = grepl("Others", Genus)) %>%
  arrange(is_other, Family, Genus) %>%   # non-Others first, by Family then Genus; Others last
  pull(Genus) %>%
  unique()

#Order the melted object
phyloseq.controls.zymo.ra.genus.filt.melt <- phyloseq.controls.zymo.ra.genus.filt.melt %>%
  mutate(Genus = factor(Genus, levels = genus_order_zymo_controls))
levels(phyloseq.controls.zymo.ra.genus.filt.melt$Genus) #OK

#Palette
set.seed(20)
zymo.genus.palette <- as.character(palette36.colors(length(unique(phyloseq.controls.zymo.ra.genus.filt.melt$Genus))))
names(zymo.genus.palette) <- unique(phyloseq.controls.zymo.ra.genus.filt.melt$Genus)
zymo.genus.palette$'Cryptococcus'<- "#1f77b4"
zymo.genus.palette$'Saccharomyces' <-"#ff7f0e"
zymo.genus.palette$'Listeria'<- "#2ca02c"
zymo.genus.palette$'Enterococcus'<- "#d62728"
zymo.genus.palette$'Staphylococcus' <- "#9467bd"
zymo.genus.palette$'Bacillus' <- "#8c564b"
zymo.genus.palette$'Limosilactobacillus' <-"#e377c2"
zymo.genus.palette$'Pseudomonas' <- "#7f7f7f"
zymo.genus.palette$'Salmonella'  <- "#bcbd22"
zymo.genus.palette$'Escherichia'<- "#17becf"
zymo.genus.palette$"G047199115" <- 'cyan'
zymo.genus.palette$'Others <0.5% RA' <- 'grey90'


#Italicized labels (markdown)
#"G047199115" in GTDB is classified as E coli on NCBI
zymo_genus_level_italic_labels <- c(
  "Bacillus" = "*Bacillus*",
  "Cryptococcus" = "*Cryptococcus*", 
  "Escherichia" = "Escherichia",
  "G047199115" = "*G047199115* (NCBI *E. Coli*)", 
  "Salmonella" = "*Salmonella*", 
  "unclassified Enterobacteriaceae" = "unclassified *Enterobacteriaceae*", 
  "Enterococcus" = "*Enterococcus*", 
  "Limosilactobacillus" = "*Limosilactobacillus*",
  "Listeria" = "*Listeria*",
  "Pseudomonas" = "*Pseudomonas*",
  "Saccharomyces" = "*Saccharomyces*",
  "Staphylococcus" = "*Staphylococcus*",
  "unclassified Bacteria" = "unclassified Bacteria", 
  "Others <0.5% RA" = "Others <0.5% RA"
)

#Plot
phyloseq.controls.zymo.ra.genus.plot <- ggplot(
  phyloseq.controls.zymo.ra.genus.filt.melt, aes(x=Sample, y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Zymo Mock Communities") +
  geom_bar(stat = "identity", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = zymo.genus.palette, 
                    breaks = levels(phyloseq.controls.zymo.ra.genus.filt.melt$Genus),
                    labels = zymo_genus_level_italic_labels) +
  guides(fill=guide_legend(title.position="top", nrow = 7))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 24),
        legend.text = element_markdown(size = 22), 
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
phyloseq.controls.zymo.ra.genus.plot

####SPECIES LEVEL########
phyloseq.controls.zymo.ra.melt <- psmelt(phyloseq.controls.zymo.ra)

#What are the top species?
top_species_zymo <- phyloseq.controls.zymo.ra.melt %>%
  group_by(Species) %>%
  summarise(
    mean_ra = mean(Abundance, na.rm = TRUE),
    sd_ra   = sd(Abundance,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  #Rank by abundance, not alphabetically
  arrange(desc(mean_ra)) %>%
  mutate(
    `Mean Relative Abundance (%) ± SD` = paste0(
      round(mean_ra, 2),
      " ± ",
      round(sd_ra, 3)))%>%
  head(n = 40)%>%
  dplyr::select(Species,`Mean Relative Abundance (%) ± SD`)
top_species_zymo 

#Plot - SPECIES 
phyloseq.controls.zymo.ra.species.filt <- merge_low_abundance_grouped_ra(phyloseq.controls.zymo.ra, 
                                                                                  variable = "SampleID",
                                                                                  level = "Species", 
                                                                                  threshold = 0.5)
phyloseq.controls.zymo.ra.species.filt #20 species over 0.5% mean RA

#Melt object
phyloseq.controls.zymo.ra.species.filt.melt <- psmelt(phyloseq.controls.zymo.ra.species.filt)

#Order species by genus
species_order_zymo_controls <- phyloseq.controls.zymo.ra.species.filt.melt %>%
  distinct(Species, Genus) %>%
  mutate(is_other = grepl("Others", Species)) %>%
  arrange(is_other, Genus, Species) %>%   # non-Others first, by Family then Genus; Others last
  pull(Species) %>%
  unique()

#Order the melted object
phyloseq.controls.zymo.ra.species.filt.melt <- phyloseq.controls.zymo.ra.species.filt.melt %>%
  mutate(Species = factor(Species, levels = species_order_zymo_controls))
levels(phyloseq.controls.zymo.ra.species.filt.melt$Species) #OK

#Color palette
#Create base colors based on genera
zymo_genus_base_colors <- zymo.genus.palette

#Make hues based on species within each genus
palette_zymo_genus_df <- phyloseq.controls.zymo.ra.species.filt.melt %>% 
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

#Italicized labels (markdown)
#"G047199115" in GTDB is classified as E coli on NCBI
zymo_species_level_italic_labels <- c(
  "Bacillus rugosus" = "*Bacillus rugosus*",
  "Bacillus spizizenii" = "*Bacillus spizizenii*", 
  "Bacillus spizizenii_A" = "*Bacillus spizizenii* A",
  "unclassified Bacillus" = "unclassified *Bacillus*", 
  "Cryptococcus neoformans" = "*Cryptococcus neoformans*", 
  "Enterococcus faecalis" = "*Enterococcus faecalis*", 
  "unclassified Escherichia" = "unclassified *Escherichia*",
  "G047199115 sp047199115" = "*G047199115* (NCBI *E. Coli*) *sp047199115*",
  "Limosilactobacillus fermentum" = "*Limosilactobacillus fermentum*",
  "Listeria monocytogenes_B" = "*Listeria monocytogenes* B",
  "unclassified Listeria" = "unclassified *Listeria*",
  "Pseudomonas aeruginosa" = "*Pseudomonas aeruginosa*",
  "Saccharomyces cerevisiae" = "*Saccharomyces cerevisiae*",
  "Salmonella enterica" = "*Salmonella enterica*",
  "unclassified Salmonella" = "unclassified *Salmonella*",
  "Staphylococcus aureus" = "*Staphylococcus aureus*",
  "unclassified Staphylococcus" = "unclassified *Staphylococcus*",
  "unclassified Bacteria" = "unclassified Bacteria", 
  "unclassified Enterobacteriaceae" = "unclassified *Enterobacteriaceae*",
  "Others <0.5% RA" = "Others <0.5% RA" 
)
zymo_species_level_italic_labels

#PLOT
phyloseq.controls.zymo.ra.species.plot <-
  ggplot(phyloseq.controls.zymo.ra.species.filt.melt, aes(x=Sample, y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", title = "Zymo Mock Communities") +
  geom_bar(stat = "summary", colour = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_x_discrete(expand = c(0.03,0,0.03,0)) +
  scale_fill_manual(values = palette_zymo_genus,
                    labels = zymo_species_level_italic_labels) +
  guides(fill=guide_legend(title.position="top", ncol = 2))+
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_markdown(size = 18),
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
phyloseq.controls.zymo.ra.species.plot

#####SUPPLEMENTARY FIGURE S2 #####
sfigure2 <- phyloseq.controls.zymo.ra.species.plot
ggsave("SupplementaryFigure2.png", 
       sfigure2, 
       device = "png", 
       width = 10, height =10, 
       dpi = 500)

#####SUPPLEMENTARY TABLE 2 - MOCK COMMUNITIES SPECIES####
stable2 <- top_species_zymo
write_xlsx(stable2, 
           "SupplementaryTable2.xlsx")


#Adding an OTU to the phyloseq object to represent unclassified reads by kraken#########
#First, calculate number of unclassified reads. Will calculate it as the trimmomatic total read counts - sample sums of the phyloseq object
#Trimmed reads
cuso4_trimmed_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Trimmed/CuSO4_clean_trimmed_read_counts.csv')
cuso4_trimmed_read_counts <- cuso4_trimmed_read_counts %>%
  select(SampleID, Num_Reads_Forward_Trimmed_Paired)%>%
  rename(Total_Counts_Paired_End_Trimmomatic_Output = Num_Reads_Forward_Trimmed_Paired)

#Sample sums from the phyloseq object
sample_sums_ntcore_gtdb_kraken_classified <- data.frame(
  SampleID = names(sample_sums(phyloseq)),
  ntcore_GTDB_classified_counts = as.numeric(sample_sums(phyloseq)),
  row.names = NULL
)

#Unclassified reads will be  = trimmed reads - sample sums (ntcore and GTDB classifications)
kraken_ntcore_gtdb_kraken_unclassified <- merge(
  cuso4_trimmed_read_counts, 
  sample_sums_ntcore_gtdb_kraken_classified, 
  by = 'SampleID')%>%
  mutate(
    Total_Unclassified_Kraken_Counts = 
      Total_Counts_Paired_End_Trimmomatic_Output - ntcore_GTDB_classified_counts, 
    Percent_Reads_Unclassified_Kraken = 
      Total_Unclassified_Kraken_Counts/Total_Counts_Paired_End_Trimmomatic_Output *100
      )

#Make an 'Unclassifieds' OTU table 
phyloseq.unclassified_otutable_df <- kraken_ntcore_gtdb_kraken_unclassified%>%
  dplyr::select(SampleID, Total_Unclassified_Kraken_Counts)%>%
  rename(Unclassified_Reads = Total_Unclassified_Kraken_Counts)%>%
  column_to_rownames(var = 'SampleID')%>%
  t(.)%>%
  data.frame()

#OTU table from the original phyloseq
phyloseq_otutable_df <- data.frame(otu_table(phyloseq))

#Add the 'unclassified' OTU to the main phyloseq object
phyloseq_unclassified_otutable <- bind_rows(
  phyloseq_otutable_df, 
  phyloseq.unclassified_otutable_df)%>% #unclassified reads OTU
  as.matrix()%>%
  otu_table(.,  taxa_are_rows = TRUE) #make into OTU table
nrow(phyloseq_unclassified_otutable) #316,779 taxa (316778 plus 1 Unclassified reads OTU)
phyloseq #316,778 taxa

#Add the 'unclassified' OTU taxonomy (all unclassified reads)
phyloseq.unclassified_taxtable_df <- data.frame(
  Phylum = "Unclassified Reads",
  Class = "Unclassified Reads",
  Order = "Unclassified Reads",
  Family = "Unclassified Reads",
  Genus = "Unclassified Reads",
  Species = "Unclassified Reads",
  check.names = FALSE,
  row.names = "Unclassified_Reads")

#Merge with the original OTU table (with only classified reads)
phyloseq_unclassified_taxtable <- bind_rows(
  data.frame(phyloseq::tax_table(phyloseq)), #original tax table, without 'unclassified reads' OTU
  phyloseq.unclassified_taxtable_df)%>% #unclassified reads tax table
  as.matrix()%>%
  phyloseq::tax_table(.)
phyloseq_unclassified_taxtable

#Now, make phyloseq object
phyloseq_with_unclassifieds_otu <- phyloseq(
  phyloseq_unclassified_otutable,
  phyloseq_unclassified_taxtable,
  sample_data(phyloseq)
)

#Now, all counts should add up to the initial trimmed read counts 
sample_sums(phyloseq_with_unclassifieds_otu)
kraken_ntcore_gtdb_kraken_unclassified%>%
  dplyr::select(SampleID, Total_Counts_Paired_End_Trimmomatic_Output) #OK!


#SAMPLES ONLY######
##Phyloseq without unclassified OTUs########
phyloseq.samples <- prune_samples(!grepl("NTC|EB|Zymo", 
                                         sample_names(phyloseq)), 
                                  phyloseq) 
phyloseq.samples <- prune_taxa(taxa_sums(phyloseq.samples) > 0, 
                               phyloseq.samples) 
phyloseq.samples #309,784 taxa and 207 samples
sort(sample_sums(phyloseq.samples))
#Taking out those with low counts (P1_0308 has 2, H21_0109 has 167, H21_0120 has 357,519 and then jumps to 9,364,439 in H21_1101)
phyloseq.samples <- prune_samples(sample_sums(phyloseq.samples) > 400000, phyloseq.samples) 
phyloseq.samples <- prune_taxa(taxa_sums(phyloseq.samples) > 0, phyloseq.samples) 
phyloseq.samples #309,747 taxa and 204 samples (dropped P1_0308, H21_0109, and H21_0120)
sort(sample_sums(phyloseq.samples)) #OK

##Phyloseq with unclassified OTUs########
phyloseq_with_unclassifieds_otu.samples <- prune_samples(!grepl("NTC|EB|Zymo", 
                                         sample_names(phyloseq_with_unclassifieds_otu)), 
                                  phyloseq_with_unclassifieds_otu) 
phyloseq_with_unclassifieds_otu.samples <- prune_taxa(taxa_sums(phyloseq_with_unclassifieds_otu.samples) > 0, 
                               phyloseq_with_unclassifieds_otu.samples) 
phyloseq_with_unclassifieds_otu.samples #309,785 taxa and 207 samples
sort(sample_sums(phyloseq_with_unclassifieds_otu.samples))

#Taking out those with low counts (P1_0308 has 5, H21_0109 has 439, H21_0120 has 687,179 and then jumps to 16,825,028 in P1_0223)
phyloseq_with_unclassifieds_otu.samples <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples) > 700000, 
                                                         phyloseq_with_unclassifieds_otu.samples) 
phyloseq_with_unclassifieds_otu.samples <- prune_taxa(taxa_sums(phyloseq_with_unclassifieds_otu.samples) > 0, 
                                                      phyloseq_with_unclassifieds_otu.samples) 
phyloseq_with_unclassifieds_otu.samples #309,748 taxa and 204 samples (dropped P1_0308, H21_0109, and H21_0120)
sort(sample_sums(phyloseq_with_unclassifieds_otu.samples)) #OK

##COMPARING SEQUENCING DEPTHS#######
cuso4_raw_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Raw/CuSO4_clean_raw_read_counts.csv')
str(cuso4_raw_read_counts)

#Keeping just those samples that I'm analyzing in phyloseq.samples (dropped P1_0308, H21_0109, and H21_0120 because of very low counts)
sampleIDs_phyloseq.samples <- phyloseq.samples@sam_data$SampleID 
length(sampleIDs_phyloseq.samples) #204 samples

###Zymo mock communities #####
cuso4_raw_read_counts_zymo <- cuso4_raw_read_counts %>%
  filter(SampleID %in% c("ZymoMock1a_S142", 
                         "ZymoMock1_S238",
                         "ZymoMock2_S78"))
cuso4_raw_read_counts_zymo
# SampleID        Num_Reads_Forward_Raw 
# ZymoMock1a_S142              80646780              
# ZymoMock1_S238               42456785              
# ZymoMock2_S78                66356203          


###Negative controls ####
cuso4_raw_read_counts_negative <- cuso4_raw_read_counts %>%
  filter(grepl('NTC|EB', SampleID))
cuso4_raw_read_counts_negative
# SampleID   Num_Reads_Forward_Raw
# EB2_S186                     495                   
# NTC10_S231                   164                   
# NTC11_S155                   128                   
# NTC1_S163                    523                  
# NTC2_S7                     1115                  
# NTC3_S180                   1277                  
# NTC4_S22                     518                   
# NTC5_S112                    123                  
# NTC6_S42                     416                   
# NTC7_S129                    299                   
# NTC8_S60                     146                   
# NTC9_S141                    116                   


#Filtering just phyloseq.samples SampleIDs, adding metadata
cuso4_raw_read_counts_samples_metadata <- cuso4_raw_read_counts %>%
  filter(SampleID %in% sampleIDs_phyloseq.samples)%>%
  dplyr::left_join(metadata, by = "SampleID")
nrow(cuso4_raw_read_counts_samples_metadata) #Ok, 204 samples

#Summary of raw reads
summary(cuso4_raw_read_counts_samples_metadata$Num_Reads_Forward_Raw)
# Min.      1st Qu.   Median    Mean      3rd Qu.  Max. 
# 17144696  36646154  47571222  55223521  66455354 237358233 

sd(cuso4_raw_read_counts_samples_metadata$Num_Reads_Forward_Raw)
#27,086,557

###Established vs Naive####
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T) 
sequencing_depth_P1vsH21
wilcox_test(cuso4_raw_read_counts_samples_metadata, Num_Reads_Forward_Raw ~ Enclosure)
# group1 group2    n1    n2 statistic       p
# H21    P1        87   117  4012        0.00980

###Established and Naive Over Time####
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
          breaks = c("1","53","65","104","141"),
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

##COMPARING TRIMMED READS#######
cuso4_trimmed_read_counts <- read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Read_counts/Trimmed/CuSO4_clean_trimmed_read_counts.csv')
str(cuso4_trimmed_read_counts)

##Trimmed read percentages
cuso4_trimmed_read_percentages <- cuso4_raw_read_counts %>%
  dplyr::left_join(metadata, by = "SampleID") %>%
  dplyr::left_join(cuso4_trimmed_read_counts, by = "SampleID")%>%
  mutate(Num_Trimmed_Reads = Num_Reads_Forward_Raw - Total_Counts_Paired_End_Trimmomatic_Output, 
         Percentage_Trimmed_Reads = Num_Trimmed_Reads/Num_Reads_Forward_Raw *100)
cuso4_trimmed_read_percentages

#Filtering just phyloseq.samples SampleIDs
cuso4_trimmed_read_counts_samples_metadata <- cuso4_trimmed_read_percentages %>%
  filter(SampleID %in% sampleIDs_phyloseq.samples)
nrow(cuso4_trimmed_read_counts_samples_metadata) #Ok, 204 samples
summary(cuso4_trimmed_read_counts_samples_metadata$Percentage_Trimmed_Reads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.106   1.531   1.842   1.821   2.079   3.764 
sd(cuso4_trimmed_read_counts_samples_metadata$Percentage_Trimmed_Reads)
#0.389698


#Filtering just zymos
cuso4_trimmed_read_counts_zymos_metadata <- cuso4_trimmed_read_percentages %>%
  filter(grepl("Zymo", SampleID))
nrow(cuso4_trimmed_read_counts_zymos_metadata) #Ok, 3 zymos
cuso4_trimmed_read_counts_zymos_metadata$Percentage_Trimmed_Reads


#Filtering just negative controls
cuso4_trimmed_read_counts_negative_metadata <- cuso4_trimmed_read_percentages %>%
  filter(grepl("NTC|EB", SampleID))
summary(cuso4_trimmed_read_counts_negative_metadata$Percentage_Trimmed_Reads) #Ok, 12 negative controls
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 24.14   57.45   66.66   65.57   82.79   93.31 
sd(cuso4_trimmed_read_counts_negative_metadata$Percentage_Trimmed_Reads)
#21.78944

###Established vs Naive####
trimmed_reads_P1vsH21<- ggplot(cuso4_trimmed_read_counts_samples_metadata, 
                               aes(x = Enclosure, y= Total_Counts_Paired_End_Trimmomatic_Output, 
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
trimmed_reads_P1vsH21
wilcox_test(cuso4_trimmed_read_counts_samples_metadata, 
            Total_Counts_Paired_End_Trimmomatic_Output ~ Enclosure)
# group1 group2    n1    n2 statistic       p
# H21    P1        87   117  4010        0.00967

###Established and Naive, Raw and Trimmed over time####
cuso4_trimmed_read_counts_samples_metadata_long <- 
  cuso4_trimmed_read_counts_samples_metadata %>%
  pivot_longer(cols = c(Total_Counts_Paired_End_Trimmomatic_Output,
                        Num_Reads_Forward_Raw),
               names_to = "Read_Status",
               values_to = "Num_Paired_Reads") %>%
  mutate(
    Read_Status = dplyr::recode(Read_Status,
                                "Total_Counts_Paired_End_Trimmomatic_Output" = "QC - Trimmomatic",
                                "Num_Reads_Forward_Raw" = "Raw (Sequencing Depth)"))%>%
  mutate(
    Read_Status = factor(Read_Status, levels = c("Raw (Sequencing Depth)",
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
  scale_y_continuous(expand = c(0.05, 0, 0.15, 0), labels = scales::comma) +
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
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
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

##COMPARING TRIMMED READS AND SEQUENCING DEPTH ON SAME PLOT#######
###Raw and Trimmed comparisons between systems (Established and Naive) ####
trimmed_and_raw_enclosure_comp <- ggplot(
  cuso4_trimmed_read_counts_samples_metadata_long,
  aes(x = Enclosure, y = Num_Paired_Reads,
      color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  facet_grid(~Read_Status, scales = "free_y") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand = c(0.05, 0, 0.15, 0), labels = scales::comma) +
  labs(y= "Paired-end Reads", color = "System", fill = "System") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22, face = "bold"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(colour = "white", size = 28, face = "bold"),
        axis.title.y = element_text(size = 28, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.text.x = element_blank(),
        axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.x = element_blank()) +
  geom_pwc(method = "wilcox_test",
           label = "Wilcoxon, p = {p.format}",
           tip.length = 0.02,
           label.size = 5,
           hide.ns = TRUE)
trimmed_and_raw_enclosure_comp

###Established and Naive, Raw vs Trimmed within same panel####
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
  scale_y_continuous(expand= c(0.01,0,0.1,0),labels = scales::comma) +
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

##COMPARING CLASSIFIED READ PERCENTAGES BY KRAKEN#######
kraken_ntcore_gtdb_kraken_unclassified

###Mock communities#########
kraken_unclassified_reads_pos_control <- kraken_ntcore_gtdb_kraken_unclassified %>%
  mutate(Percent_Reads_Classified_Kraken = 
           (ntcore_GTDB_classified_counts/Total_Counts_Paired_End_Trimmomatic_Output) *100)%>%
  filter(grepl("Zymo", SampleID))
nrow(kraken_unclassified_reads_pos_control) #Ok, 3 mock communities
kraken_unclassified_reads_pos_control$Percent_Reads_Classified_Kraken
#97.03094, 97.03889, 97.02522

###Samples########
kraken_unclassified_reads_samples_metadata <- kraken_ntcore_gtdb_kraken_unclassified %>%
  mutate(Percent_Reads_Classified_Kraken = 
           (ntcore_GTDB_classified_counts/Total_Counts_Paired_End_Trimmomatic_Output) *100)%>%
  merge(metadata, by = 'SampleID')%>%
  filter(SampleID %in% sampleIDs_phyloseq.samples)
nrow(kraken_unclassified_reads_samples_metadata) #204 samples

#Descriptive stats
summary(kraken_unclassified_reads_samples_metadata$Percent_Reads_Classified_Kraken)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 35.32   50.70   67.50   64.29   78.13   86.66 
sd(kraken_unclassified_reads_samples_metadata$Percent_Reads_Classified_Kraken) #14.81898

#Descriptive stats per group
kraken_unclassified_reads_samples_metadata %>%
  group_by(Enclosure)%>%
  summarise(mean_percentage_classified_reads = mean(Percent_Reads_Classified_Kraken),
            sd_percentage_classified_reads = sd(Percent_Reads_Classified_Kraken), 
            min_percentage_classified_reads = min(Percent_Reads_Classified_Kraken), 
            max_percentage_classified_reads = max(Percent_Reads_Classified_Kraken))
# Enclosure mean_percentage_classified_reads sd_percentage_classified_reads  min_percentage_classified_reads   max_percentage_classified_reads
# H21                                   51.0                           9.28                      35.3                          82.1
# P1                                    74.2                           9.41                      49.0                          86.7


####Kraken2 Classified Percentages Established vs Naive####
kraken2_classified_read_percentages_P1vsH21<- ggplot(kraken_unclassified_reads_samples_metadata, 
                                                     aes(x = Enclosure, 
                                                         y= Percent_Reads_Classified_Kraken, 
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
kraken2_classified_read_percentages_P1vsH21
wilcox_test(kraken_unclassified_reads_samples_metadata, Percent_Reads_Classified_Kraken ~ Enclosure) 
# group1 group2    n1    n2 statistic        p
# H21    P1        87   117       649 1.79e-26


####Kraken2 Classified Percentages Established and Naive over time#### 
kraken2_classified_read_percentages_P1andH21_overtime <- ggplot(kraken_unclassified_reads_samples_metadata,
                                                                aes(x = factor(Date_num), 
                                                                    y= Percent_Reads_Classified_Kraken, 
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
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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


##COMPARING SAMPLE SUMS OTUs (CLASSIFIED READ COUNTS FROM KRAKEN)#######
####ALL TAXA#####
kraken_unclassified_reads_samples_metadata$ntcore_GTDB_classified_counts

#Descriptive stats
summary(kraken_unclassified_reads_samples_metadata$ntcore_GTDB_classified_counts)
# Min.      1st Qu.    Median      Mean   3rd Qu.      Max. 
# 10985229  21083821  29012127  35360498  45662805 137049923 
sd(kraken_unclassified_reads_samples_metadata$ntcore_GTDB_classified_counts) #20,007,809

#Descriptive stats per group
kraken_unclassified_reads_samples_metadata %>%
  group_by(Enclosure)%>%
  summarise(mean_OTU_counts = mean(ntcore_GTDB_classified_counts),
            sd_OTU_counts = sd(ntcore_GTDB_classified_counts), 
            min_OTU_counts = min(ntcore_GTDB_classified_counts), 
            max_OTU_counts = max(ntcore_GTDB_classified_counts))%>%
  as.data.frame() %>%
  print(digits = 12) 
# Enclosure mean_OTU_counts   sd_OTU_counts    min_OTU_counts   max_OTU_counts
# H21      25778185.8046      16410250.3527       10985229      137049923
# P1       42485806.2479      19514680.9041       13769517      109551194

#Samplesums of Zymo mock communities
sample_sums(phyloseq.controls)
# ZymoMock1a_S142  ZymoMock1_S238  ZymoMock2_S78 
#      76530152        40440500    62900637 

#Samplesums of negative controls
sample_sums(phyloseq.controls)
# EB2_S186      NTC10_S231      NTC11_S155       NTC1_S163         NTC2_S7       NTC3_S180        NTC4_S22       NTC5_S112        NTC6_S42 
# 38              31              51              24              74             321              73              36             220 
# NTC7_S129        NTC8_S60       NTC9_S141
# 82              44              70 

#####Established vs Naive####
all_taxa_samplesums_P1vsH21<- ggplot(kraken_unclassified_reads_samples_metadata, 
                                             aes(x = Enclosure, y= ntcore_GTDB_classified_counts, 
                                                 color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Classified Reads", color = "System", fill = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0), labels = scales::comma) +
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
all_taxa_samplesums_P1vsH21

##Stats 
wilcox_test(kraken_unclassified_reads_samples_metadata, ntcore_GTDB_classified_counts~Enclosure) 
# group1 group2    n1    n2 statistic        p
# H21    P1        87   117      1986 9.98e-14

#####Established and Naive over time####
all_taxa_samplesums_P1andH21_overtime<- ggplot(kraken_unclassified_reads_samples_metadata, 
                                                       aes(x = factor(Date_num), 
                                                           y= ntcore_GTDB_classified_counts, 
                                                           color = Enclosure)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Classified Reads", x = "Day", color = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.05,0,0.1,0), labels = scales::comma) +
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
          breaks = c("1","53","65","104","141"),
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
all_taxa_samplesums_P1andH21_overtime

#####H21 - Naive system####
phyloseq.samples_H21 <- subset_samples(phyloseq.samples, Enclosure == "H21")
phyloseq.samples_H21 <- prune_taxa(taxa_sums(phyloseq.samples_H21) > 0, 
                                            phyloseq.samples_H21)
phyloseq.samples_H21 #287,870 taxa and 87 samples
range(phyloseq.samples_H21@sam_data$Collection_Date)#OK, "2023-10-09" through "2024-02-20"

#####P1 - Established system####
phyloseq.samples_P1 <- subset_samples(phyloseq.samples, Enclosure == "P1")
phyloseq.samples_P1 <- prune_taxa(taxa_sums(phyloseq.samples_P1) > 0, 
                                           phyloseq.samples_P1)
phyloseq.samples_P1 #264,693 taxa and 117 samples
range(phyloseq.samples_P1@sam_data$Collection_Date)#OK, "2023-11-14" through "2024-04-02"


####NITRIFYING TAXA####
phyloseq.samples.nitrifiers <- subset_taxa(phyloseq.samples, 
                              #AOB
                              Family == "Nitrosomonadaceae" | # AOB; Includes genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                              Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                              #AOA
                                Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). InGTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes genera Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                              #NOB
                                #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                Family == "NS-4" | # NOB (mostly uncultured genomes)
                                Family == "UBA8639" | # NOB
                                #In the order Rhizobiales and family Xanthobacteraceae
                                Genus == "Nitrobacter" | #NOB; most representative one
                                #Order Burkholderiales, family Gallionellaceae
                                Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                #Order Nitrococcales, family Nitrococcaceae
                                Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                              # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                              Genus == "Nitrolancea"  #NOB; new one
                              ) 

phyloseq.samples.nitrifiers #1690 taxa and 204 samples
phyloseq.samples.nitrifiers <- subset_samples(phyloseq.samples.nitrifiers, 
                             sample_sums(phyloseq.samples.nitrifiers) > 0)
phyloseq.samples.nitrifiers #1690 taxa and 204 samples 


#####H21 - Naive system####
phyloseq.samples.nitrifiers_H21 <- subset_samples(phyloseq.samples.nitrifiers, Enclosure == "H21")
phyloseq.samples.nitrifiers_H21 <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers_H21) > 0, 
                                              phyloseq.samples.nitrifiers_H21)
phyloseq.samples.nitrifiers_H21 # 1686 taxa and 87 samples

#####P1 - Established system####
phyloseq.samples.nitrifiers_P1 <- subset_samples(phyloseq.samples.nitrifiers, Enclosure == "P1")
phyloseq.samples.nitrifiers_P1 <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers_P1) > 0, 
                                             phyloseq.samples.nitrifiers_P1)
phyloseq.samples.nitrifiers_P1 # 1685 taxa and 87 samples


#####Sample sums
sample.sums.nit <- sample_sums(phyloseq.samples.nitrifiers) #making a sample sums object
nitrifiers.samplesums.df <- cbind(phyloseq.samples.nitrifiers@sam_data, 
                                  sample.sums.nit) #combining sample sums with metadata
nitrifiers.samplesums.df
nitrifiers.samplesums.df$SampleID <- rownames(nitrifiers.samplesums.df) ##making a sampleID column

#####Established vs Naive####
nitrifier_bacteria_archaea_samplesums_P1vsH21<- ggplot(nitrifiers.samplesums.df, 
                                                       aes(x = Enclosure, y= sample.sums.nit, 
                                                           color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Nitrifier Classified Reads", color = "System", fill = "System",) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0), labels = scales::comma) +
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
nitrifier_bacteria_archaea_samplesums_P1vsH21
##Stats 
wilcox_test(nitrifiers.samplesums.df, sample.sums.nit~Enclosure) 
# group1 group2    n1    n2 statistic         p
# H21    P1        87   117      7032 0.00000321

#####Established and Naive over time####
nitrifiers_samplesums_P1andH21_overtime<- ggplot(nitrifiers.samplesums.df, 
                                                       aes(x = factor(Date_num), 
                                                           y= sample.sums.nit, 
                                                           color = Enclosure)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Nitrifier Classified Reads", x = "Day", color = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.05,0,0.1,0), labels = scales::comma) +
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
          breaks = c("1","53","65","104","141"),
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
nitrifiers_samplesums_P1andH21_overtime


####BACTERIAL - ARCHAEAL TAXA####
phyloseq.samples.bacteria_1 <- subset_taxa(phyloseq.samples, Domain %in% c("Bacteria", "Archaea"))
phyloseq.samples.bacteria_1 # 219,873 taxa and 204 samples
phyloseq.samples.bacteria <- prune_samples(sample_sums(phyloseq.samples.bacteria_1) > 0, phyloseq.samples.bacteria_1)
phyloseq.samples.bacteria #219,873 taxa and 204 samples

# some QC checks of the "classified" reads per samples
min(sample_sums(phyloseq.samples.bacteria)) # 10,652,712
max(sample_sums(phyloseq.samples.bacteria)) # 135,834,095  (H21_0119) 
mean(sample_sums(phyloseq.samples.bacteria)) # 35,056,339
median(sample_sums(phyloseq.samples.bacteria)) # 28,864,604
sort(sample_sums(phyloseq.samples.bacteria))

#####H21 - Naive system####
phyloseq.samples.bacteria.H21 <- subset_samples(phyloseq.samples.bacteria, Enclosure == "H21")
phyloseq.samples.bacteria.H21 <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.H21) > 0, 
                                            phyloseq.samples.bacteria.H21)
phyloseq.samples.bacteria.H21 # 218,414 taxa and 87 samples

#####P1 - Established system####
phyloseq.samples.bacteria.P1 <- subset_samples(phyloseq.samples.bacteria, Enclosure == "P1")
phyloseq.samples.bacteria.P1 <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.P1) > 0, 
                                            phyloseq.samples.bacteria.P1)
phyloseq.samples.bacteria.P1 # 217,728 taxa and 117 samples

#####Sample sums
sample.sums.bacteria <- sample_sums(phyloseq.samples.bacteria) #making a sample sums object
bacteria.samplesums.df <- cbind(phyloseq.samples.bacteria@sam_data, 
                                sample.sums.bacteria) #combining sample sums with metadata
bacteria.samplesums.df

#####Established vs Naive####
bacteria_archaea_samplesums_P1vsH21<- ggplot(bacteria.samplesums.df, 
                                             aes(x = Enclosure, y= sample.sums.bacteria, 
                                                 color = Enclosure, fill = Enclosure)) +
  theme_bw() +
  labs(y= "Bacterial - Archaeal Classified Reads", color = "System", fill = "System",) +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_color_manual(values = enclosure.palette, labels = c("H21" = "Naive", "P1" = "Established"))+
  scale_y_continuous(expand= c(0.05,0,0.1,0), labels = scales::comma) +
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T)
bacteria_archaea_samplesums_P1vsH21
##Stats 
wilcox_test(bacteria.samplesums.df, sample.sums.bacteria~Enclosure) 
# group1 group2    n1    n2 statistic        p
# H21    P1        87   117      1955 5.67e-14

#####Established and Naive over time####
bacteria_archaea_samplesums_P1andH21_overtime<- ggplot(bacteria.samplesums.df, 
                                                 aes(x = factor(Date_num), 
                                                     y= sample.sums.bacteria, 
                                                     color = Enclosure)) +
  theme_bw() +
  facet_grid(~Enclosure, 
             scales = "free", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  labs(y= "Bacterial - Archaeal Classified Reads", x = "Day", color = "System") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8) +
  scale_y_continuous(expand= c(0.05,0,0.1,0), labels = scales::comma) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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


#TAX GLOMMING AND RELATIVE ABUNDANCES####
##Phyloseq without unclassified reads - all taxa########
phyloseq.samples #309,747 taxa and 204 samples
###Phylum#####
phyloseq.samples.phylum <- tax_glom(phyloseq.samples, taxrank = "Phylum", NArm = F) 
phyloseq.samples.phylum #2754 Phyla and 204 samples
#Relative abundance
phyloseq.samples.phylum.ra <- transform_sample_counts(phyloseq.samples.phylum,
                                                      function(x) x/sum(x)*100) ##Relative abundance 
###Class#####
phyloseq.samples.class <- tax_glom(phyloseq.samples, taxrank = "Class", NArm = F) 
phyloseq.samples.class #3932 Classes and 204 samples
#Relative abundance
phyloseq.samples.class.ra <- transform_sample_counts(phyloseq.samples.class,
                                                      function(x) x/sum(x)*100) ##Relative abundance 
###Order#####
phyloseq.samples.order <- tax_glom(phyloseq.samples, taxrank = "Order", NArm = F) 
phyloseq.samples.order #11,708 orders and 204 samples
#Relative abundance
phyloseq.samples.order.ra <- transform_sample_counts(phyloseq.samples.order,
                                                     function(x) x/sum(x)*100) ##Relative abundance 
###Family#####
phyloseq.samples.family <- tax_glom(phyloseq.samples, taxrank = "Family", NArm = F) 
phyloseq.samples.family #23,426 families and 204 samples
#Relative abundance
phyloseq.samples.family.ra <- transform_sample_counts(phyloseq.samples.family,
                                                     function(x) x/sum(x)*100) ##Relative abundance 
###Genus (will take a while)####
phyloseq.samples.genus <- tax_glom(phyloseq.samples, taxrank = "Genus", NArm = F) 
phyloseq.samples.genus #89899 genera and 204 samples
#Relative abundance
phyloseq.samples.genus.ra <- transform_sample_counts(phyloseq.samples.genus,
                                                     function(x) x/sum(x)*100) ##Relative abundance
phyloseq.samples.genus.ra #89,899 genera and 204 samples 

###Species#####
phyloseq.samples.species <- phyloseq.samples #Lowest rank is species
#Relative abundance
phyloseq.samples.species.ra <- transform_sample_counts(phyloseq.samples.species,
                                                     function(x) x/sum(x)*100) ##Relative abundance 


##Phyloseq without unclassified reads - just bacterial and archaeal taxa########
phyloseq.samples.bacteria #219,873 taxa and 204 samples

###Phylum#####
phyloseq.samples.bacteria.phylum <- tax_glom(phyloseq.samples.bacteria, taxrank = "Phylum", NArm = F) 
phyloseq.samples.bacteria.phylum #199 phyla and 204 samples
#Relative abundance
phyloseq.samples.bacteria.phylum.ra <- transform_sample_counts(phyloseq.samples.bacteria.phylum,
                                                      function(x) x/sum(x)*100) ##Relative abundance from normalized data
###Class#####
phyloseq.samples.bacteria.class <- tax_glom(phyloseq.samples.bacteria, taxrank = "Class", NArm = F) 
phyloseq.samples.bacteria.class # 696 classes and 204 samples 
#Relative abundance
phyloseq.samples.bacteria.class.ra <- transform_sample_counts(phyloseq.samples.bacteria.class,
                                                     function(x) x/sum(x)*100) ##Relative abundance from normalized data
phyloseq.samples.bacteria.class.ra
###Order#####
phyloseq.samples.bacteria.order <- tax_glom(phyloseq.samples.bacteria, taxrank = "Order", NArm = F) 
phyloseq.samples.bacteria.order #2583 orders and 204 samples 
#Relative abundance
phyloseq.samples.bacteria.order.ra <- transform_sample_counts(phyloseq.samples.bacteria.order,
                                                     function(x) x/sum(x)*100) ##Relative abundance from normalized data
phyloseq.samples.bacteria.order.ra
###Family#####
phyloseq.samples.bacteria.family <- tax_glom(phyloseq.samples.bacteria, taxrank = "Family", NArm = F) 
phyloseq.samples.bacteria.family #7596 families and 204 samples
#Relative abundance
phyloseq.samples.bacteria.family.ra <- transform_sample_counts(phyloseq.samples.bacteria.family,
                                                      function(x) x/sum(x)*100) ##Relative abundance from normalized data
###Genus (will take a while)#####
phyloseq.samples.bacteria.genus <- tax_glom(phyloseq.samples.bacteria, taxrank = "Genus", NArm = F) 
phyloseq.samples.bacteria.genus #41,284 genera and 204 samples 
#Relative abundance
phyloseq.samples.bacteria.genus.ra <- transform_sample_counts(phyloseq.samples.bacteria.genus,
                                                     function(x) x/sum(x)*100) ##Relative abundance 

###Species#####
phyloseq.samples.bacteria.species <- phyloseq.samples.bacteria #Lowest rank is species
phyloseq.samples.bacteria.species #219873 species and 204 samples
#Relative abundance
phyloseq.samples.bacteria.species.ra <- transform_sample_counts(phyloseq.samples.bacteria.species,
                                                       function(x) x/sum(x)*100) ##Relative abundance 

##Phyloseq with unclassified reads########
phyloseq_with_unclassifieds_otu.samples #309,748 taxa and 204 samples

###Phylum#####
phyloseq_with_unclassifieds_otu.samples.phylum <- tax_glom(phyloseq_with_unclassifieds_otu.samples, taxrank = "Phylum", NArm = F) 
phyloseq_with_unclassifieds_otu.samples.phylum #2,755 phyla and 204 samples
#Relative abundance
phyloseq_with_unclassifieds_otu.samples.phylum.ra <- transform_sample_counts(phyloseq_with_unclassifieds_otu.samples.phylum,
                                                      function(x) x/sum(x)*100) ##Relative abundance from normalized data
###Class#####
phyloseq_with_unclassifieds_otu.samples.class <- tax_glom(phyloseq_with_unclassifieds_otu.samples, taxrank = "Class", NArm = F) 
phyloseq_with_unclassifieds_otu.samples.class #3933 classes and 204 samples
#Relative abundance
phyloseq_with_unclassifieds_otu.samples.class.ra <- transform_sample_counts(phyloseq_with_unclassifieds_otu.samples.class,
                                                     function(x) x/sum(x)*100) ##Relative abundance from normalized data
###Order#####
phyloseq_with_unclassifieds_otu.samples.order <- tax_glom(phyloseq_with_unclassifieds_otu.samples, taxrank = "Order", NArm = F) 
phyloseq_with_unclassifieds_otu.samples.order #11709 orders and 204 samples
#Relative abundance
phyloseq_with_unclassifieds_otu.samples.order.ra <- transform_sample_counts(phyloseq_with_unclassifieds_otu.samples.order,
                                                     function(x) x/sum(x)*100) ##Relative abundance from normalized data

###Family#####
phyloseq_with_unclassifieds_otu.samples.family <- tax_glom(phyloseq_with_unclassifieds_otu.samples, taxrank = "Family", NArm = F) 
phyloseq_with_unclassifieds_otu.samples.family #23,427 taxa and 204 samples
#Relative abundance
phyloseq_with_unclassifieds_otu.samples.family.ra <- transform_sample_counts(phyloseq_with_unclassifieds_otu.samples.family,
                                                      function(x) x/sum(x)*100) ##Relative abundance from normalized data
###Genus (will take a while)#####
phyloseq_with_unclassifieds_otu.samples.genus <- tax_glom(phyloseq_with_unclassifieds_otu.samples, taxrank = "Genus", NArm = F) 
phyloseq_with_unclassifieds_otu.samples.genus #89,900 genera and 204 samples
#Relative abundance
phyloseq_with_unclassifieds_otu.samples.genus.ra <- transform_sample_counts(phyloseq_with_unclassifieds_otu.samples.genus,
                                                     function(x) x/sum(x)*100) ##Relative abundance from normalized data
phyloseq_with_unclassifieds_otu.samples.genus.ra #89,900 genera and 204 samples

###Species#####
phyloseq_with_unclassifieds_otu.samples.species <- phyloseq_with_unclassifieds_otu.samples #Lowest rank is species
#Relative abundance
phyloseq_with_unclassifieds_otu.samples.species.ra <- transform_sample_counts(phyloseq_with_unclassifieds_otu.samples.species,
                                                                function(x) x/sum(x)*100) ##Relative abundance 


##CLASSIFICATION PERCENTAGES AT DIFFERENT TAXONOMIC LEVELS - BACTERIAL AND ARCHAEAL TAXA ####
###PHYLUM######
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.phylum.ra)[, "Phylum"])) #199 phyla (so No duplicates)

Unclassified_phylum_abundance <- phyloseq.samples.bacteria.phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##12% abundance by Unclassified Phyla

Classified_phylum_abundance <- phyloseq.samples.bacteria.phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##88% abundance by Classified Phyla

##Checking on excel
write.csv(phyloseq.samples.bacteria.phylum.ra@otu_table, "phylum_otus.csv")
write.csv(phyloseq.samples.bacteria.phylum.ra@tax_table, "phylum_taxa.csv")  

#Just classified Phyla
phyloseq.samples.bacteria.phylum.classified.ra <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.phylum.ra)[, "Phylum"]),
  phyloseq.samples.bacteria.phylum.ra)
phyloseq.samples.bacteria.phylum.classified.ra ##197 classified Phyla

#Just classified Phyla
phyloseq.samples.bacteria.phylum.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.phylum.ra)[, "Phylum"]),
  phyloseq.samples.bacteria.phylum.ra)
phyloseq.samples.bacteria.phylum.unclassified.ra ##2 unclassified Phyla

###CLASS#####
#Are there duplicates? 
phyloseq.samples.bacteria.class.ra
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.class.ra)[, "Class"])) #696 classes (no duplicates)

Unclassified_class_abundance <- phyloseq.samples.bacteria.class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##14.5% abundance by Unclassified Classes

Classified_class_abundance <- phyloseq.samples.bacteria.class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##85.5% abundance by Classified Phyla

##Checking on excel
write.csv(phyloseq.samples.bacteria.class.ra@otu_table, "class_otus.csv")
write.csv(phyloseq.samples.bacteria.class.ra@tax_table, "class_taxa.csv")  

#Just unclassified Classes
phyloseq.samples.bacteria.class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.class.ra)[, "Class"]),
  phyloseq.samples.bacteria.class.ra)
phyloseq.samples.bacteria.class.unclassified.ra ##51 classified Classes

#Just classified Classes
phyloseq.samples.bacteria.class.classified.ra <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.class.ra)[, "Class"]),
  phyloseq.samples.bacteria.class.ra)
phyloseq.samples.bacteria.class.classified.ra #645 classified Classes

###ORDER######
#Are there duplicates? 
phyloseq.samples.bacteria.order.ra
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.order.ra)[, "Order"])) #2583 orders (so No duplicates)

Unclassified_order_abundance <- phyloseq.samples.bacteria.order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unorderified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##20.5% abundance by Unclassified orders

Classified_order_abundance <- phyloseq.samples.bacteria.order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified", Order, ignore.case = TRUE)) %>%  ##filter out unorderified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##79.5% abundance by Classified orders

##Checking on excel
write.csv(phyloseq.samples.bacteria.order.ra@otu_table, "order_otus.csv")
write.csv(phyloseq.samples.bacteria.order.ra@tax_table, "order_taxa.csv")  

#Just unclassified Orders
phyloseq.samples.bacteria.order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.order.ra)[, "Order"]),
  phyloseq.samples.bacteria.order.ra)
phyloseq.samples.bacteria.order.unclassified.ra ##232 unclassified orders

#Just classified Orders
phyloseq.samples.bacteria.order.classified.ra <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.order.ra)[, "Order"]),
  phyloseq.samples.bacteria.order.ra)
phyloseq.samples.bacteria.order.classified.ra ##2351 classified orders

###FAMILY######
#Are there duplicates? 
phyloseq.samples.bacteria.family.ra #7596
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.family.ra)[, "Family"])) #7575 taxa (there are duplicates; 7596 total)
family_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.samples.bacteria.family.ra)[, "Family"])
unique(family_taxa_vec[duplicated(family_taxa_vec)]) #Duplicates:
# "unclassified J058"         "unclassified UBA8108"      "unclassified UBA9160"      "unclassified UBA796"       "unclassified UBA1018"     
# "unclassified UBA2968"      "unclassified UBA2361"      "unclassified RBG-16-71-46" "unclassified WGA-4E"       "unclassified UBA10199"    
# "unclassified GWC2-55-46"   "unclassified UBA11872"     "unclassified SZUA-567"     "unclassified UBA6077"      "unclassified MSB-5A5"     
# "unclassified 4484-113"     "unclassified SAR324"       "unclassified UBA8248"      "unclassified UBA5829"      "unclassified UBA7883"     
# "unclassified SBBH01"    

Unclassified_family_abundance <- phyloseq.samples.bacteria.family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##22.3% abundance by Unclassified Families

Classified_family_abundance <- phyloseq.samples.bacteria.family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##77.7% abundance by Classified Families

#Checking on excel
write.csv(phyloseq.samples.bacteria.family.ra@otu_table, "family_otus.csv")
write.csv(phyloseq.samples.bacteria.family.ra@tax_table, "family_taxa.csv") 

#How many unclassified?
phyloseq.samples.bacteria.family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.family.ra)[, "Family"]),
  phyloseq.samples.bacteria.family.ra)
phyloseq.samples.bacteria.family.unclassified.ra #849 unclassified families
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.family.unclassified.ra)[, "Family"])) ##828 unclassified families (unique - without duplicates)

#Just classified Families
phyloseq.samples.bacteria.family.classified.ra <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.family.ra)[, "Family"]),
  phyloseq.samples.bacteria.family.ra)
phyloseq.samples.bacteria.family.classified.ra #6747 classified families
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.family.classified.ra)[, "Family"]))#6747 classified families (unique - without duplicates)

###GENUS ######
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.genus.ra)[, "Genus"])) #41031 taxa (there are duplicates; total 41284)
genus_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.samples.bacteria.genus.ra)[, "Genus"])
unique(genus_taxa_vec[duplicated(genus_taxa_vec)]) #237 duplicated 

Unclassified_genus_abundance <- phyloseq.samples.bacteria.genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance)) #Sum across OTUs
Unclassified_genus_abundance ##30.7% abundance by unclassified genera

Classified_genus_abundance <- phyloseq.samples.bacteria.genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##69.3% abundance by Classified Genera

#Checking on excel
write.csv(phyloseq.samples.bacteria.genus.ra@otu_table, "genus_otus.csv")
write.csv(phyloseq.samples.bacteria.genus.ra@tax_table, "genus_taxa.csv") 

#How many unclassified?
phyloseq.samples.bacteria.genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.genus.ra)[, "Genus"]),
  phyloseq.samples.bacteria.genus.ra)
phyloseq.samples.bacteria.genus.unclassified.ra #3781 unclassified genera
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.genus.unclassified.ra)[, "Genus"])) ##3528 unclassified genera (unique - without duplicates)

#Keep just classified Genera
phyloseq.samples.bacteria.genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.genus.ra)[, "Genus"]),
  phyloseq.samples.bacteria.genus.ra)
phyloseq.samples.bacteria.genus.classified.ra #37503 classified genera
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.genus.classified.ra)[, "Genus"])) ##37503 classified genera (unique - without duplicates)


###SPECIES######
phyloseq.samples.bacteria.species.ra ##219,873  Species- OTUs

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.species.ra)[, "Species"])) #218,522 species (there are duplicates)
species_taxa_vec <- as.character(phyloseq::tax_table(phyloseq.samples.bacteria.species.ra)[, "Species"])
unique(species_taxa_vec[duplicated(species_taxa_vec)]) #duplicated unique ones

#Extract OTU and tax table
species_level_otu <- as(phyloseq::otu_table(phyloseq.samples.bacteria.species.ra), "matrix")
species_level_tax <- as(phyloseq::tax_table(phyloseq.samples.bacteria.species.ra), "matrix")

# Identify unclassified and classified species
unclassified_idx_species <- grepl("unclassified", species_level_tax[, "Species"], ignore.case = TRUE)
classified_idx_species <- !grepl("unclassified|unknown", species_level_tax[, "Species"], ignore.case = TRUE)

# Subset OTUs
otu_unclassified_species <- species_level_otu[unclassified_idx_species, , drop = FALSE]
otu_classified_species <- species_level_otu[classified_idx_species, , drop = FALSE]

# Mean abundance per OTU (row means)
otu_unclassified_species_means <- rowMeans(otu_unclassified_species)
otu_classified_species_means <- rowMeans(otu_classified_species)

# Sum across OTUs
unclassified_species_sum <- sum(otu_unclassified_species_means) #35.20% abundance by unclassified species
classified_species_sum <- sum(otu_classified_species_means) #64.8% abundance by classified species

#Checking on excel
write.csv(phyloseq.samples.bacteria.species.ra@otu_table, "species_otus.csv")
write.csv(phyloseq.samples.bacteria.species.ra@tax_table, "species_taxa.csv") 

#How many unclassified?
phyloseq.samples.bacteria.species.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.species.ra)[, "Species"]),
  phyloseq.samples.bacteria.species.ra)
phyloseq.samples.bacteria.species.unclassified.ra #20,273 unclassified species
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.species.unclassified.ra)[, "Species"])) ##18,922 unclassified species (unique - without duplicates)

#Keep just classified Genera
phyloseq.samples.bacteria.species.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq.samples.bacteria.species.ra)[, "Species"]),
  phyloseq.samples.bacteria.species.ra)
phyloseq.samples.bacteria.species.classified.ra #199,600 classified species
length(unique(phyloseq::tax_table(phyloseq.samples.bacteria.species.classified.ra)[, "Species"])) ##199,600 classified species (unique - without duplicates)



##CLASSIFICATION PERCENTAGES AT DIFFERENT TAXONOMIC LEVELS - BACTERIAL AND ARCHAEAL TAXA WITHIN ALL COMMUNITIES (INCLUDING UNCLASSIFIED READS) ####
###PHYLUM######
phyloseq_with_unclassifieds_otu.samples.phylum.ra #Relative abundances with all taxa and unclassified reads
#Pick just bacterial and archaeal taxa
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.phylum.ra, 
                                                                     Domain %in% c("Bacteria", "Archaea"))
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc # 199 phyla and 204 samples
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc) > 0, 
                                                                          phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc #199 phyla and 204 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc)[, "Phylum"])) #199 phyla (so No duplicates)

Unclassified_phylum_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance_withunclassifieds ##6.41% abundance by Unclassified Phyla

Classified_phylum_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc %>%
  psmelt()%>%
  filter(!grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance_withunclassifieds ##57.3% abundance by Classified Phyla

#6.41% + 57.3% = 63.71%. The rest are unclassified reads, eukaryotic, and viral taxa

#Just classified Phyla
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc_classified <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc)[, "Phylum"]),
  phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc_classified ##197 classified Phyla

#Just unclassified Phyla
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc_unclassified <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc)[, "Phylum"]),
  phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.phylum.ra_bacarc_unclassified ##2 unclassified Phyla 

###CLASS######
phyloseq_with_unclassifieds_otu.samples.class.ra #Relative abundances with all taxa and unclassified reads
#Pick just bacterial and archaeal taxa
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.class.ra, 
                                                                        Domain %in% c("Bacteria", "Archaea"))
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc # 696 Classes and 204 samples
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc) > 0, 
                                                                          phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc #696 Classes and 204 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc)[, "Class"])) #696 Classes (so No duplicates)

Unclassified_class_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified Classes
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance_withunclassifieds ##7.74% abundance by Unclassified Classes

Classified_class_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc %>%
  psmelt()%>%
  filter(!grepl("unclassified", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance_withunclassifieds ##56.0% abundance by Classified Classes

#7.74% + 56.0% = 63.74%. The rest are unclassified reads, eukaryotic, and viral taxa

#Just classified Classes
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc_classified <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc)[, "Class"]),
  phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc_classified ##645 classified Classes

#Just unclassified Classes
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc_unclassified <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc)[, "Class"]),
  phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.class.ra_bacarc_unclassified ##51 unclassified Classes 


###ORDER######
phyloseq_with_unclassifieds_otu.samples.order.ra #Relative abundances with all taxa and unclassified reads
#Pick just bacterial and archaeal taxa
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.order.ra, 
                                                                        Domain %in% c("Bacteria", "Archaea"))
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc # 2583 orders and 204 samples
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc) > 0, 
                                                                          phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc #2583 orders and 204 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc)[, "Order"])) #2583 orders (so No duplicates)

Unclassified_order_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified orders
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance_withunclassifieds ##11.2% abundance by Unclassified Orders

Classified_order_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc %>%
  psmelt()%>%
  filter(!grepl("unclassified", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance_withunclassifieds ##52.5% abundance by Classified Orders

#11.21% + 52.5% = 63.7%. The rest are unclassified reads, eukaryotic, and viral taxa

#Just classified Orders
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc_classified <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc)[, "Order"]),
  phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc_classified ##2351 classified Orders

#Just unclassified Orders
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc_unclassified <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc)[, "Order"]),
  phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.order.ra_bacarc_unclassified ##232 unclassified Orders 


###FAMILY######
phyloseq_with_unclassifieds_otu.samples.family.ra #Relative abundances with all taxa and unclassified reads
#Pick just bacterial and archaeal taxa
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.family.ra, 
                                                                        Domain %in% c("Bacteria", "Archaea"))
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc # 7596 Families and 204 samples
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc) > 0, 
                                                                          phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc #7596 Families and 204 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)[, "Family"])) #7575 Families (so there are duplicates)
family_taxa_vec_unc <- as.character(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)[, "Family"])
unique(family_taxa_vec_unc[duplicated(family_taxa_vec_unc)]) #Duplicates:
# "unclassified J058"         "unclassified UBA8108"      "unclassified UBA9160"      "unclassified UBA796"       "unclassified UBA1018"     
# "unclassified UBA2968"      "unclassified UBA2361"      "unclassified RBG-16-71-46" "unclassified WGA-4E"       "unclassified UBA10199"    
# "unclassified GWC2-55-46"   "unclassified UBA11872"     "unclassified SZUA-567"     "unclassified UBA6077"      "unclassified MSB-5A5"     
# "unclassified 4484-113"     "unclassified SAR324"       "unclassified UBA8248"      "unclassified UBA5829"      "unclassified UBA7883"     
# "unclassified SBBH01"  

Unclassified_family_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unclassified Families
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance_withunclassifieds ##12.1% abundance by Unclassified Families

Classified_family_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc %>%
  psmelt()%>%
  filter(!grepl("unclassified", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance_withunclassifieds ##51.6% abundance by Classified Families

#12.1% + 51.6% = 63.71%. The rest are unclassified reads, eukaryotic, and viral taxa

#Just classified Families
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc_classified <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)[, "Family"]),
  phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc_classified ##6747 classified Families

#Just unclassified Families
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc_unclassified <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)[, "Family"]),
  phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.family.ra_bacarc_unclassified ##849 unclassified Families 


###GENUS######
phyloseq_with_unclassifieds_otu.samples.genus.ra #Relative abundances with all taxa and unclassified reads
#Pick just bacterial and archaeal taxa
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.genus.ra, 
                                                                        Domain %in% c("Bacteria", "Archaea"))
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc # 41284 Genera and 204 samples
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc) > 0, 
                                                                          phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc #41284 Genera and 204 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)[, "Genus"])) #41031 Genera (so there are duplicates)
genus_taxa_vec_unc <- as.character(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)[, "Genus"])
unique(genus_taxa_vec_unc[duplicated(genus_taxa_vec_unc)]) #237 duplicates   
  

Unclassified_genus_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unclassified Genera
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_genus_abundance_withunclassifieds ##% abundance by Unclassified Genera

Classified_genus_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc %>%
  psmelt()%>%
  filter(!grepl("unclassified", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance_withunclassifieds ##% abundance by Classified Genera

#% + % = 63.8%. The rest are unclassified reads, eukaryotic, and viral taxa

#Just classified Genera
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc_classified <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)[, "Genus"]),
  phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc_classified ##37503 classified Genera

#Just unclassified Genera
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc_unclassified <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)[, "Genus"]),
  phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.genus.ra_bacarc_unclassified ##3781 unclassified Genera 

###SPECIES######
phyloseq_with_unclassifieds_otu.samples.species.ra #Relative abundances with all taxa and unclassified reads
#Pick just bacterial and archaeal taxa
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.species.ra, 
                                                                       Domain %in% c("Bacteria", "Archaea"))
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc # 219873 Species and 204 samples
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc <- prune_samples(sample_sums(phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc) > 0, 
                                                                         phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc #219873 Species and 204 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)[, "Species"])) #218522 Species (so there are duplicates)
species_taxa_vec_unc <- as.character(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)[, "Species"])
unique(species_taxa_vec_unc[duplicated(species_taxa_vec_unc)]) #Duplicates

Unclassified_species_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc %>%
  psmelt()%>%
  filter(grepl("unclassified", Species, ignore.case = TRUE)) %>%  # Filter unclassified Species
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_species_abundance_withunclassifieds ##% abundance by Unclassified Species

Classified_species_abundance_withunclassifieds <- phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc %>%
  psmelt()%>%
  filter(!grepl("unclassified", Species, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_species_abundance_withunclassifieds ##% abundance by Classified Species

#% + % = 63.71%. The rest are unclassified reads, eukaryotic, and viral taxa

#Just classified Species
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc_classified <- prune_taxa(
  !grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)[, "Species"]),
  phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc_classified ## 199600 classified Species

#Just unclassified Species
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc_unclassified <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)[, "Species"]),
  phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc)
phyloseq_with_unclassifieds_otu.samples.species.ra_bacarc_unclassified ##20273 unclassified Species 


###SUPPLEMENTARY TABLE 2####
stable2 <- data.frame(
  "Taxonomic level" = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  # Total taxa at each level
  "Number of Taxa" = c(
    length(taxa_names(phyloseq.samples.bacteria.phylum.ra)),
    length(taxa_names(phyloseq.samples.bacteria.class.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.order.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.family.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.genus.ra)),
    length(taxa_names(phyloseq.samples.bacteria.species.ra))),

  
  # Unclassified taxa
  "Number of Unclassified Taxa" = c(
    length(taxa_names(phyloseq.samples.bacteria.phylum.unclassified.ra)),
    length(taxa_names(phyloseq.samples.bacteria.class.unclassified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.order.unclassified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.family.unclassified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.genus.unclassified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.species.unclassified.ra))),
  
  # Classified taxa
  "Number of Classified Taxa" = c(
    length(taxa_names(phyloseq.samples.bacteria.phylum.classified.ra)),
    length(taxa_names(phyloseq.samples.bacteria.class.classified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.order.classified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.family.classified.ra)), 
    length(taxa_names(phyloseq.samples.bacteria.genus.classified.ra)),
    length(taxa_names(phyloseq.samples.bacteria.species.classified.ra))),
  
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
stable2
#Make into excel file
write_xlsx(stable2, 
          "SupplementaryTable2.xlsx")



#ALPHA DIVERSITY ######
##BACTERIAL AND ARCHAEAL TAXA#####
alpha_div1 <- phyloseq::estimate_richness(phyloseq.samples.bacteria, 
                                          measures = c("Observed", "Shannon")) # richness, diversity
alpha_div2 <- microbiome::evenness(phyloseq.samples.bacteria, index = "pielou", 
                                   zeroes = FALSE, #Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
                                   detection = 0) ##evenness

#Combine alpha div metrics with metadata
alpha_div <- cbind(alpha_div1, alpha_div2)
alpha_div

#Metadata and div metrics
alpha_div_meta <- cbind(phyloseq.samples.bacteria@sam_data, 
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

###Naive (H21) vs Established (P1)#####
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_P1vsH21

####Alpha diversity indeces and Water quality over time#####
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

#####On which days do the phases change?#######
metadata_phyloseq.samples.bacteria <- data.frame(phyloseq.samples.bacteria@sam_data)
#When do phases switch?
#Naive
metadata_phyloseq.samples.bacteria %>%
  filter(Enclosure == "H21") %>%
  group_by(Date_num_phase_naive) %>%
  summarise(
    min_naive = min(Date_num_naive, na.rm = TRUE),
    max_naive = max(Date_num_naive, na.rm = TRUE)
  )

#Established
metadata_phyloseq.samples.bacteria %>%
  filter(Enclosure == "P1") %>%
  group_by(Date_num_phase_established) %>%
  summarise(
    min_naive = min(Date_num_established, na.rm = TRUE),
    max_naive = max(Date_num_established, na.rm = TRUE)
  )

#####Plot - Just  Shannon, Copper and Ammonia levels####
alpha_div_wq_date_num_factor <- ggplot(alpha_div_wq_time_long%>%
                                filter(Index %in% c("Copper_level_mg_L",
                                                    "Shannon",
                                                    "Ammonia_mg_L")),
                              aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
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
  ggh4x::facetted_pos_scales(
    x = list(
      
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
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

#####Plot - Other Water quality measures#####
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
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
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
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive"))) +
  theme_bw() +
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

##NITRIFYING TAXA#####
phyloseq.samples.nitrifiers #1690 taxa and 204 samples
alpha_div1_nit <- phyloseq::estimate_richness(phyloseq.samples.nitrifiers, 
                                              measures = c("Observed", "Shannon")) # richness, diversity
#Evenness based only on taxa actually present in each sample, so zeroes set to FALSE.  Keeps the focus on the taxa actually observed.
alpha_div2_nit <- microbiome::evenness(phyloseq.samples.nitrifiers, index = "pielou", 
                                   zeroes = TRUE, 
                                   detection = 0) ##evenness

# combine metrics with metadata
alpha_div_nit <- cbind(alpha_div1_nit, alpha_div2_nit)
alpha_div_nit
alpha_div_nit_meta <- cbind(phyloseq.samples.nitrifiers@sam_data, 
                        alpha_div_nit) 
alpha_div_nit_meta


#Pivot to long format 
alpha_div_nit_meta_long <- 
  alpha_div_nit_meta %>%
  pivot_longer(cols = c(Observed, Shannon, pielou),  
               names_to = "alpha_div_metric", 
               values_to = "alpha_div_value") 
alpha_div_nit_meta_long

##Factoring alpha div metrics
alpha_div_nit_meta_long$alpha_div_metric<- factor(alpha_div_nit_meta_long$alpha_div_metric, levels = c("Observed","pielou", "Shannon"))

###Naive (H21) vs Established (P1)#####
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
            label = "Wilcoxon, p = {p.format}",
            step.increase = 0.1,
            size = 0.5,
            label.size = 5,
            tip.length = 0.02,
            hide.ns = T) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))
alpha_div_nit_P1vsH21

###Alpha div indices over time####
alpha_div_nit_time <- ggplot(alpha_div_nit_meta_long, 
                         aes(x = factor(Date_num), y= alpha_div_value)) +
  theme_bw() +
  labs(title= "ALPHA DIVERSITY") +
  facet_nested(alpha_div_metric ~Enclosure,
               scales = "free",
               labeller = as_labeller(c("Observed" = "RICHNESS\n(OBSERVED)",
                                        "Shannon" = "DIVERSITY\n(SHANNON)",
                                        "pielou" = "EVENNESS",
                                        "P1" = "Established",
                                        "H21" = "Naive"))) +
  geom_point(size = 3, shape = 18)+
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  theme(
    legend.position = 'none',
    strip.background = element_rect(fill = "black"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", linewidth= 1),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
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
alpha_div_nit_time

###Alpha diversity indices and Water quality over time#####
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

####When do phases change?#####
metadata_phyloseq.samples.nitrifiers <- data.frame(phyloseq.samples.nitrifiers@sam_data)
#Naive
metadata_phyloseq.samples.nitrifiers %>%
  filter(Enclosure == "H21") %>%
  group_by(Date_num_phase_naive) %>%
  summarise(
    min_naive = min(Date_num_naive, na.rm = TRUE),
    max_naive = max(Date_num_naive, na.rm = TRUE)
  )

#Established
metadata_phyloseq.samples.nitrifiers %>%
  filter(Enclosure == "P1") %>%
  group_by(Date_num_phase_established) %>%
  summarise(
    min_naive = min(Date_num_established, na.rm = TRUE),
    max_naive = max(Date_num_established, na.rm = TRUE)
  )

###Plot - Just  Shannon, Copper and Ammonia levels#####
alpha_div_nit_wq_date_num_factor <- ggplot(alpha_div_nit_wq_time_long%>%
                                             filter(Index %in% c("Copper_level_mg_L",
                                                                 "Shannon",
                                                                 "Ammonia_mg_L"
                                             )),
                                           aes(x = factor(Date_num), y = Index_value, color = Copper_keep)) +
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
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
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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
        panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
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


###Plot - Shannon with Other Water quality measures#####
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
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
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
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
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


###Plot - Shannon with Other Water quality measures#####
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
  geom_point(size = 3, shape = 18)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
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
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
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


#RELATIVE ABUNDANCE######
##ALL TAXA######
###INCLUDING ALL TAXA AND UNCLASSIFIED READS####
#### ORDER LEVEL #####
phyloseq_with_unclassifieds_otu.samples.order.ra #11,709 orders and 204 samples

#Grouping the low abundance orders into one category
phyloseq_with_unclassifieds_otu.samples.order.filt <- merge_low_abundance_grouped_ra(phyloseq_with_unclassifieds_otu.samples.order.ra, 
                                                                       "Enclosure", 
                                                                       level = "Order", 
                                                                       threshold = 0.5)
phyloseq_with_unclassifieds_otu.samples.order.filt #16 orders over 0.5% mean RA
#Melt the filtered ps object to long format
phyloseq_with_unclassifieds_otu.samples.order.filt.melt <- psmelt(phyloseq_with_unclassifieds_otu.samples.order.filt)

##Factoring the Order column so that they are ordered by class, and 'Others' and 'Unclassified Reads' are last
order_levels_with_unclassifieds_df <- phyloseq_with_unclassifieds_otu.samples.order.filt.melt %>%
  mutate(last_rank = case_when(
    grepl("Others", Order, ignore.case = TRUE) ~ 2,
    grepl("Reads",   Order, ignore.case = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  arrange(last_rank, Class, Order)
  
order_levels_with_unclassifieds <- order_levels_with_unclassifieds_df%>%
  distinct(Class, Order, last_rank) %>%
  arrange(last_rank, Class, Order) %>%
  pull(Order) %>%
  unique()

#Apply Order levels
phyloseq_with_unclassifieds_otu.samples.order.filt.melt$Order <- 
  factor(phyloseq_with_unclassifieds_otu.samples.order.filt.melt$Order, 
         levels = order_levels_with_unclassifieds)
levels(phyloseq_with_unclassifieds_otu.samples.order.filt.melt$Order) ##ok

##Create color palette
order_filt_palette_with_unclassifieds <- unname(alphabet2())

#Create a named vector for the palette
order_filt_names_with_unclassifieds <- 
  unique(phyloseq_with_unclassifieds_otu.samples.order.filt.melt$Order)
order_named_palette_with_unclassifieds <- setNames((order_filt_palette_with_unclassifieds)[1:length(order_filt_names_with_unclassifieds)], 
                                                   order_filt_names_with_unclassifieds)
order_named_palette_with_unclassifieds$'Others <0.5% RA' <- "grey95"
order_named_palette_with_unclassifieds$'Caulobacterales' <-  "orange4"
order_named_palette_with_unclassifieds$'Rhizobiales' <- "lightgoldenrod3"
order_named_palette_with_unclassifieds$'Sphingomonadales' <-  'sienna1'
order_named_palette_with_unclassifieds$'Rhodobacterales' <- "#E3B199"
order_named_palette_with_unclassifieds$'unclassified Alphaproteobacteria' <- "gold"
order_named_palette_with_unclassifieds$'Flavobacteriales' <-  'darkslategray'
order_named_palette_with_unclassifieds$'unclassified Bacteroidia' <- "olivedrab"
order_named_palette_with_unclassifieds$'unclassified Gammaproteobacteria' <- "#4B0082"
order_named_palette_with_unclassifieds$'Pseudomonadales' <- "#DEA0FD"
order_named_palette_with_unclassifieds$'Enterobacterales' <- "mediumorchid4"
order_named_palette_with_unclassifieds$'Burkholderiales' <- 'palevioletred'
order_named_palette_with_unclassifieds$'Nitrososphaerales' <- "#951942"
order_named_palette_with_unclassifieds$'unclassified Bacteria' <- "navyblue"
order_named_palette_with_unclassifieds$'unclassified Pseudomonadota' <- "dodgerblue"
order_named_palette_with_unclassifieds$'Unclassified Reads' <- "grey50"


##### Plot RA at the order level (over time) #####
RA_order_alltaxa_unclassifieds_plot_datenum <- ggplot(phyloseq_with_unclassifieds_otu.samples.order.filt.melt,
                                                   aes(x=factor(Date_num), y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = order_named_palette_with_unclassifieds,
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_order_alltaxa_unclassifieds_plot_datenum

#### FAMILY LEVEL #####
phyloseq_with_unclassifieds_otu.samples.family.ra #23,427 families and 204 samples 

#Group the low abundance families
phyloseq_with_unclassifieds_otu.samples.family.filt <- merge_low_abundance_grouped_ra(phyloseq_with_unclassifieds_otu.samples.family.ra, 
                                                                        "Enclosure", 
                                                                        level = "Family", 
                                                                        threshold = 0.3)
phyloseq_with_unclassifieds_otu.samples.family.filt #23 families over 0.3% mean RA

#Melt the filtered ps object to long format
phyloseq_with_unclassifieds_otu.samples.family.filt.melt <- psmelt(phyloseq_with_unclassifieds_otu.samples.family.filt)

##Factoring the Family column so that they are ordered by Phylum, Class, Order, and 'Others' and 'Unclassified Reads' are last
family_levels_with_unclassifieds_df <- phyloseq_with_unclassifieds_otu.samples.family.filt.melt %>%
  mutate(last_rank = case_when(
    grepl("Others", Family, ignore.case = TRUE) ~ 2,
    grepl("Reads", Family, ignore.case = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  arrange(last_rank, Class, Order, Family)

family_levels_with_unclassifieds <- family_levels_with_unclassifieds_df%>%
  distinct(Phylum, Class, Order, Family, last_rank) %>%
  arrange(last_rank, Phylum, Class, Order, Family) %>%
  pull(Family) %>%
  unique()

#Apply Family levels
phyloseq_with_unclassifieds_otu.samples.family.filt.melt$Family <- 
  factor(phyloseq_with_unclassifieds_otu.samples.family.filt.melt$Family, 
         levels = family_levels_with_unclassifieds)
levels(phyloseq_with_unclassifieds_otu.samples.family.filt.melt$Family) ##ok

#Have to add some orders that are not in the main palette but need them for the family plot: 
order_named_palette_with_unclassifieds$'Mycobacteriales' <- "firebrick1"
order_named_palette_with_unclassifieds$'Nevskiales' <- 'slateblue'
order_named_palette_with_unclassifieds$'Pirellulales' <- 'tan3'
order_named_palette_with_unclassifieds$'Rhodospirillales' <- 'chocolate3'
order_named_palette_with_unclassifieds$'unclassified Actinomycetes' <- 'turquoise'
order_named_palette_with_unclassifieds$'unclassified Actinomycetota' <- 'turquoise4'
order_named_palette_with_unclassifieds$'Others <0.3% RA' <- 'grey95'

##Create color palette - based on families within the same order
palette_family_level_with_unclassifieds_df <- phyloseq_with_unclassifieds_otu.samples.family.filt.melt %>% 
  arrange(Order, Family) %>%   # ensure consistent shading order
  group_by(Order) %>%
  mutate(
    base_color = order_named_palette_with_unclassifieds[Order],
    shade = seq(0.0005, 0.9, length.out = n()),  # avoid extremes
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_family_level_with_unclassifieds_df

#Set up final palette
family_named_palette_with_unclassifieds <- setNames(
  palette_family_level_with_unclassifieds_df$color,
  palette_family_level_with_unclassifieds_df$Family)
family_named_palette_with_unclassifieds
family_named_palette_with_unclassifieds$'Nitrosopumilaceae' <- "#951942"

#If too many families for the legend, filter how many to indlude:
top_families_with_unclassifieds <- top_taxa_legend(phyloseq_with_unclassifieds_otu.samples.family.filt.melt,
                                taxlevel = "Family", n = 15)
top_families_with_unclassifieds
#Keep the factor levels
top_families_with_unclassifieds <- family_levels_with_unclassifieds[family_levels_with_unclassifieds %in% top_families_with_unclassifieds]

##### Plot RA at the order level (over time) #####
RA_family_alltaxa_unclassifieds_plot_datenum <- ggplot(phyloseq_with_unclassifieds_otu.samples.family.filt.melt,
                                                    aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = family_named_palette_with_unclassifieds,
                    breaks = top_families_with_unclassifieds,
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_alltaxa_unclassifieds_plot_datenum


###INCLUDING ALL TAXA####
#### ORDER LEVEL #####
phyloseq.samples.order.ra

#Grouping the low abundance orders into one category
phyloseq.samples.order.filt <- merge_low_abundance_grouped_ra(phyloseq.samples.order.ra, 
                                                                                     "Enclosure", 
                                                                                     level = "Order", 
                                                                                     threshold = 0.5)
phyloseq.samples.order.filt #26 orders over 0.5% mean RA

#Melt the filtered ps object to long format
phyloseq.samples.order.filt.melt <- psmelt(phyloseq.samples.order.filt)

##Factoring the Order column so that they are ordered by class, and 'Others' are last
order_levels_alltaxa_df <- phyloseq.samples.order.filt.melt %>%
  mutate(last_rank = case_when(
    grepl("Others", Order, ignore.case = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  arrange(last_rank, Phylum, Class, Order)

order_levels_alltaxa <- order_levels_alltaxa_df%>%
  distinct(Phylum, Class, Order, last_rank) %>%
  arrange(last_rank, Phylum, Class, Order) %>%
  pull(Order) %>%
  unique()

#Apply Order levels
phyloseq.samples.order.filt.melt$Order <- 
  factor(phyloseq.samples.order.filt.melt$Order, 
         levels = order_levels_alltaxa)
levels(phyloseq.samples.order.filt.melt$Order) ##ok

##Create color palette
order_filt_palette_alltaxa <- unname(alphabet2())

#Create a named vector for the palette
order_filt_names_alltaxa <- 
  unique(phyloseq.samples.order.filt.melt$Order)
order_named_palette_alltaxa <- setNames((order_filt_palette_alltaxa)[1:length(order_filt_names_alltaxa)], 
                                        order_filt_names_alltaxa)
order_named_palette_alltaxa$'Others <0.5% RA' <- "grey95"
order_named_palette_alltaxa$'Caulobacterales' <-  "orange4"
order_named_palette_alltaxa$'Rhizobiales' <- "lightgoldenrod3"
order_named_palette_alltaxa$'Sphingomonadales' <-  'sienna1'
order_named_palette_alltaxa$'Rhodobacterales' <- "#E3B199"
order_named_palette_alltaxa$'unclassified Alphaproteobacteria' <- "gold"
order_named_palette_alltaxa$'Flavobacteriales' <-  'darkslategray'
order_named_palette_alltaxa$'unclassified Bacteroidia' <- "olivedrab"
order_named_palette_alltaxa$'unclassified Gammaproteobacteria' <- "#4B0082"
order_named_palette_alltaxa$'Pseudomonadales' <- "#DEA0FD"
order_named_palette_alltaxa$'Enterobacterales' <- "mediumorchid4"
order_named_palette_alltaxa$'Burkholderiales' <- 'palevioletred'
order_named_palette_alltaxa$'Nitrososphaerales' <- "#951942"
order_named_palette_alltaxa$'unclassified Bacteria' <- "navyblue"
order_named_palette_alltaxa$'unclassified Pseudomonadota' <- "dodgerblue"
order_named_palette_alltaxa$'Unclassified Reads' <- "grey50"
order_named_palette_alltaxa$'Mycobacteriales' <- "firebrick1"
order_named_palette_alltaxa$'Nevskiales' <- 'slateblue'
order_named_palette_alltaxa$'Pirellulales' <- 'tan3'
order_named_palette_alltaxa$'Rhodospirillales' <- 'chocolate3'
order_named_palette_alltaxa$'unclassified Actinomycetes' <- 'turquoise'
order_named_palette_alltaxa$'unclassified Actinomycetota' <- 'turquoise4'
order_named_palette_alltaxa$'Chlamydiales' <- 'greenyellow'
order_named_palette_alltaxa$'Planctomycetales' <- 'gold4'
order_named_palette_alltaxa$'Others <0.3% RA' <- 'grey95'

# If too many taxa for legend, apply function to keep only top ones
top_orders_alltaxa <- top_taxa_legend(phyloseq.samples.order.filt.melt,
                              taxlevel = "Order", n = 20)
top_orders_alltaxa
#Keep factor levels 
top_orders_alltaxa <- order_levels_alltaxa[order_levels_alltaxa %in% top_orders_alltaxa]

##### Plot RA at the order level (over time) #####
RA_order_alltaxa_plot_datenum <- ggplot(phyloseq.samples.order.filt.melt,
                                                   aes(x=factor(Date_num), y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = order_named_palette_alltaxa,
                    breaks = top_orders_alltaxa,
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_order_alltaxa_plot_datenum

#### FAMILY LEVEL #####
phyloseq.samples.family.ra #23,426 families and 216 samples 

#Merge low abundance families
phyloseq.samples.family.filt <- merge_low_abundance_grouped_ra(phyloseq.samples.family.ra, 
                                                                        "Enclosure", 
                                                                        level = "Family", 
                                                                        threshold = 0.5)
phyloseq.samples.family.filt #25 families over 0.5% mean RA

#Melt the filtered ps object to long format
phyloseq.samples.family.filt.melt <- psmelt(phyloseq.samples.family.filt)

##Factoring the Family column so that they are ordered by Phylum, Class, Order, and 'Others' are last
family_levels_alltaxa_df <- phyloseq.samples.family.filt.melt %>%
  mutate(last_rank = case_when(
    grepl("Others", Family, ignore.case = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  arrange(last_rank, Class, Order, Family)

family_levels_alltaxa <- family_levels_alltaxa_df%>%
  distinct(Phylum, Class, Order, Family, last_rank) %>%
  arrange(last_rank, Phylum, Class, Order, Family) %>%
  pull(Family) %>%
  unique()

#Apply Family levels
phyloseq.samples.family.filt.melt$Family <- 
  factor(phyloseq.samples.family.filt.melt$Family, 
         levels = family_levels_alltaxa)
levels(phyloseq.samples.family.filt.melt$Family) ##ok

##Create color palette - based on families within the same order
palette_family_level_alltaxa_df <- phyloseq.samples.family.filt.melt %>% 
  arrange(Phylum, Class, Order, Family) %>%   # ensure consistent shading order
  group_by(Order) %>%
  mutate(
    base_color = order_named_palette_alltaxa[Order],
    shade = seq(0.02, 0.5, length.out = n()),  # avoid extremes
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_family_level_alltaxa_df

#Set up final palette
family_named_palette_alltaxa <- setNames(
  palette_family_level_alltaxa_df$color,
  palette_family_level_alltaxa_df$Family)
family_named_palette_alltaxa
family_named_palette_alltaxa$'Nitrosopumilaceae' <- "#951942"
family_named_palette_alltaxa$'Others <0.5% RA' <- "grey95"

##If too many taxa for the legend, apply the function to obtain top ones
top_families_alltaxa <- top_taxa_legend(phyloseq.samples.family.filt.melt,
                                taxlevel = "Family", 
                                n = 17)
top_families_alltaxa
#Keep the factor levels 
top_families_alltaxa <- family_levels_alltaxa[family_levels_alltaxa %in% top_families_alltaxa]

##### Plot RA at the family level (over time) #####
RA_family_alltaxa_plot_datenum <- ggplot(phyloseq.samples.family.filt.melt,
                                                    aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color =    "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = family_named_palette_alltaxa,
                    breaks = top_families_alltaxa,
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_alltaxa_plot_datenum

###ONLY ARCHAEAL AND BACTERIAL TAXA####
#### ORDER LEVEL #####
phyloseq.samples.bacteria.order.ra #2583 taxa and 204 samples 

#Grouping the low abundance orders into one category
phyloseq.samples.bacteria.order.filt <- merge_low_abundance_grouped_ra(phyloseq.samples.bacteria.order.ra, 
                                                                        "Enclosure", 
                                                                        level = "Order", 
                                                                        threshold = 0.5)
phyloseq.samples.bacteria.order.filt #25 orders over 0.5% mean RA

#Melt the filtered ps object to long format
phyloseq.samples.bacteria.order.filt.melt <- psmelt(phyloseq.samples.bacteria.order.filt)

##Factoring the Order column so that they are ordered by Phylum and Class, and 'Others'  are last
order_levels_bacarc_df <- phyloseq.samples.bacteria.order.filt.melt %>%
  mutate(last_rank = case_when(
    grepl("Others", Order, ignore.case = TRUE) ~ 1,
    TRUE ~ 0
  )) %>%
  arrange(last_rank, Phylum, Class, Order)

order_levels_bacarc <- order_levels_bacarc_df%>%
  distinct(Phylum, Class, Order, last_rank) %>%
  arrange(last_rank, Phylum, Class, Order) %>%
  pull(Order) %>%
  unique()

#Apply Order levels
phyloseq.samples.bacteria.order.filt.melt$Order <- 
  factor(phyloseq.samples.bacteria.order.filt.melt$Order, 
         levels = order_levels_bacarc)
levels(phyloseq.samples.bacteria.order.filt.melt$Order) ##ok

##Create color palette
order_filt_palette_bacarc <- unname(alphabet2())

#Create a named vector for the palette
order_filt_names_bacarc <- 
  unique(phyloseq.samples.bacteria.order.filt.melt$Order)
order_named_palette_bacarc <- setNames((order_filt_palette_bacarc)[1:length(order_filt_names_bacarc)], 
                                        order_filt_names_bacarc)
order_named_palette_bacarc$'Others <0.5% RA' <- "grey95"
order_named_palette_bacarc$'Caulobacterales' <-  "orange4"
order_named_palette_bacarc$'Rhizobiales' <- "lightgoldenrod3"
order_named_palette_bacarc$'Sphingomonadales' <-  'sienna1'
order_named_palette_bacarc$'Rhodobacterales' <- "#E3B199"
order_named_palette_bacarc$'unclassified Alphaproteobacteria' <- "gold"
order_named_palette_bacarc$'Flavobacteriales' <-  'darkslategray'
order_named_palette_bacarc$'unclassified Bacteroidia' <- "olivedrab"
order_named_palette_bacarc$'unclassified Gammaproteobacteria' <- "#4B0082"
order_named_palette_bacarc$'Pseudomonadales' <- "#DEA0FD"
order_named_palette_bacarc$'Enterobacterales' <- "mediumorchid4"
order_named_palette_bacarc$'Burkholderiales' <- 'palevioletred'
order_named_palette_bacarc$'Nitrososphaerales' <- "#951942"
order_named_palette_bacarc$'unclassified Bacteria' <- "navyblue"
order_named_palette_bacarc$'unclassified Pseudomonadota' <- "dodgerblue"
order_named_palette_bacarc$'Unclassified Reads' <- "grey50"
order_named_palette_bacarc$'Mycobacteriales' <- "firebrick1"
order_named_palette_bacarc$'Nevskiales' <- 'slateblue'
order_named_palette_bacarc$'Pirellulales' <- 'tan3'
order_named_palette_bacarc$'Rhodospirillales' <- 'chocolate3'
order_named_palette_bacarc$'unclassified Actinomycetes' <- 'turquoise'
order_named_palette_bacarc$'unclassified Actinomycetota' <- 'turquoise4'
order_named_palette_bacarc$'Chlamydiales' <- 'greenyellow'
order_named_palette_bacarc$'Planctomycetales' <- 'gold4'
order_named_palette_bacarc$'Others <0.3% RA' <- 'grey95'

# If too many taxa for legend, apply function to keep only top ones
top_orders_bacarc <- top_taxa_legend(phyloseq.samples.bacteria.order.filt.melt,
                                      taxlevel = "Order", n = 20)
top_orders_bacarc
#Keep factor levels in legend
top_orders_bacarc <- order_levels_bacarc[order_levels_bacarc %in% top_orders_bacarc]


##### Plot RA at the order level (over time) #####
RA_order_bacteriaarchaea_plot_datenum <- ggplot(phyloseq.samples.bacteria.order.filt.melt,
                                     aes(x=factor(Date_num), y= Abundance, fill = Order)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = order_named_palette_bacarc,
                    breaks = top_orders_bacarc,
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
        panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
        axis.line.y = element_line(linewidth = 0.7, colour = "black"),
        axis.text.x = element_text(colour = "black", size = 20,
                                   vjust = 0.5, hjust = 0.5),
        axis.title = element_text(colour = "black", size = 22),
        axis.text.y = element_text(colour = "black", size = 20),
        axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_order_bacteriaarchaea_plot_datenum

#### FAMILY LEVEL#####
phyloseq.samples.bacteria.family.ra #7596 families and 204 samples 

phyloseq.samples.bacteria.family.filt <- merge_low_abundance_grouped_ra(phyloseq.samples.bacteria.family.ra, 
                                                                             "Enclosure", 
                                                                             level = "Family", 
                                                                             threshold = 0.5)
phyloseq.samples.bacteria.family.filt #28 families over 0.5% mean RA

#Melt the filtered ps object to long format
phyloseq.samples.bacteria.family.filt.melt <- psmelt(phyloseq.samples.bacteria.family.filt)

##Factoring the Family column so that they are ordered by Phylum, Class, Order; 'Others'  are last
family_levels_bacarc_df <- phyloseq.samples.bacteria.family.filt.melt %>%
  mutate(last_rank = case_when(
    grepl("Others", Family, ignore.case = TRUE) ~1,
    TRUE ~ 0
  )) %>%
  arrange(last_rank, Phylum, Class, Order, Family)

family_levels_bacarc <- family_levels_bacarc_df%>%
  distinct(Phylum, Class, Order, Family, last_rank) %>%
  arrange(last_rank, Phylum, Class, Order, Family) %>%
  pull(Family) %>%
  unique()

#Apply Family levels
phyloseq.samples.bacteria.family.filt.melt$Family <- 
  factor(phyloseq.samples.bacteria.family.filt.melt$Family, 
         levels = family_levels_bacarc)
levels(phyloseq.samples.bacteria.family.filt.melt$Family) ##ok

##Create color palette - based on families within the same order
palette_family_level_bacarc_df <- phyloseq.samples.bacteria.family.filt.melt %>% 
  arrange(Phylum, Class, Order, Family) %>%   # ensure consistent shading order
  group_by(Order) %>%
  mutate(
    base_color = order_named_palette_bacarc[Order],
    shade = seq(0.02, 0.5, length.out = n()),  # avoid extremes
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_family_level_bacarc_df

#Set up final palette
family_named_palette_bacarc <- setNames(
  palette_family_level_bacarc_df$color,
  palette_family_level_bacarc_df$Family)
family_named_palette_bacarc
family_named_palette_bacarc$'Nitrosopumilaceae' <- "#951942"
family_named_palette_bacarc$'Others <0.5% RA' <- "grey95"


##Create color palette - based on families within the same order
palette_family_level_bacarc_df <- phyloseq.samples.bacteria.family.filt.melt %>% 
  arrange(Order, Family) %>%   # ensure consistent shading order
  group_by(Order) %>%
  mutate(
    base_color = order_named_palette_bacarc[Order],
    shade = seq(0.02, 0.5, length.out = n()),  # avoid extremes
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_family_level_bacarc_df

#Set up final palette
family_named_palette_bacarc <- setNames(
  palette_family_level_bacarc_df$color,
  palette_family_level_bacarc_df$Family)
family_named_palette_bacarc
family_named_palette_bacarc$'Nitrosopumilaceae' <- "#951942"
family_named_palette_bacarc$'Others <0.5% RA' <- "grey95"

##If too many taxa for the legend, apply the function to obtain top ones
top_families_bacarc <- top_taxa_legend(phyloseq.samples.bacteria.family.filt.melt,
                                       taxlevel = "Family", 
                                       n = 20)
top_families_bacarc
#Keep the factor levels 
top_families_bacarc <- family_levels_bacarc[family_levels_bacarc %in% top_families_bacarc]


##### Plot RA at the family level (over time) #####
RA_family_bacteriaarchaea_plot_datenum <- ggplot(phyloseq.samples.bacteria.family.filt.melt,
                                             aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = family_named_palette_bacarc,
                    breaks = top_families_bacarc,
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_family_bacteriaarchaea_plot_datenum


#####SUPPLEMENTARY TABLE 5 #####
#Using relative abundances from the object that has all taxa, plus unclassified reads
stable5 <-  phyloseq_with_unclassifieds_otu.samples.family.filt.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Family, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  mutate(`Mean Relative Abundance (%) ± SD` = 
           paste0(
             mean_abun,
             " ± ", 
             sd_abun
           ), 
         `Min - Max Relative Abundance (%)` = paste0(
           min_abun,
           " - ",
           max_abun
         ))%>%
  rename(Phase = Date_num_phase_abbrv)%>%
  select(Family, System, Phase, `Mean Relative Abundance (%) ± SD`, `Min - Max Relative Abundance (%)`)
stable5
write_xlsx(stable5, 
           "SupplementaryTable5.xlsx")

# ###RA OF ONLY RHODOBACTERACEAE OVER TIME ##############
# phyloseq.samples.bacteria.family.ra #6755 families, 216 samples
# 
# #Out of this overall communities object, select only Nitrosopumilaceae
# phyloseq.samples.bacteria.family.ra.rhodobacteraceae <- subset_taxa(phyloseq.samples.bacteria.family.ra, 
#                                                                      Family == "Rhodobacteraceae") 
# phyloseq.samples.bacteria.family.ra.rhodobacteraceae <- subset_samples(phyloseq.samples.bacteria.family.ra.rhodobacteraceae, 
#                                                                         sample_sums(phyloseq.samples.bacteria.family.ra.rhodobacteraceae) > 0)
# phyloseq.samples.bacteria.family.ra.rhodobacteraceae #Rhodobacteraceae (1 taxa) in 216 samples 
# 
# #Melt to long format
# phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt <- psmelt(phyloseq.samples.bacteria.family.ra.rhodobacteraceae)
# 
# #Going to make a single column for abbreviated phases, so I can color by this column 
# phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt <-  phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt%>%
#   mutate(System = ifelse(Enclosure == "H21", "Naive", "Established"), 
#          System = factor(System, levels = c("Naive", "Established")))%>%
#   mutate(Date_num_colors_systems = case_when(
#     System == "Naive" & !is.na(Date_num_phase_naive_abbrv) ~
#       paste0(Date_num_phase_naive_abbrv, "_", System),
#     
#     System == "Established" & !is.na(Date_num_phase_established_abbrv) ~
#       paste0(Date_num_phase_established_abbrv, "_", System),
#     
#     TRUE ~ NA_character_
#   ))%>%
#   mutate(Date_num_colors_systems = factor(Date_num_colors_systems, 
#                                           levels = c("L_Naive", 
#                                                      "T1_Naive",
#                                                      "E1_Naive", 
#                                                      "T2_Naive", 
#                                                      "T3_Naive", 
#                                                      "E2_Naive",
#                                                      "P_Naive", 
#                                                      "L_Established", 
#                                                      "T1_Established", 
#                                                      "E_Established", 
#                                                      "P_Established"
#                                           )))
# 
# labels_rhodobacteraceae <- phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt %>%
#   distinct(Date_num_colors_systems, Date_num_phase_abbrv) %>%
#   mutate(Date_num_phase_abbrv = as.character(Date_num_phase_abbrv)) %>%
#   tibble::deframe()
# 
# #Get wilcox stats
# #Have to add another variable for each system to get their respective stats 
# phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt_2 <- phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt%>%
#   mutate(Date_num_phase_naive_abbrv_stats = paste0(Date_num_phase_naive_abbrv, "_Naive"), 
#          Date_num_phase_naive_abbrv_stats = factor(Date_num_phase_naive_abbrv_stats, 
#                                                    levels = c(
#                                                      "L_Naive", 
#                                                      "T1_Naive",
#                                                      "E1_Naive", 
#                                                      "T2_Naive", 
#                                                      "T3_Naive", 
#                                                      "E2_Naive",
#                                                      "P_Naive"))
#   )%>%
#   mutate(Date_num_phase_established_abbrv_stats = paste0(Date_num_phase_established_abbrv, "_Established"), 
#          Date_num_phase_established_abbrv_stats = factor(Date_num_phase_established_abbrv_stats, 
#                                                          levels = c(
#                                                            "L_Established", 
#                                                            "T1_Established", 
#                                                            "E_Established", 
#                                                            "P_Established")))
# 
# #Get stats for naive system - only comparisons in order of time (L vs T1, T1 vs E1, E1 vs T2, T2 vs T3, T3 vs E2, E2 vs P)
# stat_wilcox_test_rhodobacteraceae_H21 <- phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt_2 %>%
#   filter(Enclosure == "H21")%>%
#   wilcox_test(
#     Abundance ~ Date_num_phase_naive_abbrv_stats,
#     comparisons = list(
#       c("L_Naive", "T1_Naive"),
#       c("T1_Naive", "E1_Naive"),
#       c("E1_Naive", "T2_Naive"),
#       c("T2_Naive", "T3_Naive"),
#       c("T3_Naive", "E2_Naive"),
#       c("E2_Naive", "P_Naive")
#     ),
#     p.adjust.method = "BH"
#   ) %>%
#   filter(p.adj <= 0.05) 
# 
# stat_wilcox_test_rhodobacteraceae_df_H21 <- stat_wilcox_test_rhodobacteraceae_H21 %>%
#   add_xy_position(x = "Date_num_phase_naive_abbrv_stats", 
#                   step.increase = 0.10, 
#                   scales = "free"
#   )%>%
#   mutate(System = "Naive")
# 
# #Get stats for established system - only comparisons in order of time (L vs T1, T1 vs E, E vs P)
# stat_wilcox_test_rhodobacteraceae_P1 <- 
#   phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt_2 %>%
#   filter(Enclosure == "P1") %>%
#   wilcox_test(
#     Abundance ~ Date_num_phase_established_abbrv_stats,
#     comparisons = list(
#       c("L_Established", "T1_Established"),
#       c("T1_Established", "E_Established"),
#       c("E_Established", "P_Established")
#     ),
#     p.adjust.method = "BH"
#   ) %>%
#   filter(p.adj <= 0.05)
# 
# stat_wilcox_test_rhodobacteraceae_df_P1 <- stat_wilcox_test_rhodobacteraceae_P1 %>%
#   add_xy_position(x = "Date_num_phase_established_abbrv_stats",
#                   step.increase = 0.30)%>%
#   mutate(System = "Established")
# 
# stat_wilcox_test_rhodobacteraceae_df_both_H21_P1 <- bind_rows(stat_wilcox_test_rhodobacteraceae_df_P1, 
#                                                                stat_wilcox_test_rhodobacteraceae_df_H21)%>%
#   mutate(System = factor(System, levels = c("Naive", "Established")))
# 
# #Plot
# rhodobacteraceae.RA.phases.boxplot <- ggplot(phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt, 
#                                               aes(x = Date_num_colors_systems, 
#                                                   y= Abundance, 
#                                                   color = Date_num_colors_systems, 
#                                                   fill = Date_num_colors_systems)) +
#   theme_bw() +
#   facet_wrap(~System, 
#              scales = "free")+
#   labs(x = "Phase", 
#        y = expression(italic(Rhodobacteraceae) ~ "RA (%)"),
#        color = "Phase", 
#        fill = "Phase") +
#   geom_jitter(size = 3, shape = 18, 
#               alpha = 0.8, width = 0.2) +
#   geom_boxplot(alpha = 0.3) +
#   scale_x_discrete(labels = labels_rhodobacteraceae)+
#   scale_fill_manual(values = phases_naive_established_palette)+
#   scale_color_manual(values = phases_naive_established_palette)+
#   scale_y_continuous(expand= c(0.03,0,0.1,0)) +
#   theme(
#     legend.position = "none",
#     # legend.text = element_text(size = 20),
#     # legend.title = element_text(size = 22, face = "bold"),
#     panel.border = element_rect(colour = "black", linewidth= 1),
#     strip.background = element_rect(fill = "black"),
#     strip.text = element_text(color = "white", face = "bold", size = 32),
#     plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
#     axis.title = element_text(size = 24, colour = "black", face = "bold"),
#     axis.text = element_text(colour = "black", size = 20),
#     axis.ticks = element_line(colour = "black", linewidth = 0.5))+
#   stat_pvalue_manual(stat_wilcox_test_rhodobacteraceae_df_both_H21_P1,
#                    label = "p.adj.signif",
#                    hide.ns = TRUE,
#                    tip.length = 0.01)
# rhodobacteraceae.RA.phases.boxplot

##NITRIFIERS######
###WITHIN ALL TAXA INCLUDING UNCLASSIFIED READS####
#### FAMILY LEVEL#######
phyloseq_with_unclassifieds_otu.samples.family.ra #23,427 families

##Which families are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
unique(phyloseq.samples.nitrifiers.melt$Family) 
# "Nitrosopumilaceae"          "Nitrospiraceae"             "Nitrospinaceae"             "Nitrosomonadaceae"          "VA-1"                      
# "UBA8639"                    "JAJDAT01"                   "unclassified Nitrospinales" "Nitrococcaceae"             "Xanthobacteraceae"         
# "Nitrosococcaceae"           "JACRGP01"                   "Thermomicrobiaceae"         "Nitrososphaeraceae"         "Gallionellaceae"           
# "UBA213"                     "NS-4"                       "Nitrosocaldaceae"    

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq_with_unclassifieds_otu.samples
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.family.ra, 
                                                                            Family == "Nitrosopumilaceae"| #AOA
                                                                            Family == "Nitrospiraceae"| #NOB
                                                                            Family == "Nitrospinaceae"| #NOB
                                                                            Family == "Nitrosomonadaceae"| #AOB
                                                                            Family == "VA-1" | #NOB
                                                                            Family == "UBA8639"| #NOB
                                                                            Family == "JAJDAT01" | #NOB
                                                                            Family == "unclassified Nitrospinales"| #NOB
                                                                            Family == "Nitrococcaceae"| #NOB
                                                                            Family == "Xanthobacteraceae" | #NOB
                                                                            Family == "Nitrosococcaceae"| #AOB
                                                                            Family == "JACRGP01"| #NOB
                                                                            Family == "Thermomicrobiaceae"| #NOB
                                                                            Family == "Nitrososphaeraceae"| #AOA
                                                                            Family == "Gallionellaceae"| #NOB
                                                                            Family == "UBA9942"| #NOB
                                                                            Family == "JACRGO01"| #NOB
                                                                            Family == "NS-4"|#NOB
                                                                            Family == "unclassified Nitrospinia"| #NOB
                                                                            Family == "CAJXCL01"| #NOB
                                                                            Family == "SZUA-224"| #NOB
                                                                            Family == "Nitrosocaldaceae"| #AOA
                                                                            Family == "JAHJSF01") #NOB
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers <- subset_samples(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers, 
                                                                 sample_sums(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers) > 0)
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers #23 nitrifying families in 204 samples 

#Melt the ps object to long format 
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt <- 
  psmelt(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers)

#Which families? 
unique(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt$Family)#ok, 23

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Chromatiaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "Candidatus Nitrosocaldaceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "Ectothiorhodospiraceae" ~ "NOB",
    Family == "Nitrobacteraceae" ~ "NOB",
    Family == "Gallionellaceae" ~ "NOB",
    Family == "Nitrospinaceae" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Family) %>%
  mutate(Family = factor(Family, levels = unique(Family)))
levels(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt$Family) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on families within each ammonia-nitrite oxidizing group
palette_nitrifiers_with_unclassifieds_family_df <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt %>% 
  distinct(Family, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Family) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    #shade = seq(-0.1, 0.1, length.out = n()),
    shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_with_unclassifieds_family_df

#Set up final palette
palette_nitrifiers_with_unclassifieds_family <- setNames(
  palette_nitrifiers_with_unclassifieds_family_df$color,
  palette_nitrifiers_with_unclassifieds_family_df$Family)
palette_nitrifiers_with_unclassifieds_family

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_families_with_unclassifieds <- top_taxa_legend(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt,
                                           n = 5)
top_nitrifying_families_with_unclassifieds
#Keep the factor levels for legend
nitifier_families_levels_with_unclassifieds <- levels(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt$Family)
top_nitrifying_families_with_unclassifieds <- 
  nitifier_families_levels_with_unclassifieds[nitifier_families_levels_with_unclassifieds %in% top_nitrifying_families_with_unclassifieds]
top_nitrifying_families_with_unclassifieds

#Plot
RA_nitrifiers_family_with_unclassifieds_plot_datenum <- ggplot(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrifiers.melt,
                                                aes(x=factor(Date_num), 
                                                    y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_with_unclassifieds_family,
                    breaks = top_nitrifying_families_with_unclassifieds
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_family_with_unclassifieds_plot_datenum

#### GENUS LEVEL#######
phyloseq_with_unclassifieds_otu.samples.genus.ra #89,900 genera

##Which genera are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
length(unique(phyloseq.samples.nitrifiers.melt$Genus)) #162 Nitrifying genera

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.genus.ra,
                                                                           #AOB
                                                                           Family == "Nitrosomonadaceae" | # AOB; Includes genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                                                                             Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                                                                             #AOA
                                                                             Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). InGTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                                                             Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes genera Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                                                             Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                                                             Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                                                                             #NOB
                                                                             #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                                                             Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                                                             Family == "NS-4" | # NOB (mostly uncultured genomes)
                                                                             Family == "UBA8639" | # NOB
                                                                             #In the order Rhizobiales and family Xanthobacteraceae
                                                                             Genus == "Nitrobacter" | #NOB; most representative one
                                                                             #Order Burkholderiales, family Gallionellaceae
                                                                             Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                                                             #Order Nitrococcales, family Nitrococcaceae
                                                                             Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                                                             #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                                                             Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                                                                             # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                                                                             Genus == "Nitrolancea")  #NOB; new one
phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers <- subset_samples(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers, 
                                                                               sample_sums(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers) > 0)
phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers #162 nitrifying genera in 204 samples 

#Melt the ps object to long format 
phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt <- 
  psmelt(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers)

#Which genera? 
unique(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt$Genus)#OK, 162

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt <- phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Genus) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))
levels(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt$Genus) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on genera within each ammonia-nitrite oxidizing group
palette_nitrifiers_with_unclassifieds_genus_df <- phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt %>% 
  distinct(Genus, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Genus) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    shade = seq(0, 0.8, length.out = n()),
    #shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_with_unclassifieds_genus_df

#Set up final palette
palette_nitrifiers_with_unclassifieds_genus <- setNames(
  palette_nitrifiers_with_unclassifieds_genus_df$color,
  palette_nitrifiers_with_unclassifieds_genus_df$Genus)
palette_nitrifiers_with_unclassifieds_genus

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_genera_with_unclassifieds <- top_taxa_legend(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt,
                                                            taxlevel = "Genus",
                                                              n = 6)
top_nitrifying_genera_with_unclassifieds
#Keep the factor levels for legend
nitifier_genera_levels_with_unclassifieds <- levels(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt$Genus)
top_nitrifying_genera_with_unclassifieds <- 
  nitifier_genera_levels_with_unclassifieds[nitifier_genera_levels_with_unclassifieds %in% top_nitrifying_genera_with_unclassifieds]
top_nitrifying_genera_with_unclassifieds

#Plot
RA_nitrifiers_genus_with_unclassifieds_plot_datenum <- ggplot(phyloseq_with_unclassifieds_otu.samples.genus.ra.nitrifiers.melt,
                                                               aes(x=factor(Date_num), 
                                                                   y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_with_unclassifieds_genus,
                    labels = function(x) str_wrap(x, width = 20), 
                    breaks = top_nitrifying_genera_with_unclassifieds
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.5),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_genus_with_unclassifieds_plot_datenum

#### SPECIES LEVEL#######
phyloseq_with_unclassifieds_otu.samples.species.ra #309748 Species

##Which Species are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
length(unique(phyloseq.samples.nitrifiers.melt$Species)) #1688 Nitrifying Species

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.species.ra,
                                                                           #AOB
                                                                           Family == "Nitrosomonadaceae" | # AOB; Includes Genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                                                                             Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                                                                             #AOA
                                                                             Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). In GTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                                                             Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes Species Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                                                             Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                                                             Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                                                                             #NOB
                                                                             #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                                                             Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                                                             Family == "NS-4" | # NOB (mostly uncultured genomes)
                                                                             Family == "UBA8639" | # NOB
                                                                             #In the order Rhizobiales and family Xanthobacteraceae
                                                                             Genus == "Nitrobacter" | #NOB; most representative one
                                                                             #Order Burkholderiales, family Gallionellaceae
                                                                             Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                                                             #Order Nitrococcales, family Nitrococcaceae
                                                                             Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                                                             #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                                                             Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                                                                             # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                                                                             Genus == "Nitrolancea")  #NOB; new one
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers <- subset_samples(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers, 
                                                                              sample_sums(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers) > 0)
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers #1690 nitrifying Species in 204 samples 

#Melt the ps object to long format 
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt <- 
  psmelt(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers)

#Which Species? 
length(unique(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt$Species))#OK, 1688

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt <- phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Species) %>%
  mutate(Species = factor(Species, levels = unique(Species)))
levels(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt$Species) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on Species within each ammonia-nitrite oxidizing group
palette_nitrifiers_with_unclassifieds_species_df <- phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt %>% 
  distinct(Species, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Species) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    shade = seq(0, 0.8, length.out = n()),
    #shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_with_unclassifieds_species_df

#Set up final palette
palette_nitrifiers_with_unclassifieds_species <- setNames(
  palette_nitrifiers_with_unclassifieds_species_df$color,
  palette_nitrifiers_with_unclassifieds_species_df$Species)
palette_nitrifiers_with_unclassifieds_species

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_Species_with_unclassifieds <- top_taxa_legend(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt,
                                                            taxlevel = "Species",
                                                            n = 15)
top_nitrifying_Species_with_unclassifieds
#Keep the factor levels for legend
nitifier_Species_levels_with_unclassifieds <- levels(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt$Species)
top_nitrifying_Species_with_unclassifieds <- 
  nitifier_Species_levels_with_unclassifieds[nitifier_Species_levels_with_unclassifieds %in% top_nitrifying_Species_with_unclassifieds]
top_nitrifying_Species_with_unclassifieds

#Plot
RA_nitrifiers_species_with_unclassifieds_plot_datenum <- ggplot(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers.melt,
                                                              aes(x=factor(Date_num), 
                                                                  y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_with_unclassifieds_species,
                    labels = function(x) str_wrap(x, width = 20), 
                    breaks = top_nitrifying_Species_with_unclassifieds
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.5),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_species_with_unclassifieds_plot_datenum

###WITHIN ALL TAXA####
#### GENUS LEVEL#######
phyloseq.samples.genus.ra #89,899 genera

##Which genera are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
length(unique(phyloseq.samples.nitrifiers.melt$Genus)) #162 Nitrifying genera

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq.samples.genus.ra.nitrifiers <- subset_taxa(phyloseq.samples.genus.ra,
                                                                           #AOB
                                                                           Family == "Nitrosomonadaceae" | # AOB; Includes genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                                                                             Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                                                                             #AOA
                                                                             Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). InGTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                                                             Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes genera Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                                                             Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                                                             Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                                                                             #NOB
                                                                             #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                                                             Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                                                             Family == "NS-4" | # NOB (mostly uncultured genomes)
                                                                             Family == "UBA8639" | # NOB
                                                                             #In the order Rhizobiales and family Xanthobacteraceae
                                                                             Genus == "Nitrobacter" | #NOB; most representative one
                                                                             #Order Burkholderiales, family Gallionellaceae
                                                                             Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                                                             #Order Nitrococcales, family Nitrococcaceae
                                                                             Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                                                             #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                                                             Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                                                                             # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                                                                             Genus == "Nitrolancea")  #NOB; new one
phyloseq.samples.genus.ra.nitrifiers <- subset_samples(phyloseq.samples.genus.ra.nitrifiers, 
                                                                              sample_sums(phyloseq.samples.genus.ra.nitrifiers) > 0)
phyloseq.samples.genus.ra.nitrifiers #162 nitrifying genera in 204 samples 

#Melt the ps object to long format 
phyloseq.samples.genus.ra.nitrifiers.melt <- 
  psmelt(phyloseq.samples.genus.ra.nitrifiers)

#Which genera? 
unique(phyloseq.samples.genus.ra.nitrifiers.melt$Genus)#OK, 162

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.samples.genus.ra.nitrifiers.melt <- phyloseq.samples.genus.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Genus) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))
levels(phyloseq.samples.genus.ra.nitrifiers.melt$Genus) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on genera within each ammonia-nitrite oxidizing group
palette_nitrifiers_all_taxa_genus_df <- phyloseq.samples.genus.ra.nitrifiers.melt %>% 
  distinct(Genus, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Genus) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    shade = seq(0, 0.8, length.out = n()),
    #shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_all_taxa_genus_df

#Set up final palette
palette_nitrifiers_all_taxa_genus <- setNames(
  palette_nitrifiers_all_taxa_genus_df$color,
  palette_nitrifiers_all_taxa_genus_df$Genus)
palette_nitrifiers_all_taxa_genus

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_genera_all_taxa <- top_taxa_legend(phyloseq.samples.genus.ra.nitrifiers.melt,
                                                            taxlevel = "Genus",
                                                            n = 10)
top_nitrifying_genera_all_taxa
#Keep the factor levels for legend
nitifier_genera_levels_all_taxa <- levels(phyloseq.samples.genus.ra.nitrifiers.melt$Genus)
top_nitrifying_genera_all_taxa <- 
  nitifier_genera_levels_all_taxa[nitifier_genera_levels_all_taxa %in% top_nitrifying_genera_all_taxa]
top_nitrifying_genera_all_taxa

#Plot
RA_nitrifiers_genus_all_taxa_plot_datenum <- ggplot(phyloseq.samples.genus.ra.nitrifiers.melt,
                                                              aes(x=factor(Date_num), 
                                                                  y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_all_taxa_genus,
                    labels = function(x) str_wrap(x, width = 20), 
                    breaks = top_nitrifying_genera_all_taxa
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.5),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_genus_all_taxa_plot_datenum

#### SPECIES LEVEL#######
phyloseq.samples.species.ra #309747 Species

##Which Species are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
length(unique(phyloseq.samples.nitrifiers.melt$Species)) #1688 Nitrifying Species

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq.samples.species.ra.nitrifiers <- subset_taxa(phyloseq.samples.species.ra,
                                                                             #AOB
                                                                             Family == "Nitrosomonadaceae" | # AOB; Includes Genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                                                                               Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                                                                               #AOA
                                                                               Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). In GTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                                                               Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes Species Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                                                               Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                                                               Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                                                                               #NOB
                                                                               #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                                                               Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                                                               Family == "NS-4" | # NOB (mostly uncultured genomes)
                                                                               Family == "UBA8639" | # NOB
                                                                               #In the order Rhizobiales and family Xanthobacteraceae
                                                                               Genus == "Nitrobacter" | #NOB; most representative one
                                                                               #Order Burkholderiales, family Gallionellaceae
                                                                               Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                                                               #Order Nitrococcales, family Nitrococcaceae
                                                                               Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                                                               #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                                                               Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                                                                               # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                                                                               Genus == "Nitrolancea")  #NOB; new one
phyloseq.samples.species.ra.nitrifiers <- subset_samples(phyloseq.samples.species.ra.nitrifiers, 
                                                                                sample_sums(phyloseq.samples.species.ra.nitrifiers) > 0)
phyloseq.samples.species.ra.nitrifiers #1690 nitrifying Species in 204 samples 

#Melt the ps object to long format 
phyloseq.samples.species.ra.nitrifiers.melt <- 
  psmelt(phyloseq.samples.species.ra.nitrifiers)

#Which Species? 
length(unique(phyloseq.samples.species.ra.nitrifiers.melt$Species))#OK, 1688

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.samples.species.ra.nitrifiers.melt <- phyloseq.samples.species.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Species) %>%
  mutate(Species = factor(Species, levels = unique(Species)))
levels(phyloseq.samples.species.ra.nitrifiers.melt$Species) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on Species within each ammonia-nitrite oxidizing group
palette_nitrifiers_all_taxa_species_df <- phyloseq.samples.species.ra.nitrifiers.melt %>% 
  distinct(Species, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Species) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    shade = seq(0, 0.8, length.out = n()),
    #shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_all_taxa_species_df

#Set up final palette
palette_nitrifiers_all_taxa_species <- setNames(
  palette_nitrifiers_all_taxa_species_df$color,
  palette_nitrifiers_all_taxa_species_df$Species)
palette_nitrifiers_all_taxa_species

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_species_all_taxa <- top_taxa_legend(phyloseq.samples.species.ra.nitrifiers.melt,
                                                             taxlevel = "Species",
                                                             n = 15)
top_nitrifying_species_all_taxa
#Keep the factor levels for legend
nitifier_species_levels_all_taxa <- levels(phyloseq.samples.species.ra.nitrifiers.melt$Species)
top_nitrifying_species_all_taxa <- 
  nitifier_species_levels_all_taxa[nitifier_species_levels_all_taxa %in% top_nitrifying_species_all_taxa]
top_nitrifying_species_all_taxa

#Plot
RA_nitrifiers_species_all_taxa_plot_datenum <- ggplot(phyloseq.samples.species.ra.nitrifiers.melt,
                                                                aes(x=factor(Date_num), 
                                                                    y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_all_taxa_species,
                    labels = function(x) str_wrap(x, width = 20), 
                    breaks = top_nitrifying_species_all_taxa
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.5),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_species_all_taxa_plot_datenum


###WITHIN BACTERIAL AND ARCHAEAL TAXA ####
#### GENUS LEVEL#######
phyloseq.samples.bacteria.genus.ra #41,284 genera

##Which genera are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
length(unique(phyloseq.samples.nitrifiers.melt$Genus)) #162 Nitrifying genera

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq.samples.bacteria.genus.ra.nitrifiers <- subset_taxa(phyloseq.samples.bacteria.genus.ra,
                                                    #AOB
                                                    Family == "Nitrosomonadaceae" | # AOB; Includes genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                                                      Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                                                      #AOA
                                                      Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). InGTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                                      Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes genera Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                                      Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                                      Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                                                      #NOB
                                                      #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                                      Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                                      Family == "NS-4" | # NOB (mostly uncultured genomes)
                                                      Family == "UBA8639" | # NOB
                                                      #In the order Rhizobiales and family Xanthobacteraceae
                                                      Genus == "Nitrobacter" | #NOB; most representative one
                                                      #Order Burkholderiales, family Gallionellaceae
                                                      Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                                      #Order Nitrococcales, family Nitrococcaceae
                                                      Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                                      #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                                      Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                                                      # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                                                      Genus == "Nitrolancea")  #NOB; new one
phyloseq.samples.bacteria.genus.ra.nitrifiers <- subset_samples(phyloseq.samples.bacteria.genus.ra.nitrifiers, 
                                                       sample_sums(phyloseq.samples.bacteria.genus.ra.nitrifiers) > 0)
phyloseq.samples.bacteria.genus.ra.nitrifiers #162 nitrifying genera in 204 samples 

#Melt the ps object to long format 
phyloseq.samples.bacteria.genus.ra.nitrifiers.melt <- 
  psmelt(phyloseq.samples.bacteria.genus.ra.nitrifiers)

#Which genera? 
unique(phyloseq.samples.bacteria.genus.ra.nitrifiers.melt$Genus)#OK, 162

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.samples.bacteria.genus.ra.nitrifiers.melt <- phyloseq.samples.bacteria.genus.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Genus) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus)))
levels(phyloseq.samples.bacteria.genus.ra.nitrifiers.melt$Genus) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on genera within each ammonia-nitrite oxidizing group
palette_nitrifiers_bacarc_genus_df <- phyloseq.samples.bacteria.genus.ra.nitrifiers.melt %>% 
  distinct(Genus, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Genus) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    shade = seq(0, 0.8, length.out = n()),
    #shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_bacarc_genus_df

#Set up final palette
palette_nitrifiers_bacarc_genus <- setNames(
  palette_nitrifiers_bacarc_genus_df$color,
  palette_nitrifiers_bacarc_genus_df$Genus)
palette_nitrifiers_bacarc_genus

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_genera_bacarc <- top_taxa_legend(phyloseq.samples.bacteria.genus.ra.nitrifiers.melt,
                                                  taxlevel = "Genus",
                                                  n = 10)
top_nitrifying_genera_bacarc
#Keep the factor levels for legend
nitifier_genera_levels_bacarc <- levels(phyloseq.samples.bacteria.genus.ra.nitrifiers.melt$Genus)
top_nitrifying_genera_bacarc <- 
  nitifier_genera_levels_bacarc[nitifier_genera_levels_bacarc %in% top_nitrifying_genera_bacarc]
top_nitrifying_genera_bacarc

#Plot
RA_nitrifiers_genus_bacarc_plot_datenum <- ggplot(phyloseq.samples.bacteria.genus.ra.nitrifiers.melt,
                                                    aes(x=factor(Date_num), 
                                                        y= Abundance, fill = Genus)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_bacarc_genus,
                    labels = function(x) str_wrap(x, width = 20), 
                    breaks = top_nitrifying_genera_bacarc
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.5),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_genus_bacarc_plot_datenum

#### SPECIES LEVEL#######
phyloseq.samples.bacteria.species.ra #309747 Species

##Which Species are nitrifiers? 
phyloseq.samples.nitrifiers.melt <- psmelt(phyloseq.samples.nitrifiers)
length(unique(phyloseq.samples.nitrifiers.melt$Species)) #1688 Nitrifying Species

#Out of the overall (all taxa and unclassified reads) phyloseq object, select only those from the list of nitrifiers
phyloseq.samples.bacteria.species.ra.nitrifiers <- subset_taxa(phyloseq.samples.bacteria.species.ra,
                                                      #AOB
                                                      Family == "Nitrosomonadaceae" | # AOB; Includes Genera Nitrosomonas, Nitrosospira, BJGV01 (some NCBI Nitrosovibrio), VFJL01 (some NCBI Nitrosovibrio)
                                                        Family == "Nitrosococcaceae" | #AOB; Includes Nitrosacidococcus, Nitrosococcus, Nitrosoglobus
                                                        #AOA
                                                        Family == "Nitrosocaldaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a). In GTDB 232, Candidatus Nitrosocaldales are all classified as Nitrosocaldaceae
                                                        Family == "Nitrosopumilaceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 a-associated; includes Species Nitrosopumilus, Nitrosotalea, Nitrosoarchaeum)
                                                        Family == "Nitrososphaeraceae" | #AOA Family within Nitrososphaerales (Thaumarchaeota - Group 1.1 b; includes Nitrososphaera)
                                                        Family == "UBA213" | #AOA, This is "Candidatus Nitrosomirales" from NCBI, in GTDB it is within Family UBA213. Though rarely present in ocean metagenomes. 
                                                        #NOB
                                                        #In the Order Nitrospirales, GTDB families Nitrospiraceae, NS-4 and UBA8639 have gene nxr (nitrite oxidizers), but UBA2166 doesnt. 
                                                        Family == "Nitrospiraceae" | # NOB/Commamox; predominant deep-sea waters, and sediments
                                                        Family == "NS-4" | # NOB (mostly uncultured genomes)
                                                        Family == "UBA8639" | # NOB
                                                        #In the order Rhizobiales and family Xanthobacteraceae
                                                        Genus == "Nitrobacter" | #NOB; most representative one
                                                        #Order Burkholderiales, family Gallionellaceae
                                                        Genus == "Nitrotoga" | #NOB Within Family == "Gallionellaceae" 
                                                        #Order Nitrococcales, family Nitrococcaceae
                                                        Genus == "Nitrococcus" | #NOB within Family Nitrococcaceae
                                                        #Within the phylum Nitrospinota in GTDB, only class is Nitrospinia. GTDB also has Nitrospinota_A phylum, but these may not be nitrite oxidizers. Also, within Nitrospinia class there are other Orders that may not be nitrite oxidizers
                                                        Order == "Nitrospinales" | # NOB; most rep Nitrospina; predominant marine NOB. Includes family Nitrosopinaceae, Nitrospina-like genommes, including Candidatus Nitromaritima and Nitrohelix 
                                                        # Phylum Chloroflexi, Class Thermomicrobia; Order Thermomicrobiales 
                                                        Genus == "Nitrolancea")  #NOB; new one
phyloseq.samples.bacteria.species.ra.nitrifiers <- subset_samples(phyloseq.samples.bacteria.species.ra.nitrifiers, 
                                                         sample_sums(phyloseq.samples.bacteria.species.ra.nitrifiers) > 0)
phyloseq.samples.bacteria.species.ra.nitrifiers #1690 nitrifying Species in 204 samples 

#Melt the ps object to long format 
phyloseq.samples.bacteria.species.ra.nitrifiers.melt <- 
  psmelt(phyloseq.samples.bacteria.species.ra.nitrifiers)

#Which Species? 
length(unique(phyloseq.samples.bacteria.species.ra.nitrifiers.melt$Species))#OK, 1688

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.samples.bacteria.species.ra.nitrifiers.melt <- phyloseq.samples.bacteria.species.ra.nitrifiers.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Species) %>%
  mutate(Species = factor(Species, levels = unique(Species)))
levels(phyloseq.samples.bacteria.species.ra.nitrifiers.melt$Species) #OK

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  
  AOB = "#1E88E5",  
  NOB = "#FFC107")

#Make hues based on Species within each ammonia-nitrite oxidizing group
palette_nitrifiers_bacarc_species_df <- phyloseq.samples.bacteria.species.ra.nitrifiers.melt %>% 
  distinct(Species, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Species) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    shade = seq(0, 0.8, length.out = n()),
    #shade = seq(0.01, 0.7, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_bacarc_species_df

#Set up final palette
palette_nitrifiers_bacarc_species <- setNames(
  palette_nitrifiers_bacarc_species_df$color,
  palette_nitrifiers_bacarc_species_df$Species)
palette_nitrifiers_bacarc_species

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_species_bacarc <- top_taxa_legend(phyloseq.samples.bacteria.species.ra.nitrifiers.melt,
                                                   taxlevel = "Species",
                                                   n = 15)
top_nitrifying_species_bacarc
#Keep the factor levels for legend
nitifier_species_levels_bacarc <- levels(phyloseq.samples.bacteria.species.ra.nitrifiers.melt$Species)
top_nitrifying_species_bacarc <- 
  nitifier_species_levels_bacarc[nitifier_species_levels_bacarc %in% top_nitrifying_species_bacarc]
top_nitrifying_species_bacarc

#Plot
RA_nitrifiers_species_bacarc_plot_datenum <- ggplot(phyloseq.samples.bacteria.species.ra.nitrifiers.melt,
                                                      aes(x=factor(Date_num), 
                                                          y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "identity", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.01,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(values = palette_nitrifiers_bacarc_species,
                    labels = function(x) str_wrap(x, width = 20), 
                    breaks = top_nitrifying_species_bacarc
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.5),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = 'black'),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_species_bacarc_plot_datenum


##NITROSOPUMILACEAE RA ACROSS PHASES ########
###WITHIN ALL TAXA INCLUDING UNCLASSIFIED READS######
phyloseq_with_unclassifieds_otu.samples.family.ra #23,427 families

#Out of this overall communities object, select only Nitrosopumilaceae
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae <- subset_taxa(phyloseq_with_unclassifieds_otu.samples.family.ra, 
                                                                    Family == "Nitrosopumilaceae") 
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae <- subset_samples(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae, 
                                                                       sample_sums(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae) > 0)
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae #Nitrosopumilaceae (1 taxa) in 204 samples 

#Melt to long format
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt <- psmelt(
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae)

#Going to make a single column for abbreviated phases, so I can color by this column 
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt <-  
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt%>%
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
                                                     "L_Established", 
                                                     "T1_Established", 
                                                     "E_Established", 
                                                     "P_Established"
                                                     )))
#Labels for the plot
labels_nitrosopumilaceae <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt %>%
  distinct(Date_num_colors_systems, Date_num_phase_abbrv) %>%
  mutate(Date_num_phase_abbrv = as.character(Date_num_phase_abbrv)) %>%
  tibble::deframe()

#Get wilcox stats
#Have to add another variable for each system to get their respective stats 
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_2 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt%>%
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
         
#Get stats for naive system - only comparisons in order of time (L vs T1, T1 vs E1, E1 vs T2, T2 vs T3, T3 vs E2)
stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_H21 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_2 %>%
  filter(Enclosure == "H21")%>%
  wilcox_test(
    Abundance ~ Date_num_phase_naive_abbrv_stats,
    comparisons = list(
      c("L_Naive", "T1_Naive"),
      c("T1_Naive", "E1_Naive"),
      c("E1_Naive", "T2_Naive"),
      c("T2_Naive", "T3_Naive"),
      c("T3_Naive", "E2_Naive")
    ),
    p.adjust.method = "BH"
  ) %>%
  filter(p.adj <= 0.05) 

#Add y position on plot for p values
stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_df_H21 <- stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_H21 %>%
  add_xy_position(x = "Date_num_phase_naive_abbrv_stats", 
                  step.increase = 0.04, 
                  scales = "free"
                  )%>%
  mutate(System = "Naive")

#Get stats for established system - only comparisons in order of time (L vs T1, T1 vs E, E vs P)
stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_P1 <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_2 %>%
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
#Add y position on plot for p values
stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_df_P1 <- stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_P1 %>%
  add_xy_position(x = "Date_num_phase_established_abbrv_stats",
                  step.increase = 0.04)%>%
  mutate(System = "Established")

#Naive and established wilcox results in same df
stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_df_both_H21_P1 <- bind_rows(stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_df_P1, 
                                                               stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_df_H21)%>%
  mutate(System = factor(System, levels = c("Naive", "Established")))

#Plot
nitrosopumilaceae.RA.phases.boxplot <- ggplot(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt, 
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
  stat_pvalue_manual(stat_wilcox_test_nitrosopumilaceae_with_unclassifieds_df_both_H21_P1,
                     label = "p.adj.signif",
                     hide.ns = TRUE,
                     tip.length = 0.01)
nitrosopumilaceae.RA.phases.boxplot


##GAM Generalized Additive Model - Nitrosopumilaceae RA (within all taxa including unclassified reads) modelled by phase####
#What do the relative abundance distributions look like
ggplot(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt, aes(x = Abundance)) +
  geom_histogram()

##Modelling each system (naive, established) separately
#Established
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt.P1 <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt %>%
  filter(Enclosure == "P1")%>%
  mutate(Date_num_phase_established_abbrv_model = factor(Date_num_phase_established_abbrv, 
                                                         levels = c("P", "L", "T1", "E"))
  )
#Naive
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt.H21 <-  
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt %>%
  filter(Enclosure == "H21")%>%
  mutate(Date_num_phase_naive_abbrv_model = factor(Date_num_phase_naive_abbrv, 
                                                   levels = c("E2", "L", "T1", "E1", "T2", "T3")))

###Naive system#####
#Model - H21
####Comparisons between phases #####
#In the model, E2 is the reference
gam_model_nitrosopumilaceae_H21_phases <- gam(Abundance ~ Date_num_phase_naive_abbrv_model,
                                             family = tw(link = "log"), #log link, so it can't model negative abundances. Also tweedie distribution
                                             data = phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt.H21)
gam_model_nitrosopumilaceae_H21_phases_summary <- (summary(gam_model_nitrosopumilaceae_H21_phases))$p.table  

#Get coefficients with SE
gam_model_nitrosopumilaceae_H21_phases_results <- data.frame(
  Term = rownames(gam_model_nitrosopumilaceae_H21_phases_summary),
  Estimate = gam_model_nitrosopumilaceae_H21_phases_summary[, "Estimate"],
  SE = gam_model_nitrosopumilaceae_H21_phases_summary[, "Std. Error"],
  p = gam_model_nitrosopumilaceae_H21_phases_summary[, "Pr(>|t|)"]
)

#Get CI (log scale)
gam_model_nitrosopumilaceae_H21_phases_results$Lower_log <- gam_model_nitrosopumilaceae_H21_phases_results$Estimate - 1.96 * gam_model_nitrosopumilaceae_H21_phases_results$SE
gam_model_nitrosopumilaceae_H21_phases_results$Upper_log <- gam_model_nitrosopumilaceae_H21_phases_results$Estimate + 1.96 * gam_model_nitrosopumilaceae_H21_phases_results$SE

#Get CI (original scale)
gam_model_nitrosopumilaceae_H21_phases_results$`Fold Change` <- exp(gam_model_nitrosopumilaceae_H21_phases_results$Estimate)
gam_model_nitrosopumilaceae_H21_phases_results$`Lower CI` <- exp(gam_model_nitrosopumilaceae_H21_phases_results$Lower_log)
gam_model_nitrosopumilaceae_H21_phases_results$`Upper CI` <- exp(gam_model_nitrosopumilaceae_H21_phases_results$Upper_log)

#Final df (coefficients) for comparisons between phases
gam_model_nitrosopumilaceae_H21_phases_results <- gam_model_nitrosopumilaceae_H21_phases_results%>%
  filter(Term != "(Intercept)")%>%
  mutate(Comparison = case_when(Term == "Date_num_phase_naive_abbrv_modelL" ~ "E2 vs L",
                                Term == "Date_num_phase_naive_abbrv_modelT1" ~ "E2 vs T1",
                                Term == "Date_num_phase_naive_abbrv_modelE1" ~ "E2 vs E1", 
                                Term == "Date_num_phase_naive_abbrv_modelT2" ~ "E2 vs T2", 
                                Term == "Date_num_phase_naive_abbrv_modelT3" ~ "E2 vs T3"
                                ))%>%
  select(Comparison, `Fold Change`, `Lower CI`, `Upper CI`, p)
gam_model_nitrosopumilaceae_H21_phases_results

####Estimates per phase #####
#First,get the phases as a df
gam_nitrosopumilaceae_H21_phases <- data.frame(
  Date_num_phase_naive_abbrv_model =
    levels(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt.H21$Date_num_phase_naive_abbrv_model)
)

#Now, get the estimate fit for each phase
prediction_GAM_nitrosopumilaceae_H21_phases <- predict(
  gam_model_nitrosopumilaceae_H21_phases,
  newdata = gam_nitrosopumilaceae_H21_phases,
  type = "link",
  se.fit = TRUE
)

#Get the estimate (link scale is log)
gam_nitrosopumilaceae_H21_phases$fit_link <- prediction_GAM_nitrosopumilaceae_H21_phases$fit

#Get the estimate on original scale (exp the fit_link)
gam_nitrosopumilaceae_H21_phases$fit <- exp(gam_nitrosopumilaceae_H21_phases$fit_link)

#Get the lower CI (exp the calculation based on log link) for the predicted value
gam_nitrosopumilaceae_H21_phases$lower_CI <-
  exp(gam_nitrosopumilaceae_H21_phases$fit_link - 1.96 * prediction_GAM_nitrosopumilaceae_H21_phases$se.fit)

#Get the upper CI (exp the calculation based on log link) for the predicted value
gam_nitrosopumilaceae_H21_phases$upper_CI <-
  exp(gam_nitrosopumilaceae_H21_phases$fit_link + 1.96 * prediction_GAM_nitrosopumilaceae_H21_phases$se.fit)
gam_nitrosopumilaceae_H21_phases

#Now, adding a column for system
gam_nitrosopumilaceae_H21_phases <- gam_nitrosopumilaceae_H21_phases%>%
  mutate(System = "Naive")


###Established system#####
#Model - P1
####Comparisons between phases #####
#In the model, P is the reference
gam_model_nitrosopumilaceae_P1_phases <- gam(Abundance ~ Date_num_phase_established_abbrv_model,
                                              family = tw(link = "log"), #log link, so it can't model negative abundances. Also tweedie distribution
                                              data = phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt.P1)
gam_model_nitrosopumilaceae_P1_phases_summary <- (summary(gam_model_nitrosopumilaceae_P1_phases))$p.table  

#######Get coefficients with SE
gam_model_nitrosopumilaceae_P1_phases_results <- data.frame(
  Term = rownames(gam_model_nitrosopumilaceae_P1_phases_summary),
  Estimate = gam_model_nitrosopumilaceae_P1_phases_summary[, "Estimate"],
  SE = gam_model_nitrosopumilaceae_P1_phases_summary[, "Std. Error"],
  p = gam_model_nitrosopumilaceae_P1_phases_summary[, "Pr(>|t|)"]
)

#Get CI (log scale)
gam_model_nitrosopumilaceae_P1_phases_results$Lower_log <- 
  gam_model_nitrosopumilaceae_P1_phases_results$Estimate - 1.96 * gam_model_nitrosopumilaceae_P1_phases_results$SE
gam_model_nitrosopumilaceae_P1_phases_results$Upper_log <- 
  gam_model_nitrosopumilaceae_P1_phases_results$Estimate + 1.96 * gam_model_nitrosopumilaceae_P1_phases_results$SE

#Get CI (original scale)
gam_model_nitrosopumilaceae_P1_phases_results$`Fold Change` <- exp(gam_model_nitrosopumilaceae_P1_phases_results$Estimate)
gam_model_nitrosopumilaceae_P1_phases_results$`Lower CI` <- exp(gam_model_nitrosopumilaceae_P1_phases_results$Lower_log)
gam_model_nitrosopumilaceae_P1_phases_results$`Upper CI` <- exp(gam_model_nitrosopumilaceae_P1_phases_results$Upper_log)

#Final df (coefficients)
gam_model_nitrosopumilaceae_P1_phases_results <- gam_model_nitrosopumilaceae_P1_phases_results%>%
  filter(Term != "(Intercept)")%>%
  mutate(Comparison = case_when(Term == "Date_num_phase_established_abbrv_modelL" ~ "P vs L",
                                Term == "Date_num_phase_established_abbrv_modelT1" ~ "P vs T1",
                                Term == "Date_num_phase_established_abbrv_modelE" ~ "P vs E"))%>%
  select(Comparison, `Fold Change`, `Lower CI`, `Upper CI`, p)
gam_model_nitrosopumilaceae_P1_phases_results

####Estimates per phase #####
#First, get a dataframe with the phases
gam_nitrosopumilaceae_P1_phases <- data.frame(
  Date_num_phase_established_abbrv_model =
    levels(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt.P1$Date_num_phase_established_abbrv_model)
)

#Get the estimate for each phase
prediction_GAM_nitrosopumilaceae_P1_phases <- predict(
  gam_model_nitrosopumilaceae_P1_phases,
  newdata = gam_nitrosopumilaceae_P1_phases,
  type = "link",
  se.fit = TRUE
)

#Get the estimate for each phase (link scale is log)
gam_nitrosopumilaceae_P1_phases$fit_link <- prediction_GAM_nitrosopumilaceae_P1_phases$fit

#Get the estimate for each phase on original scale (exp the fit_link)
gam_nitrosopumilaceae_P1_phases$fit <- exp(gam_nitrosopumilaceae_P1_phases$fit_link)

#Get the upper CI (exp the calculation based on log link)  for the predicted value
gam_nitrosopumilaceae_P1_phases$lower_CI <-
  exp(gam_nitrosopumilaceae_P1_phases$fit_link - 1.96 * prediction_GAM_nitrosopumilaceae_P1_phases$se.fit)

#Get the lower CI (exp the calculation based on log link)  for the predicted value
gam_nitrosopumilaceae_P1_phases$upper_CI <-
  exp(gam_nitrosopumilaceae_P1_phases$fit_link + 1.96 * prediction_GAM_nitrosopumilaceae_P1_phases$se.fit)
gam_nitrosopumilaceae_P1_phases

#Now add a column for system
gam_nitrosopumilaceae_P1_phases <- gam_nitrosopumilaceae_P1_phases%>%
  mutate(System = "Established")

###Plot GAM Phase comparisons as boxplot#######
#Get GAM coefficients for both naive and established systems
gam_model_nitrosopumilaceae_phases_results <- bind_rows(
  gam_model_nitrosopumilaceae_P1_phases_results%>%mutate(System = 'Established'), 
  gam_model_nitrosopumilaceae_H21_phases_results%>%mutate(System = 'Naive')
  )%>%
  mutate(Phase = word(Comparison, sep = " ", -1), 
         Date_num_phase_abbrv_model = paste0(Phase, "_", System))

#Add annotations according to p values 
gam_model_nitrosopumilaceae_phases_results$label_signif <- case_when(
  gam_model_nitrosopumilaceae_phases_results$p < 0.001 ~ "***",
  gam_model_nitrosopumilaceae_phases_results$p < 0.01 ~ "**",
  gam_model_nitrosopumilaceae_phases_results$p < 0.05 ~ "*",
  TRUE ~ "ns"
)

#Get GAM estimates for both naive and established system per phase
gam_nitrosopumilaceae_estimates_phases <- bind_rows(gam_nitrosopumilaceae_P1_phases, 
                                                    gam_nitrosopumilaceae_H21_phases)%>%
  mutate(Date_num_phase_abbrv_model = coalesce(Date_num_phase_established_abbrv_model, Date_num_phase_naive_abbrv_model))%>%
  mutate(Date_num_phase_abbrv_model = ifelse(System == "Established", 
                                             paste0(Date_num_phase_abbrv_model, "_", "Established"),
                                             paste0(Date_num_phase_abbrv_model, "_", "Naive")))

phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_2 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt%>%
  mutate(Date_num_phase_abbrv_model = ifelse(Enclosure == "P1", 
                                           paste0(Date_num_phase_abbrv, "_", "Established"),
                                           paste0(Date_num_phase_abbrv, "_", "Naive")))

#Join with Nitrosopumilaceae melted DF
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_model <- 
  merge(
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_2, 
  gam_nitrosopumilaceae_estimates_phases, 
  by = c("Date_num_phase_abbrv_model", "System")) %>%
  mutate(Date_num_phase_abbrv_model = factor(Date_num_phase_abbrv_model, 
                                             levels = c(
                                               "L_Naive", 
                                               "T1_Naive",
                                               "E1_Naive", 
                                               "T2_Naive", 
                                               "T3_Naive", 
                                               "E2_Naive",
                                               "L_Established", 
                                               "T1_Established", 
                                               "E_Established",
                                               "P_Established")))


#Plot
nitrosopumilaceae.RA.phases.GAM_model_plot <- ggplot(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_model, 
                                              aes(x = Date_num_phase_abbrv_model, 
                                                  y= Abundance, 
                                                  color = Date_num_phase_abbrv_model, 
                                                  fill = Date_num_phase_abbrv_model)) +
  theme_bw() +
  facet_wrap(~System, 
             scales = "free")+
  labs(x = "Phase", 
       y = expression(italic(Nitrosopumilaceae) ~ "RA (%)"),
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.1, 
               linewidth = 0.2) +
  geom_point(aes(
      x = Date_num_phase_abbrv_model,
      y = fit,
      color = Date_num_phase_abbrv_model,
      group = System),
    size = 3,
    shape = 16
  ) +
  geom_errorbar(
    aes(
      x = Date_num_phase_abbrv_model,
      ymin = lower_CI,
      ymax = upper_CI,
      color = Date_num_phase_abbrv_model,
      group = System
    ),
    width = 0.2,
    linewidth = 0.6
  )+
  geom_text(
    data = gam_model_nitrosopumilaceae_phases_results,
    aes(
      x = Date_num_phase_abbrv_model,
      #group = System,
      y = `Upper CI` * 6,   # pushes stars above CI
      label = label_signif
    ),
    size = 6,
    color = "black",
    fontface = "bold"
  )+
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
    axis.ticks = element_line(colour = "black", linewidth = 0.5)) 
nitrosopumilaceae.RA.phases.GAM_model_plot


##GAM Generalized Additive Model - Nitrosopumilaceae RA (within all taxa including unclassified reads) modeled with spline for each phase for copper levels####

#First, doing imputations to fill Copper levels where they were not measured 
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt %>%
  arrange(Enclosure, Date_num)

#How many samples are missing copper levels: 
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
filter(is.na(Copper_level_mg_L)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure         Date_num_phase                      n
# H21 (naive)       Low Copper Levels (Day 1-27)        9
# H21 (naive)       Transition Period 2 (Day 52-81)     3
# P1 (established)    Low Copper Levels (Day 1-53)        1
# P1 (established)    Tx Copper Exposure (Day 66-104)     1

#pH?
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  filter(is.na(pH_spu)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure Date_num_phase                           n
# H21       First Tx Copper Exposure (Day 39-51)     1
# P1        Low Copper Levels (Day 1-53)             1

#Salinity?
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  filter(is.na(Salinity_ppt)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure Date_num_phase                              n
# H21       First Tx Copper Exposure (Day 39-51)        1
# H21       Transition Period 3 (Day 86-108)            1
# H21       Second Tx Copper Exposure (Day 109-135)     1
# P1        Low Copper Levels (Day 1-53)                1

#Temperature?
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  filter(is.na(Temperature_F)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure Date_num_phase                           n
# H21       First Tx Copper Exposure (Day 39-51)     1
# P1        Low Copper Levels (Day 1-53)             3
# P1        Tx Copper Exposure (Day 66-104)          2
# P1        Post-Treatment (Day 105-141)             2

#Add a column that will have the measured copper values as well as the imputed ones
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$Copper_level_mg_L_imp <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$Copper_level_mg_L

#Using na.approx() to fill copper missing values, fills them by drawing a straight line between known points and estimating values along that line.
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  group_by(Enclosure) %>%
  mutate(
    Copper_level_mg_L_imp = zoo::na.approx(
      Copper_level_mg_L_imp,
      x = Date_num,
      na.rm = FALSE), 
      Copper_imputation = case_when(
        is.na(Copper_level_mg_L) & !is.na(Copper_level_mg_L_imp) ~ "yes",
        TRUE ~ "no"
      ) #Fill in copper with linear interpolation. If copper level imputed, 'yes' on Copper_imputation
    ) %>%
      ungroup()

#Using na.approx() to fill Temperature, salinity, and pH
#Add a column for pH
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$pH_spu_imp <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$pH_spu
#Add a column for Temperature
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$Temperature_F_imp <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$Temperature_F
#Add a column for Salinity
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$Salinity_ppt_imp <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date$Salinity_ppt

#Now impute those (Temperature, salinity, and pH)
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  group_by(Enclosure) %>%
  mutate(
    pH_spu_imp = zoo::na.approx(
      pH_spu_imp,
      x = Date_num,
      na.rm = FALSE), 
    Temperature_F_imp = zoo::na.approx(
      Temperature_F_imp,
      x = Date_num,
      na.rm = FALSE), 
    Salinity_ppt_imp = zoo::na.approx(
      Salinity_ppt_imp,
      x = Date_num,
      na.rm = FALSE)
  ) %>%
  ungroup()

#How does it look now? - Copper
ggplot(phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date, 
       aes(x = Date_num, 
            y = Copper_level_mg_L_imp))+
  facet_grid(~Enclosure)+
  geom_point() #OK

#Fitting models for each system separately
#Established
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  filter(Enclosure == "P1")%>%
  mutate(Date_num_phase_established_abbrv_model = factor(Date_num_phase_established_abbrv, 
                                                         levels = c("P", "L", "T1", "E"))
  )%>%
  arrange(Date_num) #Make sure they are arranged by date

#Naive
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 <-  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date %>%
  filter(Enclosure == "H21")%>%
  mutate(Date_num_phase_naive_abbrv_model = factor(Date_num_phase_naive_abbrv, 
                                                   levels = c("E2", "L", "T1", "E1", "T2", "T3")))%>%
  arrange(Date_num) #Make sure they are arranged by date

###Naive system#####
#Model - H21
#Is the copper–abundance relationship (for each phase) different from zero?”
gam_model_nitrosopumilaceae_H21_phases_spline <- gam(
  Abundance ~ 
    Date_num_phase_naive_abbrv_model +
    s(Copper_level_mg_L_imp, by = Date_num_phase_naive_abbrv_model, k = 3)+
    s(pH_spu_imp, k = 5) +
    s(Salinity_ppt_imp, k = 5) +
    s(Temperature_F_imp, k = 5),
  family = tw(link = "log"), #log link, so it can't model negative abundances. Also tweedie distribution
  data = phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21,
  method = "REML"
)
summary(gam_model_nitrosopumilaceae_H21_phases_spline)
gam.check(gam_model_nitrosopumilaceae_H21_phases_spline)


#Check residual autocorrelation (Did the model adequately remove temporal structure?)
#First, get residuals
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21$resid_gam <-
  residuals(gam_model_nitrosopumilaceae_H21_phases_spline)

#Now, check for autocorrelation in the whole time series in t he naive system
autocorrelation_table_whole_series_gam_h21_nitrosopumilaceae <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
  arrange(Date_num_naive) %>%   # ensure ordering
  summarize(
    System = "Naive",
    n = n(),
    lag1_acf = acf(resid_gam, plot = FALSE, lag.max = 1)$acf[2],
    # ljung_p_val = Box.test(resid_gam, lag = 1, type = "Ljung-Box")$p.value
  )
autocorrelation_table_whole_series_gam_h21_nitrosopumilaceae #Correlation coefficient r = 0.07 (n = 87). Almost negligible correlation betweeen each resifual and the one before

# #Now, check for autocorrelation per phase
# autocorrelation_table_per_phase_gam_h21_nitrosopumilaceae <- 
#   phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
#   arrange(Date_num_naive) %>%                            # chronological within each group
#   group_by(Date_num_phase_naive_abbrv_model) %>%
#   summarize(
#     n = n(),
#     lag1_acf = acf(resid_gam, plot = FALSE, lag.max = 1)$acf[2],
#     ljung_p_val = Box.test(resid_gam, lag = 1, type = "Ljung-Box")$p.value
#   )%>%
#   rename(Date_num_phase_abbrv_model = Date_num_phase_naive_abbrv_model)
# autocorrelation_table_per_phase_gam_h21_nitrosopumilaceae
# #Put them together in a table
# autocorrelation_table_gam_h21_nitrosopumilaceae <- bind_rows(autocorrelation_table_whole_series_gam_h21_nitrosopumilaceae,
#                                                              autocorrelation_table_per_phase_gam_h21_nitrosopumilaceae,
#                                                              )
# autocorrelation_table_gam_h21_nitrosopumilaceae

#Build a prediction grid for plotting a GAM’s fitted copper effect, separately for each phase 
gam_model_nitrosopumilaceae_H21_phases_spline_data_copper <-
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
  group_by(Date_num_phase_naive_abbrv_model) %>%
  summarise(
    min = min(Copper_level_mg_L_imp, na.rm = TRUE),
    max = max(Copper_level_mg_L_imp, na.rm = TRUE),
    .groups = "drop"
  ) %>% #Finds the lowest and highest imputed copper concentration (mg/L) observed within each phase
  mutate(
    Copper_level_mg_L_imp = purrr::map2(min, max,
                                        ~ seq(.x, .y, length.out = 100))
  ) %>% #For each phase, get 100 evenly spaced copper values from that phase’s min to its max
  tidyr::unnest(Copper_level_mg_L_imp) %>% #Expand the list-column so each phase has 100 rows
  select(-min, -max)

#Get means of other environmental covariates included in the model by phase
gam_model_nitrosopumilaceae_H21_phases_spline_data_means_covariates <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
  group_by(Date_num_phase_naive_abbrv_model) %>%
  summarise(
    pH_spu_imp = mean(pH_spu_imp, na.rm = TRUE),
    Salinity_ppt_imp = mean(Salinity_ppt_imp, na.rm = TRUE),
    Temperature_F_imp = mean(Temperature_F_imp, na.rm = TRUE),
    .groups = "drop"
  )

#Now, join copper values with the covariates. Each of the 100 copper values in a phase then gets that phase’s mean pH, salinity, and temperature.
gam_model_nitrosopumilaceae_H21_phases_spline_data_to_predict <- 
  gam_model_nitrosopumilaceae_H21_phases_spline_data_copper %>%
  left_join(gam_model_nitrosopumilaceae_H21_phases_spline_data_means_covariates, 
            by = "Date_num_phase_naive_abbrv_model")

#Predict values 
gam_model_nitrosopumilaceae_H21_phases_spline_predicted <- 
  predict(gam_model_nitrosopumilaceae_H21_phases_spline, 
          newdata = gam_model_nitrosopumilaceae_H21_phases_spline_data_to_predict,
          se.fit = TRUE,
          type = "link") #Predictions are on the link scale (log)

#Join with prediction data, calculate lower and upper confidence intervals 
gam_predicted_H21_nitrosopumilaceae_copper_spline_df <- gam_model_nitrosopumilaceae_H21_phases_spline_data_to_predict %>%
  dplyr::mutate(fit = exp(gam_model_nitrosopumilaceae_H21_phases_spline_predicted$fit),
                lower = exp(gam_model_nitrosopumilaceae_H21_phases_spline_predicted$fit - 1.96 * 
                              gam_model_nitrosopumilaceae_H21_phases_spline_predicted$se.fit),
                upper = exp(gam_model_nitrosopumilaceae_H21_phases_spline_predicted$fit + 1.96 * 
                              gam_model_nitrosopumilaceae_H21_phases_spline_predicted$se.fit))%>%
  mutate(System = "Naive", 
         Date_num_phase_abbrv_model = Date_num_phase_naive_abbrv_model, 
         Date_num_phase_abbrv_plot = 
           paste0(Date_num_phase_naive_abbrv_model, 
                  "_", 
                  System))

####Emmeans compare E1 and E2#######
#Which range of copper levels do they share?
E1_copper_levels_H21 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
  filter(Date_num_phase_naive_abbrv_model == "E1") %>%
  pull(Copper_level_mg_L_imp)

E2_copper_levels_H21 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
  filter(Date_num_phase_naive_abbrv_model == "E2") %>%
  pull(Copper_level_mg_L_imp)

# Get the copper level ranges
range(E1_copper_levels_H21) #0.149 to 0.174
range(E2_copper_levels_H21) #0.145 to 0.210

#Get a grid of values for copper levels intersecting in E1 and E2
copper_grid_E1E2_naive <- seq(min(E1_copper_levels_H21), max(E1_copper_levels_H21), length.out = 50)  

#Get emmeans at levels shared between the two
gam_model_nitrosopumilaceae_H21_phases_spline_emmeans_E1_E2 <- emmeans(
  gam_model_nitrosopumilaceae_H21_phases_spline,
  ~ Date_num_phase_naive_abbrv_model,
  at = list(
    Date_num_phase_naive_abbrv_model = c("E1", "E2"),
    Copper_level_mg_L_imp = copper_grid_E1E2_naive #range of copper therapeutic levels shared by E1 and E2
  )
)

# E2 vs E1 at each copper value, on the log scale
E2_vs_E1_nitrosopumilaceae_H21_contrast_df <- 
  data.frame(contrast(gam_model_nitrosopumilaceae_H21_phases_spline_emmeans_E1_E2,
         method = "revpairwise", #Just so it reads as E2 - E1
         infer = TRUE, 
         type = 'link'))
E2_vs_E1_nitrosopumilaceae_H21_contrast_df

# #Plot the difference 
# ggplot(E2_vs_E1_nitrosopumilaceae_H21_contrast_df, 
#        aes(Copper_level_mg_L_imp, estimate)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL), alpha = 0.3) +
#   geom_line() +
#   labs(x = "Copper (mg/L)", y = "E2 minus E1 (log abundance)") +
#   theme_bw()



###Established system#####
#Model - P1
#Is the copper–abundance relationship (for each phase) different from zero?”
#“Does the copper–abundance relationship change across experimental phases?”
gam_model_nitrosopumilaceae_P1_phases_spline <- gam(
  Abundance ~ 
    Date_num_phase_established_abbrv_model +
    s(Copper_level_mg_L_imp, by = Date_num_phase_established_abbrv_model, k = 3)+
    s(pH_spu_imp, k = 5) +
    s(Salinity_ppt_imp, k = 5) +
    s(Temperature_F_imp, k = 5),
  family = tw(link = "log"), #log link, so it can't model negative abundances. Also tweedie distribution
  data = phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1,
  method = "REML"
)
summary(gam_model_nitrosopumilaceae_P1_phases_spline)

#Check residual autocorrelation (Did the model adequately remove temporal structure?)
#First, get residuals
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1$resid_gam <-
  residuals(gam_model_nitrosopumilaceae_P1_phases_spline)

#Now, check for autocorrelation in the whole time series
autocorrelation_table_whole_series_gam_p1_nitrosopumilaceae <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1 %>%
  arrange(Date_num_established) %>%   # ensure ordering by date number
  summarize(
    System = "Established",
    n = n(),
    lag1_acf = acf(resid_gam, plot = FALSE, lag.max = 1)$acf[2],
    #ljung_p_val = Box.test(resid_gam, lag = 1, type = "Ljung-Box")$p.value
    )
autocorrelation_table_whole_series_gam_p1_nitrosopumilaceae #Correlation coefficient r = 0.165 (n = 117). Very weak positive lag-1 autocorrelation. 
# 
# #Now, check for autocorrelation per phase
# autocorrelation_table_per_phase_gam_p1_nitrosopumilaceae <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1 %>%
#   filter(Date_num_phase_established_abbrv_model != "P") %>% #not inclusing P since it doesnt have enough samples
#   group_by(Date_num_phase_established_abbrv_model) %>%
#   summarize(
#     n = n(),
#     lag1_acf = acf(resid_gam, plot = FALSE, lag.max = 1)$acf[2],
#     ljung_p_val = Box.test(resid_gam, lag = 1, type = "Ljung-Box")$p.value
#   )%>%
#   rename(Date_num_phase_abbrv_model = Date_num_phase_established_abbrv_model)
# autocorrelation_table_per_phase_gam_p1_nitrosopumilaceae
# 
# #Put them together in a table
# autocorrelation_table_gam_p1_nitrosopumilaceae <- bind_rows(autocorrelation_table_whole_series_gam_p1_nitrosopumilaceae, 
#                                                              autocorrelation_table_per_phase_gam_p1_nitrosopumilaceae, 
#                                                              )
# autocorrelation_table_gam_p1_nitrosopumilaceae

#Build a prediction grid for plotting a GAM’s fitted copper effect, separately for each phase 
gam_model_nitrosopumilaceae_P1_phases_spline_data_copper <-
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1 %>%
  group_by(Date_num_phase_established_abbrv_model) %>%
  summarise(
    min = min(Copper_level_mg_L_imp, na.rm = TRUE),
    max = max(Copper_level_mg_L_imp, na.rm = TRUE),
    .groups = "drop"
  ) %>% #Finds the lowest and highest imputed copper concentration (mg/L) observed within each phase
  mutate(
    Copper_level_mg_L_imp = purrr::map2(min, max,
                                        ~ seq(.x, .y, length.out = 100))
  ) %>% #For each phase, get 100 evenly spaced copper values from that phase’s min to its max
  tidyr::unnest(Copper_level_mg_L_imp) %>% #Expand the lsit column so each phase has 100 rows
  select(-min, -max)

#Get means of other environmental covariates included in the model by phase
gam_model_nitrosopumilaceae_P1_phases_spline_data_means_covariates <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1 %>%
  group_by(Date_num_phase_established_abbrv_model) %>%
  summarise(
    pH_spu_imp = mean(pH_spu_imp, na.rm = TRUE),
    Salinity_ppt_imp = mean(Salinity_ppt_imp, na.rm = TRUE),
    Temperature_F_imp = mean(Temperature_F_imp, na.rm = TRUE),
    .groups = "drop"
  )

#Now, join copper values with the covariates. Each of the 100 copper values in a phase then gets that phase’s mean pH, salinity, and temperature.
gam_model_nitrosopumilaceae_P1_phases_spline_data_to_predict <- 
  gam_model_nitrosopumilaceae_P1_phases_spline_data_copper %>%
  left_join(gam_model_nitrosopumilaceae_P1_phases_spline_data_means_covariates, 
            by = "Date_num_phase_established_abbrv_model")

#Predict values 
gam_model_nitrosopumilaceae_P1_phases_spline_predicted <- predict(gam_model_nitrosopumilaceae_P1_phases_spline, 
                                                                   newdata = gam_model_nitrosopumilaceae_P1_phases_spline_data_to_predict,
                                                                   se.fit = TRUE,
                                                                   type = "link") #Predictions on the link scale (log)

#Join with prediction data
gam_predicted_P1_nitrosopumilaceae_copper_spline_df <- gam_model_nitrosopumilaceae_P1_phases_spline_data_to_predict %>%
  dplyr::mutate(fit = exp(gam_model_nitrosopumilaceae_P1_phases_spline_predicted$fit),
                lower = exp(gam_model_nitrosopumilaceae_P1_phases_spline_predicted$fit - 1.96 * 
                              gam_model_nitrosopumilaceae_P1_phases_spline_predicted$se.fit),
                upper = exp(gam_model_nitrosopumilaceae_P1_phases_spline_predicted$fit + 1.96 * 
                              gam_model_nitrosopumilaceae_P1_phases_spline_predicted$se.fit))%>%
  mutate(System = "Established", 
         Date_num_phase_abbrv_model = Date_num_phase_established_abbrv_model, 
         Date_num_phase_abbrv_plot = 
           paste0(Date_num_phase_established_abbrv_model, 
                  "_", 
                  System))


###PLOT FOR BOTH SYSTEMS#########
#Join GAM model data from both systems
gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_df <- bind_rows(gam_predicted_P1_nitrosopumilaceae_copper_spline_df, 
                                                                     gam_predicted_H21_nitrosopumilaceae_copper_spline_df)
gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_df

#Just adding these two "dummy" variables to get the plot to have 2 rows, top Naive and Bottom established
gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_df <-  
  gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_df%>%
  add_row(Date_num_phase_abbrv_model = "XX", 
          System = "Established", fit = 0, 
          Copper_level_mg_L_imp = 0, 
          Date_num_phase_abbrv_plot = "P_Established")%>%
  add_row(Date_num_phase_abbrv_model = "WW", 
          System = "Established", fit = 0, 
          Copper_level_mg_L_imp = 0, 
          Date_num_phase_abbrv_plot = "P_Established")%>%
  mutate(System = factor(System , levels = c("Naive", "Established")), 
         Date_num_phase_abbrv_plot = factor(Date_num_phase_abbrv_plot, 
                                            levels = c("L_Naive",
                                                       "L_Established",
                                                       "T1_Naive", 
                                                       "T1_Established",
                                                       "E1_Naive", 
                                                       "T2_Naive", 
                                                       "T3_Naive", 
                                                       "E2_Naive", 
                                                       "E_Established", 
                                                       "P_Established"
                                            )),
         Date_num_phase_abbrv_model = factor(Date_num_phase_abbrv_model, 
                                             levels = c("L", 
                                                        "T1", 
                                                        "E1", 
                                                        "T2", 
                                                        "T3", 
                                                        "E2", 
                                                        "E", 
                                                        "P", 
                                                        "WW", 
                                                        "XX"
                                             )))

#Also need this data frame to have the original relative abundance values (not predicted ones)
phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date_plot <- 
  phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date%>%
  mutate(System = ifelse(Enclosure == "H21", "Naive", "Established"), 
         Date_num_phase_abbrv_plot = 
           paste0(Date_num_phase_abbrv, 
                  "_", 
                  System), 
         Date_num_phase_abbrv_model = Date_num_phase_abbrv)%>%
  add_row(Date_num_phase_abbrv_model = "XX", 
          System = "Established", 
          Copper_level_mg_L = 0, 
          Date_num_phase_abbrv_plot = "P_Established")%>%
  add_row(Date_num_phase_abbrv_model = "WW", 
          System = "Established", 
          Copper_level_mg_L = 0, 
          Date_num_phase_abbrv_plot = "P_Established")%>%
  mutate(System = factor(System , levels = c("Naive", "Established")), 
         Date_num_phase_abbrv_model = factor(Date_num_phase_abbrv_model, 
                                             levels = c("L", 
                                                        "T1", 
                                                        "E1", 
                                                        "T2", 
                                                        "T3", 
                                                        "E2", 
                                                        "E", 
                                                        "P", 
                                                        "WW", 
                                                        "XX"
                                            )))

#PLOT IT 
gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_plot <- 
  ggplot(gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_df, 
       aes(x = Copper_level_mg_L_imp, 
           y = fit, 
           color = Date_num_phase_abbrv_plot, 
           fill = Date_num_phase_abbrv_plot)) +
  geom_line() +
  facet_nested_wrap(
    . ~ System + Date_num_phase_abbrv_model,
             scales = "free",
             nrow = 2)+
  geom_point(
    data = phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date_plot,
    aes(
      x = Copper_level_mg_L_imp,
      y = Abundance, 
      color = Date_num_phase_abbrv_plot,
      shape = Copper_imputation),
    inherit.aes = FALSE,
    alpha = 0.7,
    size = 5
  )+
  geom_ribbon(aes(ymin = lower, ymax = upper, 
                  fill = Date_num_phase_abbrv_plot),
              alpha = 0.2, color = NA)+
  #Dashed lines around the ribbon
  geom_line(
    aes(y = lower),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(
    aes(y = upper),
    linetype = "dashed",
    linewidth = 0.5
  )+
  labs(x = "Copper levels (mg/L)", y = expression(italic(Nitrosopumilaceae) ~ "RA (%)"))+
  scale_color_manual(values = phases_naive_established_palette)+
  scale_fill_manual(values = phases_naive_established_palette)+
  scale_x_continuous(labels = function(x) sprintf("%.3f", x))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", linewidth= 1),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(size = 36, color = "white", face = 'bold'),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
    axis.title= element_text(size = 32, colour = "black"),
    axis.text = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.5))
gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_plot

###FIGURE 3#######
figure3 <- gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_plot
ggsave("Figure3.png", 
       gam_predicted_P1_H21_nitrosopumilaceae_copper_spline_plot, 
       device = "png", 
       height = 12, 
       width = 36, 
       dpi = 500
       )



###SUPPLEMENTARY TABLE 8 & 9######
#Naive system stats
gam_model_nitrosopumilaceae_H21_phases_spline_summary <- summary(gam_model_nitrosopumilaceae_H21_phases_spline)
gam_model_nitrosopumilaceae_H21_phases_spline_summary$dev.expl #Model explains 81.3% of deviance
gam_model_nitrosopumilaceae_H21_phases_spline_summary$r.sq# Model explains 77.8% of the variance 

#Make a dataframe summary of F and P values (smooth terms)
gam_model_nitrosopumilaceae_H21_phases_spline_summary_table <- data.frame(gam_model_nitrosopumilaceae_H21_phases_spline_summary$s.table)
gam_model_nitrosopumilaceae_H21_phases_spline_summary_table

#Make a dataframe summary of F and P values (parametric terms)
gam_model_nitrosopumilaceae_H21_phases_parametric_summary_table <- data.frame(gam_model_nitrosopumilaceae_H21_phases_spline_summary$p.table, 
                                                                              check.names = F)
gam_model_nitrosopumilaceae_H21_phases_parametric_summary_table

#Number of samples per phase
n_by_phase_H21 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.H21 %>%
  count(Date_num_phase_abbrv, name = "No. Samples")

#Stable 8 (naive system) with parametric terms 
stable8.1 <- 
  gam_model_nitrosopumilaceae_H21_phases_parametric_summary_table %>%
  mutate(
    System = "Naive",
    Phase = case_when(
      str_detect(rownames(.), "modelL")  ~ "L",
      str_detect(rownames(.), "modelT1") ~ "T1",
      str_detect(rownames(.), "modelE1") ~ "E1",
      str_detect(rownames(.), "modelT2") ~ "T2",
      str_detect(rownames(.), "modelT3") ~ "T3",
      TRUE ~ NA_character_
    ),
    Contrast = case_when(
      str_detect(rownames(.), "modelL")  ~ "L vs E2",
      str_detect(rownames(.), "modelT1") ~ "T1 vs E2",
      str_detect(rownames(.), "modelE1") ~ "E1 vs E2",
      str_detect(rownames(.), "modelT2") ~ "T2 vs E2",
      str_detect(rownames(.), "modelT3") ~ "T3 vs E2",
      str_detect(rownames(.), "Intercept") ~ "Intercept",
      TRUE ~ NA_character_
    ),
    Effect = "Phase",
    SE = round(`Std. Error`, 3),
    `t value` = round(`t value`, 3),
    `p value` = format(`Pr(>|t|)`, scientific = TRUE, digits = 3),
  )%>%
  filter(!Contrast == "Intercept")%>%
  left_join(n_by_phase_H21, by = c("Phase" = "Date_num_phase_abbrv"))%>%
  #Contrast samples
  mutate(`No. Samples` = paste(`No. Samples`, "-", (n_by_phase_H21%>%filter(Date_num_phase_abbrv == "E2")%>%pull(`No. Samples`))))%>%
  select(System, Contrast, Effect, `No. Samples`, Estimate, SE, `t value`, `p value`)
stable8.1

#Stable 9 (naive system) with smooth terms 
stable9.1 <- 
  gam_model_nitrosopumilaceae_H21_phases_spline_summary_table %>%
  mutate(
    System = "Naive",
    Phase = case_when(
      str_detect(rownames(.), "modelL")  ~ "L",
      str_detect(rownames(.), "modelT1") ~ "T1",
      str_detect(rownames(.), "modelE1") ~ "E1",
      str_detect(rownames(.), "modelT2") ~ "T2",
      str_detect(rownames(.), "modelT3") ~ "T3",
      str_detect(rownames(.), "modelE2") ~ "E2",
      TRUE ~ NA_character_
    ),
    Effect = case_when(
      str_detect(rownames(.), "pH")  ~ "pH",
      str_detect(rownames(.), "Salinity") ~ "Salinity",
      str_detect(rownames(.), "Temperature") ~ "Temperature",
      str_detect(rownames(.), "model")  ~ "Copper"
    ),
    EDF = round(edf, 3),
    `Reference EDF` = round(Ref.df, 3),
    `F Statistic` = round(`F`, 3),
    `p value` = format(p.value, scientific = TRUE, digits = 3),
  )%>%
  left_join(n_by_phase_H21, by = c("Phase" = "Date_num_phase_abbrv"))%>%
  select(System, Phase, Effect, `No. Samples`, EDF, `F Statistic`, `p value`)
stable9.1

#Established system stats
gam_model_nitrosopumilaceae_P1_phases_spline_summary <- summary(gam_model_nitrosopumilaceae_P1_phases_spline)
gam_model_nitrosopumilaceae_P1_phases_spline_summary$dev.expl #Model explains 72.9% of deviance
gam_model_nitrosopumilaceae_P1_phases_spline_summary$r.sq #Model explains 67.3 of the variance 

#Smooth terms
gam_model_nitrosopumilaceae_P1_phases_spline_summary_table <- data.frame(gam_model_nitrosopumilaceae_P1_phases_spline_summary$s.table)
gam_model_nitrosopumilaceae_P1_phases_spline_summary_table

#Parametric terms
gam_model_nitrosopumilaceae_P1_phases_parametric_summary_table <- data.frame(gam_model_nitrosopumilaceae_P1_phases_spline_summary$p.table,
                                                                             check.names = F)
gam_model_nitrosopumilaceae_P1_phases_parametric_summary_table


#Number of samples per phase
n_by_phase_P1 <- phyloseq_with_unclassifieds_otu.samples.family.ra.nitrosopumilaceae.melt_arranged_date.P1 %>%
  count(Date_num_phase_abbrv, name = "No. Samples")


#Stable 8 (established system) with parametric terms 
stable8.2 <- 
  gam_model_nitrosopumilaceae_P1_phases_parametric_summary_table %>%
  mutate(
    System = "Established",
    Phase = case_when(
      str_detect(rownames(.), "modelL")  ~ "L",
      str_detect(rownames(.), "modelT1") ~ "T1",
      str_detect(rownames(.), "modelE") ~ "E",
      TRUE ~ NA_character_
    ),
    Contrast = case_when(
      str_detect(rownames(.), "modelL")  ~ "L vs P",
      str_detect(rownames(.), "modelT1") ~ "T1 vs P",
      str_detect(rownames(.), "modelE") ~ "E vs P",
      str_detect(rownames(.), "Intercept") ~ "Intercept",
      TRUE ~ NA_character_
    ),
    Effect = "Phase",
    SE = round(`Std. Error`, 3),
    `t value` = round(`t value`, 3),
    `p value` = format(`Pr(>|t|)`, scientific = TRUE, digits = 3),
  )%>%
  filter(!Contrast == "Intercept")%>%
  left_join(n_by_phase_P1, by = c("Phase" = "Date_num_phase_abbrv"))%>%
  #Contrast samples
  mutate(`No. Samples` = paste(`No. Samples`, "-", (n_by_phase_P1%>%filter(Date_num_phase_abbrv == "P")%>%pull(`No. Samples`))))%>%
  select(System, Contrast, Effect, `No. Samples`, Estimate, SE, `t value`, `p value`)
stable8.2

#Stable 9.2 (established system) with smooth terms 
stable9.2 <- 
  gam_model_nitrosopumilaceae_P1_phases_spline_summary_table %>%
  mutate(
    System = "Established",
    Phase = case_when(
      str_detect(rownames(.), "modelL")  ~ "L",
      str_detect(rownames(.), "modelT1") ~ "T1",
      str_detect(rownames(.), "modelE") ~ "E",
      str_detect(rownames(.), "modelP") ~ "P",
      TRUE ~ NA_character_
    ),
    Effect = case_when(
      str_detect(rownames(.), "pH")  ~ "pH",
      str_detect(rownames(.), "Salinity") ~ "Salinity",
      str_detect(rownames(.), "Temperature") ~ "Temperature",
      str_detect(rownames(.), "model")  ~ "Copper"
    ),
    EDF = round(edf, 3),
    `Reference EDF` = round(Ref.df, 3),
    `F Statistic` = round(`F`, 3),
    `p value` = format(p.value, scientific = TRUE, digits = 3)
  )%>%  
  left_join(n_by_phase_P1, by = c("Phase" = "Date_num_phase_abbrv"))%>%
  select(System, Phase, Effect, `No. Samples`, EDF, `F Statistic`, `p value`)
stable9.2

#Join into one stable 8 (for parametric terms)
stable8 <- bind_rows(stable8.1, 
                     stable8.2)
stable8
write_xlsx(stable8, "SupplementaryTable8.xlsx")

#Join into one stable 8 (for parametric terms)
stable9 <- bind_rows(stable9.1, 
                     stable9.2)
stable9
write_xlsx(stable9, "SupplementaryTable9.xlsx")


#RELATIVE ABUNDANCES OF JUST THE NITRIFYING COMMUNITIES#####
any(sample_sums(phyloseq.samples.nitrifiers)== 0) ## no samples with 0 OTUs

###FAMILY#######
phyloseq.samples.nitrifiers.family <- tax_glom(phyloseq.samples.nitrifiers, taxrank = "Family", NArm = F)
phyloseq.samples.nitrifiers.family #18 families and 204 samples

#Relative abundance
phyloseq.samples.nitrifiers.family.ra <- transform_sample_counts(phyloseq.samples.nitrifiers.family, 
                                                          function(x) x/sum(x)*100) ##Relative abundance from normalized data

phyloseq.samples.nitrifiers.family.ra #18 families and 204 samples (RA)

##Melt ps object
phyloseq.samples.nitrifiers.family.ra.melt <- psmelt(phyloseq.samples.nitrifiers.family.ra)

####SUPPLEMENTARY TABLE 7 #####
stable7 <- phyloseq.samples.nitrifiers.family.ra.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established")) %>%
  group_by(System, Family, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    min_abun  = min(Abundance, na.rm = TRUE),
    max_abun  = max(Abundance, na.rm = TRUE),
    .groups = "drop_last"
  ) %>%
  arrange(System, Date_num_phase_abbrv, desc(mean_abun)) %>%
  mutate(
    `Mean Relative Abundance (%) ± SD` = paste0(
      format(mean_abun, scientific = TRUE, digits = 3),
      " ± ",
      format(sd_abun, scientific = TRUE, digits = 3)
    ),
    `Min - Max Relative Abundance (%)` = paste0(
      format(min_abun, scientific = TRUE, digits = 3),
      " - ",
      format(max_abun, scientific = TRUE, digits = 3)
    )
  ) %>%
  rename(Phase = Date_num_phase_abbrv) %>%
  select(
    Family,
    System,
    Phase,
    `Mean Relative Abundance (%) ± SD`,
    `Min - Max Relative Abundance (%)`
  )
stable7
write_xlsx(stable7, 
           "SupplementaryTable7.xlsx")

##Which are the top most abundant taxa by family group? 
phyloseq.samples.nitrifiers.family.ra.melt_most_abundant <- phyloseq.samples.nitrifiers.family.ra.melt %>%
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
phyloseq.samples.nitrifiers.family.ra.melt_most_abundant

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
unique(phyloseq.samples.nitrifiers.family.ra.melt$Family)
# "Nitrosopumilaceae"          "Nitrospiraceae"             "Nitrosomonadaceae"          "UBA8639"                    "Nitrospinaceae"            
# "VA-1"                       "JAJDAT01"                   "Nitrococcaceae"             "Nitrososphaeraceae"         "unclassified Nitrospinales"
# "Xanthobacteraceae"          "Nitrosococcaceae"           "Gallionellaceae"            "Thermomicrobiaceae"         "JACRGP01"                  
# "NS-4"                       "UBA213"                     "Nitrosocaldaceae"       
phyloseq.samples.nitrifiers.family.ra.melt <- phyloseq.samples.nitrifiers.family.ra.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Family == "Xanthobacteraceae" ~ "NOB",
    Family == "Gallionellaceae" ~ "NOB",
    Family == "Nitrococcaceae" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Family == "unclassified Nitrospinales" ~ "NOB",
    Family == "Thermomicrobiaceae" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB"))) %>%
  arrange(Nitrifying_group, Phylum, Class, Order, Family) %>%
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
palette_nitrifiers_only_family_df <- phyloseq.samples.nitrifiers.family.ra.melt %>% 
  distinct(Family, Nitrifying_group) %>%
  group_by(Nitrifying_group) %>%
  arrange(Family) %>%   
  mutate(
    base_color = nitrifier_base_colors[Nitrifying_group],
    # shade = seq(-0.1, 0.1, length.out = n()),
    shade = seq(0.01, 0.8, length.out = n()),
    color = darken(base_color, amount = shade))%>%
  ungroup()
palette_nitrifiers_only_family_df

#Set up final palette
palette_nitrifiers_only_family <- setNames(
  palette_nitrifiers_only_family_df$color,
  palette_nitrifiers_only_family_df$Family)
palette_nitrifiers_only_family

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_families_only <- top_taxa_legend(phyloseq.samples.nitrifiers.family.ra.melt,
                                                taxlevel = "Family",
                                                              n = 5)
top_nitrifying_families_only
#Keep the factor levels for legend
nitifier_families_levels_only <- levels(phyloseq.samples.nitrifiers.family.ra.melt$Family)
top_nitrifying_families_only <- 
  nitifier_families_levels_only[nitifier_families_levels_only %in% top_nitrifying_families_only]
top_nitrifying_families_only

#Plot
RA_nitrifiers_only_family.plot <- ggplot(phyloseq.samples.nitrifiers.family.ra.melt,
                                                    aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
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
  scale_fill_manual(values = palette_nitrifiers_only_family, 
                    breaks = top_nitrifying_families_only) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.07, 0.6),  # x, y inside plot
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
    panel.grid.major.x = element_line(linetype = "dashed", color = "black"),
    axis.line.y = element_line(linewidth = 0.7, colour = "black"),
    axis.text.x = element_text(colour = "black", size = 20,
                               vjust = 0.5, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 22),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
RA_nitrifiers_only_family.plot

###SPECIES#######
phyloseq.samples.nitrifiers.species <- tax_glom(phyloseq.samples.nitrifiers, taxrank = "Species", NArm = F)
phyloseq.samples.nitrifiers.species #1690 species and 204 samples

#Relative abundance
phyloseq.samples.nitrifiers.species.ra <- transform_sample_counts(phyloseq.samples.nitrifiers.species, 
                                                                 function(x) x/sum(x)*100) ##Relative abundance from normalized data

phyloseq.samples.nitrifiers.species.ra #1690 species and 204 samples (RA)

#Merge low abun species
phyloseq.samples.nitrifiers.ra.species.filt <- merge_low_abundance_grouped_ra(phyloseq.samples.nitrifiers.species.ra, 
                                                                              "Enclosure", 
                                                                              level = "Species", threshold = 0.6)
phyloseq.samples.nitrifiers.ra.species.filt #26 Species over 0.5% mean RA

#Melt ps object
phyloseq.samples.nitrifiers.ra.species.filt.melt <- psmelt(phyloseq.samples.nitrifiers.ra.species.filt)

##Which are the top most abundant taxa by group? 
top_nitrifier_only_species_list <- phyloseq.samples.nitrifiers.ra.species.filt.melt %>%
  group_by(Enclosure, Species) %>%
  summarise(
    mean_abun = mean(Abundance, na.rm = TRUE),
    sd_abun   = sd(Abundance, na.rm = TRUE),
    .groups = "drop_last") %>%
  arrange(Enclosure,  desc(mean_abun))
top_nitrifier_only_species_list

##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq.samples.nitrifiers.ra.species.filt.melt <- phyloseq.samples.nitrifiers.ra.species.filt.melt %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB")))%>%
  arrange(Nitrifying_group, Phylum, Class, Order, Family, Species) %>%
  mutate(Species = factor(Species, levels = unique(Species)))

####SUPPLEMENTARY TABLE 9#####
##Which are the top most abundant taxa by group? 
stable9 <-  phyloseq.samples.nitrifiers.ra.species.filt.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Species, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 3),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 3), 
    max_abun = round(max(Abundance, na.rm = TRUE), 3),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  mutate(
    `Mean Relative Abundance Within the Nitrifying Community (%) ± SD` = paste0(
      mean_abun,
      " ± ",
      sd_abun
    ),
    `Min - Max Relative Abundance Within the Nitrifying Community (%)` = paste0(
      min_abun,
      " - ",
      max_abun
    ))%>%
  rename(Phase = Date_num_phase_abbrv)%>%
  select(System, Species, Phase, 
         `Mean Relative Abundance Within the Nitrifying Community (%) ± SD`,
         `Min - Max Relative Abundance Within the Nitrifying Community (%)`
  )
stable9
write_xlsx(stable9, 
           "SupplementaryTable9.xlsx")

#If too many taxa for legend, apply the function to obtain top ones
top_nitrifying_species_only <- top_taxa_legend(phyloseq.samples.nitrifiers.ra.species.filt.melt,
                                                taxlevel = "Species",
                                                n = 10)
top_nitrifying_species_only
#Keep the factor levels for legend
nitifier_species_levels_only <- levels(phyloseq.samples.nitrifiers.ra.species.filt.melt$Species)
top_nitrifying_species_only <- 
  nitifier_species_levels_only[nitifier_species_levels_only %in% top_nitrifying_species_only]
top_nitrifying_species_only

#Color palette
#Create base colors based on ammonia-nitrate oxidizing groups
nitrifier_base_colors <- c(
  AOA = "#D81B60",  # bright pink/red
  AOB = "#1E88E5",  # strong blue
  NOB = "#FFC107"   # vivid amber/yellow
)

#Make hues based on species within each ammonia-nitrite oxidizing group
# palette_nitrifiers_only_species_df <- phyloseq.samples.nitrifiers.ra.species.filt.melt %>% 
#   distinct(Species, Nitrifying_group) %>%
#   group_by(Nitrifying_group) %>%
#   arrange(Species) %>%   
#   mutate(
#     base_color = nitrifier_base_colors[Nitrifying_group],
#     # shade = seq(-0.1, 0.1, length.out = n()),
#     shade = seq(0.01, 0.8, length.out = n()),
#     color = darken(base_color, amount = shade))%>%
#   ungroup()
# palette_nitrifiers_only_species_df

#Make hues based on species within each ammonia-nitrite oxidizing group
palette_nitrifiers_only_species_df <-
  phyloseq.samples.nitrifiers.ra.species.filt.melt %>%
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
          "#99003A",
          "#D5A0C3",
          "#A64D79",
          "#E8C4CE",
          "#B36A83",
          "pink",
          "#65001F",
          "#7D3150"
          # "#FCE4EC",  # very light pink
          # "#880E4F",  # wine
          # "#4A148C",   # deep purple
          # #"#D50000",  # red (clear shift toward red)
          # #"#E53935",
          # "#880E4F",  # deep wine/magenta (dark endpoint)
          # "#FF4081", # hot pink
          # "#B71C1C",
          # "#EF9A9A"
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
          "lightgoldenrod",
          "#FFD54F", 
          "#FFB300", 
          "goldenrod4",   
          "darkgoldenrod" 
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
palette_nitrifiers_only_species$'Others <0.6% RA' <- "grey95"
palette_nitrifiers_only_species$'unclassified Nitrosopumilaceae' <- "#420017"


#Plot
RA_nitrifiers_only_species.plot <- ggplot(phyloseq.samples.nitrifiers.ra.species.filt.melt, 
                                                     aes(x=factor(Date_num), y= Abundance, fill = Species)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_col(color = "black")+
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  ggh4x::facetted_pos_scales(
    x = list(
      Enclosure == "H21" ~
        scale_x_discrete(
          breaks = c("1", "27", "38", "51", "81", "108", "135"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        ),
      
      Enclosure == "P1" ~
        scale_x_discrete(
          breaks = c("1","53","65","104","141"),
          expand = expansion(mult = c(0.03, 0.03)),
          drop = TRUE
        )))+
  scale_fill_manual(
    values = palette_nitrifiers_only_species,
    breaks = top_nitrifying_species_only,
    labels = function(x) str_wrap(x, width = 20),
    drop = FALSE
  )+
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.position = c(1.085, 0.5),  # x, y inside plot
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
RA_nitrifiers_only_species.plot


#SPECIES DIFFERENTIAL ABUNDANCE BETWEEN PHASES ##### 
#First, doing imputations to fill Copper levels where they were not measured 
phyloseq.samples_metadata <- 
  data.frame(sample_data(phyloseq.samples))%>%
  arrange(Enclosure, Date_num)

#How many samples are missing copper levels: 
phyloseq.samples_metadata %>%
  filter(is.na(Copper_level_mg_L)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure         Date_num_phase                      n
# H21 (naive)       Low Copper Levels (Day 1-27)        9
# H21 (naive)       Transition Period 2 (Day 52-81)     3
# P1 (established)    Low Copper Levels (Day 1-53)        1
# P1 (established)    Tx Copper Exposure (Day 66-104)     1

#pH?
phyloseq.samples_metadata %>%
  filter(is.na(pH_spu)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure Date_num_phase                           n
# H21       First Tx Copper Exposure (Day 39-51)     1
# P1        Low Copper Levels (Day 1-53)             1

#Salinity?
phyloseq.samples_metadata %>%
  filter(is.na(Salinity_ppt)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure Date_num_phase                              n
# H21       First Tx Copper Exposure (Day 39-51)        1
# H21       Transition Period 3 (Day 86-108)            1
# H21       Second Tx Copper Exposure (Day 109-135)     1
# P1        Low Copper Levels (Day 1-53)                1

#Temperature?
phyloseq.samples_metadata %>%
  filter(is.na(Temperature_F)) %>%
  group_by(Enclosure, Date_num_phase) %>%
  summarise(n = n(), .groups = "drop")
# Enclosure Date_num_phase                           n
# H21       First Tx Copper Exposure (Day 39-51)     1
# P1        Low Copper Levels (Day 1-53)             3
# P1        Tx Copper Exposure (Day 66-104)          2
# P1        Post-Treatment (Day 105-141)             2

#Add a column that will have the measured copper values as well as the imputed ones
phyloseq.samples_metadata$Copper_level_mg_L_imp <- 
  phyloseq.samples_metadata$Copper_level_mg_L

#Using na.approx() to fill copper missing values, fills them by drawing a straight line between known points and estimating values along that line.
phyloseq.samples_metadata <- 
  phyloseq.samples_metadata %>%
  group_by(Enclosure) %>%
  mutate(
    Copper_level_mg_L_imp = zoo::na.approx(
      Copper_level_mg_L_imp,
      x = Date_num,
      na.rm = FALSE), 
    Copper_imputation = case_when(
      is.na(Copper_level_mg_L) & !is.na(Copper_level_mg_L_imp) ~ "yes",
      TRUE ~ "no"
    ) #Fill in copper with linear interpolation. If copper level imputed, 'yes' on Copper_imputation
  ) %>%
  ungroup()

#Using na.approx() to fill Temperature, salinity, and pH
#Add a column for pH
phyloseq.samples_metadata$pH_spu_imp <- 
  phyloseq.samples_metadata$pH_spu
#Add a column for Temperature
phyloseq.samples_metadata$Temperature_F_imp <- 
  phyloseq.samples_metadata$Temperature_F
#Add a column for Salinity
phyloseq.samples_metadata$Salinity_ppt_imp <- 
  phyloseq.samples_metadata$Salinity_ppt

#Now impute those (Temperature, salinity, and pH)
phyloseq.samples_metadata <- 
  phyloseq.samples_metadata %>%
  group_by(Enclosure) %>%
  mutate(
    pH_spu_imp = zoo::na.approx(
      pH_spu_imp,
      x = Date_num,
      na.rm = FALSE), 
    Temperature_F_imp = zoo::na.approx(
      Temperature_F_imp,
      x = Date_num,
      na.rm = FALSE), 
    Salinity_ppt_imp = zoo::na.approx(
      Salinity_ppt_imp,
      x = Date_num,
      na.rm = FALSE)
  ) %>%
  ungroup()


##NAIVE SYSTEM####
#For metadata, will use the imputed values (ph, Copper, Salinity, Temperature)
ancombc2_naive_species_metadata <- phyloseq.samples_metadata%>%
  filter(Enclosure == "H21")%>%
  column_to_rownames(var = "SampleID")

#Only naive system nitrifying taxa raw counts
phyloseq.samples.nitrifiers.H21 <- subset_samples(phyloseq.samples.nitrifiers, Enclosure == "H21")
phyloseq.samples.nitrifiers.H21 <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers.H21) > 0, 
                                              phyloseq.samples.nitrifiers.H21)
phyloseq.samples.nitrifiers.H21 #1686 taxa and 87 samples
ancombc2_naive_nitrifiers_otu_table <- data.frame(otu_table(phyloseq.samples.nitrifiers.H21))
nrow(ancombc2_naive_nitrifiers_otu_table) #1686 taxa


###ANCOMBC#######
##Untransformed (raw) counts in naive samples of nitrifying taxa 
nrow(ancombc2_naive_nitrifiers_otu_table) #1686 taxa

##Reorder the Date_num_phase_abbrv as factor, in naive system E2 as reference
ancombc2_naive_species_metadata$Date_num_phase_abbrv <- 
  factor(ancombc2_naive_species_metadata$Date_num_phase_abbrv, levels = c("E2",
                                                                      "L", 
                                                                      "T1", 
                                                                      "E1",
                                                                      "T2", 
                                                                      "T3"))
##running ancombc2
ancombc_output_H21_nitrifiers <-ancombc2(
  data= ancombc2_naive_nitrifiers_otu_table, 
  meta_data = ancombc2_naive_species_metadata,
  assay_name = "counts", 
  tax_level = "Species",
  fix_formula = "Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp",
  prv_cut = 0.05, 
  lib_cut = 0, 
  p_adj_method = "BH",
  dunnet = TRUE, 
  group= "Date_num_phase_abbrv", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, #default significance
  verbose = TRUE)

## extract results from comparisons 
res_H21_nitrifiers <- ancombc_output_H21_nitrifiers$res_dunn

#Pivot longer the results (lfc, se, p val, q val, passed sensitivity, diff)
ancom_H21_phase_nitrifiers <- res_H21_nitrifiers %>%
  mutate(
    across(
      starts_with("lfc"),
      ~ round(.x, 2)
    )
  ) %>%
  pivot_longer(
    cols = matches("^(lfc|se|passed_ss|p|q|diff_robust)_"),
    names_to = c(".value", "group"),
    names_pattern = "(lfc|passed_ss|se|p|q|diff_robust)_(.*)"
  ) %>%
  filter(group != "(Intercept)")%>%
  arrange(taxon)
ancom_H21_phase_nitrifiers

##rework our group names so they're shorter and more manageable 
ancom_H21_phase_nitrifiers <- ancom_H21_phase_nitrifiers %>%
  mutate(group= case_when(
    group == "Date_num_phase_abbrvL" ~ "L vs E2",
    group == "Date_num_phase_abbrvT1" ~ "T1 vs E2",
    group == "Date_num_phase_abbrvE1" ~ "E1 vs E2",
    group == "Date_num_phase_abbrvT2" ~ "T2 vs E2",
    group == "Date_num_phase_abbrvT3" ~ "T3 vs E2",
    TRUE ~ group ##keeps original name for groups not specified 
  ))
ancom_H21_phase_nitrifiers$group ##Now the group names are shorter and more manageable

##Obtaining the confidence intervals for the log fold change
ancom_H21_phase_nitrifiers_2 <- ancom_H21_phase_nitrifiers %>%
  mutate(lower.ci = lfc - 1.96*se, 
         upper.ci = lfc + 1.96*se)

#Taxonomy of species in the naive system
taxonomy.H21 <- data.frame(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples_H21))%>%
  rownames_to_column(var = "OTU")

##Final fixup to make compatible with plotting. Only keeping nitrifiers. 
ancom_H21_phase_nitrifiers_nitirfiers <- ancom_H21_phase_nitrifiers_2 %>%
  filter (passed_ss == 1)%>%##Only want those that passed sensitivity testing
  # #Will only include the classified taxa
  rename(stderr = se, ##Renaming
         pval = p,  
         qval = q)%>%
  mutate(direction = ifelse(lfc > 0, "elevated", "depleted"),
         # Exponentiate the values first
         lfc_exp = exp(lfc),
         lower.ci_exp = exp(lower.ci),
         upper.ci_exp = exp(upper.ci)) %>%
  mutate(
    ##Getting log 2 so I can compare with maaslin
    coef =  log2(lfc_exp),
    lower.ci =  log2(lower.ci_exp),
    upper.ci =  log2(upper.ci_exp),
    plot = "Log2 Fold change with 95%CI", 
    test = "ANCOM-BC",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  filter(DA == "q ≤ 0.05")%>%
  left_join(taxonomy.H21, by = c("taxon" = 'OTU'))%>%
  select(Species, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
length(unique(ancom_H21_phase_nitrifiers_nitirfiers$Species)) #225 taxa

#How many DA species per comparison? 
ancom_H21_phase_nitrifiers_nitirfiers %>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_nitrifiers = n())
# group             n_DA_species
# E1 vs E2                    85
# L vs E2                     39
# T1 vs E2                    30
# T2 vs E2                    36
# T3 vs E2                   174


###MaAsLin3#######
#OTU table (RA) with unclassified reads 
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers #Nitrifier RA within the overall community (including unclassified reads)
#Select only naive system
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21 <- 
  subset_samples(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers, 
                 Enclosure == "H21")
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21 <- 
  prune_taxa(taxa_sums(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21) > 0, 
             phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21)
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21 #  taxa and 87 samples

#For otu table:
data_maaslin_H21 <- data.frame(t(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21@otu_table))

##Sample metadata. Will use imputed copper, salinity, ph, and temperature (and factor Date_num_phase_abbrv with E2 as reference)
metadata_maaslin_H21 <- phyloseq.samples_metadata%>%
  filter(Enclosure == "H21")%>%
  column_to_rownames(var = "SampleID")%>%
  select(c("Date_num_phase_abbrv", "pH_spu_imp", "Temperature_F_imp", "Salinity_ppt_imp"))%>%
  mutate(Date_num_phase_abbrv = 
           factor(Date_num_phase_abbrv, 
                  levels = c("E2",
                             "L", 
                             "T1", 
                             "E1",
                             "T2", 
                             "T3")
                  )
         )#making Date_num_phase_abbrv a factor since I'll be adding it as a fixed effect on the MaAslin model

##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_H21_nitrifiers <-  maaslin3(
  input_data = data_maaslin_H21, 
  input_metadata = metadata_maaslin_H21, 
  correction = 'BH',
  output = "MaAsLin3_H21",  
  fixed_effect = c("Date_num_phase_abbrv", "pH_spu_imp", "Temperature_F_imp", "Salinity_ppt_imp"), 
  min_prevalence= 0.05,
  median_comparison_abundance = T, #default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, #None
  normalization="NONE", #data is already RA
  transform="LOG", #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
length(unique(maaslin_H21_nitrifiers$fit_data_abundance$results$feature)) #Abundance results from MaAslin - 1678

#Calculate confidence intervals
maaslin_H21_nitrifiers_results <- maaslin_H21_nitrifiers$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr)%>%
  filter(!grepl("imp", name)) #not including temperature, salinity, pH results, only Date_num_phase
maaslin_H21_nitrifiers_results
length(unique(maaslin_H21_nitrifiers_results$feature)) #1678 taxa

##Taxonomy of H21 OTUs 
input_taxonomy.H21 <- data.frame(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_H21@tax_table) %>%
  rownames_to_column(var = "feature")

##Final edits to put together for plot, keep only nitrifiers
maaslin_H21_nitrifiers_results_2 <- maaslin_H21_nitrifiers_results%>%
  mutate(group = case_when(value == "L" ~ "L vs E2", 
                            value == "T1" ~ "T1 vs E2", 
                            value == "E1" ~ "E1 vs E2", 
                            value == "T2" ~ "T2 vs E2", 
                            value == "T3" ~ "T3 vs E2"))%>%
  rename(pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  filter(DA == "q ≤ 0.05")%>%
  left_join(input_taxonomy.H21, by = "feature")%>%
  select(Species, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) 
maaslin_H21_nitrifiers_results_2

##ANCOM and MaAslin together
#Get the intersect per group 
DA_intersect_nitrifiers_H21 <-
  inner_join(ancom_H21_phase_nitrifiers_nitirfiers,
             maaslin_H21_nitrifiers_results_2,
             by = c("group", "Species"))

#Species that are differentially abundant per group AND agreed on by both MaAsLin and ANCOM”
DA_H21_nitrifiers_plot_MaAslinANCOM.data <-
  bind_rows(
    ancom_H21_phase_nitrifiers_nitirfiers,
    maaslin_H21_nitrifiers_results_2
  ) %>%
  semi_join(DA_intersect_nitrifiers_H21, by = c("group", "Species")) %>%
  mutate(group = factor(group, levels = c(
    "L vs E2",
    "T1 vs E2",
    "E1 vs E2",
    "T2 vs E2",
    "T3 vs E2"
  )))

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_H21_nitrifiers <- maaslin_H21_nitrifiers$transformed_data %>% #transformed data is TSS transformed otu counts from MaAslin output
  rownames_to_column(var = "SampleID")%>%
  pivot_longer(cols = -SampleID, names_to = "OTU", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy.H21, by = c("OTU" = "feature"))%>%
  filter(Species %in% intersect(maaslin_H21_nitrifiers_results_2$Species,
                                ancom_H21_phase_nitrifiers_nitirfiers$Species))%>%
  left_join(metadata_maaslin_H21%>%rownames_to_column(var = "SampleID"), by = "SampleID")%>%
  select(Species, SampleID, logvalue, Date_num_phase_abbrv, plot)

# Non-E2 samples
RA_MaaslinAncom_H21_nitrifiers_nonE2 <-
  RA_MaaslinAncom_H21_nitrifiers %>%
  mutate(group = case_when(
    Date_num_phase_abbrv == "L"  ~ "L vs E2",
    Date_num_phase_abbrv == "T1" ~ "T1 vs E2",
    Date_num_phase_abbrv == "E1" ~ "E1 vs E2",
    Date_num_phase_abbrv == "T2" ~ "T2 vs E2",
    Date_num_phase_abbrv == "T3" ~ "T3 vs E2"
    )) %>%
  filter(!is.na(group))

# Duplicate E2 samples into every comparison
RA_MaaslinAncom_H21_nitrifiers_E2 <-
  RA_MaaslinAncom_H21_nitrifiers %>%
  filter(Date_num_phase_abbrv == "E2") %>%
  tidyr::crossing(group = c(
    "L vs E2", 
    "T1 vs E2",
    "E1 vs E2",
    "T2 vs E2",
    "T3 vs E2"
  ))

# Join data
RA_MaaslinAncom_H21_nitrifiers <-
  bind_rows(RA_MaaslinAncom_H21_nitrifiers_nonE2, RA_MaaslinAncom_H21_nitrifiers_E2)

#Bias-corrected abundances (ANCOMBC2)
nitrifiers_H21_nitrifiers_log_corr_abn <- ancombc_output_H21_nitrifiers$bias_correct_log_table%>%
  data.frame()%>% ##make into data frame
  rownames_to_column("OTU")%>%
  left_join(input_taxonomy.H21, by = c("OTU" = "feature"))%>%
  filter(Species %in% intersect(maaslin_H21_nitrifiers_results_2$Species,
                                ancom_H21_phase_nitrifiers_nitirfiers$Species))%>%
  select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, -OTU) %>%
  pivot_longer(cols = -Species, names_to = "SampleID", values_to = "logvalue")%>% ##Pivot to longer data format
  left_join(metadata_maaslin_H21%>%rownames_to_column(var="SampleID"), by = "SampleID")%>%
  select(Species, SampleID, logvalue, Date_num_phase_abbrv)%>%
  mutate(plot = 'log(Bias-corrected abundances)')
nitrifiers_H21_nitrifiers_log_corr_abn

# Non-E2 samples
nitrifiers_H21_species_log_corr_abn_nonE2 <-
  nitrifiers_H21_species_log_corr_abn %>%
  mutate(group = case_when(
    Date_num_phase_abbrv == "L"  ~ "L vs E2",
    Date_num_phase_abbrv == "T1" ~ "T1 vs E2",
    Date_num_phase_abbrv == "E1" ~ "E1 vs E2",
    Date_num_phase_abbrv == "T2" ~ "T2 vs E2",
    Date_num_phase_abbrv == "T3" ~ "T3 vs E2",
  )) %>%
  filter(!is.na(group))

# Duplicate E2 samples into every comparison
nitrifiers_H21_species_log_corr_abn_E2 <-
  nitrifiers_H21_species_log_corr_abn %>%
  filter(Date_num_phase_abbrv == "E2") %>%
  tidyr::crossing(group = c(
    "L vs E2", 
    "T1 vs E2",
    "E1 vs E2",
    "T2 vs E2",
    "T3 vs E2"
  ))

# Binding datasets
nitrifiers_H21_species_log_corr_abn <-
  bind_rows(nitrifiers_H21_species_log_corr_abn_nonE2, 
            nitrifiers_H21_species_log_corr_abn_E2)

#Put together objects to plot DA
#First, only keeping RA for those differentially abundant taxa 
#Maaslin
RA_MaaslinAncom_H21_nitrifiers_differential  <-
  RA_MaaslinAncom_H21_nitrifiers %>%
  semi_join(
    DA_H21_nitrifiers_plot_MaAslinANCOM.data %>% select(group, Species),
    by = c("group", "Species")
  )
#Ancombc
nitrifiers_H21_species_log_corr_abn_differential <- 
  nitrifiers_H21_species_log_corr_abn %>%
  semi_join(
    DA_H21_nitrifiers_plot_MaAslinANCOM.data %>% select(group, Species),
    by = c("group", "Species")
  )

#Put RA, bias corrected abundances, DA coefficients dataframes together
DA_H21_phase_plot_together <- bind_rows(DA_H21_nitrifiers_plot_MaAslinANCOM.data, 
                                                   nitrifiers_H21_species_log_corr_abn_differential, 
                                                   RA_MaaslinAncom_H21_nitrifiers_differential) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)
DA_H21_phase_plot_together

#How many DA nitrifying species per comparison? (without RA filtering)
DA_H21_phase_plot_together %>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# E1 vs E2           78
# L vs E2            39
# T1 vs E2           30
# T2 vs E2           35
# T3 vs E2          135

#DA species - all nitrifying taxa (without RA filtering)
DA_H21_all_nitrifying_taxa_phase_comparisons <- DA_H21_phase_plot_together %>%
  filter(!is.na(coef)) %>%
  group_by(Species, group) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(n_groups = n_distinct(group)) %>%   # how many comparisons each species appears in
  ungroup() %>%
  arrange(desc(n_groups), Species, group)
DA_H21_all_nitrifying_taxa_phase_comparisons
unique((DA_H21_all_nitrifying_taxa_phase_comparisons %>%filter(group == "T2 vs E2"))$Species) #OK!

#How many belonging to Nitrosopumilus (without RA filtering)
DA_H21_phase_plot_together %>%
  filter(grepl("Nitrosopumilus", Species))%>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# L vs E2            27
# T1 vs E2           2
# E1 vs E2           2
# T2 vs E2           25
# T3 vs E2           75

#Which Nitrosopumilus species are they? (without RA filtering)
DA_H21_Nitrosopumilus_phase_comparisons <- DA_H21_phase_plot_together %>%
  filter(!is.na(coef)) %>%
  filter(grepl("Nitrosopumilus", Species)) %>%
  group_by(Species, group) %>%
  #slice(1) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(n_groups = n_distinct(group)) %>%   # how many comparisons each species appears in
  ungroup() %>%
  arrange(desc(n_groups), Species, group)
DA_H21_Nitrosopumilus_phase_comparisons

#Which Nitrosopumilus species are DA against E2 across all phases? (without RA filtering)
DA_H21_all_phases_comparisons_E2 <- DA_H21_phase_plot_together %>%
  filter(!is.na(coef)) %>%
  filter(grepl("Nitrosopumilus", Species)) %>%
  group_by(Species) %>%
  mutate(n_species = n()) %>%   # total rows per species
  ungroup() %>%
  filter(n_species == 10) %>% #10 rows means they were DA across all comparisons (5 groups, 2 tests)
  arrange(desc(n_species), Species)
unique(DA_H21_all_phases_comparisons_E2$Species)#none

####Filtering out those species that make up less than 0.1% RA of the nitrifying community########
any(sample_sums(phyloseq.samples.nitrifiers_H21)== 0) ## no samples with 0 OTUs
phyloseq.samples.nitrifiers_H21.ra <- transform_sample_counts(phyloseq.samples.nitrifiers_H21, 
                                                              function(x) x/sum(x)*100)
phyloseq.samples.nitrifiers_H21.ra #1686 taxa and 87 samples
sample_sums(phyloseq.samples.nitrifiers_H21.ra) #RA to 100

phyloseq.samples.nitrifiers_H21.ra.filt.DA <- filter_taxa(phyloseq.samples.nitrifiers_H21.ra, function(x) mean(x) > 0.1, TRUE) 
phyloseq.samples.nitrifiers_H21.ra.filt.DA ##52 OTUs with mean RA > 0.1% across samples (naive system)

#Filter them from the DA plot object
DA_H21_phase_plot_together_01_RA <- DA_H21_phase_plot_together %>%
  filter(Species %in% phyloseq::tax_table(phyloseq.samples.nitrifiers_H21.ra.filt.DA)[,"Species"])

#How many DA nitrifiers per comparison? 
DA_H21_phase_plot_together_01_RA %>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# L vs E2            18
# T1 vs E2           16
# E1 vs E2           11
# T2 vs E2           8
# T3 vs E2           15

#How many belonging to Nitrosopumilus
DA_H21_phase_plot_together_01_RA %>%
  filter(grepl("Nitrosopumilus", Species))%>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# L vs E2            8
# T1 vs E2            2
# E1 vs E2           2
# T2 vs E2           8
# T3 vs E2           8

#Which Nitrosopumilus species are they?
DA_H21_Nitrosopumilus_phase_comparisons <- DA_H21_phase_plot_together_01_RA %>%
  filter(!is.na(coef)) %>%
  filter(grepl("Nitrosopumilus", Species)) %>%
  group_by(Species, group) %>%
  #slice(1) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(n_groups = n_distinct(group)) %>%   # how many comparisons each species appears in
  ungroup() %>%
  arrange(desc(n_groups), Species, group)
DA_H21_Nitrosopumilus_phase_comparisons


####SUPPPLEMENTARY TABLE 11.1 ######
stable11.1 <- DA_H21_all_nitrifying_taxa_phase_comparisons%>%
  mutate(`Effect Size (95% CI)` = 
           paste0(round(coef, 2), " (", 
                  round(lower.ci,2),
                  ", ",
                  round(upper.ci,2),
                  ")"), 
         `P value` = format(pval, scientific = TRUE, digits = 3), 
         `Q value` = format(qval, scientific = TRUE, digits = 3), 
         Comparison = group,
         Test = test,
         System = "Naive"
         )%>%
  select(System, Species, Comparison,  Test, `Effect Size (95% CI)`, `P value`, `Q value`)
stable11.1

###PLOTTING DA - NITRIFIERS#####
####Only E1 VS E2 ########
##Ordering how I want the "Family" tax level to show up on the plot 
input_taxonomy.H21 ##dataframe object for Taxonomy of H21 taxa

#Add nitrifier group
input_taxonomy.H21 <- input_taxonomy.H21 %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB")))

#Taxonomy plot data
taxonomy_plot_data_H21_E1E2 <- DA_H21_phase_plot_together_01_RA %>%
  left_join(input_taxonomy.H21, by= "Species")%>%
  filter(group == "E1 vs E2")%>%
  distinct(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  arrange(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  mutate(Species = factor(Species, levels = rev(Species)))%>% ##Since I arranged by family, this is the order I want the genera to show up
  dplyr::group_by(Family) %>%
  dplyr::mutate(label_Family = ifelse(row_number() == 1, Family, "")) %>%  # Create 'label_Family' with only the first occurrence of each Family
  ungroup()
taxonomy_plot_data_H21_E1E2

##Factor "Species" level by the order I want (taxonomy_plot_data_H21_E1E2$Species)
DA_H21_phase_plot_together_01_RA_E1E2 <- DA_H21_phase_plot_together_01_RA %>%
  filter(group == "E1 vs E2")
DA_H21_phase_plot_together_01_RA_E1E2$Species <-  factor(DA_H21_phase_plot_together_01_RA_E1E2$Species, 
                                                        levels = rev(taxonomy_plot_data_H21_E1E2$Species))

# Create the updated taxonomy plot
taxonomy_plot_H21_E1E2 <- ggplot(taxonomy_plot_data_H21_E1E2) +
  geom_text(aes(x = 0, y = Species, label = label_Family), hjust = 0, size = 7.85, family = "sans") +  # Move text left by adjusting x
  labs(title = "Family") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_DNA$Genus)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 28, vjust = -0.5, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_H21_E1E2

#To plot DA, factor levels
DA_H21_phase_plot_together_01_RA_E1E2$plot <- factor(DA_H21_phase_plot_together_01_RA_E1E2$plot, 
                                      levels = c("Log2 Fold change with 95%CI", 
                                                 "log(Bias-corrected abundances)", 
                                                 "log2(Relative abundances)"))
DA_H21_phase_plot_together_01_RA_E1E2$group <- factor(DA_H21_phase_plot_together_01_RA_E1E2$group, 
                                                             levels = c("L vs E2", 
                                                                        "T1 vs E2",
                                                                        "E1 vs E2",
                                                                        "T2 vs E2",
                                                                        "T3 vs E2"
                                                                        ))
DA_H21_phase_plot_together_01_RA_E1E2$Date_num_phase_abbrv <- factor(DA_H21_phase_plot_together_01_RA_E1E2$Date_num_phase_abbrv, 
                                                      levels = c(
                                                                 "L",
                                                                 "T1", 
                                                                 "E1", 
                                                                 "T2", 
                                                                 "T3",
                                                                 "E2"))

DA_H21_phase_plot_together_01_RA_E1E2$test <- factor(DA_H21_phase_plot_together_01_RA_E1E2$test, 
                                                     levels = c("ANCOM-BC",
                                                                "MaAsLin3"))

#DA plot - Just E1 vs E2
DA_H21_plot_MaAslinANCOM_E1E2 <-
  ggplot(data=DA_H21_phase_plot_together_01_RA_E1E2%>%filter(grepl("abundances", plot)),
         aes(y=Species, 
             x=logvalue, 
             fill = Date_num_phase_abbrv, 
             color = Date_num_phase_abbrv
         )) +
  geom_boxplot(orientation = "y", notch=FALSE, outlier.size=0.5, 
               size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75, orientation = "y")) +
  guides(color=guide_legend(order=1, title="Phase", title.position="top",
                            nrow = 1,
                            byrow = TRUE))+
  scale_color_manual(values = naive_phase_palette)+
  scale_fill_manual(values= naive_phase_palette)+
  new_scale_color()+
  geom_errorbar(
    inherit.aes = FALSE,
    data = DA_H21_phase_plot_together_01_RA_E1E2 %>%
      filter(plot == "Log2 Fold change with 95%CI"),
    aes(
      y = Species,
      xmin = lower.ci,
      xmax = upper.ci,
      color = direction,
      linetype = test),
    orientation = "y",
    #position = position_identity(),
    position = position_dodge(0.75, orientation = "y"),
    width=0.8,
    size = 0.75
  )+
  geom_point(
    inherit.aes = FALSE,
    data = DA_H21_phase_plot_together_01_RA_E1E2 %>%
      filter(plot == "Log2 Fold change with 95%CI"),
    aes(
      y = Species,
      x = coef,
      color = direction,
      shape = test),
    position = position_dodge(0.75, orientation = "y"),
    size = 1.75
  )+
  geom_vline(data=DA_H21_phase_plot_together_01_RA_E1E2%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(xintercept=0),
             size=0.5, 
             linetype='dashed', 
             alpha=0.5) +
  facet_nested(group~plot, 
               scales='free',
               #space='free',
               strip=strip_nested(text_y=list(element_text(angle=0))),
               labeller=labeller(group=label_wrap_gen(width=10),
                                 sub_group=label_wrap_gen(width=10))) +
  scale_y_discrete(position='left') +
  scale_x_continuous(position='top') +  # Default to break_labels otherwise 
  scale_color_manual(values=c(
    "depleted" = "red", 
    "elevated" = "blue")) +
  scale_linetype_manual(values = c(
    "MaAsLin3" = "solid",
    "ANCOM-BC" = "22"  
  ))+
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Phase", title.position="top",
                           nrow = 1,
                           byrow = TRUE),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 4)),
         linetype=guide_legend(title="Fold change source", title.position="top", override.aes = list(linewidth = 1.2)),
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 4))) +
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=22),
        legend.key.width = unit(1.5, "cm"),
        axis.title.x=element_blank(), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.text.x=element_text(size=24),
        axis.title.y=element_text(size= 28, angle=0, vjust= 1.03, face = "bold"), 
        axis.text.y=element_text(size=22, vjust = 0.5, face = 'italic'),
        strip.text.x=element_text(size=18, color = "white", face = "bold"),
        strip.text.y=element_text(size=26, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", 
        panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_H21_plot_MaAslinANCOM_E1E2

##Putting together DA plot with taxonomy (family) plot
combined_plot_H21_plot_MaAslinANCOM_E1E2 <- plot_grid(
  taxonomy_plot_H21_E1E2+ 
    theme(plot.margin = unit(c(0.28, -5, 0.2, 0.5), "cm")),  # Change margins
  DA_H21_plot_MaAslinANCOM_E1E2 + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm")
  ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  labels = c("A", " "), 
  label_size = 40, 
  axis = "tb")
combined_plot_H21_plot_MaAslinANCOM_E1E2

####Other phase comparisons######
#If doing Family:
#Create taxonomy plot data
taxonomy_plot_data_H21_others <- DA_H21_phase_plot_together_01_RA %>%
  filter(group != "E1 vs E2") %>%
  droplevels() %>%  # drop the unused "E1 vs E2" level if group is a factor
  left_join(input_taxonomy.H21, by = "Species") %>%
  distinct(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  arrange(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  mutate(
    Species_key = paste(group, Species, sep = "___"),
    Species_key = factor(Species_key, levels = rev(unique(Species_key)))
  ) %>%
  group_by(group, Family) %>%
  mutate(label_Family = ifelse(row_number() == 1, as.character(Family), "")) %>%
  ungroup()

#Order the species in the original DA plot data frame
DA_H21_phase_plot_together_01_RA_others <- DA_H21_phase_plot_together_01_RA %>%
  filter(group != "E1 vs E2") %>%
  droplevels() %>%
  mutate(
    Species_key = factor(paste(group, Species, sep = "___"),
                         levels = levels(taxonomy_plot_data_H21_others$Species_key))
  )

# Create the updated taxonomy plot
taxonomy_plot_H21_others <- ggplot(taxonomy_plot_data_H21_others) +
  geom_text(aes(x = 0, y = Species_key, label = label_Family),
            hjust = 0, size = 7.85, family = "sans") +
  labs(title = "Family") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 28, vjust = -0.5,
                              face = "bold", family = "sans"),
    strip.text = element_blank(),        # hide group labels here (show them on the main plot instead)
    panel.spacing.y = unit(0.5, "lines") # must match the main plot's spacing
  ) +
  coord_cartesian(xlim = c(0.05, 1))     # keeps the text pushed left
taxonomy_plot_H21_others

#To plot DA, factor levels
DA_H21_phase_plot_together_01_RA_others$plot <- factor(DA_H21_phase_plot_together_01_RA_others$plot, 
                                                          levels = c("Log2 Fold change with 95%CI", 
                                                                     "log(Bias-corrected abundances)", 
                                                                     "log2(Relative abundances)"))
DA_H21_phase_plot_together_01_RA_others$group <- factor(DA_H21_phase_plot_together_01_RA_others$group, 
                                                           levels = c("L vs E2", 
                                                                      "T1 vs E2",
                                                                      "E1 vs E2",
                                                                      "T2 vs E2",
                                                                      "T3 vs E2"
                                                           ))
DA_H21_phase_plot_together_01_RA_others$Date_num_phase_abbrv <- factor(DA_H21_phase_plot_together_01_RA_others$Date_num_phase_abbrv, 
                                                                          levels = c(
                                                                            "L",
                                                                            "T1", 
                                                                            "E1", 
                                                                            "T2", 
                                                                            "T3",
                                                                            "E2"))

DA_H21_phase_plot_together_01_RA_others$test <- factor(DA_H21_phase_plot_together_01_RA_others$test, 
                                                          levels = c("ANCOM-BC",
                                                                     "MaAsLin3"))

#DA plot - Other phase comparisons
DA_H21_plot_MaAslinANCOM_others <-
  ggplot(data = DA_H21_phase_plot_together_01_RA_others %>% filter(grepl("abundances", plot)),
         aes(y = Species_key, x = logvalue,
             fill = Date_num_phase_abbrv, color = Date_num_phase_abbrv)) +
  geom_boxplot(orientation = "y", notch = FALSE, outlier.size = 0.5, size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order=1, title="Phase", title.position="top",
                            nrow = 1,
                            byrow = TRUE))+
  scale_color_manual(values = naive_phase_palette)+
  scale_fill_manual(values= naive_phase_palette)+
  new_scale_color()+
  geom_errorbar(
    inherit.aes = FALSE,
    data = DA_H21_phase_plot_together_01_RA_others %>% filter(plot == "Log2 Fold change with 95%CI"),
    aes(y = Species_key, xmin = lower.ci, xmax = upper.ci, color = direction, linetype = test),
    orientation = "y", position = position_dodge(0.75), width = 0.8, linewidth = 0.75) +
  geom_point(
    inherit.aes = FALSE,
    data = DA_H21_phase_plot_together_01_RA_others %>% filter(plot == "Log2 Fold change with 95%CI"),
    aes(y = Species_key, x = coef, color = direction, shape = test),
    position = position_dodge(0.75), size = 1.75) +
  geom_vline(data = DA_H21_phase_plot_together_01_RA_others %>% filter(plot == "Log2 Fold change with 95%CI"),
             aes(xintercept = 0), linewidth = 0.5, linetype = "dashed", alpha = 0.5) +
  facet_nested(group ~ plot,
               scales = "free", space = "free_y",
               strip = strip_nested(text_y = list(element_text(angle = 0))),
               labeller = labeller(group = label_wrap_gen(width = 10))) +
  scale_y_discrete(position = "left", labels = function(x) sub(".*___", "", x)) +
  scale_x_continuous(position = "top") +
  scale_color_manual(values=c(
    "depleted" = "red", 
    "elevated" = "blue")) +
  scale_linetype_manual(values = c(
    "MaAsLin3" = "solid",
    "ANCOM-BC" = "22"  
  ))+
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Phase", title.position="top",
                           nrow = 1,
                           byrow = TRUE),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 4)),
         linetype=guide_legend(title="Fold change source", title.position="top", override.aes = list(linewidth = 1.2)),
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 4))) +
  theme_bw()+
  labs(y = "Species") +
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=22),
        legend.key.width = unit(1.5, "cm"),
        axis.title.x=element_blank(), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size= 28, angle=0, vjust= 1.03, face = "bold"), 
        axis.text.y=element_text(size=22, vjust = 0.5, face = 'italic'),
        strip.text.x=element_text(size=18, color = "white", face = "bold"),
        strip.text.y=element_text(size=26, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", 
        panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_H21_plot_MaAslinANCOM_others

##Putting together DA plot with taxonomy (family) plot
combined_plot_H21_plot_MaAslinANCOM_others <- plot_grid(
  taxonomy_plot_H21_others+ 
    theme(plot.margin = unit(c(0.28, 0, 0.2, 0.5), "cm")),  # Change margins
  DA_H21_plot_MaAslinANCOM_others + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm")
  ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb")
combined_plot_H21_plot_MaAslinANCOM_others

##ESTABLISHED SYSTEM ########
#For metadata, will use the imputed values (ph, Copper, Salinity, Temperature)
ancombc2_established_species_metadata <- phyloseq.samples_metadata%>%
  filter(Enclosure == "P1")%>%
  column_to_rownames(var = "SampleID")

#Only established system nitrifying taxa raw counts
phyloseq.samples.nitrifiers.P1 <- subset_samples(phyloseq.samples.nitrifiers, Enclosure == "P1")
phyloseq.samples.nitrifiers.P1 <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers.P1) > 0, 
                                              phyloseq.samples.nitrifiers.P1)
phyloseq.samples.nitrifiers.P1 #1685 taxa and 117 samples
ancombc2_established_nitrifiers_otu_table <- data.frame(otu_table(phyloseq.samples.nitrifiers.P1))
nrow(ancombc2_established_nitrifiers_otu_table) #1685 taxa


###ANCOMBC#######
##Untransformed (raw) counts in established samples of nitrifying taxa 
nrow(ancombc2_established_nitrifiers_otu_table) #1685 taxa

##Reorder the Date_num_phase_abbrv as factor, in established system E2 as reference
ancombc2_established_species_metadata$Date_num_phase_abbrv <- 
  factor(ancombc2_established_species_metadata$Date_num_phase_abbrv, levels = c("P",
                                                                          "L", 
                                                                          "T1", 
                                                                          "E"))
##running ancombc2
ancombc_output_P1_nitrifiers <-ancombc2(
  data= ancombc2_established_nitrifiers_otu_table, 
  meta_data = ancombc2_established_species_metadata,
  assay_name = "counts", 
  tax_level = "Species",
  fix_formula = "Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp",
  prv_cut = 0.05, 
  lib_cut = 0, 
  p_adj_method = "BH",
  dunnet = TRUE, 
  group= "Date_num_phase_abbrv", 
  struc_zero = TRUE, 
  neg_lb = TRUE,
  alpha = 0.05, #default significance
  verbose = TRUE)

## extract results from comparisons 
res_P1_nitrifiers <- ancombc_output_P1_nitrifiers$res_dunn

#Pivot longer the results (lfc, se, p val, q val, passed sensitivity, diff)
ancom_P1_phase_nitrifiers <- res_P1_nitrifiers %>%
  mutate(
    across(
      starts_with("lfc"),
      ~ round(.x, 2)
    )
  ) %>%
  pivot_longer(
    cols = matches("^(lfc|se|passed_ss|p|q|diff_robust)_"),
    names_to = c(".value", "group"),
    names_pattern = "(lfc|passed_ss|se|p|q|diff_robust)_(.*)"
  ) %>%
  filter(group != "(Intercept)")%>%
  arrange(taxon)
ancom_P1_phase_nitrifiers

##rework our group names so they're shorter and more manageable 
ancom_P1_phase_nitrifiers <- ancom_P1_phase_nitrifiers %>%
  mutate(group= case_when(
    group == "Date_num_phase_abbrvL" ~ "L vs P",
    group == "Date_num_phase_abbrvT1" ~ "T1 vs P",
    group == "Date_num_phase_abbrvE" ~ "E vs P",
    TRUE ~ group ##keeps original name for groups not specified 
  ))
ancom_P1_phase_nitrifiers$group ##Now the group names are shorter and more manageable

##Obtaining the confidence intervals for the log fold change
ancom_P1_phase_nitrifiers_2 <- ancom_P1_phase_nitrifiers %>%
  mutate(lower.ci = lfc - 1.96*se, 
         upper.ci = lfc + 1.96*se)

#Taxonomy of taxa in the established system
taxonomy.P1 <- data.frame(phyloseq::tax_table(phyloseq_with_unclassifieds_otu.samples_P1))%>%
  rownames_to_column(var = "OTU")

##Final fixup to make compatible with plotting. Only keeping nitrifiers. 
ancom_P1_phase_nitrifiers_nitirfiers <- ancom_P1_phase_nitrifiers_2 %>%
  left_join(taxonomy.P1, by = c("taxon" = 'OTU'))%>%
  filter (passed_ss == 1)%>%##Only want those that passed sensitivity testing
  # #Will only include the classified taxa
  rename(stderr = se, ##Renaming
         pval = p,  
         qval = q)%>%
  mutate(direction = ifelse(lfc > 0, "elevated", "depleted"),
         # Exponentiate the values first
         lfc_exp = exp(lfc),
         lower.ci_exp = exp(lower.ci),
         upper.ci_exp = exp(upper.ci)) %>%
  mutate(
    ##Getting log 2 so I can compare with maaslin
    coef =  log2(lfc_exp),
    lower.ci =  log2(lower.ci_exp),
    upper.ci =  log2(upper.ci_exp),
    plot = "Log2 Fold change with 95%CI", 
    test = "ANCOM-BC",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  filter(DA == "q ≤ 0.05")%>%
  select(Species, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA)##Only need these columns
length(unique(ancom_P1_phase_nitrifiers_nitirfiers$Species)) #298 taxa

#How many DA species per comparison? 
ancom_P1_phase_nitrifiers_nitirfiers %>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_nitrifiers = n())
# group             n_DA_nitrifiers
# E vs P            210
# L vs P            193
# T1 vs P           144


###MaAsLin3#######
#OTU table (RA) with unclassified reads 
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers #Nitrifier RA within the overall community (including unclassified reads)
#Select only established system
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1 <- 
  subset_samples(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers, 
                 Enclosure == "P1")
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1 <- 
  prune_taxa(taxa_sums(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1) > 0, 
             phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1)
phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1 #1685 taxa and 117 samples

#For otu table:
data_maaslin_P1 <- data.frame(t(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1@otu_table))

##Sample metadata. Will use imputed copper, salinity, ph, and temperature (and factor Date_num_phase_abbrv with E2 as reference)
metadata_maaslin_P1 <- phyloseq.samples_metadata%>%
  filter(Enclosure == "P1")%>%
  column_to_rownames(var = "SampleID")%>%
  select(c("Date_num_phase_abbrv", "pH_spu_imp", "Temperature_F_imp", "Salinity_ppt_imp"))%>%
  mutate(Date_num_phase_abbrv = 
           factor(Date_num_phase_abbrv, 
                  levels = c("P",
                             "L", 
                             "T1",
                             "E")
           )
  )#making Date_num_phase_abbrv a factor since I'll be adding it as a fixed effect on the MaAslin model

##MaAslin to peform differential abundance analysis using linear regression with log2 transformed relative abundances
maaslin_P1_nitrifiers <-  maaslin3(
  input_data = data_maaslin_P1, 
  input_metadata = metadata_maaslin_P1, 
  correction = 'BH',
  output = "MaAsLin3_P1",  
  fixed_effect = c("Date_num_phase_abbrv", "pH_spu_imp", "Temperature_F_imp", "Salinity_ppt_imp"), 
  min_prevalence= 0.05,
  median_comparison_abundance = T, #default
  median_comparison_prevalence = FALSE, #default  
  min_abundance = 0, #None
  normalization="NONE", #data is already RA
  transform="LOG", #(default LOG, base 2)
  max_significance=0.1, ##max q value
  standardize=F)
length(unique(maaslin_P1_nitrifiers$fit_data_abundance$results$feature)) #Abundance results from MaAslin - 1666 


#Calculate confidence intervals
maaslin_P1_nitrifiers_results <- maaslin_P1_nitrifiers$fit_data_abundance$results %>%
  mutate(lower.ci = coef - 1.96*stderr,
         upper.ci = coef + 1.96*stderr)%>%
  filter(!grepl("imp", name)) #not including temperature, salinity, pH results, only Date_num_phase
maaslin_P1_nitrifiers_results
length(unique(maaslin_P1_nitrifiers_results$feature)) #1666 taxa

##Taxonomy of P1 OTUs 
input_taxonomy.P1 <- data.frame(phyloseq_with_unclassifieds_otu.samples.species.ra.nitrifiers_P1@tax_table) %>%
  rownames_to_column(var = "feature")

##Final edits to put together for plot.Keep only nitrifiers
maaslin_P1_nitrifiers_results_2 <- maaslin_P1_nitrifiers_results%>%
  mutate(group = case_when(value == "L" ~ "L vs P", 
                           value == "T1" ~ "T1 vs P", 
                           value == "E" ~ "E vs P"))%>%
  rename(pval = pval_individual,
         qval = qval_individual)%>% ##Renaming
  mutate(direction = ifelse(coef > 0, "elevated", "depleted"))%>%
  mutate(
    plot = "Log2 Fold change with 95%CI", 
    test = "MaAsLin3",
    DA = case_when(
      qval <= 0.05 ~ "q ≤ 0.05",
      qval <= 0.1 ~ "q ≤ 0.1",
      TRUE ~ "Not significant")
  )%>%
  filter(DA == "q ≤ 0.05")%>%
  left_join(input_taxonomy.P1, by = "feature")%>%
  select(Species, group, coef, stderr, pval, qval, lower.ci, upper.ci, plot, test, direction, DA) 
maaslin_P1_nitrifiers_results_2

##ANCOM and MaAslin together
#Get the intersect per group 
DA_intersect_nitrifiers_P1 <-
  inner_join(ancom_P1_phase_nitrifiers_nitirfiers,
             maaslin_P1_nitrifiers_results_2,
             by = c("group", "Species"))
#Species that are differentially abundant per group AND agreed on by both MaAsLin and ANCOM”
DA_P1_nitrifiers_plot_MaAslinANCOM.data <-
  bind_rows(
    ancom_P1_phase_nitrifiers_nitirfiers,
    maaslin_P1_nitrifiers_results_2
  ) %>%
  semi_join(DA_intersect_nitrifiers_P1, by = c("group", "Species")) %>%
  mutate(group = factor(group, levels = c(
    "L vs P",
    "T1 vs P",
    "E vs P"
  )))

##BOX PLOTS OF RELATIVE ABUNDANCES (FOR MAASLIN) AND BIAS ADJUSTED ABUNDANCES (FOR ANCOM)
##Calculate the RA of the Genera included in the  MaAslin and ANCOM analyses for each sample
RA_MaaslinAncom_P1_nitrifiers <- maaslin_P1_nitrifiers$transformed_data %>% #transformed data is TSS transformed otu counts from MaAslin output
  rownames_to_column(var = "SampleID")%>%
  pivot_longer(cols = -SampleID, names_to = "OTU", values_to = "logvalue") %>%
  mutate(plot = 'log2(Relative abundances)') %>%
  left_join(input_taxonomy.P1, by = c("OTU" = "feature"))%>%
  filter(Species %in% intersect(maaslin_P1_nitrifiers_results_2$Species,
                                ancom_P1_phase_nitrifiers$Species))%>%
  left_join(metadata_maaslin_P1%>%rownames_to_column(var = "SampleID"), by = "SampleID")%>%
  select(Species, SampleID, logvalue, Date_num_phase_abbrv, plot)

# Non-P samples
RA_MaaslinAncom_P1_nitrifiers_nonP <-
  RA_MaaslinAncom_P1_nitrifiers %>%
  mutate(group = case_when(
    Date_num_phase_abbrv == "L"  ~ "L vs P",
    Date_num_phase_abbrv == "T1" ~ "T1 vs P",
    Date_num_phase_abbrv == "E" ~ "E vs P"
  )) %>%
  filter(!is.na(group))

# Duplicate P samples into every comparison
RA_MaaslinAncom_P1_nitrifiers_P <-
  RA_MaaslinAncom_P1_nitrifiers %>%
  filter(Date_num_phase_abbrv == "P") %>%
  tidyr::crossing(group = c(
    "L vs P", 
    "T1 vs P",
    "E vs P"
  ))

# Join data
RA_MaaslinAncom_P1_nitrifiers <-
  bind_rows(RA_MaaslinAncom_P1_nitrifiers_nonP, RA_MaaslinAncom_P1_nitrifiers_P)

#Bias-corrected abundances (ANCOMBC2)
nitrifiers_P1_species_log_corr_abn <- ancombc_output_P1_nitrifiers$bias_correct_log_table%>%
  data.frame()%>% ##make into data frame
  rownames_to_column("OTU")%>%
  left_join(input_taxonomy.P1, by = c("OTU" = "feature"))%>%
  filter(Species %in% intersect(maaslin_P1_nitrifiers_results_2$Species,
                                ancom_P1_phase_nitrifiers$Species))%>%
  select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, -OTU) %>%
  pivot_longer(cols = -Species, names_to = "SampleID", values_to = "logvalue")%>% ##Pivot to longer data format
  left_join(metadata_maaslin_P1%>%rownames_to_column(var="SampleID"), by = "SampleID")%>%
  select(Species, SampleID, logvalue, Date_num_phase_abbrv)%>%
  mutate(plot = 'log(Bias-corrected abundances)')
nitrifiers_P1_species_log_corr_abn

# Non-P samples
nitrifiers_P1_species_log_corr_abn_nonP <-
  nitrifiers_P1_species_log_corr_abn %>%
  mutate(group = case_when(
    Date_num_phase_abbrv == "L"  ~ "L vs P",
    Date_num_phase_abbrv == "T1" ~ "T1 vs P",
    Date_num_phase_abbrv == "E" ~ "E vs P"
  )) %>%
  filter(!is.na(group))

# Duplicate P samples into every comparison
nitrifiers_P1_species_log_corr_abn_P <-
  nitrifiers_P1_species_log_corr_abn %>%
  filter(Date_num_phase_abbrv == "P") %>%
  tidyr::crossing(group = c(
    "L vs P", 
    "T1 vs P",
    "E vs P"
  ))

# Binding datasets
nitrifiers_P1_species_log_corr_abn <-
  bind_rows(nitrifiers_P1_species_log_corr_abn_nonP, 
            nitrifiers_P1_species_log_corr_abn_P)

#Put together objects to plot DA
#First, only keeping RA for those differentially abundant taxa 
#Maaslin
RA_MaaslinAncom_P1_nitrifiers_differential  <-
  RA_MaaslinAncom_P1_nitrifiers %>%
  semi_join(
    DA_P1_nitrifiers_plot_MaAslinANCOM.data %>% select(group, Species),
    by = c("group", "Species")
  )
#Ancombc
nitrifiers_P1_species_log_corr_abn_differential <- 
  nitrifiers_P1_species_log_corr_abn %>%
  semi_join(
    DA_P1_nitrifiers_plot_MaAslinANCOM.data %>% select(group, Species),
    by = c("group", "Species")
  )

#Put RA, bias corrected abundances, DA coefficients dataframes together
DA_P1_phase_plot_together <- bind_rows(DA_P1_nitrifiers_plot_MaAslinANCOM.data, 
                                                   nitrifiers_P1_species_log_corr_abn_differential, 
                                                   RA_MaaslinAncom_P1_nitrifiers_differential) ##This object will have DA values, bias-corrected abundances from ANCOM, and relative abundances)

#How many DA nitrifying species per comparison? (without RA filtering)
DA_P1_phase_plot_together %>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# E vs P           197
# L vs P           152
# T1 vs P          110

#DA species - all nitrifying taxa (without RA filtering)
DA_P1_all_nitrifying_taxa_phase_comparisons <- DA_P1_phase_plot_together %>%
  filter(!is.na(coef)) %>%
  group_by(Species, group) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(n_groups = n_distinct(group)) %>%   # how many comparisons each species appears in
  ungroup() %>%
  arrange(desc(n_groups), Species, group)
DA_P1_all_nitrifying_taxa_phase_comparisons
unique((DA_P1_all_nitrifying_taxa_phase_comparisons %>%filter(group == "E vs P"))$Species) #OK!

#How many belonging to Nitrosopumilus (without RA filtering)
DA_P1_phase_plot_together %>%
  filter(grepl("Nitrosopumilus", Species))%>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# E vs P           125
# L vs P            51
# T1 vs P           69

#Which Nitrosopumilus species are they? (without RA filtering)
DA_P1_Nitrosopumilus_phase_comparisons <- DA_P1_phase_plot_together %>%
  filter(!is.na(coef)) %>%
  filter(grepl("Nitrosopumilus", Species)) %>%
  group_by(Species, group) %>%
  #slice(1) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(n_groups = n_distinct(group)) %>%   # how many comparisons each species appears in
  ungroup() %>%
  arrange(desc(n_groups), Species, group)
DA_P1_Nitrosopumilus_phase_comparisons

#Which Nitrosopumilus species are DA against P across all phases? (without RA filtering)
DA_P1_all_phases_comparisons_P <- DA_P1_phase_plot_together %>%
  filter(!is.na(coef)) %>%
  filter(grepl("Nitrosopumilus", Species)) %>%
  group_by(Species) %>%
  mutate(n_species = n()) %>%   # total rows per species
  ungroup() %>%
  filter(n_species == 10) %>% #10 rows means they were DA across all comparisons (5 groups, 2 tests)
  arrange(desc(n_species), Species)
unique(DA_P1_all_phases_comparisons_P$Species)#none


####Filtering out those species that make up less than 0.1% RA of the nitrifying community########
any(sample_sums(phyloseq.samples.nitrifiers_P1)== 0) ## no samples with 0 OTUs
phyloseq.samples.nitrifiers_P1.ra <- transform_sample_counts(phyloseq.samples.nitrifiers_P1, 
                                                              function(x) x/sum(x)*100)
phyloseq.samples.nitrifiers_P1.ra #1685 taxa and 117 samples
sample_sums(phyloseq.samples.nitrifiers_P1.ra) #RA to 100

phyloseq.samples.nitrifiers_P1.ra.filt.DA <- filter_taxa(phyloseq.samples.nitrifiers_P1.ra, function(x) mean(x) > 0.1, TRUE) 
phyloseq.samples.nitrifiers_P1.ra.filt.DA ##48 OTUs with mean RA > 0.1% across samples (established system)


#Filter them from the DA plot object
DA_P1_phase_plot_together_01_RA <- DA_P1_phase_plot_together %>%
  filter(Species %in% phyloseq::tax_table(phyloseq.samples.nitrifiers_P1.ra.filt.DA)[,"Species"])
#How many DA nitrifying species per comparison? 
DA_P1_phase_plot_together_01_RA %>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group   n_DA_species
# E vs P            29
# L vs P            30
# T1 vs P           27

#How many belonging to Nitrosopumilus
DA_P1_phase_plot_together_01_RA %>%
  filter(grepl("Nitrosopumilus", Species))%>%
  distinct(group, Species) %>%
  group_by(group) %>%
  summarise(n_DA_species = n())
# group    n_DA_species
# E vs P            20
# L vs P            16
# T1 vs P           17

#Which Nitrosopumilus species are they?
DA_P1_Nitrosopumilus_phase_comparisons <- DA_P1_phase_plot_together_01_RA %>%
  filter(!is.na(coef)) %>%
  filter(grepl("Nitrosopumilus", Species)) %>%
  group_by(Species, group) %>%
  #slice(1) %>%
  ungroup() %>%
  group_by(Species) %>%
  mutate(n_groups = n_distinct(group)) %>%   # how many comparisons each species appears in
  ungroup() %>%
  arrange(desc(n_groups), Species, group)
DA_P1_Nitrosopumilus_phase_comparisons

####SUPPPLEMENTARY TABLE 11.2 ######
stable11.2 <- DA_P1_all_nitrifying_taxa_phase_comparisons%>%
  mutate(`Effect Size (95% CI)` = 
           paste0(round(coef, 2), " (", 
                  round(lower.ci,2),
                  ", ",
                  round(upper.ci,2),
                  ")"), 
         `P value` = format(pval, scientific = TRUE, digits = 3), 
         `Q value` = format(qval, scientific = TRUE, digits = 3), 
         Comparison = group,
         Test = test,
         System = "established"
  )%>%
  select(System, Species, Comparison,  Test, `Effect Size (95% CI)`, `P value`, `Q value`)
stable11.2


###PLOTTING DA - NITRIFIERS#####
####Only E VS P ########
##Ordering how I want the "Family" tax level to show up on the plot 
input_taxonomy.P1 ##dataframe object for Taxonomy of P1 taxa

#Add nitrifier group
input_taxonomy.P1 <- input_taxonomy.P1 %>%
  mutate(Nitrifying_group = case_when(
    Family == "Nitrosomonadaceae" ~ "AOB",
    Family == "Nitrosococcaceae" ~ "AOB",
    Family == "Nitrosocaldaceae" ~ "AOA",
    Family == "Nitrosopumilaceae" ~ "AOA",
    Family == "Nitrososphaeraceae" ~ "AOA",
    Family == "UBA213" ~ "AOA",
    Family == "Nitrospiraceae" ~ "NOB",
    Family == "NS-4" ~ "NOB",
    Family == "UBA8639" ~ "NOB",
    Genus == "Nitrobacter" ~ "NOB",
    Genus == "Nitrotoga" ~ "NOB",
    Genus == "Nitrococcus" ~ "NOB",
    Order == "Nitrospinales" ~ "NOB",
    Genus == "Nitrolancea" ~ "NOB",
    TRUE ~ NA_character_))%>%
  mutate(Nitrifying_group = factor(Nitrifying_group, levels = c("AOA", "AOB", "NOB")))

#Taxonomy plot data
taxonomy_plot_data_P1_EP <- DA_P1_phase_plot_together_01_RA %>%
  left_join(input_taxonomy.P1, by= "Species")%>%
  filter(group == "E vs P")%>%
  distinct(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  arrange(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  mutate(Species = factor(Species, levels = rev(Species)))%>% ##Since I arranged by family, this is the order I want the genera to show up
  dplyr::group_by(Family) %>%
  dplyr::mutate(label_Family = ifelse(row_number() == 1, Family, "")) %>%  # Create 'label_Family' with only the first occurrence of each Family
  ungroup()
taxonomy_plot_data_P1_EP

##Factor "Species" level by the order I want (taxonomy_plot_data_P1_EP$Species)
DA_P1_phase_plot_together_01_RA_EP <- DA_P1_phase_plot_together_01_RA %>%
  filter(group == "E vs P")
DA_P1_phase_plot_together_01_RA_EP$Species <-  factor(DA_P1_phase_plot_together_01_RA_EP$Species, 
                                                              levels = rev(taxonomy_plot_data_P1_EP$Species))

# Create the updated taxonomy plot
taxonomy_plot_P1_EP <- ggplot(taxonomy_plot_data_P1_EP) +
  geom_text(aes(x = 0, y = Species, label = label_Family), hjust = 0, size = 7.85, family = "sans") +  # Move text left by adjusting x
  labs(title = "Family") +
  theme_void() +
  #scale_y_discrete(limits =  rev(taxonomy_plot_data_DNA$Genus)) +  # Ensure y-axis matches the taxon order
  theme(plot.title = element_text(hjust = 0, size = 28, vjust = -0.5, face = "bold", family = "sans"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()
  )+
  coord_cartesian(xlim = c(0.05, 1)) #had to add this to move the geom_text more to the left
taxonomy_plot_P1_EP

#To plot DA, factor levels
DA_P1_phase_plot_together_01_RA_EP$plot <- factor(DA_P1_phase_plot_together_01_RA_EP$plot, 
                                                          levels = c("Log2 Fold change with 95%CI", 
                                                                     "log(Bias-corrected abundances)", 
                                                                     "log2(Relative abundances)"))
DA_P1_phase_plot_together_01_RA_EP$group <- factor(DA_P1_phase_plot_together_01_RA_EP$group, 
                                                           levels = c("L vs P", 
                                                                      "T1 vs P",
                                                                      "E vs P",
                                                                      "T2 vs P",
                                                                      "T3 vs P"
                                                           ))
DA_P1_phase_plot_together_01_RA_EP$Date_num_phase_abbrv <- factor(DA_P1_phase_plot_together_01_RA_EP$Date_num_phase_abbrv, 
                                                                          levels = c(
                                                                            "L",
                                                                            "T1", 
                                                                            "E", 
                                                                            "T2", 
                                                                            "T3",
                                                                            "P"))

DA_P1_phase_plot_together_01_RA_EP$test <- factor(DA_P1_phase_plot_together_01_RA_EP$test, 
                                                          levels = c("ANCOM-BC",
                                                                     "MaAsLin3"))

#DA plot - Just E vs P
DA_P1_plot_MaAslinANCOM_EP <-
  ggplot(data=DA_P1_phase_plot_together_01_RA_EP%>%filter(grepl("abundances", plot)),
         aes(y=Species, 
             x=logvalue, 
             fill = Date_num_phase_abbrv, 
             color = Date_num_phase_abbrv
         )) +
  geom_boxplot(orientation = "y", notch=FALSE, outlier.size=0.5, 
               size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75, orientation = "y")) +
  guides(color=guide_legend(order=1, title="Phase", title.position="top",
                            nrow = 1,
                            byrow = TRUE))+
  scale_color_manual(values = established_phase_palette)+
  scale_fill_manual(values= established_phase_palette)+
  new_scale_color()+
  geom_errorbar(
    inherit.aes = FALSE,
    data = DA_P1_phase_plot_together_01_RA_EP %>%
      filter(plot == "Log2 Fold change with 95%CI"),
    aes(
      y = Species,
      xmin = lower.ci,
      xmax = upper.ci,
      color = direction,
      linetype = test),
    orientation = "y",
    #position = position_identity(),
    position = position_dodge(0.75, orientation = "y"),
    width=0.8,
    size = 0.75
  )+
  geom_point(
    inherit.aes = FALSE,
    data = DA_P1_phase_plot_together_01_RA_EP %>%
      filter(plot == "Log2 Fold change with 95%CI"),
    aes(
      y = Species,
      x = coef,
      color = direction,
      shape = test),
    position = position_dodge(0.75, orientation = "y"),
    size = 1.75
  )+
  geom_vline(data=DA_P1_phase_plot_together_01_RA_EP%>%filter(plot == 'Log2 Fold change with 95%CI'),
             aes(xintercept=0),
             size=0.5, 
             linetype='dashed', 
             alpha=0.5) +
  facet_nested(group~plot, 
               scales='free',
               #space='free',
               strip=strip_nested(text_y=list(element_text(angle=0))),
               labeller = labeller(
                 group = as_labeller(c("Log2 Fold change with 95%CI" = "Log2 Fold change with 95%CI", 
                                       "log(Bias-corrected abundance)" = "log(Bias-corrected abundance)",
                                       "log2(Relative abundances)" = "log2(Relative abundances)",
                                       "E vs P" = "\u00A0\u00A0E vs P\u00A0" #adding padding
                                       ), 
                 default = label_wrap_gen(width = 10)))
  )+
  scale_y_discrete(position='left') +
  scale_x_continuous(position='top') +  # Default to break_labels otherwise 
  scale_color_manual(values=c(
    "depleted" = "red", 
    "elevated" = "blue")) +
  scale_linetype_manual(values = c(
    "MaAsLin3" = "solid",
    "ANCOM-BC" = "22"  
  ))+
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Phase", title.position="top",
                           nrow = 1,
                           byrow = TRUE),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 4)),
         linetype=guide_legend(title="Fold change source", title.position="top", override.aes = list(linewidth = 1.2)),
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 4))) +
  theme_bw()+
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=22),
        legend.key.width = unit(1.5, "cm"),
        axis.title.x=element_blank(), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size= 28, angle=0, vjust= 1.03, face = "bold"), 
        axis.text.y=element_text(size=22, vjust = 0.5, face = 'italic'),
        strip.text.x=element_text(size=18, color = "white", face = "bold"),
        strip.text.y=element_text(size=26, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", 
        panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_P1_plot_MaAslinANCOM_EP

##Putting together DA plot with taxonomy (family) plot
combined_plot_P1_plot_MaAslinANCOM_EP <- plot_grid(
  taxonomy_plot_P1_EP+ 
    theme(plot.margin = unit(c(0.28, -5, 0.2, 0.5), "cm")),
  DA_P1_plot_MaAslinANCOM_EP + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -1.5), "cm")
  ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  labels = c("B", " "), 
  label_size = 40, 
  axis = "tb")
combined_plot_P1_plot_MaAslinANCOM_EP


####Other phase comparisons######
#If doing Family:
#Create taxonomy plot data
taxonomy_plot_data_P1_others <- DA_P1_phase_plot_together_01_RA %>%
  filter(group != "E vs P") %>%
  droplevels() %>%  # drop the unused "E vs P" level if group is a factor
  left_join(input_taxonomy.P1, by = "Species") %>%
  distinct(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  arrange(group, Nitrifying_group, Phylum, Order, Family, Genus, Species) %>%
  mutate(
    Species_key = paste(group, Species, sep = "___"),
    Species_key = factor(Species_key, levels = rev(unique(Species_key)))
  ) %>%
  group_by(group, Family) %>%
  mutate(label_Family = ifelse(row_number() == 1, as.character(Family), "")) %>%
  ungroup()

#Order the species in the original DA plot data frame
DA_P1_phase_plot_together_01_RA_others <- DA_P1_phase_plot_together_01_RA %>%
  filter(group != "E vs P") %>%
  droplevels() %>%
  mutate(
    Species_key = factor(paste(group, Species, sep = "___"),
                         levels = levels(taxonomy_plot_data_P1_others$Species_key))
  )

# Create the updated taxonomy plot
taxonomy_plot_P1_others <- ggplot(taxonomy_plot_data_P1_others) +
  geom_text(aes(x = 0, y = Species_key, label = label_Family),
            hjust = 0, size = 7.85, family = "sans") +
  labs(title = "Family") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 28, vjust = -0.5,
                              face = "bold", family = "sans"),
    strip.text = element_blank(),        # hide group labels here (show them on the main plot instead)
    panel.spacing.y = unit(0.5, "lines") # must match the main plot's spacing
  ) +
  coord_cartesian(xlim = c(0.05, 1))     # keeps the text pushed left
taxonomy_plot_P1_others

#To plot DA, factor levels
DA_P1_phase_plot_together_01_RA_others$plot <- factor(DA_P1_phase_plot_together_01_RA_others$plot, 
                                                            levels = c("Log2 Fold change with 95%CI", 
                                                                       "log(Bias-corrected abundances)", 
                                                                       "log2(Relative abundances)"))
DA_P1_phase_plot_together_01_RA_others$group <- factor(DA_P1_phase_plot_together_01_RA_others$group, 
                                                             levels = c("L vs P", 
                                                                        "T1 vs P",
                                                                        "E vs P",
                                                                        "T2 vs P",
                                                                        "T3 vs P"
                                                             ))
DA_P1_phase_plot_together_01_RA_others$Date_num_phase_abbrv <- factor(DA_P1_phase_plot_together_01_RA_others$Date_num_phase_abbrv, 
                                                                            levels = c(
                                                                              "L",
                                                                              "T1", 
                                                                              "E", 
                                                                              "T2", 
                                                                              "T3",
                                                                              "P"))

DA_P1_phase_plot_together_01_RA_others$test <- factor(DA_P1_phase_plot_together_01_RA_others$test, 
                                                            levels = c("ANCOM-BC",
                                                                       "MaAsLin3"))

#DA plot - Other phase comparisons
DA_P1_plot_MaAslinANCOM_others <-
  ggplot(data = DA_P1_phase_plot_together_01_RA_others %>% filter(grepl("abundances", plot)),
         aes(y = Species_key, x = logvalue,
             fill = Date_num_phase_abbrv, color = Date_num_phase_abbrv)) +
  geom_boxplot(orientation = "y", notch = FALSE, outlier.size = 0.5, size = 0.5, alpha = 0.3) +
  geom_point(size = 2, shape = 18, position = position_dodge(width = 0.75)) +
  guides(color=guide_legend(order=1, title="Phase", title.position="top",
                            nrow = 1,
                            byrow = TRUE))+
  scale_color_manual(values = established_phase_palette)+
  scale_fill_manual(values= established_phase_palette)+
  new_scale_color()+
  geom_errorbar(
    inherit.aes = FALSE,
    data = DA_P1_phase_plot_together_01_RA_others %>% filter(plot == "Log2 Fold change with 95%CI"),
    aes(y = Species_key, xmin = lower.ci, xmax = upper.ci, color = direction, linetype = test),
    orientation = "y", position = position_dodge(0.75), width = 0.8, linewidth = 0.75) +
  geom_point(
    inherit.aes = FALSE,
    data = DA_P1_phase_plot_together_01_RA_others %>% filter(plot == "Log2 Fold change with 95%CI"),
    aes(y = Species_key, x = coef, color = direction, shape = test),
    position = position_dodge(0.75), size = 1.75) +
  geom_vline(data = DA_P1_phase_plot_together_01_RA_others %>% filter(plot == "Log2 Fold change with 95%CI"),
             aes(xintercept = 0), linewidth = 0.5, linetype = "dashed", alpha = 0.5) +
  facet_nested(group ~ plot,
               scales = "free", 
               space = "free_y",
               strip = strip_nested(text_y = list(element_text(angle = 0))),
               labeller=labeller(group=label_wrap_gen(width=10),
                                 sub_group=label_wrap_gen(width=10))) +
  scale_y_discrete(position = "left", labels = function(x) sub(".*___", "", x)) +
  scale_x_continuous(position = "top") +
  scale_color_manual(values=c(
    "depleted" = "red", 
    "elevated" = "blue")) +
  scale_linetype_manual(values = c(
    "MaAsLin3" = "solid",
    "ANCOM-BC" = "22"  
  ))+
  scale_shape_manual(values=c(16, 15)) +
  guides(fill=guide_legend(order=1, title="Phase", title.position="top",
                           nrow = 1,
                           byrow = TRUE),
         color=guide_legend(order=2, title="Fold change direction", title.position="top", override.aes = list(size = 4)),
         linetype=guide_legend(title="Fold change source", title.position="top", override.aes = list(linewidth = 1.2)),
         shape=guide_legend(title="Fold change source", title.position="top", override.aes = list(size = 4))) +
  theme_bw()+
  labs(y = "Species") +
  theme(legend.position="top", legend.key=element_blank(),
        legend.title=element_text(size=22), 
        legend.text=element_text(size=22),
        legend.key.width = unit(1.5, "cm"),
        axis.title.x=element_blank(), 
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10),  # top, right, bottom, left
        axis.text.x=element_text(size=16),
        axis.title.y=element_text(size=28, angle=0, vjust= 1.03, face = "bold"), 
        axis.text.y=element_text(size=26, vjust = 0.5, face = 'italic'),
        strip.text.x=element_text(size=16, color = "white", face = "bold"),
        strip.text.y=element_text(size=28, color = "white", face = "bold"),
        strip.background=element_rect(fill='black'
                                      , color='white'),
        strip.placement="outside", 
        panel.spacing.y=unit(0.5, "lines"),
        panel.grid.minor = element_blank())
DA_P1_plot_MaAslinANCOM_others

##Putting together DA plot with taxonomy (family) plot
combined_plot_P1_plot_MaAslinANCOM_others <- plot_grid(
  taxonomy_plot_P1_others+ 
    theme(plot.margin = unit(c(0.28, 0, 0.2, 0.5), "cm")),  # Change margins
  DA_P1_plot_MaAslinANCOM_others + theme(plot.margin = unit(c(0.2, 0.5, 0.2, -0.5), "cm")
  ),
  ncol = 2, 
  rel_widths = c(0.4,2),  # Adjust the width ratio as needed
  align = "h", 
  axis = "tb")
combined_plot_P1_plot_MaAslinANCOM_others

###FIGURE 5#######
figure5 <- plot_grid(combined_plot_H21_plot_MaAslinANCOM_E1E2, 
                     combined_plot_P1_plot_MaAslinANCOM_EP, 
                     align = "v", 
                     ncol = 1)
figure5
ggsave("Figure5.png", 
       figure5, 
       device = "png", 
       dpi = 500, 
       width = 22.5, 
       height = 24, 
       bg = "white")


###SUPPLEMENTARY FIGURES 5 & 6#######
#Supp figure 5 other phase comparisons in the naive system
sfigure5 <- combined_plot_H21_plot_MaAslinANCOM_others
ggsave("SupplementaryFigure5.png", 
       sfigure5, 
       device = "png",
       dpi = 500, 
       width = 22.5, 
       height = 24,
       bg ="white")

#Supp figure 6 other phase comparisons in the established system
sfigure6 <- combined_plot_P1_plot_MaAslinANCOM_others
ggsave("SupplementaryFigure6.png", 
       sfigure6, 
       device = "png", 
       dpi = 500, 
       width = 22.5, 
       height = 24, 
       bg ="white")

#BETA DIV####
##OVERALL COMMUNITIES - DISTANCES BETWEEN CONSECUTIVE DAYS#########
###DISTANCES IN NAIVE (H21) SYSTEM####
phyloseq.samples_H21 #This phyloseq object only has samples from the naive system (with classified taxa)
#First, need to order samples according to Collection Date
str(phyloseq.samples_H21@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for naive system (H21)
phyloseq.samples_H21_metadata <- data.frame(sample_data(phyloseq.samples_H21))

# Order by Date
phyloseq.samples_H21_metadata <- phyloseq.samples_H21_metadata[order(phyloseq.samples_H21_metadata$Collection_Date),
                                                                                 ]
phyloseq.samples_H21_metadata$SampleID #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0220 (Feb 20, 2024)

#Extract OTU table for naive system (H21), then reorder according to date as in metadata
phyloseq.samples_H21_otu <- data.frame(phyloseq.samples_H21@otu_table)
phyloseq.samples_H21_otu <- phyloseq.samples_H21_otu %>%
  select(phyloseq.samples_H21_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
phyloseq.samples_H21.ordered <- phyloseq.samples_H21
colnames(phyloseq.samples_H21.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
phyloseq.samples_H21.ordered@otu_table <- otu_table(phyloseq.samples_H21_otu, taxa_are_rows = T)
phyloseq.samples_H21.ordered@sam_data <- sample_data(phyloseq.samples_H21_metadata) 
colnames(phyloseq.samples_H21.ordered@otu_table) #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0220 (Feb 20, 2024)


####BRAY CURTIS#######
#First, normalize counts (Relative abundance)
phyloseq.samples_H21.ordered_RA <- transform_sample_counts(phyloseq.samples_H21.ordered, 
                                                                  function(x) x/sum(x)*100)
sample_sums(phyloseq.samples_H21.ordered_RA) #OK

#Calculate BC distances
phyloseq.samples_H21.ordered_RA.bray <- vegdist(t(phyloseq.samples_H21.ordered_RA@otu_table), method = "bray")
phyloseq.samples_H21.ordered_RA.bray

#Make into square matrix
phyloseq.samples_H21.ordered_RA.bray_matrix <- as.matrix(phyloseq.samples_H21.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_H21_BC_dist <- expand.grid(
  i = 1:nrow(phyloseq.samples_H21.ordered_RA.bray_matrix),
  j = 1:nrow(phyloseq.samples_H21.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_H21_BC_dist) #7569 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_H21_BC_dist <- all_pairs_comparisons_H21_BC_dist[all_pairs_comparisons_H21_BC_dist$i < all_pairs_comparisons_H21_BC_dist$j, ]
nrow(all_pairs_comparisons_H21_BC_dist) #3741 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_H21_BC_dist$BC_dist <- mapply(function(i, j) {
  phyloseq.samples_H21.ordered_RA.bray_matrix[i, j]
  }, all_pairs_comparisons_H21_BC_dist$i, all_pairs_comparisons_H21_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_H21_BC_dist$time_diff <- as.numeric(
  phyloseq.samples_H21_metadata$Collection_Date[all_pairs_comparisons_H21_BC_dist$j] - phyloseq.samples_H21_metadata$Collection_Date[all_pairs_comparisons_H21_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_H21_BC_dist <- all_pairs_comparisons_H21_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_H21 <- rownames(phyloseq.samples_H21.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_H21_BC_dist <- all_pairs_comparisons_H21_BC_dist%>%
  mutate(sample_i = sample_ids_H21[all_pairs_comparisons_H21_BC_dist$i], 
         sample_j = sample_ids_H21[all_pairs_comparisons_H21_BC_dist$j])
all_pairs_comparisons_H21_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_H21_BC <- all_pairs_comparisons_H21_BC_dist %>%
  filter(j == i + 1)
first_diff_df_H21_BC$sample_j <- factor(first_diff_df_H21_BC$sample_j, levels = first_diff_df_H21_BC$sample_j)

#Add some more metadata for plotting 
first_diff_df_H21_BC <- first_diff_df_H21_BC %>%
  mutate(
    Collection_Date_i = phyloseq.samples_H21_metadata$Collection_Date[i],
    Collection_Date_j = phyloseq.samples_H21_metadata$Collection_Date[j],
    Date_num_i = phyloseq.samples_H21_metadata$Date_num[i],
    Date_num_j = phyloseq.samples_H21_metadata$Date_num[j], 
    Enclosure = "H21"
    )


###DISTANCES IN ESTABLISHED (P1) SYSTEM####
phyloseq.samples_P1 #This phyloseq object only has samples from the established system
#First, need to order samples according to Collection Date
str(phyloseq.samples_P1@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for established system (P1)
phyloseq.samples_P1_metadata <- data.frame(sample_data(phyloseq.samples_P1))

# Order by Date
phyloseq.samples_P1_metadata <- phyloseq.samples_P1_metadata[order(phyloseq.samples_P1_metadata$Collection_Date), ]
phyloseq.samples_P1_metadata$SampleID #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0402 (April 02, 2024)

#Extract OTU table for established system (P1), then reorder according to date as in metadata
phyloseq.samples_P1_otu <- data.frame(phyloseq.samples_P1@otu_table)
phyloseq.samples_P1_otu <- phyloseq.samples_P1_otu %>%
  select(phyloseq.samples_P1_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
phyloseq.samples_P1.ordered <- phyloseq.samples_P1
colnames(phyloseq.samples_P1.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
phyloseq.samples_P1.ordered@otu_table <- otu_table(phyloseq.samples_P1_otu, taxa_are_rows = T)
phyloseq.samples_P1.ordered@sam_data <- sample_data(phyloseq.samples_P1_metadata) 
colnames(phyloseq.samples_P1.ordered@otu_table) #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0402 (April 02, 2024)


####BRAY CURTIS#######
#Normalize counts (Relative abundance)
phyloseq.samples_P1.ordered_RA <- transform_sample_counts(phyloseq.samples_P1.ordered, 
                                                                          function(x) x/sum(x)*100)
sample_sums(phyloseq.samples_P1.ordered_RA)

#Calculate BC distances
phyloseq.samples_P1.ordered_RA.bray <- vegdist(t(phyloseq.samples_P1.ordered_RA@otu_table), method = "bray")
phyloseq.samples_P1.ordered_RA.bray

#Make into square matrix
phyloseq.samples_P1.ordered_RA.bray_matrix <- as.matrix(phyloseq.samples_P1.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_P1_BC_dist <- expand.grid(
  i = 1:nrow(phyloseq.samples_P1.ordered_RA.bray_matrix),
  j = 1:nrow(phyloseq.samples_P1.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_P1_BC_dist) #13689 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_P1_BC_dist <- all_pairs_comparisons_P1_BC_dist[all_pairs_comparisons_P1_BC_dist$i < all_pairs_comparisons_P1_BC_dist$j, ]
nrow(all_pairs_comparisons_P1_BC_dist) #6786 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_P1_BC_dist$BC_dist <- mapply(function(i, j) {
  phyloseq.samples_P1.ordered_RA.bray_matrix[i, j]
}, all_pairs_comparisons_P1_BC_dist$i, all_pairs_comparisons_P1_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_P1_BC_dist$time_diff <- as.numeric(
  phyloseq.samples_P1_metadata$Collection_Date[all_pairs_comparisons_P1_BC_dist$j] - phyloseq.samples_P1_metadata$Collection_Date[all_pairs_comparisons_P1_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_P1_BC_dist <- all_pairs_comparisons_P1_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_P1 <- rownames(phyloseq.samples_P1.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_P1_BC_dist <- all_pairs_comparisons_P1_BC_dist%>%
  mutate(sample_i = sample_ids_P1[all_pairs_comparisons_P1_BC_dist$i], 
         sample_j = sample_ids_P1[all_pairs_comparisons_P1_BC_dist$j])
all_pairs_comparisons_P1_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_P1_BC <- all_pairs_comparisons_P1_BC_dist %>%
  filter(j == i + 1)
first_diff_df_P1_BC$sample_j <- factor(first_diff_df_P1_BC$sample_j, levels = first_diff_df_P1_BC$sample_j)

#Add some more metadata for plotting 
first_diff_df_P1_BC <- first_diff_df_P1_BC %>%
  mutate(
    Collection_Date_i = phyloseq.samples_P1_metadata$Collection_Date[i],
    Collection_Date_j = phyloseq.samples_P1_metadata$Collection_Date[j],
    Date_num_i = phyloseq.samples_P1_metadata$Date_num[i],
    Date_num_j = phyloseq.samples_P1_metadata$Date_num[j], 
    Enclosure = "P1"
  )


###JOINING DISTANCES IN ESTABLISHED (P1) AND NAIVE (H21) SYSTEMS####
first_diff_df_all_BC <- bind_rows(first_diff_df_H21_BC, 
                               first_diff_df_P1_BC)
first_diff_df_all_BC

##Add a row for date_num_j to be "1" just so it is plotted, but it will have no values
first_diff_df_all_BC <- first_diff_df_all_BC %>%
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
first_diff_df_all_BC

#### PLOT#######
first_diff_BC_H21_P1_plot <- ggplot(first_diff_df_all_BC, aes(x = factor(Date_num_j))) +
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


##NITRIFYING COMMUNITIES - DISTANCES BETWEEN CONSECUTIVE DAYS#########
phyloseq.samples.nitrifiers #nitrifying taxa

###DISTANCES IN NAIVE (H21) SYSTEM####
phyloseq.samples.nitrifiers_H21 # 1686 taxa and 87 samples

#First, need to order samples according to Collection Date
str(phyloseq.samples.nitrifiers_H21@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for naive system (H21)
phyloseq.samples.nitrifiers_H21_metadata <- data.frame(sample_data(phyloseq.samples.nitrifiers_H21))

# Order by Date
phyloseq.samples.nitrifiers_H21_metadata <- phyloseq.samples.nitrifiers_H21_metadata[order(phyloseq.samples.nitrifiers_H21_metadata$Collection_Date), ]
phyloseq.samples.nitrifiers_H21_metadata$SampleID #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0220 (February 20, 2024)

#Extract OTU table for naive system (H21), then reorder according to date as in metadata
phyloseq.samples.nitrifiers_H21.ra_OTU <- data.frame(phyloseq.samples.nitrifiers_H21@otu_table)
phyloseq.samples.nitrifiers_H21.ra_OTU <- phyloseq.samples.nitrifiers_H21.ra_OTU %>%
  select(phyloseq.samples.nitrifiers_H21_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
phyloseq.samples.nitrifiers_H21.ordered <- phyloseq.samples.nitrifiers_H21
colnames(phyloseq.samples.nitrifiers_H21.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
phyloseq.samples.nitrifiers_H21.ordered@otu_table <- otu_table(phyloseq.samples.nitrifiers_H21.ra_OTU, taxa_are_rows = T)
phyloseq.samples.nitrifiers_H21.ordered@sam_data <- sample_data(phyloseq.samples.nitrifiers_H21_metadata) 
colnames(phyloseq.samples.nitrifiers_H21.ordered@otu_table) #Ok, now it starts H21_1009 (Oct 09, 2023) and ends H21_0220 (February 20, 2024)

####BRAY CURTIS######
#Now, normalize counts (Relative abundance)
phyloseq.samples.nitrifiers_H21.ordered_RA <- transform_sample_counts(phyloseq.samples.nitrifiers_H21.ordered, 
                                                                          function(x) x/sum(x)*100)
sample_sums(phyloseq.samples.nitrifiers_H21.ordered_RA)

#Calculate BC distances
phyloseq.samples.nitrifiers_H21.ordered_RA.bray <- vegdist(t(phyloseq.samples.nitrifiers_H21.ordered_RA@otu_table), method = "bray")
phyloseq.samples.nitrifiers_H21.ordered_RA.bray

#Make into square matrix
phyloseq.samples.nitrifiers_H21.ordered_RA.bray_matrix <- as.matrix(phyloseq.samples.nitrifiers_H21.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_H21_nit_BC_dist <- expand.grid(
  i = 1:nrow(phyloseq.samples.nitrifiers_H21.ordered_RA.bray_matrix),
  j = 1:nrow(phyloseq.samples.nitrifiers_H21.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_H21_nit_BC_dist) #8464 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_H21_nit_BC_dist <- all_pairs_comparisons_H21_nit_BC_dist[all_pairs_comparisons_H21_nit_BC_dist$i < all_pairs_comparisons_H21_nit_BC_dist$j, ]
nrow(all_pairs_comparisons_H21_nit_BC_dist) #4186 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_H21_nit_BC_dist$BC_dist <- mapply(function(i, j) {
  phyloseq.samples.nitrifiers_H21.ordered_RA.bray_matrix[i, j]
}, all_pairs_comparisons_H21_nit_BC_dist$i, all_pairs_comparisons_H21_nit_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_H21_nit_BC_dist$time_diff <- as.numeric(
  phyloseq.samples.nitrifiers_H21_metadata$Collection_Date[all_pairs_comparisons_H21_nit_BC_dist$j] - phyloseq.samples.nitrifiers_H21_metadata$Collection_Date[all_pairs_comparisons_H21_nit_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_H21_nit_BC_dist <- all_pairs_comparisons_H21_nit_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_H21_nit <- rownames(phyloseq.samples.nitrifiers_H21.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_H21_nit_BC_dist <- all_pairs_comparisons_H21_nit_BC_dist%>%
  mutate(sample_i = sample_ids_H21_nit[all_pairs_comparisons_H21_nit_BC_dist$i], 
         sample_j = sample_ids_H21_nit[all_pairs_comparisons_H21_nit_BC_dist$j])
all_pairs_comparisons_H21_nit_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_H21_nit_BC <- all_pairs_comparisons_H21_nit_BC_dist %>%
  filter(j == i + 1)
first_diff_df_H21_nit_BC$sample_j <- factor(first_diff_df_H21_nit_BC$sample_j, levels = first_diff_df_H21_nit_BC$sample_j)

#Add some more metadata for plotting 
first_diff_df_H21_nit_BC <- first_diff_df_H21_nit_BC %>%
  mutate(
    Collection_Date_i = phyloseq.samples.nitrifiers_H21_metadata$Collection_Date[i],
    Collection_Date_j = phyloseq.samples.nitrifiers_H21_metadata$Collection_Date[j],
    Date_num_i = phyloseq.samples.nitrifiers_H21_metadata$Date_num[i],
    Date_num_j = phyloseq.samples.nitrifiers_H21_metadata$Date_num[j], 
    Enclosure = "H21"
  )


###DISTANCES IN ESTABLISHED (P1) SYSTEM####
phyloseq.samples.nitrifiers_P1 #1685 taxa and 117 samples from the established system

#First, need to order samples according to Collection Date
str(phyloseq.samples.nitrifiers_P1@sam_data$Collection_Date) #Collection Date is an as.Date format 

# Extract metadata for established system (P1)
phyloseq.samples.nitrifiers_P1_metadata <- data.frame(sample_data(phyloseq.samples.nitrifiers_P1))

# Order by Date
phyloseq.samples.nitrifiers_P1_metadata <- phyloseq.samples.nitrifiers_P1_metadata[order(phyloseq.samples.nitrifiers_P1_metadata$Collection_Date), ]
phyloseq.samples.nitrifiers_P1_metadata$SampleID #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0402 (April 02, 2024)

#Extract OTU table for established system (P1), then reorder according to date as in metadata
phyloseq.samples.nitrifiers_P1_otu <- data.frame(phyloseq.samples.nitrifiers_P1@otu_table)
phyloseq.samples.nitrifiers_P1_otu <- phyloseq.samples.nitrifiers_P1_otu %>%
  select(phyloseq.samples.nitrifiers_P1_metadata$SampleID)%>%
  as.matrix() ##turn it back into a matrix so it is compatible with otu_table function from phyloseq

#Make a new phyloseq object, where samples are actually ordered by Date
phyloseq.samples.nitrifiers_P1.ordered <- phyloseq.samples.nitrifiers_P1
colnames(phyloseq.samples.nitrifiers_P1.ordered@otu_table) #Not ordered yet
#Replace the OTU table and metadata with ordered samples
phyloseq.samples.nitrifiers_P1.ordered@otu_table <- otu_table(phyloseq.samples.nitrifiers_P1_otu, taxa_are_rows = T)
phyloseq.samples.nitrifiers_P1.ordered@sam_data <- sample_data(phyloseq.samples.nitrifiers_P1_metadata) 
colnames(phyloseq.samples.nitrifiers_P1.ordered@otu_table) #Ok, now it starts P1_1114 (Nov 14, 2023) and ends P1_0402 (April 02, 2024)

####BRAY CURTIS######
#Now, normalize counts (Relative abundance)
phyloseq.samples.nitrifiers_P1.ordered_RA <- transform_sample_counts(phyloseq.samples.nitrifiers_P1.ordered, 
                                                                         function(x) x/sum(x)*100)
sample_sums(phyloseq.samples.nitrifiers_P1.ordered_RA)

#Calculate BC distances
phyloseq.samples.nitrifiers_P1.ordered_RA.bray <- vegdist(t(phyloseq.samples.nitrifiers_P1.ordered_RA@otu_table), method = "bray")
phyloseq.samples.nitrifiers_P1.ordered_RA.bray

#Make into square matrix
phyloseq.samples.nitrifiers_P1.ordered_RA.bray_matrix <- as.matrix(phyloseq.samples.nitrifiers_P1.ordered_RA.bray)

#Get all possible combinations of sample indices
all_pairs_comparisons_P1_nit_BC_dist <- expand.grid(
  i = 1:nrow(phyloseq.samples.nitrifiers_P1.ordered_RA.bray_matrix),
  j = 1:nrow(phyloseq.samples.nitrifiers_P1.ordered_RA.bray_matrix))
nrow(all_pairs_comparisons_P1_nit_BC_dist) #13689 comparisons (including duplicates, ex 1,2 and 2,1)

#Now, keep only unique pairwise comparisons
all_pairs_comparisons_P1_nit_BC_dist <- all_pairs_comparisons_P1_nit_BC_dist[all_pairs_comparisons_P1_nit_BC_dist$i < all_pairs_comparisons_P1_nit_BC_dist$j, ]
nrow(all_pairs_comparisons_P1_nit_BC_dist) #6786 unique comparisons

#Now, get BC distance for each pairwise comparison (i, j). Loops over row index i and column index j. 
all_pairs_comparisons_P1_nit_BC_dist$BC_dist <- mapply(function(i, j) {
  phyloseq.samples.nitrifiers_P1.ordered_RA.bray_matrix[i, j]
}, all_pairs_comparisons_P1_nit_BC_dist$i, all_pairs_comparisons_P1_nit_BC_dist$j)

#Now, compute time difference between sample pairs
all_pairs_comparisons_P1_nit_BC_dist$time_diff <- as.numeric(
  phyloseq.samples.nitrifiers_P1_metadata$Collection_Date[all_pairs_comparisons_P1_nit_BC_dist$j] - phyloseq.samples.nitrifiers_P1_metadata$Collection_Date[all_pairs_comparisons_P1_nit_BC_dist$i]
)

#Get a rate of change
all_pairs_comparisons_P1_nit_BC_dist <- all_pairs_comparisons_P1_nit_BC_dist%>%
  mutate(bc_rate_of_change = BC_dist / time_diff)

#Add samples 
#First, get sampleIDs from orderes matrix
sample_ids_P1_nit <- rownames(phyloseq.samples.nitrifiers_P1.ordered_RA.bray_matrix)
#Now, add them to the dataframe
all_pairs_comparisons_P1_nit_BC_dist <- all_pairs_comparisons_P1_nit_BC_dist%>%
  mutate(sample_i = sample_ids_P1_nit[all_pairs_comparisons_P1_nit_BC_dist$i], 
         sample_j = sample_ids_P1_nit[all_pairs_comparisons_P1_nit_BC_dist$j])
all_pairs_comparisons_P1_nit_BC_dist

#Okay, finally, get just the first differential (how much the community composition changes between consecutive index samples)
first_diff_df_P1_BC_nit <- all_pairs_comparisons_P1_nit_BC_dist %>%
  filter(j == i + 1)
first_diff_df_P1_BC_nit$sample_j <- factor(first_diff_df_P1_BC_nit$sample_j, 
                                        levels = first_diff_df_P1_BC_nit$sample_j)

#Add some more metadata for plotting 
first_diff_df_P1_BC_nit <- first_diff_df_P1_BC_nit %>%
  mutate(
    Collection_Date_i = phyloseq.samples.nitrifiers_P1_metadata$Collection_Date[i],
    Collection_Date_j = phyloseq.samples.nitrifiers_P1_metadata$Collection_Date[j],
    Date_num_i = phyloseq.samples.nitrifiers_P1_metadata$Date_num[i],
    Date_num_j = phyloseq.samples.nitrifiers_P1_metadata$Date_num[j], 
    Enclosure = "P1"
  )


###JOINING DISTANCES IN ESTABLISHED (P1) AND NAIVE (H21) SYSTEMS####
first_diff_df_nit_BC <- bind_rows(first_diff_df_H21_nit_BC, 
                               first_diff_df_P1_BC_nit)
first_diff_df_nit_BC

##Add a row for date_num_j to be "1" just so it is plotted, but it will have no values
first_diff_df_nit_BC <- first_diff_df_nit_BC %>%
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
first_diff_df_nit_BC

#### PLOT#######
first_diff_nit_BC_H21_P1_plot <- ggplot(first_diff_df_nit_BC, aes(x = factor(Date_num_j))) +
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


####OVERALL COMMUNITIES#####
#####NAIVE SYSTEM#####
#Pull the pairwise BC distances. Vector of day-to-day change n community composition
first_diff_df_H21_BC <- first_diff_df_H21_BC %>%
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
changepoint_overall_communities_h21_days <- first_diff_df_H21_BC %>%
  filter(i %in% changepoint_overall_communities_h21_indeces) %>% #i is the index
  select(i, Date_num_i, Date_num_j) #Date_num_j is the day+1 change
changepoint_overall_communities_h21_days

gam_fit_H21_overrall <- gamm(BC_dist ~ s(Date_num_i),
                             correlation = corAR1(form = ~ Date_num_i),
                             data = first_diff_df_H21_BC)
plot(gam_fit_H21_overrall$gam, shade = TRUE)

##BC DISTANCES TO MEAN OF NITRIFYING COMMUNITY AFTER FIRST EXPOSURE TO COPPER########
###Using just the nitrifier community composition######
####NAIVE SYSTEM#####
#Do phases differ in degree of convergence to E2?##

#Only naive system samples
phyloseq.samples.nitrifiers_H21 #1686 taxa and 87 samples from the naive system
phyloseq.samples.nitrifiers_H21.metadata <- data.frame(phyloseq.samples.nitrifiers_H21@sam_data) #metadata

#Get RA
phyloseq.samples.nitrifiers_H21.ra #1686 taxa and 87 samples

#Get the mean community composition of nitrifiers during E2
phyloseq.samples.nitrifiers_H21.ra_phase_E2 <- subset_samples(phyloseq.samples.nitrifiers_H21.ra,
                                                                      Date_num_phase_abbrv == "E2")
# phyloseq.samples.nitrifiers_H21.ra_phase_E2 <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers_H21.ra_phase_E2) > 0,
#                                                           phyloseq.samples.nitrifiers_H21.ra_phase_E2) #Not pruning taxa, want to keep all 1686 even if they are 0 in E2
phyloseq.samples.nitrifiers_H21.ra_phase_E2 #1686 taxa and 19 samples

#Get mean RA composition of E2
H21_E2_mean_nitrifier_composition_RA <- rowMeans(phyloseq.samples.nitrifiers_H21.ra_phase_E2@otu_table)

#Add it to the RA otu table as a named row
phyloseq.samples.nitrifiers_H21.ra_plus_ref_E2 <- rbind(
  as.matrix(t(phyloseq.samples.nitrifiers_H21.ra@otu_table)),
  E2_reference = H21_E2_mean_nitrifier_composition_RA
)

#Bray-Curtis, then pull out distances to the E2 mean
BC_nitrifiers_H21_distance <- as.matrix(
  vegdist(phyloseq.samples.nitrifiers_H21.ra_plus_ref_E2, method = "bray"))

BC_nitrifiers_H21_distance_to_E2 <-
  BC_nitrifiers_H21_distance[sample_names(phyloseq.samples.nitrifiers_H21.ra), "E2_reference"]

#Make into dataframe
BC_nitrifiers_H21_distance_to_E2_df <- data.frame(
  SampleID = names(BC_nitrifiers_H21_distance_to_E2),
  BC_to_E2 = as.numeric(BC_nitrifiers_H21_distance_to_E2)
)
#Merge with metadata
BC_nitrifiers_H21_distance_to_E2_metadata_df <- left_join(BC_nitrifiers_H21_distance_to_E2_df, 
                                                      phyloseq.samples.nitrifiers_H21.metadata, 
                                                      by = "SampleID")
BC_nitrifiers_H21_distance_to_E2_metadata_df

#Omnibus test (Kruskal)
BC_nitrifiers_H21_distance_to_E2_metadata_df %>%
  filter(Date_num_phase_abbrv != "E2") %>%
  kruskal_test(BC_to_E2 ~ Date_num_phase_abbrv)
# .y.          n statistic    df     p method        
# BC_to_E2    68      6.62     4 0.157 Kruskal-Wallis

#Wilcox Test
BC_nitrifiers_H21_distance_to_E2_metadata_df %>%
  pairwise_wilcox_test(
    BC_to_E2 ~ Date_num_phase_abbrv,
    p.adjust.method = "BH"
  ) 
# wilcox_effsize(BC_to_E2 ~ Date_num_phase_abbrv)

#Reorder for plot
phase_order_for_plot_h21_bc <- c("L", "T1", "E1", "T2", "T3")   # chronological, no E2

#Reorder factor on dataframe
BC_nitrifiers_H21_distance_to_E2_metadata_df_plot <- 
  BC_nitrifiers_H21_distance_to_E2_metadata_df %>%
  filter(Date_num_phase_abbrv != "E2") %>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = phase_order_for_plot_h21_bc))

#Keep only comparisons in chronological order
BC_nitrifiers_H21_distance_to_E2_metadata_stat_test <- 
  BC_nitrifiers_H21_distance_to_E2_metadata_df_plot %>%
  pairwise_wilcox_test(
    BC_to_E2 ~ Date_num_phase_abbrv,
    p.adjust.method = "BH",
    comparisons = list(
      c("L", "T1"),
      c("T1", "E1"),
      c("E1", "T2"),
      c("T2", "T3")
    )
  ) %>%
  add_xy_position(x = "Date_num_phase_abbrv")
BC_nitrifiers_H21_distance_to_E2_metadata_stat_test

#Plot
BC_nitrifiers_distance_to_E2_H21_boxplot <- ggplot(BC_nitrifiers_H21_distance_to_E2_metadata_df_plot,
                                                           aes(x = Date_num_phase_abbrv, 
                                                               y= BC_to_E2, 
                                                               color = Date_num_phase_abbrv, 
                                                               fill = Date_num_phase_abbrv)) +
  theme_bw() +
  labs(x = "Phase", 
       y = "BC Distance",
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3, 
               outlier.shape = NA) +
  scale_fill_manual(values = naive_phase_palette)+
  scale_color_manual(values = naive_phase_palette)+
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
  stat_pvalue_manual(BC_nitrifiers_H21_distance_to_E2_metadata_stat_test, 
                     label = "p.adj.signif", 
                     label.size = 10,
                     hide.ns = T)
BC_nitrifiers_distance_to_E2_H21_boxplot

####ESTABLISHED SYSTEM######
#Only established system samples
phyloseq.samples.nitrifiers_P1 #1685 taxa and 117 samples from the established system
phyloseq.samples.nitrifiers_P1.metadata <- data.frame(phyloseq.samples.nitrifiers_P1@sam_data) #metadata

#Get RA
phyloseq.samples.nitrifiers_P1.ra #1685 taxa and 117 samples 

#Get the mean community composition of nitrifiers during P
phyloseq.samples.nitrifiers_P1.ra_phase_P <- subset_samples(phyloseq.samples.nitrifiers_P1.ra,
                                                              Date_num_phase_abbrv == "P")
# phyloseq.samples.nitrifiers_P1.ra_phase_P <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers_P1.ra_phase_P) > 0,
#                                                           phyloseq.samples.nitrifiers_P1.ra_phase_P) #Not pruning taxa, want to keep all 1685 even if they are 0 in P
phyloseq.samples.nitrifiers_P1.ra_phase_P #1685 taxa and 24 samples
P1_P_mean_nitrifier_composition_RA <- rowMeans(phyloseq.samples.nitrifiers_P1.ra_phase_P@otu_table)

#Add the mean to the RA otu table as a named row
phyloseq.samples.nitrifiers_P1.ra_plus_ref_P <- rbind(
  as.matrix(t(phyloseq.samples.nitrifiers_P1.ra@otu_table)),
  P_reference = P1_P_mean_nitrifier_composition_RA
)
dim(phyloseq.samples.nitrifiers_P1.ra_plus_ref_P) #should be 118 x 1685

#Bray-Curtis, then pull out distances to the P mean
BC_nitrifiers_P1_distance <- as.matrix(
  vegdist(phyloseq.samples.nitrifiers_P1.ra_plus_ref_P, method = "bray"))

BC_nitrifiers_P1_distance_to_P <-
  BC_nitrifiers_P1_distance[sample_names(phyloseq.samples.nitrifiers_P1.ra), "P_reference"]
BC_nitrifiers_P1_distance_to_P

#Make into dataframe
BC_nitrifiers_P1_distance_to_P_df <- data.frame(
  SampleID = names(BC_nitrifiers_P1_distance_to_P),
  BC_to_P = as.numeric(BC_nitrifiers_P1_distance_to_P)
)
#Merge with metadata
BC_nitrifiers_P1_distance_to_P_metadata_df <- merge(BC_nitrifiers_P1_distance_to_P_df, 
                                                      phyloseq.samples.nitrifiers_P1.metadata, 
                                                      by = "SampleID")
BC_nitrifiers_P1_distance_to_P_metadata_df

#Omnibus test
BC_nitrifiers_P1_distance_to_P_metadata_df %>%
  filter(Date_num_phase_abbrv != "P") %>%
  kruskal_test(BC_to_P ~ Date_num_phase_abbrv) 
# .y.         n statistic    df        p method        
# BC_to_P    93      14.0     2 0.000923 Kruskal-Wallis

#Wilcox Test
BC_nitrifiers_P1_distance_to_P_metadata_df %>%
  pairwise_wilcox_test(
    BC_to_P ~ Date_num_phase_abbrv,
    p.adjust.method = "BH"
  ) 
# wilcox_effsize(BC_to_P ~ Date_num_phase_abbrv)

#Reorder for plot
phase_order_for_plot_p1_bc <- c("L", "T1", "E")   # chronological, no P

#Reorder factor on dataframe
BC_nitrifiers_P1_distance_to_P_metadata_df_plot <- 
  BC_nitrifiers_P1_distance_to_P_metadata_df %>%
  filter(Date_num_phase_abbrv != "P") %>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = phase_order_for_plot_p1_bc))

#Keep only comparisons in chronological order
BC_nitrifiers_P1_distance_to_P_metadata_stat_test <- 
  BC_nitrifiers_P1_distance_to_P_metadata_df %>%
  pairwise_wilcox_test(
    BC_to_P ~ Date_num_phase_abbrv,
    p.adjust.method = "BH",
    comparisons = list(
      c("L", "T1"),
      c("T1", "E")
    )
  ) %>%
  add_xy_position(x = "Date_num_phase_abbrv")
BC_nitrifiers_P1_distance_to_P_metadata_stat_test


#Plot
BC_nitrifiers_distance_to_P_P1_boxplot <- ggplot(BC_nitrifiers_P1_distance_to_P_metadata_df_plot,
                                                   aes(x = Date_num_phase_abbrv, 
                                                       y= BC_to_P, 
                                                       color = Date_num_phase_abbrv, 
                                                       fill = Date_num_phase_abbrv)) +
  theme_bw() +
  labs(x = "Phase", 
       y = "BC Distance",
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3, 
               outlier.shape = NA) +
  scale_fill_manual(values = established_phase_palette)+
  scale_color_manual(values = established_phase_palette)+
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
  stat_pvalue_manual(BC_nitrifiers_P1_distance_to_P_metadata_stat_test, 
                     label = "p.adj.signif", 
                     label.size = 10,
                     hide.ns = T)
BC_nitrifiers_distance_to_P_P1_boxplot


###Using RA of nitrifiers within the overall community#####
####NAIVE SYSTEM#####
#Do phases differ in degree of convergence to E2?##
phyloseq.samples.species.ra.nitrifiers #nitrifying taxa, RA within the overall community
phyloseq.samples.species.ra.nitrifiers.metadata <- data.frame(phyloseq.samples.species.ra.nitrifiers@sam_data) #metadata

#Only naive system samples
phyloseq.samples.species.ra.nitrifiers_H21 <- subset_samples(phyloseq.samples.species.ra.nitrifiers,
                                                                      Enclosure == "H21")
phyloseq.samples.species.ra.nitrifiers_H21 <- prune_taxa(taxa_sums(phyloseq.samples.species.ra.nitrifiers_H21) > 0,
                                                         phyloseq.samples.species.ra.nitrifiers_H21) 
#Get the mean community composition of nitrifiers during E2
phyloseq.samples.species.ra.nitrifiers_H21_phase_E2 <- subset_samples(phyloseq.samples.species.ra.nitrifiers_H21,
                                                 Date_num_phase_abbrv == "E2")
phyloseq.samples.species.ra.nitrifiers_H21_phase_E2 #1686 taxa and 19 samples
# phyloseq.samples.species.ra.nitrifiers_H21_phase_E2 <- prune_taxa(taxa_sums(phyloseq.samples.species.ra.nitrifiers_H21_phase_E2) > 0,
#                                              phyloseq.samples.species.ra.nitrifiers_H21_phase_E2) #Not pruning the taxa sums that are O in E2, will use them for the distance calculations
H21_E2_mean_nitrifier_overall_composition_RA <- rowMeans(phyloseq.samples.species.ra.nitrifiers_H21_phase_E2@otu_table)

#Add the mean to the RA otu table as a named row
phyloseq.samples.species.ra.nitrifiers_H21_plus_ref_E2 <- rbind(
  as.matrix(t(phyloseq.samples.species.ra.nitrifiers_H21@otu_table)),
  E2_reference = H21_E2_mean_nitrifier_overall_composition_RA
)
dim(phyloseq.samples.species.ra.nitrifiers_H21_plus_ref_E2) #should be 88 rows x 1686 

#Bray-Curtis, then pull out distances to the E2 mean
BC_nitrifiers_overall_H21_distance <- as.matrix(
  vegdist(phyloseq.samples.species.ra.nitrifiers_H21_plus_ref_E2, method = "bray"))

BC_nitrifiers_overall_H21_distance_to_E2 <-
  BC_nitrifiers_overall_H21_distance[sample_names(phyloseq.samples.species.ra.nitrifiers_H21), "E2_reference"]
BC_nitrifiers_overall_H21_distance_to_E2

#Make into dataframe
BC_nitrifiers_overall_H21_distance_to_E2_df <- data.frame(
  SampleID = names(BC_nitrifiers_overall_H21_distance_to_E2),
  BC_to_E2 = as.numeric(BC_nitrifiers_overall_H21_distance_to_E2)
)
#Merge with metadata
BC_nitrifiers_overall_H21_distance_to_E2_metadata_df <- left_join(BC_nitrifiers_overall_H21_distance_to_E2_df, 
                                                  phyloseq.samples.species.ra.nitrifiers.metadata, 
                                                  by = "SampleID")
BC_nitrifiers_overall_H21_distance_to_E2_metadata_df

#Omnibus test (Kruskal)
BC_nitrifiers_overall_H21_distance_to_E2_metadata_df %>%
  filter(Date_num_phase_abbrv != "E2") %>%
  kruskal_test(BC_to_E2 ~ Date_num_phase_abbrv)
# .y.          n statistic    df        p method        
# BC_to_E2    68      19.3     4 0.000682 Kruskal-Wallis

#Wilcox Test
BC_nitrifiers_overall_H21_distance_to_E2_metadata_df %>%
  pairwise_wilcox_test(
    BC_to_E2 ~ Date_num_phase_abbrv,
    p.adjust.method = "BH"
  ) 
 # wilcox_effsize(BC_to_E2 ~ Date_num_phase_abbrv)

#Reorder for plot
phase_order_for_plot_h21_bc <- c("L", "T1", "E1", "T2", "T3")   # chronological, no E2

#Reorder factor on dataframe
BC_nitrifiers_overall_H21_distance_to_E2_metadata_df_plot <- 
  BC_nitrifiers_overall_H21_distance_to_E2_metadata_df %>%
  filter(Date_num_phase_abbrv != "E2") %>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = phase_order_for_plot_h21_bc))

#Keep only comparisons in chronological order
BC_nitrifiers_overall_H21_distance_to_E2_metadata_stat_test <- 
  BC_nitrifiers_overall_H21_distance_to_E2_metadata_df_plot %>%
  pairwise_wilcox_test(
    BC_to_E2 ~ Date_num_phase_abbrv,
    p.adjust.method = "BH",
    comparisons = list(
      c("L", "T1"),
      c("T1", "E1"),
      c("E1", "T2"),
      c("T2", "T3")
    )
  ) %>%
  add_xy_position(x = "Date_num_phase_abbrv")
BC_nitrifiers_overall_H21_distance_to_E2_metadata_stat_test

#Plot
BC_nitrifiers_distance_to_E2_overall_H21_boxplot <- ggplot(BC_nitrifiers_overall_H21_distance_to_E2_metadata_df_plot,
                                             aes(x = Date_num_phase_abbrv, 
                                                 y= BC_to_E2, 
                                                 color = Date_num_phase_abbrv, 
                                                 fill = Date_num_phase_abbrv)) +
  theme_bw() +
  labs(x = "Phase", 
       y = "BC Distance",
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3, 
               outlier.shape = NA) +
  scale_fill_manual(values = naive_phase_palette)+
  scale_color_manual(values = naive_phase_palette)+
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
  stat_pvalue_manual(BC_nitrifiers_overall_H21_distance_to_E2_metadata_stat_test, 
                     label = "p.adj.signif", 
                     label.size = 10,
                     hide.ns = T)
BC_nitrifiers_distance_to_E2_overall_H21_boxplot

####ESTABLISHED SYSTEM#####
#Do phases differ in degree of convergence to P?##
phyloseq.samples.species.ra.nitrifiers #nitrifying taxa, RA within tne overall community
phyloseq.samples.species.ra.nitrifiers.metadata <- data.frame(phyloseq.samples.species.ra.nitrifiers@sam_data) #metadata

#Only established system samples
phyloseq.samples.species.ra.nitrifiers.P1 <- subset_samples(phyloseq.samples.species.ra.nitrifiers,
                                                             Enclosure == "P1")
phyloseq.samples.species.ra.nitrifiers.P1 <- prune_taxa(taxa_sums(phyloseq.samples.species.ra.nitrifiers.P1) > 0,
                                                        phyloseq.samples.species.ra.nitrifiers.P1)

#Get the mean community composition of nitrifiers during P
phyloseq.samples.species.ra.nitrifiers_P1_phase_P <- subset_samples(phyloseq.samples.species.ra.nitrifiers.P1,
                                                                               Date_num_phase_abbrv == "P")
# phyloseq.samples.species.ra.nitrifiers_P1_phase_P <- prune_taxa(taxa_sums(phyloseq.samples.species.ra.nitrifiers_P1_phase_P) > 0,
#                                                                            phyloseq.samples.species.ra.nitrifiers_P1_phase_P) #Not pruning the taxa sums that are O in P, including them in the distance calculation
phyloseq.samples.species.ra.nitrifiers_P1_phase_P #1685 taxa and 24 samples
P1_P_mean_nitrifier_overall_composition_RA <- rowMeans(phyloseq.samples.species.ra.nitrifiers_P1_phase_P@otu_table)

#Add the mean to the RA otu table as a named row
phyloseq.samples.species.ra.nitrifiers.P1_plus_ref_P <- rbind(
  as.matrix(t(phyloseq.samples.species.ra.nitrifiers.P1@otu_table)),
  P_reference = P1_P_mean_nitrifier_overall_composition_RA
)
dim(phyloseq.samples.species.ra.nitrifiers.P1_plus_ref_P) #should be 118 x 1685

#Bray-Curtis, then pull out distances to the P mean
BC_nitrifiers_P1_distance_overall <- as.matrix(
  vegdist(phyloseq.samples.species.ra.nitrifiers.P1_plus_ref_P, method = "bray"))

BC_nitrifiers_P1_distance_to_P_overall <-
  BC_nitrifiers_P1_distance_overall[sample_names(phyloseq.samples.species.ra.nitrifiers.P1), "P_reference"]
BC_nitrifiers_P1_distance_to_P_overall

#Make into dataframe
BC_nitrifiers_P1_distance_to_P_overall_df <- data.frame(
  SampleID = names(BC_nitrifiers_P1_distance_to_P_overall),
  BC_to_P = as.numeric(BC_nitrifiers_P1_distance_to_P_overall)
)
#Merge with metadata
BC_nitrifiers_P1_distance_to_P_overall_metadata_df <- merge(BC_nitrifiers_P1_distance_to_P_overall_df, 
                                                            phyloseq.samples.species.ra.nitrifiers.metadata, 
                                                      by = "SampleID")
BC_nitrifiers_P1_distance_to_P_overall_metadata_df


#Omnibus test
BC_nitrifiers_P1_distance_to_P_overall_metadata_df %>%
  filter(Date_num_phase_abbrv != "P") %>%
  kruskal_test(BC_to_P ~ Date_num_phase_abbrv) 
# .y        n statistic    df            p method        
# BC_to_P    93      33.8     2 0.0000000447 Kruskal-Wallis

#Wilcox Test
BC_nitrifiers_P1_distance_to_P_overall_metadata_df %>%
  pairwise_wilcox_test(
    BC_to_P ~ Date_num_phase_abbrv,
    p.adjust.method = "BH"
  ) 
# wilcox_effsize(BC_to_P ~ Date_num_phase_abbrv)

#Reorder for plot
phase_order_for_plot_p1_bc <- c("L", "T1", "E")   # chronological, no P

#Reorder factor on dataframe
BC_nitrifiers_P1_distance_to_P_overall_metadata_df_plot <- 
  BC_nitrifiers_P1_distance_to_P_overall_metadata_df %>%
  filter(Date_num_phase_abbrv != "P") %>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = phase_order_for_plot_p1_bc))

#Keep only comparisons in chronological order
BC_nitrifiers_P1_distance_to_P_overall_metadata_stat_test <- 
  BC_nitrifiers_P1_distance_to_P_overall_metadata_df_plot %>%
  pairwise_wilcox_test(
    BC_to_P ~ Date_num_phase_abbrv,
    p.adjust.method = "BH",
    comparisons = list(
      c("L", "T1"),
      c("T1", "E")
    )
  ) %>%
  add_xy_position(x = "Date_num_phase_abbrv")
BC_nitrifiers_P1_distance_to_P_overall_metadata_stat_test

#Plot
BC_nitrifiers_distance_to_P_overall_P1_boxplot <- ggplot(BC_nitrifiers_P1_distance_to_P_overall_metadata_df_plot,
                                                   aes(x = Date_num_phase_abbrv, 
                                                       y= BC_to_P, 
                                                       color = Date_num_phase_abbrv, 
                                                       fill = Date_num_phase_abbrv)) +
  theme_bw() +
  labs(x = "Phase", 
       y = "BC Distance",
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3,
               outlier.shape = NA) +
  scale_fill_manual(values = established_phase_palette)+
  scale_color_manual(values = established_phase_palette)+
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
  stat_pvalue_manual(BC_nitrifiers_P1_distance_to_P_overall_metadata_stat_test, 
                     label = "p.adj.signif", 
                     label.size = 10,
                     hide.ns = T)
BC_nitrifiers_distance_to_P_overall_P1_boxplot

##AITCHISONS DISTANCES TO MEAN OF NITRIFYING COMMUNITY AFTER FIRST EXPOSURE TO COPPER########
###Using just the nitrifier community composition######
####NAIVE SYSTEM#####
#Do phases differ in degree of convergence to E2?##

#Only naive system samples
phyloseq.samples.nitrifiers_H21 #1686 taxa and 87 samples from the naive system
phyloseq.samples.nitrifiers_H21.metadata <- data.frame(phyloseq.samples.nitrifiers_H21@sam_data) #metadata

#Get the geometric mean composition of nitrifiers during E2
phyloseq.samples.nitrifiers_H21_phase_E2 <- subset_samples(phyloseq.samples.nitrifiers_H21,
                                                           Date_num_phase_abbrv == "E2")
phyloseq.samples.nitrifiers_H21_phase_E2 #1686 taxa and 19 samples

#Need pseudocount for Aitchisons
pseudocount <- 0.5


#Reference community = geometric mean of the E2 samples (the mean in log space)
#Aitchison distance is Euclidean distance after a CLR (log-ratio) transform,
#so the "center" of a group of samples is the mean of their log values.
#Subtract the pseudocount so vegdist's own +pseudocount returns exactly the mean log value
H21_E2_mean_nitrifier_composition <- exp(rowMeans(log(as.matrix(phyloseq.samples.nitrifiers_H21_phase_E2@otu_table) + pseudocount))
) - pseudocount

#Add it to the otu table (ALL 87 samples + the E2 mean as a named row)
phyloseq.samples.nitrifiers_H21.metadata_plus_ref_E2 <- rbind(
  as.matrix(t(phyloseq.samples.nitrifiers_H21@otu_table)),
  E2_reference = H21_E2_mean_nitrifier_composition
)
dim(phyloseq.samples.nitrifiers_H21.metadata_plus_ref_E2) #OK, 88 x 1686

#Calculate overall aitchisons, then get only those to the mean E2 composition
Aitchison_nitrifiers_H21_distance <- as.matrix(
  vegdist(phyloseq.samples.nitrifiers_H21.metadata_plus_ref_E2,
          method = "aitchison",
          pseudocount = pseudocount))
Aitchison_nitrifiers_H21_distance_to_E2 <-
  Aitchison_nitrifiers_H21_distance[sample_names(phyloseq.samples.nitrifiers_H21), "E2_reference"]
Aitchison_nitrifiers_H21_distance_to_E2

#Make into dataframe
Aitchison_nitrifiers_H21_distance_to_E2_df <- data.frame(
  SampleID = names(Aitchison_nitrifiers_H21_distance_to_E2),
  AITCH_to_E2 = as.numeric(Aitchison_nitrifiers_H21_distance_to_E2)
)

#Merge with metadata
AITCH_nitrifiers_H21_distance_to_E2_metadata_df <- left_join(Aitchison_nitrifiers_H21_distance_to_E2_df, 
                                                          phyloseq.samples.nitrifiers_H21.metadata, 
                                                          by = "SampleID")
AITCH_nitrifiers_H21_distance_to_E2_metadata_df

#Omnibus test (Kruskal)
AITCH_nitrifiers_H21_distance_to_E2_metadata_df %>%
  filter(Date_num_phase_abbrv != "E2") %>%
  kruskal_test(AITCH_to_E2 ~ Date_num_phase_abbrv)
#  .y.             n statistic    df          p method        
# AITCH_to_E2    68      29.0     4 0.00000785 Kruskal-Wallis

#Wilcox Test
AITCH_nitrifiers_H21_distance_to_E2_metadata_df %>%
  pairwise_wilcox_test(
    AITCH_to_E2 ~ Date_num_phase_abbrv,
    p.adjust.method = "BH"
  ) 
# wilcox_effsize(AITCH_to_E2 ~ Date_num_phase_abbrv)

#Reorder for plot
phase_order_for_plot_h21_bc <- c("L", "T1", "E1", "T2", "T3")   # chronological, no E2

#Reorder factor on dataframe
AITCH_nitrifiers_H21_distance_to_E2_metadata_df_plot <- 
  AITCH_nitrifiers_H21_distance_to_E2_metadata_df %>%
  filter(Date_num_phase_abbrv != "E2") %>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = phase_order_for_plot_h21_bc))

#Keep only comparisons in chronological order
AITCH_nitrifiers_H21_distance_to_E2_metadata_stat_test <- 
  AITCH_nitrifiers_H21_distance_to_E2_metadata_df_plot %>%
  pairwise_wilcox_test(
    AITCH_to_E2 ~ Date_num_phase_abbrv,
    p.adjust.method = "BH",
    comparisons = list(
      c("L", "T1"),
      c("T1", "E1"),
      c("E1", "T2"),
      c("T2", "T3")
    )
  ) %>%
  add_xy_position(x = "Date_num_phase_abbrv")
AITCH_nitrifiers_H21_distance_to_E2_metadata_stat_test

#Plot
AITCH_nitrifiers_distance_to_E2_H21_boxplot <- ggplot(AITCH_nitrifiers_H21_distance_to_E2_metadata_df_plot,
                                                   aes(x = Date_num_phase_abbrv, 
                                                       y= AITCH_to_E2, 
                                                       color = Date_num_phase_abbrv, 
                                                       fill = Date_num_phase_abbrv)) +
  theme_bw() +
  labs(x = "Phase", 
       y = "Aitchison Distance",
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3, 
               outlier.shape = NA) +
  scale_fill_manual(values = naive_phase_palette)+
  scale_color_manual(values = naive_phase_palette)+
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
  stat_pvalue_manual(AITCH_nitrifiers_H21_distance_to_E2_metadata_stat_test, 
                     label = "p.adj.signif", 
                     label.size = 10,
                     hide.ns = T)
AITCH_nitrifiers_distance_to_E2_H21_boxplot

####ESTABLISHED SYSTEM######
#Do phases differ in degree of convergence to P?##

#Only established system samples
phyloseq.samples.nitrifiers_P1 #1685 taxa and 117 samples from the established system
phyloseq.samples.nitrifiers_P1.metadata <- data.frame(phyloseq.samples.nitrifiers_P1@sam_data) #metadata

#Get the geometric mean reference of nitrifiers during P
phyloseq.samples.nitrifiers_P1_phase_P <- subset_samples(phyloseq.samples.nitrifiers_P1,
                                                           Date_num_phase_abbrv == "P")
phyloseq.samples.nitrifiers_P1_phase_P #1685 taxa and 24 samples

#Need pseudocount for Aitchisons
pseudocount <- 0.5

#Reference community = geometric mean of the E2 samples (the mean in log space)
#Aitchison distance is Euclidean distance after a CLR (log-ratio) transform,
#so the "center" of a group of samples is the mean of their log values.
P1_P_geomean_nitrifier_composition <- exp(rowMeans(log(as.matrix(phyloseq.samples.nitrifiers_P1_phase_P@otu_table) + pseudocount))
) - pseudocount

#Add it to the otu table (ALL 117 samples + the P mean as a named row)
phyloseq.samples.nitrifiers_P1.metadata_plus_ref_P <- rbind(
  as.matrix(t(phyloseq.samples.nitrifiers_P1@otu_table)),
  P_reference = P1_P_geomean_nitrifier_composition
)
dim(phyloseq.samples.nitrifiers_P1.metadata_plus_ref_P) #OK, 118 x 1685


#Calculate overall aitchisons, then get only those to the mean P composition
Aitchison_nitrifiers_P1_distance <- as.matrix(
  vegdist(phyloseq.samples.nitrifiers_P1.metadata_plus_ref_P,
          method = "aitchison",
          pseudocount = pseudocount))
Aitchison_nitrifiers_P1_distance_to_P <-
  Aitchison_nitrifiers_P1_distance[sample_names(phyloseq.samples.nitrifiers_P1), "P_reference"]
Aitchison_nitrifiers_P1_distance_to_P

#Make into dataframe
Aitchison_nitrifiers_P1_distance_to_P_df <- data.frame(
  SampleID = names(Aitchison_nitrifiers_P1_distance_to_P),
  AITCH_to_P = as.numeric(Aitchison_nitrifiers_P1_distance_to_P)
)

#Merge with metadata
AITCH_nitrifiers_P1_distance_to_P_metadata_df <- left_join(Aitchison_nitrifiers_P1_distance_to_P_df, 
                                                             phyloseq.samples.nitrifiers_P1.metadata, 
                                                             by = "SampleID")
AITCH_nitrifiers_P1_distance_to_P_metadata_df

#Omnibus test (P excluded, see below)
AITCH_nitrifiers_P1_distance_to_P_metadata_df %>%
  filter(Date_num_phase_abbrv != "P") %>%
  kruskal_test(AITCH_to_P ~ Date_num_phase_abbrv) #p 0.000000216

#Wilcox Test
AITCH_nitrifiers_P1_distance_to_P_metadata_df %>%
  pairwise_wilcox_test(
    AITCH_to_P ~ Date_num_phase_abbrv,
    p.adjust.method = "BH"
  ) 
# wilcox_effsize(AITCH_to_P ~ Date_num_phase_abbrv)

#Reorder for plot
phase_order_for_plot_p1_bc <- c("L", "T1", "E")   # chronological, no P

#Reorder factor on dataframe
AITCH_nitrifiers_P1_distance_to_P_metadata_df_plot <- 
  AITCH_nitrifiers_P1_distance_to_P_metadata_df %>%
  filter(Date_num_phase_abbrv != "P") %>%
  mutate(Date_num_phase_abbrv = factor(Date_num_phase_abbrv, 
                                       levels = phase_order_for_plot_p1_bc))

#Keep only comparisons in chronological order
AITCH_nitrifiers_P1_distance_to_P_metadata_stat_test <- 
  AITCH_nitrifiers_P1_distance_to_P_metadata_df_plot %>%
  pairwise_wilcox_test(
    AITCH_to_P ~ Date_num_phase_abbrv,
    p.adjust.method = "BH",
    comparisons = list(
      c("L", "T1"),
      c("T1", "E")
    )
  ) %>%
  add_xy_position(x = "Date_num_phase_abbrv")
AITCH_nitrifiers_P1_distance_to_P_metadata_stat_test

#Plot
AITCH_nitrifiers_distance_to_P_P1_boxplot <- ggplot(AITCH_nitrifiers_P1_distance_to_P_metadata_df_plot,
                                                      aes(x = Date_num_phase_abbrv, 
                                                          y= AITCH_to_P, 
                                                          color = Date_num_phase_abbrv, 
                                                          fill = Date_num_phase_abbrv)) +
  theme_bw() +
  labs(x = "Phase", 
       y = "Aitchison Distance",
       color = "Phase", 
       fill = "Phase") +
  geom_jitter(size = 3, shape = 18, 
              alpha = 0.8, width = 0.2) +
  geom_boxplot(alpha = 0.3, 
               outlier.shape = NA) +
  scale_fill_manual(values = established_phase_palette)+
  scale_color_manual(values = established_phase_palette)+
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
  stat_pvalue_manual(AITCH_nitrifiers_P1_distance_to_P_metadata_stat_test, 
                     label = "p.adj.signif", 
                     label.size = 10,
                     hide.ns = T)
AITCH_nitrifiers_distance_to_P_P1_boxplot


#ORDINATION####
##NAIVE SYSTEM######
###BACTERIAL AND ARCHAEAL COMMUNITY######
phyloseq.samples.bacteria.H21 #218414 taxa and 87 samples
####BC DISTANCES######
#First, normalize counts (Relative abundance)
phyloseq.samples.bacteria.H21.ra <- transform_sample_counts(phyloseq.samples.bacteria.H21, 
                                                           function(x) x/sum(x)*100)
sample_sums(phyloseq.samples.bacteria.H21.ra) #OK

#Calculate BC distances
phyloseq.samples.bacteria.H21.ra.bray <- vegdist(t(phyloseq.samples.bacteria.H21.ra@otu_table), method = "bray")
phyloseq.samples.bacteria.H21.ra.bray

#Ordination
set.seed(87)
phyloseq.samples.bacteria.H21.ra.bray.ord <- metaMDS(phyloseq.samples.bacteria.H21.ra.bray, 
                                                 k=2, 
                                                 try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)

#For metadata, will use the imputed values (ph, Copper, Salinity, Temperature)
phyloseq.samples.H21_metadata <- phyloseq.samples_metadata%>%
  filter(Enclosure == "H21")

#####Centroids by Phase#####
#Simple ordination plot
phyloseq.samples.bacteria.H21.ra.bray.ord_plot <- ordiplot(phyloseq.samples.bacteria.H21.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.samples.bacteria.H21.ra.bray.ord_plot_scrs <- scores(phyloseq.samples.bacteria.H21.ra.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.samples.bacteria.H21.ra.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.samples.bacteria.H21.ra.bray.ord_plot_scrs),
                                                Copper_level_mg_L = phyloseq.samples.H21_metadata$Copper_level_mg_L, 
                                                SampleID = phyloseq.samples.H21_metadata$SampleID, 
                                                Date_num = phyloseq.samples.H21_metadata$Date_num,
                                                Collection_Date = phyloseq.samples.H21_metadata$Collection_Date, 
                                                Ammonia_level = phyloseq.samples.H21_metadata$Ammonia_mg_L, 
                                                Date_num_phase = phyloseq.samples.H21_metadata$Date_num_phase, 
                                                Date_num_phase_abbrv = phyloseq.samples.H21_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.samples.bacteria.H21.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Date_num_phase_abbrv, 
                                                     data = phyloseq.samples.bacteria.H21.ra.bray.ord_plot_scrs, 
                                                     FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.samples.bacteria.H21.ra.bray.ord_plot_segs <- merge(phyloseq.samples.bacteria.H21.ra.bray.ord_plot_scrs, 
                             setNames(phyloseq.samples.bacteria.H21.ra.dates.bray.cent.phase, 
                                      c("Date_num_phase_abbrv","cMDS1","cMDS2")),
                             by = 'Date_num_phase_abbrv', 
                             sort = F)

######PERMANOVA - Just phase####
phyloseq.samples.H21_metadata$Date_num_phase_established <- 
  factor(phyloseq.samples.H21_metadata$Date_num_phase_abbrv, levels = c("E2", "L", "T1", 'E1', "T2", 'T3'))
set.seed(98)
phase_H21_BC_adonis  <- adonis2(phyloseq.samples.bacteria.H21.ra.bray ~ 
                                  Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                phyloseq.samples.H21_metadata, 
                                #strata = phyloseq.samples.H21_metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_H21_BC_adonis #11.6% of the variation in community composition is explained by phase
#                      Df SumOfSqs      R2      F Pr(>F)   
# Date_num_phase_abbrv  5   1.2934 0.11618 2.6890 0.0010 ***
# pH_spu_imp            1   0.3083 0.02770 3.2052 0.0223 *  
# Temperature_F_imp     1   0.1898 0.01705 1.9734 0.0919 .  
# Salinity_ppt_imp      1   0.2927 0.02629 3.0422 0.0274 *  

#PERMANOVA DF
phase_H21_BC_adonis_df <- data.frame(phase_H21_BC_adonis, 
                                     check.names = F, check.rows = F)%>%
  rownames_to_column(var = "Fixed Effect")%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))%>%
  mutate(Predictor = case_when(
    `Fixed Effect` == "Date_num_phase_abbrv" ~ "Phase", 
    `Fixed Effect` == "pH_spu_imp" ~ "pH", 
    `Fixed Effect` == "Temperature_F_imp" ~ "Temperature", 
    `Fixed Effect` == "Salinity_ppt_imp" ~ "Salinity"
  ))%>%
  mutate(Comparison = "Global", 
         System = "Naive", 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3),
         `R2` = round(`R2`, 3), 
         `F` = round(`F`, 2), 
         SumOfSqs = round(SumOfSqs, 2)
  )%>%
  select(System, Comparison, Predictor, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_H21_BC_adonis_df

#Pairwise permanova 
set.seed(98)
phase_H21_BC_adonis_pairwise  <- pairwise.adonis2(phyloseq.samples.bacteria.H21.ra.bray ~ 
                                                    Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                                  phyloseq.samples.H21_metadata,
                                                  by = "margin",
                                                  nperm = 9999)
phase_H21_BC_adonis_pairwise 


#Comparisons against reference E2
phase_H21_BC_adonis_pairwise_E2 <- bind_rows(
  data.frame(phase_H21_BC_adonis_pairwise$L_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs L"), 
  data.frame(phase_H21_BC_adonis_pairwise$T1_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs T1"), 
  data.frame(phase_H21_BC_adonis_pairwise$E1_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs E1"), 
  data.frame(phase_H21_BC_adonis_pairwise$T2_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs T2"), 
  data.frame(phase_H21_BC_adonis_pairwise$T3_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs T3")
)
phase_H21_BC_adonis_pairwise_E2 <- phase_H21_BC_adonis_pairwise_E2%>%
  mutate(SumOfSqs = round(SumOfSqs, 2), 
         R2 = round(R2, 3), 
         `F` = round(`F`, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3), 
         System = "Naive"
  )%>%
  select(System,Comparison, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_H21_BC_adonis_pairwise_E2 #All phases compared against E2 are different in compostion 

##PERMDISPS
# Run the betadisper function, average distance to centroid
bray.h21.phase.disp <- betadisper(phyloseq.samples.bacteria.H21.ra.bray, 
                                  phyloseq.samples.H21_metadata$Date_num_phase_abbrv)
bray.h21.phase.disp 

##Then test by permuting
set.seed(98)
bray.h21.phase.permdisp <- permutest(bray.h21.phase.disp, 
                                     pairwise = TRUE,
                                     permutations = 9999)
bray.h21.phase.permdisp ##S, P = 6e-04

#Pairwise permuted results:
bray.h21.phase.permdisp.pairwise.permuted <- bray.h21.phase.permdisp$pairwise$permuted
bray.h21.phase.permdisp.pairwise.permuted

#Make into dataframe
H21_permdisp_pairwise_results_table <- 
  enframe(bray.h21.phase.permdisp.pairwise.permuted,
          name = "Comparison", value = "p") %>%
  mutate(Comparison = case_when(
    Comparison == "L-E2" ~ "E2 vs L",
    Comparison == "T1-E2" ~ "E2 vs T1",
    Comparison == "E1-E2" ~ "E2 vs E1",
    Comparison == "T2-E2" ~ "E2 vs T2",
    Comparison == "T3-E2" ~ "E2 vs T3",
  ))%>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%   # all 15 pairs, adjusted with BH
  filter(grepl("E2", Comparison)) %>% # then keep E2 only
  mutate(P = format(p, scientific = TRUE, digits = 3),
         `P adjusted` = format(p_adj, scientific = TRUE, digits = 3)) %>%
  select(Comparison, P, `P adjusted`)
H21_permdisp_pairwise_results_table

# Extract R2 and p-values (for phase effect)
R2_adonis_h21_phase <- phase_H21_BC_adonis$R2[1] 
pvalue_adonis_h21_phase <- phase_H21_BC_adonis$`Pr(>F)`[1]

#PLOT
H21_phase_BC_beta_div <- ggplot(phyloseq.samples.bacteria.H21.ra.bray.ord_plot_segs) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2",
       color = "Naive System Phases") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = Date_num_phase_abbrv),
               show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10, show.legend = F) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = naive_phase_palette)+
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
  annotate("text", x = -0.5, y = -0.7,
           label = "PERMANOVA\nPhase",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = -0.5, y = -0.7,
           label = paste("R² = ", round(R2_adonis_h21_phase* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_h21_phase, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
H21_phase_BC_beta_div

###NITRIFYING COMMUNITY######
phyloseq.samples.nitrifiers_H21 #1686 taxa and 87 samples
####BC DISTANCES######
#First, normalize counts (Relative abundance)
phyloseq.samples.nitrifiers_H21.ra <- transform_sample_counts(phyloseq.samples.nitrifiers_H21, 
                                                            function(x) x/sum(x)*100)
sample_sums(phyloseq.samples.nitrifiers_H21.ra) #OK

#Calculate BC distances
phyloseq.samples.nitrifiers_H21.ra.bray <- vegdist(t(phyloseq.samples.nitrifiers_H21.ra@otu_table), method = "bray")
phyloseq.samples.nitrifiers_H21.ra.bray

#Ordination
set.seed(87)
phyloseq.samples.nitrifiers_H21.ra.bray.ord <- metaMDS(phyloseq.samples.nitrifiers_H21.ra.bray, 
                                                     k=2, 
                                                     try = 50, 
                                                     trymax = 1000,
                                                     autotransform = F)

#For metadata, will use the imputed values (ph, Copper, Salinity, Temperature)
phyloseq.samples.H21_metadata <- phyloseq.samples_metadata%>%
  filter(Enclosure == "H21")

#####Centroids by Phase#####
#Simple ordination plot
phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot <- ordiplot(phyloseq.samples.nitrifiers_H21.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_scrs <- scores(phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_scrs),
                                                             Copper_level_mg_L = phyloseq.samples.H21_metadata$Copper_level_mg_L, 
                                                             SampleID = phyloseq.samples.H21_metadata$SampleID, 
                                                             Date_num = phyloseq.samples.H21_metadata$Date_num,
                                                             Collection_Date = phyloseq.samples.H21_metadata$Collection_Date, 
                                                             Ammonia_level = phyloseq.samples.H21_metadata$Ammonia_mg_L, 
                                                             Date_num_phase = phyloseq.samples.H21_metadata$Date_num_phase, 
                                                             Date_num_phase_abbrv = phyloseq.samples.H21_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.samples.nitrifiers_H21.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Date_num_phase_abbrv, 
                                                                    data = phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_scrs, 
                                                                    FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_segs <- merge(phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_scrs, 
                                                             setNames(phyloseq.samples.nitrifiers_H21.ra.dates.bray.cent.phase, 
                                                                      c("Date_num_phase_abbrv","cMDS1","cMDS2")),
                                                             by = 'Date_num_phase_abbrv', 
                                                             sort = F)

######PERMANOVA - Just phase####
phyloseq.samples.H21_metadata$Date_num_phase_established <- 
  factor(phyloseq.samples.H21_metadata$Date_num_phase_abbrv, levels = c("E2", "L", "T1", 'E1', "T2", 'T3'))
set.seed(98)
phase_H21_nitrifiers_BC_adonis  <- adonis2(phyloseq.samples.nitrifiers_H21.ra.bray ~ 
                                  Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                phyloseq.samples.H21_metadata, 
                                #strata = phyloseq.samples.H21_metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_H21_nitrifiers_BC_adonis 
#22.1% of the variation in community composition is explained by phase
#                      Df SumOfSqs      R2      F Pr(>F)   
# Date_num_phase_abbrv  5   2.7451 0.22117 6.0054 0.0001 ***
# pH_spu_imp            1   0.4890 0.03940 5.3487 0.0039 ** 
# Temperature_F_imp     1   0.2254 0.01816 2.4652 0.0637 .  
# Salinity_ppt_imp      1   0.2811 0.02265 3.0750 0.0330 *   

#PERMANOVA DF
phase_H21_nitrifiers_BC_adonis_df <- data.frame(phase_H21_nitrifiers_BC_adonis, 
                                     check.names = F, check.rows = F)%>%
  rownames_to_column(var = "Fixed Effect")%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))%>%
  mutate(Predictor = case_when(
    `Fixed Effect` == "Date_num_phase_abbrv" ~ "Phase", 
    `Fixed Effect` == "pH_spu_imp" ~ "pH", 
    `Fixed Effect` == "Temperature_F_imp" ~ "Temperature", 
    `Fixed Effect` == "Salinity_ppt_imp" ~ "Salinity"
  ))%>%
  mutate(Comparison = "Global", 
         System = "Naive", 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3),
         `R2` = round(`R2`, 3), 
         `F` = round(`F`, 2), 
         SumOfSqs = round(SumOfSqs, 2)
  )%>%
  select(System, Comparison, Predictor, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_H21_nitrifiers_BC_adonis_df

#Pairwise permanova 
set.seed(98)
phase_H21_nitrifiers_BC_adonis_pairwise  <- pairwise.adonis2(phyloseq.samples.nitrifiers_H21.ra.bray ~ 
                                                    Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                                  phyloseq.samples.H21_metadata,
                                                  by = "margin",
                                                  nperm = 9999)
phase_H21_nitrifiers_BC_adonis_pairwise 

#Comparisons against reference E2
phase_H21_nitrifiers_BC_adonis_pairwise_E2 <- bind_rows(
  data.frame(phase_H21_nitrifiers_BC_adonis_pairwise$L_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs L"), 
  data.frame(phase_H21_nitrifiers_BC_adonis_pairwise$T1_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs T1"), 
  data.frame(phase_H21_nitrifiers_BC_adonis_pairwise$E1_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs E1"), 
  data.frame(phase_H21_nitrifiers_BC_adonis_pairwise$T2_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs T2"), 
  data.frame(phase_H21_nitrifiers_BC_adonis_pairwise$T3_vs_E2, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "E2 vs T3")
)
phase_H21_nitrifiers_BC_adonis_pairwise_E2
phase_H21_nitrifiers_BC_adonis_pairwise_E2 <- phase_H21_nitrifiers_BC_adonis_pairwise_E2%>%
  mutate(SumOfSqs = round(SumOfSqs, 2), 
         R2 = round(R2, 3), 
         `F` = round(`F`, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3), 
         System = "Naive"
  )%>%
  select(System,Comparison, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_H21_nitrifiers_BC_adonis_pairwise_E2 #All phases compared against E2 are different in compostion 

##PERMDISPS
# Run the betadisper function, average distance to centroid
bray.h21.phase.nitrifiers.disp <- betadisper(phyloseq.samples.nitrifiers_H21.ra.bray, 
                                  phyloseq.samples.H21_metadata$Date_num_phase_abbrv)
bray.h21.phase.nitrifiers.disp 

##Then test by permuting
set.seed(98)
bray.h21.phase.nitrifiers.permdisp <- permutest(bray.h21.phase.nitrifiers.disp, 
                                     pairwise = TRUE,
                                     permutations = 9999)
bray.h21.phase.nitrifiers.permdisp 

#Pairwise permuted results:
bray.h21.phase.nitrifiers.permdisp.pairwise.permuted <- bray.h21.phase.nitrifiers.permdisp$pairwise$permuted
bray.h21.phase.nitrifiers.permdisp.pairwise.permuted

#Make into dataframe
H21_nitrifiers_permdisp_pairwise_results_table <- 
  enframe(bray.h21.phase.nitrifiers.permdisp.pairwise.permuted,
          name = "Comparison", value = "p") %>%
  mutate(Comparison = case_when(
    Comparison == "L-E2" ~ "E2 vs L",
    Comparison == "T1-E2" ~ "E2 vs T1",
    Comparison == "E1-E2" ~ "E2 vs E1",
    Comparison == "T2-E2" ~ "E2 vs T2",
    Comparison == "T3-E2" ~ "E2 vs T3",
  ))%>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%   # all 15 pairs, adjusted with BH
  filter(grepl("E2", Comparison)) %>% # then keep E2 only
  mutate(P = format(p, scientific = TRUE, digits = 3),
         `P adjusted` = format(p_adj, scientific = TRUE, digits = 3)) %>%
  select(Comparison, P, `P adjusted`)
H21_nitrifiers_permdisp_pairwise_results_table #Significant for all
# Comparison P       `P adjusted`
# E2 vs L    1.0e-04 3.00e-04    
# E2 vs T1   1.0e-04 3.00e-04    
# E2 vs E1   2.4e-03 4.50e-03    
# E2 vs T2   2.0e-03 4.29e-03    
# E2 vs T3   1.0e-04 3.00e-04 

# Extract R2 and p-values (for phase effect)
R2_adonis_h21_phase.nitrifiers <- phase_H21_nitrifiers_BC_adonis$R2[1] 
pvalue_adonis_h21_phase.nitrifiers <- phase_H21_nitrifiers_BC_adonis$`Pr(>F)`[1]

#PLOT
H21_phase_nitrifiers_BC_beta_div <- ggplot(phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_segs) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2",
       color = "Naive System Phases") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = Date_num_phase_abbrv),
               show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10, show.legend = F) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = naive_phase_palette)+
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
           label = paste("R² = ", round(R2_adonis_h21_phase.nitrifiers* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_h21_phase.nitrifiers, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
H21_phase_nitrifiers_BC_beta_div

##ESTABLISHED SYSTEM######
###BACTERIAL AND ARCHAEAL COMMUNITY######
phyloseq.samples.bacteria.P1 #217728 taxa and 87 samples
####BC DISTANCES######
#First, normalize counts (Relative abundance)
phyloseq.samples.bacteria.P1.ra <- transform_sample_counts(phyloseq.samples.bacteria.P1, 
                                                            function(x) x/sum(x)*100)
sample_sums(phyloseq.samples.bacteria.P1.ra) #OK

#Calculate BC distances
phyloseq.samples.bacteria.P1.ra.bray <- vegdist(t(phyloseq.samples.bacteria.P1.ra@otu_table), method = "bray")
phyloseq.samples.bacteria.P1.ra.bray

#Ordination
set.seed(87)
phyloseq.samples.bacteria.P1.ra.bray.ord <- metaMDS(phyloseq.samples.bacteria.P1.ra.bray, 
                                                     k=2, 
                                                     try = 50, 
                                                     trymax = 1000,
                                                     autotransform = F)

#For metadata, will use the imputed values (ph, Copper, Salinity, Temperature)
phyloseq.samples.P1_metadata <- phyloseq.samples_metadata%>%
  filter(Enclosure == "P1")

#####Centroids by Phase#####
#Simple ordination plot
phyloseq.samples.bacteria.P1.ra.bray.ord_plot <- ordiplot(phyloseq.samples.bacteria.P1.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.samples.bacteria.P1.ra.bray.ord_plot_scrs <- scores(phyloseq.samples.bacteria.P1.ra.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.samples.bacteria.P1.ra.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.samples.bacteria.P1.ra.bray.ord_plot_scrs),
                                                             Copper_level_mg_L = phyloseq.samples.P1_metadata$Copper_level_mg_L, 
                                                             SampleID = phyloseq.samples.P1_metadata$SampleID, 
                                                             Date_num = phyloseq.samples.P1_metadata$Date_num,
                                                             Collection_Date = phyloseq.samples.P1_metadata$Collection_Date, 
                                                             Ammonia_level = phyloseq.samples.P1_metadata$Ammonia_mg_L, 
                                                             Date_num_phase = phyloseq.samples.P1_metadata$Date_num_phase, 
                                                             Date_num_phase_abbrv = phyloseq.samples.P1_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.samples.bacteria.P1.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Date_num_phase_abbrv, 
                                                                    data = phyloseq.samples.bacteria.P1.ra.bray.ord_plot_scrs, 
                                                                    FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.samples.bacteria.P1.ra.bray.ord_plot_segs <- merge(phyloseq.samples.bacteria.P1.ra.bray.ord_plot_scrs, 
                                                             setNames(phyloseq.samples.bacteria.P1.ra.dates.bray.cent.phase, 
                                                                      c("Date_num_phase_abbrv","cMDS1","cMDS2")),
                                                             by = 'Date_num_phase_abbrv', 
                                                             sort = F)

######PERMANOVA - Just phase####
phyloseq.samples.P1_metadata$Date_num_phase_established <- 
  factor(phyloseq.samples.P1_metadata$Date_num_phase_abbrv, levels = c("P", "L", "T1", 'E'))
set.seed(98)
phase_P1_BC_adonis  <- adonis2(phyloseq.samples.bacteria.P1.ra.bray ~ 
                                  Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                phyloseq.samples.P1_metadata, 
                                #strata = phyloseq.samples.P1_metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_P1_BC_adonis 
#                        Df SumOfSqs      R2      F Pr(>F)   
# Date_num_phase_abbrv   3   1.3859 0.15499 7.9629 0.0002 ***
# pH_spu_imp             1   0.3042 0.03402 5.2434 0.0149 *  
# Temperature_F_imp      1   0.1730 0.01934 2.9812 0.0652 .  
# Salinity_ppt_imp       1   0.0786 0.00879 1.3549 0.2259    

#PERMANOVA DF
phase_P1_BC_adonis_df <- data.frame(phase_P1_BC_adonis, 
                                     check.names = F, check.rows = F)%>%
  rownames_to_column(var = "Fixed Effect")%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))%>%
  mutate(Predictor = case_when(
    `Fixed Effect` == "Date_num_phase_abbrv" ~ "Phase", 
    `Fixed Effect` == "pH_spu_imp" ~ "pH", 
    `Fixed Effect` == "Temperature_F_imp" ~ "Temperature", 
    `Fixed Effect` == "Salinity_ppt_imp" ~ "Salinity"
  ))%>%
  mutate(Comparison = "Global", 
         System = "Established", 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3),
         `R2` = round(`R2`, 3), 
         `F` = round(`F`, 2), 
         SumOfSqs = round(SumOfSqs, 2)
  )%>%
  select(System, Comparison, Predictor, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_P1_BC_adonis_df

#Pairwise permanova 
set.seed(98)
phase_P1_BC_adonis_pairwise  <- pairwise.adonis2(phyloseq.samples.bacteria.P1.ra.bray ~ 
                                                    Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                                  phyloseq.samples.P1_metadata,
                                                  by = "margin",
                                                  nperm = 9999)
phase_P1_BC_adonis_pairwise 


#Comparisons against reference P
phase_P1_BC_adonis_pairwise_P <- bind_rows(
  data.frame(phase_P1_BC_adonis_pairwise$L_vs_P, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "P vs L"), 
  data.frame(phase_P1_BC_adonis_pairwise$T1_vs_P, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "P vs T1"), 
  data.frame(phase_P1_BC_adonis_pairwise$E_vs_P, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "P vs E")
)
phase_P1_BC_adonis_pairwise_P <- phase_P1_BC_adonis_pairwise_P%>%
  mutate(SumOfSqs = round(SumOfSqs, 2), 
         R2 = round(R2, 3), 
         `F` = round(`F`, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3), 
         System = "Established"
  )%>%
  select(System,Comparison, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_P1_BC_adonis_pairwise_P #All phases compared against P are different in compostion 

##PERMDISPS
# Run the betadisper function, average distance to centroid
bray.p1.phase.disp <- betadisper(phyloseq.samples.bacteria.P1.ra.bray, 
                                  phyloseq.samples.P1_metadata$Date_num_phase_abbrv)
bray.p1.phase.disp 

##Then test by permuting
set.seed(98)
bray.p1.phase.permdisp <- permutest(bray.p1.phase.disp, 
                                     pairwise = TRUE,
                                     permutations = 9999)
bray.p1.phase.permdisp 

#Pairwise permuted results:
bray.p1.phase.permdisp.pairwise.permuted <- bray.p1.phase.permdisp$pairwise$permuted
bray.p1.phase.permdisp.pairwise.permuted

#Make into dataframe
P1_permdisp_pairwise_results_table <- 
  enframe(bray.p1.phase.permdisp.pairwise.permuted,
          name = "Comparison", value = "p") %>%
  mutate(Comparison = case_when(
    Comparison == "L-P" ~ "P vs L",
    Comparison == "T1-P" ~ "P vs T1",
    Comparison == "E-P" ~ "P vs E"
  ))%>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%   # all 15 pairs, adjusted with BH
  filter(grepl("P", Comparison)) %>% # then keep P only
  mutate(P = format(p, scientific = TRUE, digits = 3),
         `P adjusted` = format(p_adj, scientific = TRUE, digits = 3)) %>%
  select(Comparison, P, `P adjusted`)
P1_permdisp_pairwise_results_table

# Extract R2 and p-values (for phase effect)
R2_adonis_p1_phase <- phase_P1_BC_adonis$R2[1] 
pvalue_adonis_p1_phase <- phase_P1_BC_adonis$`Pr(>F)`[1]

#PLOT
P1_phase_BC_beta_div <- ggplot(phyloseq.samples.bacteria.P1.ra.bray.ord_plot_segs) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2",
       color = "Established System Phases") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = Date_num_phase_abbrv),
               show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10, show.legend = F) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = established_phase_palette)+
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
  annotate("text", x = 0.5, y = -0.5,
           label = "PERMANOVA\nPhase",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = 0.5, y = -0.45,
           label = paste("R² = ", round(R2_adonis_p1_phase* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_p1_phase, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
P1_phase_BC_beta_div

###NITRIFYING COMMUNITY######
phyloseq.samples.nitrifiers_P1 #1685 taxa and 117 samples
####BC DISTANCES######
#First, normalize counts (Relative abundance)
phyloseq.samples.nitrifiers_P1.ra <- transform_sample_counts(phyloseq.samples.nitrifiers_P1, 
                                                              function(x) x/sum(x)*100)
sample_sums(phyloseq.samples.nitrifiers_P1.ra) #OK

#Calculate BC distances
phyloseq.samples.nitrifiers_P1.ra.bray <- vegdist(t(phyloseq.samples.nitrifiers_P1.ra@otu_table), method = "bray")
phyloseq.samples.nitrifiers_P1.ra.bray

#Ordination
set.seed(87)
phyloseq.samples.nitrifiers_P1.ra.bray.ord <- metaMDS(phyloseq.samples.nitrifiers_P1.ra.bray, 
                                                       k=2, 
                                                       try = 50, 
                                                       trymax = 1000,
                                                       autotransform = F)

#For metadata, will use the imputed values (ph, Copper, Salinity, Temperature)
phyloseq.samples.P1_metadata <- phyloseq.samples_metadata%>%
  filter(Enclosure == "P1")

#####Centroids by Phase#####
#Simple ordination plot
phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot <- ordiplot(phyloseq.samples.nitrifiers_P1.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_scrs <- scores(phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_scrs),
                                                               Copper_level_mg_L = phyloseq.samples.P1_metadata$Copper_level_mg_L, 
                                                               SampleID = phyloseq.samples.P1_metadata$SampleID, 
                                                               Date_num = phyloseq.samples.P1_metadata$Date_num,
                                                               Collection_Date = phyloseq.samples.P1_metadata$Collection_Date, 
                                                               Ammonia_level = phyloseq.samples.P1_metadata$Ammonia_mg_L, 
                                                               Date_num_phase = phyloseq.samples.P1_metadata$Date_num_phase, 
                                                               Date_num_phase_abbrv = phyloseq.samples.P1_metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.samples.nitrifiers_P1.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Date_num_phase_abbrv, 
                                                                      data = phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_scrs, 
                                                                      FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_segs <- merge(phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_scrs, 
                                                               setNames(phyloseq.samples.nitrifiers_P1.ra.dates.bray.cent.phase, 
                                                                        c("Date_num_phase_abbrv","cMDS1","cMDS2")),
                                                               by = 'Date_num_phase_abbrv', 
                                                               sort = F)

######PERMANOVA - Just phase####
phyloseq.samples.P1_metadata$Date_num_phase_established <- 
  factor(phyloseq.samples.P1_metadata$Date_num_phase_abbrv, levels = c("P", "L", "T1", 'E', "T2", 'T3'))
set.seed(98)
phase_P1_nitrifiers_BC_adonis  <- adonis2(phyloseq.samples.nitrifiers_P1.ra.bray ~ 
                                             Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                           phyloseq.samples.P1_metadata, 
                                           #strata = phyloseq.samples.P1_metadata$Collection_Month,
                                           by = "margin",
                                           permutations = 9999)
phase_P1_nitrifiers_BC_adonis 
# Df                        SumOfSqs   R2       F    Pr(>F)    
# Date_num_phase_abbrv   3   3.9628 0.41037 39.7139 0.0001 ***
# pH_spu_imp             1   0.3358 0.03477 10.0946 0.0001 ***
# Temperature_F_imp      1   0.4312 0.04466 12.9648 0.0001 ***
# Salinity_ppt_imp       1   0.0886 0.00918  2.6646 0.0543 .  

#PERMANOVA DF
phase_P1_nitrifiers_BC_adonis_df <- data.frame(phase_P1_nitrifiers_BC_adonis, 
                                                check.names = F, check.rows = F)%>%
  rownames_to_column(var = "Fixed Effect")%>%
  filter(!`Fixed Effect` %in% c("Residual", "Total"))%>%
  mutate(Predictor = case_when(
    `Fixed Effect` == "Date_num_phase_abbrv" ~ "Phase", 
    `Fixed Effect` == "pH_spu_imp" ~ "pH", 
    `Fixed Effect` == "Temperature_F_imp" ~ "Temperature", 
    `Fixed Effect` == "Salinity_ppt_imp" ~ "Salinity"
  ))%>%
  mutate(Comparison = "Global", 
         System = "Established", 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3),
         `R2` = round(`R2`, 3), 
         `F` = round(`F`, 2), 
         SumOfSqs = round(SumOfSqs, 2)
  )%>%
  select(System, Comparison, Predictor, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_P1_nitrifiers_BC_adonis_df

#Pairwise permanova 
set.seed(98)
phase_P1_nitrifiers_BC_adonis_pairwise  <- pairwise.adonis2(phyloseq.samples.nitrifiers_P1.ra.bray ~ 
                                                               Date_num_phase_abbrv + pH_spu_imp + Temperature_F_imp + Salinity_ppt_imp, 
                                                             phyloseq.samples.P1_metadata,
                                                             by = "margin",
                                                             nperm = 9999)
phase_P1_nitrifiers_BC_adonis_pairwise 

#Comparisons against reference P
phase_P1_nitrifiers_BC_adonis_pairwise_P <- bind_rows(
  data.frame(phase_P1_nitrifiers_BC_adonis_pairwise$L_vs_P, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "P vs L"), 
  data.frame(phase_P1_nitrifiers_BC_adonis_pairwise$T1_vs_P, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "P vs T1"), 
  data.frame(phase_P1_nitrifiers_BC_adonis_pairwise$E_vs_P, check.names = F, check.rows = F)%>%
    rownames_to_column(var = "Fixed Effect")%>%
    filter(`Fixed Effect` == "Date_num_phase_abbrv")%>%
    mutate(Comparison = "P vs E")
)
phase_P1_nitrifiers_BC_adonis_pairwise_P
phase_P1_nitrifiers_BC_adonis_pairwise_P <- phase_P1_nitrifiers_BC_adonis_pairwise_P%>%
  mutate(SumOfSqs = round(SumOfSqs, 2), 
         R2 = round(R2, 3), 
         `F` = round(`F`, 2), 
         `Pr(>F)` = format(`Pr(>F)`, scientific = TRUE, digits = 3), 
         System = "Established"
  )%>%
  select(System,Comparison, Df, SumOfSqs, `F`, R2, `Pr(>F)`)
phase_P1_nitrifiers_BC_adonis_pairwise_P #All phases compared against P are different in compostion 

##PERMDISPS
# Run the betadisper function, average distance to centroid
bray.p1.phase.nitrifiers.disp <- betadisper(phyloseq.samples.nitrifiers_P1.ra.bray, 
                                             phyloseq.samples.P1_metadata$Date_num_phase_abbrv)
bray.p1.phase.nitrifiers.disp 

##Then test by permuting
set.seed(98)
bray.p1.phase.nitrifiers.permdisp <- permutest(bray.p1.phase.nitrifiers.disp, 
                                                pairwise = TRUE,
                                                permutations = 9999)
bray.p1.phase.nitrifiers.permdisp 

#Pairwise permuted results:
bray.p1.phase.nitrifiers.permdisp.pairwise.permuted <- bray.p1.phase.nitrifiers.permdisp$pairwise$permuted
bray.p1.phase.nitrifiers.permdisp.pairwise.permuted

#Make into dataframe
P1_nitrifiers_permdisp_pairwise_results_table <- 
  enframe(bray.p1.phase.nitrifiers.permdisp.pairwise.permuted,
          name = "Comparison", value = "p") %>%
  mutate(Comparison = case_when(
    Comparison == "L-P" ~ "P vs L",
    Comparison == "T1-P" ~ "P vs T1",
    Comparison == "E-P" ~ "P vs E"
  ))%>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%   # all 15 pairs, adjusted with BH
  filter(grepl("P", Comparison)) %>% # then keep P only
  mutate(P = format(p, scientific = TRUE, digits = 3),
         `P adjusted` = format(p_adj, scientific = TRUE, digits = 3)) %>%
  select(Comparison, P, `P adjusted`)
P1_nitrifiers_permdisp_pairwise_results_table #Significant for all
# Comparison P        `P adjusted`
# P vs L     1.00e-04 2.00e-04    
# P vs T1    3.88e-01 3.88e-01    
# P vs E     1.00e-04 2.00e-04 

# Extract R2 and p-values (for phase effect)
R2_adonis_p1_phase.nitrifiers <- phase_P1_nitrifiers_BC_adonis$R2[1] 
pvalue_adonis_p1_phase.nitrifiers <- phase_P1_nitrifiers_BC_adonis$`Pr(>F)`[1]

#PLOT
P1_phase_nitrifiers_BC_beta_div <- ggplot(phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_segs) +
  theme_bw() +
  labs(x="NMDS1", y="NMDS2",
       color = "Established System Phases") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = Date_num_phase_abbrv),
               show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv), size = 10, show.legend = F) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = established_phase_palette)+
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
  annotate("text", x = 0.5, y = 0.42,
           label = "PERMANOVA\nPhase",
           hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  annotate("text", x = 0.5, y = 0.5,
           label = paste("R² = ", round(R2_adonis_p1_phase.nitrifiers* 100, 1), "%",
                         "\np = ", round(pvalue_adonis_p1_phase.nitrifiers, 4)),
           hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
P1_phase_nitrifiers_BC_beta_div

##FACET BY SYSTEM - NITRIFYING COMMUNITIES#######
phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_segs #Established system 
phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_segs #naive system

#Join the naive and established df
phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs <- 
  bind_rows(phyloseq.samples.nitrifiers_P1.ra.bray.ord_plot_segs %>% 
              mutate(System = "Established", 
                     Date_num_phase_abbrv_plot = paste0(Date_num_phase_abbrv, "_", System)), 
            phyloseq.samples.nitrifiers_H21.ra.bray.ord_plot_segs %>% 
              mutate(System = "Naive", 
                     Date_num_phase_abbrv_plot = paste0(Date_num_phase_abbrv, "_", System)))
phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs
phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs$System <- 
  factor(phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs$System, 
         levels = c("Naive", "Established"))

#PLOT
P1_H21_nitrifiers_phase_BC_beta_div <- ggplot(phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs) + 
  theme_bw() +
  labs(x="NMDS1", y="NMDS2", 
       title= "NITRIFYING COMMUNITY") +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  facet_wrap(~System, scales = "free")+
  geom_segment(aes(x=cMDS1, y=cMDS2,
                   xend= MDS1, yend = MDS2,
                   color = Date_num_phase_abbrv_plot),
               show.legend = F)+
  geom_point(aes(x=MDS1, y=MDS2, colour = Date_num_phase_abbrv_plot), size = 5, alpha = 0.8) + # individuals
  geom_point(aes(x=cMDS1, y= cMDS2, colour = Date_num_phase_abbrv_plot), size = 10, show.legend = F) + # centroids +
  geom_text(aes (x= cMDS1, y = cMDS2,
                 label= Date_num_phase_abbrv), colour= "white", size = 6, fontface = "bold") +
  scale_color_manual(values = phases_naive_established_palette)+
  theme(legend.position = "none",
        legend.title = element_text(colour = "black", 
                                    size = 22,
                                    face = "bold"),
        legend.text = element_text(colour = "black", size = 22),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(color = "white", size = 32, face = "bold"),
        plot.margin = margin(t = 1, r = 1, b = 1, l = 1),  # top, right, bottom, left
        plot.title = element_text(size = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 26, colour = "black"))+
  guides(
    shape = guide_legend(override.aes = list(size = 7)),
    color = guide_legend(nrow = 1)) +
  #Established system annotations of Permanova results
  geom_text(data = phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs%>%filter(System == "Naive"),
            aes(x = -0.8, y = 0.38, label = "PERMANOVA\nPhase"), 
            hjust = 0.5, vjust = -0.5, 
            size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  geom_text(data = phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs%>%filter(System == "Naive"),
            aes(x = -0.8, y = 0.40,
                label = paste("R² = ", round(R2_adonis_h21_phase.nitrifiers* 100, 1), "%",
                              "\np = ", round(R2_adonis_h21_phase.nitrifiers, 4))
            ),
            hjust = 0.5, vjust = 1.1, size = 8, 
            colour = "black")+# Annotate R² and p-values
  #Established system annotations of Permanova results
  geom_text(data = phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs%>%filter(System == "Established"), 
            aes(x = 0.35, y = 0.6,label = "PERMANOVA\nPhase"),
            hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
  geom_text(data = phyloseq.samples.nitrifiers_P1_H21.ra.bray.ord_plot_segs%>%filter(System == "Established"), 
            aes(x = 0.35, y = 0.65,
                label = paste("R² = ", round(R2_adonis_h21_phase.nitrifiers* 100, 1), "%",
                              "\np = ", round(R2_adonis_h21_phase.nitrifiers, 4))),
            hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
P1_H21_nitrifiers_phase_BC_beta_div


####FIGURE 4####
P1_H21_nitrifiers_phase_BC_beta_div #beta diversity, colored by phase
BC_nitrifiers_distance_to_E2_H21_boxplot #distance to E2 mean community structure for naive system
BC_nitrifiers_distance_to_P_P1.boxplot #distance to P mean community structure for established system

figure_phase_beta_div_distances_H21_P1 <- cowplot::plot_grid(
  P1_H21_nitrifiers_phase_BC_beta_div,
  cowplot::plot_grid(
    BC_nitrifiers_distance_to_E2_H21.boxplot,
    BC_nitrifiers_distance_to_P_P1.boxplot,
    ncol = 2
  ),
  labels = "AUTO", 
  label_size = 40,
  ncol = 1,
  rel_heights = c(1.05, 0.95)
)

figure4 <- figure_phase_beta_div_distances_H21_P1
ggsave("Figure4.png", 
       figure4, 
       height = 16, 
       width = 18, 
       dpi = 500)

####Supplementary table 10 ######
#PERMANOVAS
phase_H21_nitrifiers_BC_adonis_df
phase_H21_nitrifiers_BC_adonis_pairwise_E2
phase_P1_nitrifiers_BC_adonis_df
phase_P1_nitrifiers_BC_adonis_pairwise_P

stable10 <- bind_rows(phase_H21_nitrifiers_BC_adonis_df, 
                      phase_H21_nitrifiers_BC_adonis_pairwise_E2, 
                      phase_P1_nitrifiers_BC_adonis_df, 
                      phase_P1_nitrifiers_BC_adonis_pairwise_P)

write_xlsx(stable10, "SupplementaryTable10.xlsx")

####Supplementary table 11#####
#PERMDISPS
H21_nitrifiers_permdisp_pairwise_results_table
P1_nitrifiers_permdisp_pairwise_results_table

stable11 <- bind_rows(
  H21_nitrifiers_permdisp_pairwise_results_table,
  P1_nitrifiers_permdisp_pairwise_results_table
)
write_xlsx(stable11, "SupplementaryTable11.xlsx")


# #Naive system
# #Days 0 -10 (2023-10-07 through 2023-10-17): "First" copper exposure, very low levels of copper (phase 1)
# #Days 11 - 28 (2023-10-18 through 2023-11-04): Down time between "first" copper exposure and start of actual first exposure to high copper levels
# #Day 29 - 56 (2023-11-05 through 2023-12-02): Actual first copper exposure (phase 2)
# #Day 57-86 (2023-12-03 through 2024-01-01): Downtime between first (phase 2) and second (phase 3) copper exposure
# #Day 87-131 (2024-01-02 through 2024-02-15): Second copper exposure (phase 3)
# #Day 132 - 147 (2024-02-16 through 2024-03-02): Downtime after phase 3 
# 
# #Established system
# #Days 0 - 50 (2023-11-14 through 2024-01-03): Phase 1, conservative copper dosing (phase 1)
# #Days 51 - 65 (2024-01-04 through 2024-01-18): Phase 2, increasing copper dosing (phase 2)
# #Days 66 - 98 (2024-01-19 through 2024-02-20): Phase 3, back to copper dosage as needed
# #Day 99 - 168 (2024-02-21 through 2024-04-30) Final, after copper treatment finished
# 
# #Days since start variable
# phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs_3 <- 
#   phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs %>%
#   mutate(Date_num = case_when(
#     Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07")),
#     Enclosure == "P1"  ~ as.numeric(Collection_Date - as.Date("2023-11-14"))
#   )) %>%
#   #Separate if wanting to do separate day scales
#   mutate(Date_num_naive = case_when(
#     Enclosure == "H21" ~ as.numeric(Collection_Date - as.Date("2023-10-07"))),
#     Date_num_established = case_when(
#       Enclosure == "P1" ~ as.numeric(Collection_Date - as.Date("2023-11-14"))))%>%
#   mutate(Date_num_phase_naive = case_when(
#     Enclosure == "H21" & Date_num >= 0   & Date_num <= 10  ~ "Phase 1 (Day 0 - 10)",
#     Enclosure == "H21" & Date_num >= 11  & Date_num <= 28  ~ "Downtime Phase 1 - 2 (Day 11 - 28)",
#     Enclosure == "H21" & Date_num >= 29  & Date_num <= 56  ~ "Phase 2 (Day 29 - 56)",
#     Enclosure == "H21" & Date_num >= 57  & Date_num <= 86  ~ "Downtime Phase 2 - 3 (Day 57 - 86)",
#     Enclosure == "H21" & Date_num >= 87  & Date_num <= 131 ~ "Phase 3 (Day 87 - 131)",
#     Enclosure == "H21" & Date_num >= 132 & Date_num <= 147 ~ "Post-Treatment completion (Day 132 - 147)"
#   )) %>%
#   mutate(Date_num_phase_naive = factor(Date_num_phase_naive, 
#                                        levels = c("Phase 1 (Day 0 - 10)", 
#                                                   "Downtime Phase 1 - 2 (Day 11 - 28)", 
#                                                   "Phase 2 (Day 29 - 56)", 
#                                                   "Downtime Phase 2 - 3 (Day 57 - 86)", 
#                                                   "Phase 3 (Day 87 - 131)", 
#                                                   "Post-Treatment completion (Day 132 - 147)")))%>%
#   mutate(Date_num_phase_established = case_when(
#     Enclosure == "P1" & Date_num >= 0   & Date_num <= 50  ~ "Phase 1 (Day 0 - 50)",
#     Enclosure == "P1" & Date_num >= 51  & Date_num <= 65  ~ "Phase 2 (Day 51 - 65)",
#     Enclosure == "P1" & Date_num >= 66  & Date_num <= 98  ~ "Phase 3 (Day 66 - 98)",
#     Enclosure == "P1" & Date_num >= 99  & Date_num <= 168 ~ "Post-Treatment completion (Day 99 - 168)"
#   ))%>%
#   mutate(Date_num_phase_established = factor(Date_num_phase_established, 
#                                        levels = c("Phase 1 (Day 0 - 50)",
#                                                   "Phase 2 (Day 51 - 65)",
#                                                   "Phase 3 (Day 66 - 98)",
#                                                   "Post-Treatment completion (Day 99 - 168)")))
# #Shape by phase 
# enclosure_BC_beta_div_3 <- ggplot(phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs_3) + 
#   theme_bw() +
#   # geom_segment(aes(x=cMDS1, y=cMDS2,
#   #                  xend= MDS1, yend = MDS2,
#   #                  color = Enclosure), show.legend = F)+
#   stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, fill= Enclosure, 
#                                      colour = Enclosure), alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
#   geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
#   geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
#   geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure, shape = Date_num_phase_naive), 
#              size = 5, alpha = 0.8) +
#   scale_shape_manual(values = c(15, 3, 16, 4, 17, 18), name = "Naive system\nPhase") +
#   guides(
#     shape = guide_legend(override.aes = list(size = 7), ncol = 1)) +
#   new_scale("shape") +
#   geom_point(aes(x=MDS1, y=MDS2, colour = Enclosure, shape = Date_num_phase_established), 
#              size = 5, alpha = 0.8) +
#   scale_shape_manual(values = c(15, 16, 17, 18), name = "Established system\nPhase") +
#   geom_point(aes(x=cMDS1, y= cMDS2, colour = Enclosure), size = 10) +
#   geom_text(aes (x= cMDS1, y = cMDS2,
#                  label= Enclosure), 
#             colour= "white", size = 6, fontface = "bold") +
#   labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
#        color = "System", 
#        fill = "System") +
#   scale_color_manual(values = enclosure.palette,
#                      labels = c("H21" = "Naive", 
#                                 "P1" = "Established"))+
#   scale_fill_manual(values = enclosure.palette, 
#                     labels = c("H21" = "Naive", 
#                                "P1" = "Established"))+
#   theme(legend.position = "bottom",
#         legend.title.position = "top",
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
#     shape = guide_legend(override.aes = list(size = 7), ncol = 1)) +
#   annotate("text", x = -1, y = -0.9,
#            label = "PERMANOVA\nSystem",
#            hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
#   annotate("text", x = -1, y = -0.85,
#            label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
#                          "\np = ", round(pvalue_adonis_enclosure, 4)),
#            hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
# enclosure_BC_beta_div_3
# 
# ##Colored by days since start
# enclosure_BC_beta_div_4 <- ggplot(phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs_3) + 
#   theme_bw() +
#   labs(x="NMDS1", y="NMDS2", title= "MICROBIOME", 
#        shape = "System", 
#        color = "Day") +
#   # geom_segment(aes(x=cMDS1, y=cMDS2,
#   #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
#   stat_ellipse(geom= "polygon", aes (x= MDS1, y = MDS2, group = Enclosure), 
#                #color = "black",
#                alpha = 0.2, lty = 2, linewidth = 1, level= 0.95)+
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   geom_point(aes(x=MDS1, y=MDS2, 
#                  colour = Date_num,
#                  shape = Enclosure), 
#              size = 5, 
#              alpha = 0.8) + # individuals
#   geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure), size = 10) + # centroids +
#   geom_text(aes (x= cMDS1, y = cMDS2,
#                  label= Enclosure),
#             colour= "white", size = 5, fontface = "bold") +
#   scale_shape_manual(values = c(15,16),
#                      labels = c("H21" = "Naive", 
#                                 "P1" = "Established"))+
#   scale_color_viridis_c(option = "viridis")+
#   theme(legend.position = "bottom",
#         legend.title = element_text(colour = "black", 
#                                     size = 22,
#                                     face = "bold"),
#         legend.direction = "vertical",
#         legend.text = element_text(colour = "black", size = 20),
#         plot.title = element_text(size = 45),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth = 1),
#         axis.ticks = element_line(colour = "black", linewidth = 0.75),
#         axis.title = element_text(size = 30),
#         axis.text = element_text(size = 26, colour = "black"))+
#   guides(
#     shape = guide_legend(override.aes = list(size = 7))) +
#   annotate("text", x = -0.9, y = -0.9,
#            label = "PERMANOVA\nSystem",
#            hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
#   annotate("text", x = -0.9, y = -0.8,
#            label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
#                          "\np = ", round(pvalue_adonis_enclosure, 4)),
#            hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
# enclosure_BC_beta_div_4
# ggsave("enclosure_BC_beta_div_4.png", 
#        enclosure_BC_beta_div_4, 
#        device = "png", 
#        dpi = 600, 
#        height = 8, 
#        width = 15)
# 
# 
# #Colored by days since start, separate scales for each system
# enclosure_BC_beta_div_5 <- ggplot(phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs_3) + 
#   theme_bw() +
#   labs(x="NMDS1", y="NMDS2", 
#        title= "MICROBIOME", 
#        shape = "System") +
#   # geom_segment(aes(x=cMDS1, y=cMDS2,
#   #                  xend= MDS1, yend = MDS2, alpha = 1), show.legend = F)+
#   stat_ellipse(
#     geom = "polygon",
#     aes(x = MDS1, y = MDS2, fill = Enclosure, color = Enclosure, group = Enclosure),
#     alpha = 0.2,
#     lty = 2,
#     linewidth = 1,
#     level = 0.95,
#     show.legend = FALSE
#   )+
#   #Problem with layering. If adding centroids here, dots go on top. 
#   geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure,
#                  color = Enclosure), size = 10,
#              show.legend = FALSE) + # centroids +
#   geom_text(aes (x= cMDS1, y = cMDS2,
#                  label= Enclosure),
#             colour= "white", size = 5, fontface = "bold") +
#   scale_fill_manual(values = c(
#     "#fc8d62",
#     "#8da0cb"
#   ))+
#   scale_color_manual(values = c(
#     "#fc8d62",
#     "#8da0cb"
#   ))+
#   new_scale_color()+
#   geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
#   geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
#   #Naive system individual dots
#   geom_point(data = phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs_3 %>%
#                filter(Enclosure == "H21"), aes(x=MDS1, y=MDS2, 
#                  colour = Date_num_naive,
#                  shape = Enclosure), 
#              size = 5) + # individuals
#   scale_color_gradient2(low = "#fad2c2", 
#                         mid = "#fc8d62", 
#                         high = "#ad431a",
#                         midpoint = 74,
#                         limits = c(0,150),
#                         breaks = c(0, 50, 100, 150))+
#   # scale_color_viridis_c(option = "viridis",
#   #                       limits = c(0,150),
#   #                       breaks = c(0, 50, 100, 150))+
#   labs(color = "Days (Naive system)")+
#   #Established system individual dots. 
#   new_scale_color()+
#   geom_point(data = phyloseq.samples.bacteria.ra.dates.bray.enclosure.segs_3 %>%
#                filter(Enclosure == "P1"), aes(x=MDS1, y=MDS2, 
#                  colour = Date_num_established,
#                  shape = Enclosure), 
#              size = 5) + # individuals
#   scale_color_gradient2(low = "#c6cfe6",
#                         mid = "#8da0cb", 
#                         high = "#012983", 
#                         midpoint = 84,
#                         limits = c(0,170),
#                         breaks = c(0, 50, 100, 170))+
#   labs(color = "Days (Established system)")+
#   # scale_color_viridis_c(option = "rocket", 
#   #                       limits = c(0,170),
#   #                       breaks = c(0, 50, 100, 170))+
#   scale_shape_manual(values = c(15,16),
#                      labels = c("H21" = "Naive", 
#                                 "P1" = "Established"))+
#   theme(legend.position = "bottom",
#         legend.title = element_text(colour = "black", 
#                                     size = 22,
#                                     face = "bold"),
#         legend.direction = "vertical",
#         legend.text = element_text(colour = "black", size = 18),
#         plot.title = element_text(size = 45),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", linewidth = 1),
#         axis.ticks = element_line(colour = "black", linewidth = 0.75),
#         axis.title = element_text(size = 30),
#         axis.text = element_text(size = 26, colour = "black"))+
#   guides(
#     color = guide_colorbar(order = 3),
#     #color = guide_legend(order = 3),
#     shape = guide_legend(override.aes = list(size = 7),
#                          order = 1)) +
#   # new_scale_color()+
#   # geom_point(aes(x=cMDS1, y= cMDS2, shape = Enclosure, 
#   #                color = Enclosure), size = 10,
#   #            show.legend = FALSE) + # centroids +
#   # geom_text(aes (x= cMDS1, y = cMDS2,
#   #                label= Enclosure),
#   #           colour= "white", size = 5, fontface = "bold") +
#   # scale_color_manual(values = c(
#   #   "#fc8d62",
#   #   "#8da0cb"
#   # ))+
#   annotate("text", x = -0.9, y = -0.9,
#            label = "PERMANOVA\nSystem",
#            hjust = 0.5, vjust = -0.5, size = 8, colour = "black", fontface = "bold") + ##annotate variable (Enclosure)
#   annotate("text", x = -0.9, y = -0.8,
#            label = paste("R² = ", round(R2_adonis_enclosure* 100, 1), "%",
#                          "\np = ", round(pvalue_adonis_enclosure, 4)),
#            hjust = 0.5, vjust = 1.1, size = 8, colour = "black")# Annotate R² and p-values
# enclosure_BC_beta_div_5
# ggsave("enclosure_BC_beta_div_5.png", 
#        enclosure_BC_beta_div_5, 
#        device = "png", 
#        dpi = 600, 
#        height = 8, 
#        width = 15)
# 
# #POSTER BIMS 2026
# ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Posters/enclosure_BC_beta_div.png", 
#        enclosure_BC_beta_div+
#          theme(plot.title = element_text(size = 50, face = "bold")), 
#        device = "png", 
#        dpi = 600, 
#        height = 8, 
#        width = 19)


###BOTH NAIVE AND ESTABLISHED SYSTEM######
####OVERALL MICROBIAL COMMUNITY#####
#All samples
phyloseq.samples.bacteria
#Relative abundances
phyloseq.samples.bacteria.ra

#Bray curtis distances
phyloseq.samples.bacteria.RA.bray <- vegdist(t(phyloseq.samples.bacteria.ra@otu_table), method = "bray")

#Ordination
phyloseq.samples.bacteria.RA.bray.ord <- metaMDS(phyloseq.samples.bacteria.RA.bray, 
                                                             k=2, try = 50, 
                                                             trymax = 1000,
                                                             autotransform = F)
#Metadata 
phyloseq.samples.bacteria.metadata <- data.frame(sample_data(phyloseq.samples.bacteria))

###Centroids by System and Phase#####
#Simple ordination plot
phyloseq.samples.bacteria.RA.bray.ord_plot <- ordiplot(phyloseq.samples.bacteria.RA.bray.ord$points)

#Now, extract coordinates
phyloseq.samples.bacteria.RA.bray.ord_plot_scrs <- scores(phyloseq.samples.bacteria.RA.bray.ord_plot, display = "sites")

#Add metadata to coordinates
phyloseq.samples.bacteria.RA.bray.ord_plot_scrs <- cbind(as.data.frame(phyloseq.samples.bacteria.RA.bray.ord_plot_scrs),
                                                         Enclosure = phyloseq.samples.bacteria.metadata$Enclosure,
                                                                     Copper_level_mg_L = phyloseq.samples.bacteria.metadata$Copper_level_mg_L, 
                                                                     SampleID = phyloseq.samples.bacteria.metadata$SampleID, 
                                                                     Date_num = phyloseq.samples.bacteria.metadata$Date_num,
                                                                     Collection_Date = phyloseq.samples.bacteria.metadata$Collection_Date, 
                                                                     Ammonia_level = phyloseq.samples.bacteria.metadata$Ammonia_mg_L, 
                                                                     Date_num_phase = phyloseq.samples.bacteria.metadata$Date_num_phase, 
                                                                     Date_num_phase_abbrv = phyloseq.samples.bacteria.metadata$Date_num_phase_abbrv)
##Calculate centroids according to Phase
phyloseq.samples.nitrifiers_H21.ra.dates.bray.cent.phase <- aggregate(cbind(MDS1,MDS2) ~ Enclosure + Date_num_phase_abbrv, 
                                                                    data = phyloseq.samples.bacteria.RA.bray.ord_plot_scrs, 
                                                                    FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.samples.bacteria.RA.bray.ord_plot_segs <- merge(phyloseq.samples.bacteria.RA.bray.ord_plot_scrs, 
                                                                     setNames(phyloseq.samples.nitrifiers_H21.ra.dates.bray.cent.phase, 
                                                                              c("Enclosure", "Date_num_phase_abbrv","cMDS1","cMDS2")),
                                                                     by = c('Enclosure', 'Date_num_phase_abbrv'), 
                                                                     sort = F)

####PERMANOVA - Enclosure + Phase####
set.seed(98)
phase_enclosure_BC_adonis  <- adonis2(phyloseq.samples.bacteria.RA.bray ~ Enclosure + Date_num_phase_abbrv, 
                                phyloseq.samples.bacteria.metadata, 
                                #strata = phyloseq.samples.bacteria.metadata$Collection_Month,
                                by = "margin",
                                permutations = 9999)
phase_enclosure_BC_adonis 

# ##PERMDISPS
# # Run the betadisper function, average distance to centroid
# bray.phase.enclosure.disp <- betadisper(phyloseq.samples.bacteria.RA.bray, 
#                                   phyloseq.samples.bacteria.metadata$Date_num_phase_abbrv)
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
phase_enclosure_BC_beta_div <- ggplot(phyloseq.samples.bacteria.RA.bray.ord_plot_segs) + 
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
phyloseq.samples.bacteria.ra.dates.bray.cent.copper_quart <- aggregate(
  cbind(MDS1, MDS2) ~ Enclosure + Copper_quartile,
  data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs,
  FUN = mean)

#Merge centroids with coordinates and metadata
phyloseq.samples.bacteria.ra.dates.bray.copper_quart.segs <- merge(phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs, 
                                                          setNames(phyloseq.samples.bacteria.ra.dates.bray.cent.copper_quart, 
                                                                   c("Enclosure", "Copper_quartile",
                                                                     "cMDS1","cMDS2")),
                                                          by = c("Enclosure", "Copper_quartile"), 
                                                          sort = F)
##PERMANOVA - Just copper###
set.seed(98)
copper_quart_BC_adonis  <- adonis2(phyloseq.samples.bacteria.ra.dates.bray ~ Copper_quartile, 
                                phyloseq.samples.bacteria.df, 
                                strata = phyloseq.samples.bacteria.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
copper_quart_BC_adonis 


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_adonis  <- adonis2(phyloseq.samples.bacteria.ra.dates.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.samples.bacteria.df, 
#                                 strata = phyloseq.samples.bacteria.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_adonis  <- adonis2(phyloseq.samples.bacteria.ra.dates.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.samples.bacteria.df , 
#                                 strata = phyloseq.samples.bacteria.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS#
#COPPER quartile#
# Run the betadisper function, average distance to centroid
bray.copper.quart.disp <- betadisper(phyloseq.samples.bacteria.ra.dates.bray,
                                  phyloseq.samples.bacteria.df$Copper_quartile)
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
copper_quart_BC_beta_div <- ggplot(phyloseq.samples.bacteria.ra.dates.bray.copper_quart.segs) + 
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
phyloseq.samples.bacteria.df <- as(phyloseq.samples.bacteria.ra@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.df 


###Continuous copper facetted by enclosure #####
#No samples without copper levels:
phyloseq.samples.bacteria.ra.copper <- subset_samples(phyloseq.samples.bacteria.ra.dates, 
                                                                !is.na(Copper_level_mg_L))
phyloseq.samples.bacteria.ra.copper <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.ra.copper)> 0, 
                                                            phyloseq.samples.bacteria.ra.copper)
phyloseq.samples.bacteria.ra.copper #33481 taxa and 204 samples

####H21 #######
phyloseq.samples.bacteria.ra.copper.H21 <- subset_samples(phyloseq.samples.bacteria.ra.copper, 
                                                                    Enclosure == "H21")
phyloseq.samples.bacteria.ra.copper.H21 <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.ra.copper.H21)> 0, 
                                                                phyloseq.samples.bacteria.ra.copper.H21)
phyloseq.samples.bacteria.ra.copper.H21 #33418 taxa and 82 samples 

##BC 
phyloseq.samples.bacteria.ra.copper.H21.bray <- vegdist(t(phyloseq.samples.bacteria.ra.copper.H21@otu_table), 
                                                                  method = "bray")

#make DF from metadata
phyloseq.samples.bacteria.copper.H21.df <- as(phyloseq.samples.bacteria.ra.copper.H21@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.copper.H21.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_H21 <- D2K(
  as.matrix(phyloseq.samples.bacteria.ra.copper.H21.bray))

# Covariates (date_num)
datenum_covariate_H21 <- phyloseq.samples.bacteria.copper.H21.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_H21 

# # Covariates (month)
# month_covariate_H21 <- model.matrix(~ Collection_Month, 
#                                         data = phyloseq.samples.bacteria.copper.H21.df)[, -1]
# 

#Run model
mirkat_H21_copper_cont <- MiRKAT(
  y = phyloseq.samples.bacteria.copper.H21.df$Copper_level_mg_L,
  X = datenum_covariate_H21,
  Ks = bray_kernel_H21,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_H21_copper_cont #R2 6.9% , p_values = 0.0001017485

#####PERMANOVA #######
set.seed(98)
adonis_H21_copper_cont  <- adonis2(phyloseq.samples.bacteria.ra.copper.H21.bray ~ Copper_level_mg_L, 
                                       phyloseq.samples.bacteria.copper.H21.df, 
                                       strata = phyloseq.samples.bacteria.copper.H21.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_H21_copper_cont # R2 13.5%, p value 2e-04 *

##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.copper.cont.disp.H21 <- betadisper(phyloseq.samples.bacteria.ra.copper.H21.bray,
                                            phyloseq.samples.bacteria.copper.H21.df$Copper_quartile)
bray.copper.cont.disp.H21
##Then test by permuting
set.seed(98)
bray.copper.cont.permdisp.H21 <- permutest(bray.copper.cont.disp.H21, permutations = 9999)
bray.copper.cont.permdisp.H21 ##S, p =0.0016

####P1 #######
phyloseq.samples.bacteria.ra.copper.P1 <- subset_samples(phyloseq.samples.bacteria.ra.copper, 
                                                                   Enclosure == "P1")
phyloseq.samples.bacteria.ra.copper.P1 <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.ra.copper.P1)> 0, 
                                                               phyloseq.samples.bacteria.ra.copper.P1)
phyloseq.samples.bacteria.ra.copper.P1 #28427 taxa and 122 samples

##BC 
phyloseq.samples.bacteria.ra.copper.P1.bray <- vegdist(t(phyloseq.samples.bacteria.ra.copper.P1@otu_table), 
                                                                 method = "bray")

#make DF from metadata
phyloseq.samples.bacteria.copper.P1.df <- as(phyloseq.samples.bacteria.ra.copper.P1@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.copper.P1.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_P1 <- D2K(
  as.matrix(phyloseq.samples.bacteria.ra.copper.P1.bray))

# Covariates (date_num)
datenum_covariate_P1 <- phyloseq.samples.bacteria.copper.P1.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_P1 

# # Covariates (month)
# month_covariate_P1 <- model.matrix(~ Collection_Month, 
#                                        data = phyloseq.samples.bacteria.copper.P1.df)[, -1]
# 
# month_covariate_P1

#Run model
mirkat_P1_copper_cont <- MiRKAT(
  y = phyloseq.samples.bacteria.copper.P1.df$Copper_level_mg_L,
  X = datenum_covariate_P1,
  Ks = bray_kernel_P1,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_P1_copper_cont #p val 1.320553e-08, R2 2.5%

#####PERMANOVA #####
set.seed(98)
adonis_P1_copper_cont  <- adonis2(phyloseq.samples.bacteria.ra.copper.P1.bray ~ Copper_level_mg_L, 
                                      phyloseq.samples.bacteria.copper.P1.df, 
                                      strata = phyloseq.samples.bacteria.copper.P1.df$Collection_Month,
                                      by = "margin",
                                      permutations = 9999)
adonis_P1_copper_cont # R2 10.4%, p value 0.021 *


####MIRKAT GLMM FOR COPPER #####
#make DF from metadata
phyloseq.samples.bacteria.copper.all.df <- as(phyloseq.samples.bacteria.ra.copper@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.copper.all.df

#Bray-Curtis distances
phyloseq.samples.bacteria.ra.copper.bray <- vegdist(t(phyloseq.samples.bacteria.ra.copper@otu_table), 
                                                              method = "bray")
phyloseq.samples.bacteria.ra.copper.bray

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_all <- D2K(
  as.matrix(phyloseq.samples.bacteria.ra.copper.bray))

#Jaccard distances
phyloseq.samples.bacteria.ra.copper.jac <- vegdist(t(phyloseq.samples.bacteria.ra.copper@otu_table), 
                                                             method = "jaccard")
phyloseq.samples.bacteria.ra.copper.jac

jac_kernel_all <- D2K(
  as.matrix(phyloseq.samples.bacteria.ra.copper.jac))

# Covariates (date_num)
datenum_covariate_all <- phyloseq.samples.bacteria.copper.all.df%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_all

#Id for repeated measures 
Enclosure_d_mirkat_all <- phyloseq.samples.bacteria.copper.all.df %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_mirkat_all

#Run model - dont think it works. just two enclosures and enclosure id is the same as the outcome 
mirkat_enclosure_glm <- GLMMMiRKAT(
  y = phyloseq.samples.bacteria.copper.all.df$Copper_level_mg_L,
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
copper_cont_BC_beta_div <- ggplot(phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs) + 
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
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.5, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.7, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "P1"),
            aes(x = 0.2, y = -1.1, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
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
copper_cont_BC_beta_div_2 <- ggplot(phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs) + 
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
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.5, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.7, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
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
phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs_2 <- phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs %>%
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
copper_cont_BC_beta_div_3 <- ggplot(phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs_2) + 
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
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "H21"),
            aes(x = 0.2, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "H21"),
            aes(x =0.2, y = -0.8, 
                label = paste("R² = ", round(R2_adonis_copper_cont_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_H21, 5))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs_2 %>%
              filter(Enclosure == "P1"),
            aes(x = -0.45, y = -0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.H21.ordered_RA.bray.ord_plot_scrs_2 %>%
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
phyloseq.samples.bacteria.nitifiers.ra <- subset_taxa(phyloseq.samples.bacteria.ra.dates, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
phyloseq.samples.bacteria.nitifiers.ra #827 taxa and 218 samples
phyloseq.samples.bacteria.nitifiers.ra <- subset_samples(phyloseq.samples.bacteria.nitifiers.ra, sample_sums(phyloseq.samples.bacteria.nitifiers.ra) > 0)
phyloseq.samples.bacteria.nitifiers.ra #827 taxa and 218 samples 

##BC 
phyloseq.samples.bacteria.nitifiers.ra.bray <- vegdist(t(phyloseq.samples.bacteria.nitifiers.ra@otu_table), method = "bray")

#make DF from metadata
phyloseq.samples.bacteria.nitifiers.df <- as(phyloseq.samples.bacteria.nitifiers.ra@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.nitifiers.df 

###ORDINATION####
set.seed(98)
phyloseq.samples.bacteria.nitifiers.ra.bray.ord <- metaMDS(phyloseq.samples.bacteria.nitifiers.ra.bray, 
                                                 k=3, try = 50, 
                                                 trymax = 1000,
                                                 autotransform = F)
###Centroids by enclosure #####
## BC
#Simple ordination plot
phyloseq.samples.bacteria.nitifiers.ra.bray.plot <- ordiplot(phyloseq.samples.bacteria.nitifiers.ra.bray.ord$points)

#Now, extract coordinates
phyloseq.samples.bacteria.nitifiers.ra.bray.scrs <- scores(phyloseq.samples.bacteria.nitifiers.ra.bray.plot, display = "sites")
#Add metadata to coordinates
phyloseq.samples.bacteria.nitifiers.ra.bray.scrs <- cbind(as.data.frame(phyloseq.samples.bacteria.nitifiers.ra.bray.scrs),
                                                Copper_level_mg_L = phyloseq.samples.bacteria.nitifiers.df$Copper_level_mg_L, 
                                                SampleID = phyloseq.samples.bacteria.nitifiers.df$SampleID, 
                                                Collection_Date = phyloseq.samples.bacteria.nitifiers.df$Collection_Date, 
                                                Enclosure = phyloseq.samples.bacteria.nitifiers.df$Enclosure, 
                                                Attempt = phyloseq.samples.bacteria.nitifiers.df$Attempt, 
                                                Copper_quartile = phyloseq.samples.bacteria.nitifiers.df$Copper_quartile,
                                                Copper_quartile.abvr = phyloseq.samples.bacteria.nitifiers.df$Copper_quartile.abvr)%>%
  mutate(Collection_Month = format(Collection_Date, "%Y-%m"))
##Calculate centroids according to Enclosure
phyloseq.samples.bacteria.nitifiers.ra.bray.cent.enclosure <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
                                                              data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs, 
                                                              FUN = mean) 
#Merge centroids with coordinates and metadata
phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs <- merge(phyloseq.samples.bacteria.nitifiers.ra.bray.scrs, 
                                                          setNames(phyloseq.samples.bacteria.nitifiers.ra.bray.cent.enclosure, 
                                                                   c("Enclosure","cMDS1","cMDS2")),
                                                          by = 'Enclosure', 
                                                          sort = F)

####MIRKAT- Just enclosure#####
#Dataframe with metadata
df_metadata_nit_mirkat <- data.frame(sample_data(phyloseq.samples.bacteria.nitifiers.ra))
#BC distance matrix row order
labs_nit <- rownames(as.matrix(phyloseq.samples.bacteria.nitifiers.ra.bray))  
#Align metadata df to distance order
df_nit_mirkat <- df_metadata_nit_mirkat[labs_nit, , drop = FALSE] 

#Going to get kernels out of different distance matrices
#Bray curtis 
phyloseq.samples.bacteria.ra.dates.bray

#Jaccard distances
phyloseq.samples.bacteria.ra.dates.jac <- vegdist(t(phyloseq.samples.bacteria.ra.dates@otu_table), method = "jaccard")
phyloseq.samples.bacteria.ra.dates.jac

#Aitchison 
##Going to take out samples from P1 from april and may 2023 and september from H21 (raw counts object)
phyloseq.samples.bacteria.dates <- subset_samples(phyloseq.samples.bacteria, 
                                                     Collection_Date > "2023-10-01")
phyloseq.samples.bacteria.dates #33490 taxa and 218 samples
#Take out taxa sums = 0
phyloseq.samples.bacteria.dates  <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.dates) > 0, 
                                               phyloseq.samples.bacteria.dates)
phyloseq.samples.bacteria.dates  #33490 taxa and 218 samples
#Calculate distance
phyloseq.samples.bacteria.dates_clr <- microbiome::transform(phyloseq.samples.bacteria.dates, "clr") #convert raw counts to clr
phyloseq.samples.bacteria.dates_clr_dist_aitch <- dist(t(otu_table(phyloseq.samples.bacteria.dates_clr)), method = "euclidean") #calculate euclidean distances

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
#Bray curtis
bray_kernel_nit <- D2K(
  as.matrix(phyloseq.samples.bacteria.nitifiers.ra.bray))
#Jaccards
jac_kernel_nit <- D2K(
  as.matrix(phyloseq.samples.bacteria.ra.dates.jac))
# Aitchison
aitch_kernel_nit <- D2K(
  as.matrix(phyloseq.samples.bacteria.dates_clr_dist_aitch))

#Modelling enclosure, have to make it dichotomus (used to be phyloseq.samples.bacteria.nitifiers.df)
Enclosure_d_nit_mirkat <- df_nit_mirkat %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_nit_mirkat

# # Covariates (month)
# month_covariate_nit <- model.matrix(~ Collection_Month, 
#                                 data = phyloseq.samples.bacteria.nitifiers.df)[, -1]
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
enclosure_BC_nit_adonis  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.bray ~ Enclosure, 
                                phyloseq.samples.bacteria.nitifiers.df, 
                                strata = phyloseq.samples.bacteria.nitifiers.df$Collection_Month,
                                by = "margin",
                                permutations = 9999)
enclosure_BC_nit_adonis #21% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_nit_adonis  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.samples.bacteria.nitifiers.df, 
#                                 strata = phyloseq.samples.bacteria.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_nit_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_nit_adonis  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.samples.bacteria.nitifiers.df , 
#                                 strata = phyloseq.samples.bacteria.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_nit_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS
#ENCLOSURE#
# Run the betadisper function, average distance to centroid
bray.enclosure.nit.disp <- betadisper(phyloseq.samples.bacteria.nitifiers.ra.bray, 
                                  phyloseq.samples.bacteria.nitifiers.df$Enclosure)
bray.enclosure.nit.disp
##Then test by permuting
set.seed(98)
bray.enclosure.nit.permdisp <- permutest(bray.enclosure.nit.disp, permutations = 9999)
bray.enclosure.nit.permdisp ##S, p = 6e-04

# #COPPER#
# # Run the betadisper function, average distance to centroid
# bray.copper.disp <- betadisper(phyloseq.samples.bacteria.nitifiers.ra.bray, 
#                                   phyloseq.samples.bacteria.nitifiers.df$Copper_level_mg_L)
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
enclosure_BC_nit_beta_div <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs) + theme_bw() +
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
enclosure_BC_nit_beta_div_2 <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs) + 
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
phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs_2 <- phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs%>%
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
enclosure_BC_nit_beta_div_3 <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs_2) + 
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
enclosure_BC_nit_beta_div_4 <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs_2) + 
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
  geom_point(data = phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs_2%>%
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
  geom_point(data = phyloseq.samples.bacteria.nitifiers.ra.bray.enclosure.segs_2 %>%
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
# fit <- envfit(phyloseq.samples.bacteria.nitifiers.ra.bray.ord, 
#               phyloseq.samples.bacteria.nitifiers.df$Copper_level_mg_L, 
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
phyloseq.samples.bacteria.nitifiers.ra.bray.cent.copper_quart <- aggregate(
  cbind(MDS1, MDS2) ~ Enclosure + Copper_quartile,
  data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs,
  FUN = mean)

#Merge centroids with coordinates and metadata
phyloseq.samples.bacteria.nitifiers.ra.bray.copper_quart.segs <- merge(phyloseq.samples.bacteria.nitifiers.ra.bray.scrs, 
                                                             setNames(phyloseq.samples.bacteria.nitifiers.ra.bray.cent.copper_quart, 
                                                                      c("Enclosure", "Copper_quartile",
                                                                        "cMDS1","cMDS2")),
                                                             by = c("Enclosure", "Copper_quartile"), 
                                                             sort = F)
##PERMANOVA - Just copper###
# set.seed(98)
# copper_quart_BC_nit_adonis  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.bray ~ Copper_quartile, 
#                                    phyloseq.samples.bacteria.nitifiers.df, 
#                                    strata = phyloseq.samples.bacteria.nitifiers.df$Collection_Month,
#                                    by = "margin",
#                                    permutations = 9999)
# copper_quart_BC_nit_adonis #30.1% of the variation is due to Enclosure, p = 1e-04


# ##PERMANOVA###
# set.seed(98)
# enclosure_copper_BC_nit_adonis  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.bray ~ Enclosure + 
#                                   Copper_level_mg_L, 
#                                 phyloseq.samples.bacteria.nitifiers.df, 
#                                 strata = phyloseq.samples.bacteria.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_BC_nit_adonis #11.6% of the variation is due to Enclosure, p = 1e-04
# #4.1% of the variation is due to copper levels, p = 1e-04
# 
# 
# #With interaction
# set.seed(98)
# enclosure_copper_interac_BC_nit_adonis  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.samples.bacteria.nitifiers.df , 
#                                 strata = phyloseq.samples.bacteria.nitifiers.df$Collection_Month,
#                                 by = "margin",
#                                 permutations = 9999)
# enclosure_copper_interac_BC_nit_adonis
# #Enclosure:Copper_level_mg_L interaction not significant (p = 0.24) 


##PERMDISPS#
#COPPER quartile#
# Run the betadisper function, average distance to centroid
bray.nit.copper.quart.disp <- betadisper(phyloseq.samples.bacteria.nitifiers.ra.bray,
                                     phyloseq.samples.bacteria.nitifiers.df$Copper_quartile)
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
copper_quart_BC_nit_beta_div <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.copper_quart.segs) + 
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
phyloseq.samples.bacteria.nitifiers.ra.copper <- subset_samples(phyloseq.samples.bacteria.nitifiers.ra, 
                                                                !is.na(Copper_level_mg_L))
phyloseq.samples.bacteria.nitifiers.ra.copper <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.nitifiers.ra.copper)> 0, 
                                                            phyloseq.samples.bacteria.nitifiers.ra.copper)
phyloseq.samples.bacteria.nitifiers.ra.copper #827 taxa and 204 samples

####H21 #######
phyloseq.samples.bacteria.nitifiers.ra.copper.H21 <- subset_samples(phyloseq.samples.bacteria.nitifiers.ra.copper, 
                                                                Enclosure == "H21")
phyloseq.samples.bacteria.nitifiers.ra.copper.H21 <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.nitifiers.ra.copper.H21)> 0, 
                                                                phyloseq.samples.bacteria.nitifiers.ra.copper.H21)
phyloseq.samples.bacteria.nitifiers.ra.copper.H21 # 824 taxa and 82 samples 

##BC 
phyloseq.samples.bacteria.nitifiers.ra.copper.H21.bray <- vegdist(t(phyloseq.samples.bacteria.nitifiers.ra.copper.H21@otu_table), 
                                                              method = "bray")

#make DF from metadata
phyloseq.samples.bacteria.nitifiers.copper.H21.df <- as(phyloseq.samples.bacteria.nitifiers.ra.copper.H21@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.nitifiers.copper.H21.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_H21 <- D2K(
  as.matrix(phyloseq.samples.bacteria.nitifiers.ra.copper.H21.bray))

#Jaccard distances
phyloseq.samples.bacteria.nitifiers.ra.copper.H21.jac <- vegdist(t(phyloseq.samples.bacteria.nitifiers.ra.copper.H21@otu_table), 
                                                                 method = "jaccard")
phyloseq.samples.bacteria.nitifiers.ra.copper.H21.jac

jac_kernel_nit_H21 <- D2K(
  as.matrix(phyloseq.samples.bacteria.nitifiers.ra.copper.H21.jac))

# # Covariates (month)
# month_covariate_nit_H21 <- model.matrix(~ Collection_Month, 
#                                     data = phyloseq.samples.bacteria.nitifiers.copper.H21.df)[, -1]
# 
# month_covariate_nit_H21

# Covariates (date_num)
datenum_covariate_H21_nit <- phyloseq.samples.bacteria.nitifiers.copper.H21.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_H21_nit

#Collection month for random intercepts
collection_month_H21_nit <- phyloseq.samples.bacteria.nitifiers.copper.H21.df%>%
  pull(Collection_Month)%>%
  as.matrix()
collection_month_H21_nit

#Run model
mirkat_nit_H21_copper_cont <- MiRKAT(
  y = phyloseq.samples.bacteria.nitifiers.copper.H21.df$Copper_level_mg_L,
  X = datenum_covariate_H21_nit,
  Ks = bray_kernel_nit_H21,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_nit_H21_copper_cont #R2 3.1% , p_values = 0.004

#####PERMANOVA #######
set.seed(98)
adonis_nit_H21_copper_cont  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.copper.H21.bray ~ Copper_level_mg_L, 
                                       phyloseq.samples.bacteria.nitifiers.copper.H21.df, 
                                       strata = phyloseq.samples.bacteria.nitifiers.copper.H21.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_nit_H21_copper_cont # R2 9.7%, p value 0.027*

##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.nit.copper.cont.disp.H21 <- betadisper(phyloseq.samples.bacteria.nitifiers.ra.copper.H21.bray,
                                           phyloseq.samples.bacteria.nitifiers.copper.H21.df$Copper_quartile)
bray.nit.copper.cont.disp.H21
##Then test by permuting
set.seed(98)
bray.nit.copper.cont.permdisp.H21 <- permutest(bray.nit.copper.cont.disp.H21, permutations = 9999)
bray.nit.copper.cont.permdisp.H21 ##NS, p = 0.65

####P1 #######
phyloseq.samples.bacteria.nitifiers.ra.copper.P1 <- subset_samples(phyloseq.samples.bacteria.nitifiers.ra.copper, 
                                                                    Enclosure == "P1")
phyloseq.samples.bacteria.nitifiers.ra.copper.P1 <- prune_taxa(taxa_sums(phyloseq.samples.bacteria.nitifiers.ra.copper.P1)> 0, 
                                                                phyloseq.samples.bacteria.nitifiers.ra.copper.P1)
phyloseq.samples.bacteria.nitifiers.ra.copper.P1 #624 taxa and 122 samples

##BC 
phyloseq.samples.bacteria.nitifiers.ra.copper.P1.bray <- vegdist(t(phyloseq.samples.bacteria.nitifiers.ra.copper.P1@otu_table), 
                                                                  method = "bray")

#make DF from metadata
phyloseq.samples.bacteria.nitifiers.copper.P1.df <- as(phyloseq.samples.bacteria.nitifiers.ra.copper.P1@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.nitifiers.copper.P1.df

#####MIRKAT#######
#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_P1 <- D2K(
  as.matrix(phyloseq.samples.bacteria.nitifiers.ra.copper.P1.bray))

# # Covariates (month)
# month_covariate_nit_P1 <- model.matrix(~ Collection_Month, 
#                                         data = phyloseq.samples.bacteria.nitifiers.copper.P1.df)[, -1]
# 
# month_covariate_nit_P1

# Covariates (date_num)
datenum_covariate_P1_nit <- phyloseq.samples.bacteria.nitifiers.copper.P1.df%>%
  mutate(Date_num = as.numeric(Collection_Date - min(Collection_Date)))%>%
  ungroup()%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_P1_nit

#Run model
mirkat_nit_P1_copper_cont <- MiRKAT(
  y = phyloseq.samples.bacteria.nitifiers.copper.P1.df$Copper_level_mg_L,
  X = datenum_covariate_P1_nit,
  Ks = bray_kernel_nit_P1,
  out_type = "C", 
  returnKRV = TRUE,
  returnR2 = TRUE)
mirkat_nit_P1_copper_cont #p val 0 , R2 20.6%

#####PERMANOVA #####
set.seed(98)
adonis_nit_P1_copper_cont  <- adonis2(phyloseq.samples.bacteria.nitifiers.ra.copper.P1.bray ~ Copper_level_mg_L, 
                                       phyloseq.samples.bacteria.nitifiers.copper.P1.df, 
                                       strata = phyloseq.samples.bacteria.nitifiers.copper.P1.df$Collection_Month,
                                       by = "margin",
                                       permutations = 9999)
adonis_nit_P1_copper_cont # R2 26.9%, p value 1e-04*
##PERMDISP#
#COPPER continuous#
# Run the betadisper function, average distance to centroid
bray.nit.copper.cont.disp.P1 <- betadisper(phyloseq.samples.bacteria.nitifiers.ra.copper.P1.bray,
                                           phyloseq.samples.bacteria.nitifiers.copper.P1.df$Copper_quartile)
bray.nit.copper.cont.disp.P1
##Then test by permuting
set.seed(98)
bray.nit.copper.cont.permdisp.P1 <- permutest(bray.nit.copper.cont.disp.P1, permutations = 9999)
bray.nit.copper.cont.permdisp.P1 ##S, p = 1e-04 

###MIRKAT GLMM FOR COPPER #####
#make DF from metadata
phyloseq.samples.bacteria.nitifiers.copper.all.df <- as(phyloseq.samples.bacteria.nitifiers.ra.copper@sam_data, "data.frame") %>%
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
phyloseq.samples.bacteria.nitifiers.copper.all.df

#Bray-Curtis distances
phyloseq.samples.bacteria.nitifiers.ra.copper.bray <- vegdist(t(phyloseq.samples.bacteria.nitifiers.ra.copper@otu_table), 
                                                              method = "bray")
phyloseq.samples.bacteria.nitifiers.ra.copper.bray

#Converts a distance matrix (matrix of pairwise distances) into a kernel matrix for microbiome data. 
bray_kernel_nit_all <- D2K(
  as.matrix(phyloseq.samples.bacteria.nitifiers.ra.copper.bray))

#Jaccard distances
phyloseq.samples.bacteria.nitifiers.ra.copper.jac <- vegdist(t(phyloseq.samples.bacteria.nitifiers.ra.copper@otu_table), 
                                                             method = "jaccard")
phyloseq.samples.bacteria.nitifiers.ra.copper.jac

jac_kernel_nit_all <- D2K(
  as.matrix(phyloseq.samples.bacteria.nitifiers.ra.copper.jac))

# Covariates (date_num)
datenum_covariate_nit_all <- phyloseq.samples.bacteria.nitifiers.copper.all.df%>%
  pull(Date_num)%>%
  as.matrix
datenum_covariate_nit_all

#Id for repeated measures 
Enclosure_d_nit_mirkat_all <- phyloseq.samples.bacteria.nitifiers.copper.all.df %>%
  mutate(Enclosure = ifelse(Enclosure == "P1", "0", "1"))%>%
  pull(Enclosure)%>%
  as.numeric()
Enclosure_d_nit_mirkat_all

#Run model - dont think it works. just two enclosures and enclosure id is the same as the outcome 
mirkat_enclosure_nit_glm <- GLMMMiRKAT(
  y = phyloseq.samples.bacteria.nitifiers.copper.all.df$Copper_level_mg_L,
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
copper_cont_BC_nit_beta_div <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.scrs) + 
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
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
    aes(x = -1, y = 0.6, label = "PERMANOVA\nCopper levels (mg/L)"),
    hjust = 0.5, vjust = -0.5,
    size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.4, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0.13, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
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
copper_cont_BC_nit_beta_div_2 <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.scrs) + 
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
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = -0.65, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = -0.88, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = -0.2, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
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
phyloseq.samples.bacteria.nitifiers.ra.bray.scrs_2 <- phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
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
copper_cont_BC_nit_beta_div_3 <- ggplot(phyloseq.samples.bacteria.nitifiers.ra.bray.scrs_2) + 
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
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.55, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "H21"),
            aes(x = -1, y = 0.3, 
                label = paste("R² = ", round(R2_adonis_copper_cont_nit_H21* 100, 1), "%",
                              "\np = ", round(pvalue_adonis_copper_cont_nit_H21, 3))),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black")+
  #P1
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
              filter(Enclosure == "P1"),
            aes(x = -1, y = 0.1, label = "PERMANOVA\nCopper levels (mg/L)"),
            hjust = 0.5, vjust = -0.5,
            size = 8, colour = "black", fontface = "bold")+
  geom_text(data = phyloseq.samples.bacteria.nitifiers.ra.bray.scrs %>%
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
phyloseq.samples.nitrifiers.copper <- subset_samples(phyloseq.samples.nitrifiers, !is.na(Copper_level_mg_L))
phyloseq.samples.nitrifiers.copper <- prune_taxa(taxa_sums(phyloseq.samples.nitrifiers.copper)> 0, 
                                                 phyloseq.samples.nitrifiers.copper)

phyloseq.samples.nitrifiers.copper.ra <- transform_sample_counts(phyloseq.samples.nitrifiers.copper, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data

#make DF from metadata
phyloseq.samples.nitrifiers.copper.df <- as(phyloseq.samples.nitrifiers.copper@sam_data, "data.frame") %>%
  mutate(
    # Convert Collection_Date to Date if it isn’t already
    Collection_Date = as.Date(Collection_Date),
    # Create a new month-year factor for strata
    Collection_Month = format(Collection_Date, "%Y-%m"),
    Collection_Month = factor(Collection_Month)  # convert to factor for PERMANOVA
  )

phyloseq.samples.nitrifiers.copper.df 

###BRAY CURTIS#####
phyloseq.samples.nitrifiers.ra.bray <- vegdist(t(phyloseq.samples.nitrifiers.ra@otu_table), method = "bray")
phyloseq.samples.nitrifiers.ra.bray

# ###ORDINATION####
# set.seed(98)
# phyloseq.samples.nitrifiers.ra.bray.ord <- metaMDS(phyloseq.samples.nitrifiers.ra.bray, k=3, try = 50, 
#                                                  trymax = 1000,
#                                                  autotransform = F)
# #### ADDING CENTROIDS FOR PLOTTING
# ## BC
# #Simple ordination plot
# phyloseq.samples.nitrifiers.ra.bray.plot <- ordiplot(phyloseq.samples.nitrifiers.ra.bray.ord$points)
# 
# #Now, extract coordinates
# phyloseq.samples.nitrifiers.ra.bray.scrs <- scores(phyloseq.samples.nitrifiers.ra.bray.plot, display = "sites")
# #Add metadata to coordinates
# phyloseq.samples.nitrifiers.ra.bray.scrs <- cbind(as.data.frame(phyloseq.samples.nitrifiers.ra.bray.scrs),
#                                                 Copper_level_mg_L = phyloseq.samples.nitrifiers.copper.df$Copper_level_mg_L, 
#                                                 SampleID = phyloseq.samples.nitrifiers.copper.df$SampleID, 
#                                                 Collection_Date = phyloseq.samples.nitrifiers.copper.df$Collection_Date,
#                                                 Collection_Month = phyloseq.samples.nitrifiers.copper.df$Collection_Month, 
#                                                 Enclosure = phyloseq.samples.nitrifiers.copper.df$Enclosure)
# ##Calculate centroids according to Enclosure
# phyloseq.samples.nitrifiers.ra.bray.cent <- aggregate(cbind(MDS1,MDS2) ~ Enclosure, 
#                                                     data = phyloseq.samples.nitrifiers.ra.bray.scrs, 
#                                                     FUN = mean) 
# #Merge centroids with coordinates and metadata
# phyloseq.samples.nitrifiers.ra.bray.segs <- merge(phyloseq.samples.nitrifiers.ra.bray.scrs, 
#                                                 setNames(phyloseq.samples.nitrifiers.ra.bray.cent, c("Enclosure","cMDS1","cMDS2")),
#                                                 by = 'Enclosure', 
#                                                 sort = F)
# 
# 
# ##PERMANOVA###
# set.seed(98)
# nit_enclosure_copper_BC_adonis  <- adonis2(phyloseq.samples.nitrifiers.ra.bray ~ Enclosure +
#                                          Copper_level_mg_L, 
#                                        phyloseq.samples.nitrifiers.copper.df, 
#                                        strata = phyloseq.samples.nitrifiers.copper.df$Collection_Month,
#                                        by = "margin",
#                                        permutations = 9999)
# nit_enclosure_copper_BC_adonis #4.8% of the variation is due to Enclosure, p = 1e-04
# #3.5% of the variation is due to copper levels, p = 8e-04
# 
# 
# #With interaction
# set.seed(98)
# nit_enclosure_copper_interac_BC_adonis  <- adonis2(phyloseq.samples.nitrifiers.ra.bray ~ Enclosure*Copper_level_mg_L, 
#                                                phyloseq.samples.nitrifiers.copper.df , 
#                                                strata = phyloseq.samples.nitrifiers.copper.df$Collection_Month,
#                                                by = "margin",
#                                                permutations = 9999)
# nit_enclosure_copper_interac_BC_adonis
# #Enclosure:Copper_level_mg_L interaction significant (p = 1e-04). 7.3% of variation
# 
# ##PERMDISPS
# #ENCLOSURE#
# # Run the betadisper function, average distance to centroid
# bray.enclosure_nit.disp <- betadisper(phyloseq.samples.nitrifiers.ra.bray, 
#                                   phyloseq.samples.nitrifiers.copper.df$Enclosure)
# bray.enclosure_nit.disp
# ##Then test by permuting
# set.seed(98)
# bray.enclosure_nit.permdisp <- permutest(bray.enclosure_nit.disp, permutations = 9999)
# bray.enclosure_nit.permdisp ##NS, p = 0.65
# 
# #COPPER#
# # Run the betadisper function, average distance to centroid
# bray.copper_nit.disp <- betadisper(phyloseq.samples.nitrifiers.ra.bray, 
#                                phyloseq.samples.nitrifiers.copper.df$Copper_level_mg_L)
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
# enclosure_nit_BC_beta_div <- ggplot(phyloseq.samples.nitrifiers.ra.bray.segs) + theme_bw() +
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
# enclosure_nit_BC_beta_div_2 <- ggplot(phyloseq.samples.nitrifiers.ra.bray.segs) + 
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
# fit <- envfit(phyloseq.samples.nitrifiers.ra.bray.ord, 
#               phyloseq.samples.nitrifiers.copper.df$Copper_level_mg_L, 
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


#********RESISTOME******#######
##Import count tables#######
ARGcounts <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplusv4_0gf_05gc/SNPconfirmed_AMR_0gf_0.5qc_analytic_matrix.csv')
colnames(ARGcounts)
rownames(ARGcounts)

##Annotations - downloaded from AMRplusplus#######
tax.table.ARG <- readr::read_csv(
  '/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/AMRplusplusv4_0gf_05gc/Annotations/megares_annotations_v4.00.csv')
#There is a multi-compount category, fixing it to "Multi-compound"

#Edit to make compatible with phyloseq
tax.table.ARG <- tax.table.ARG%>%
  column_to_rownames(var = "header")%>% #make "header' rownames so it matches with the OTU table rownames
  rename_with(~ str_to_title(.))%>% #want the annotation cols with the first letter capitalized
  as.matrix() ##make into matrix so it is compatible with tax_table function from phyloseq



##OTU table ########
otu_table_ARG <- ARGcounts%>%
  column_to_rownames(var = "gene_accession")
  
#Rename the columns to not have '_alignment_sorted', then pick samples to keep
otu_table_ARG <- otu_table_ARG%>%
  rename_with(~ str_remove(.x, "_alignment_sorted$"))%>%
  select(-all_of(dropping_samples))

ncol(otu_table_ARG) #208 = 205 samples + 3 mocks 

#Same number of samples in the microbiome count table as the resistome count table? 
ncol(otu_table_ARG) #208 samples (205 samples + 3 mocks). There are 205 (instead of 207 as the microbiome) since P1_0308 and  "H21_0109" did not make it through with any alignment counts
length(sample_names(phyloseq.samples.bacteria)) #204 samples

setdiff(sample_names(OTU), colnames(otu_table_ARG)) #"H21_0109", "P1_0308" are in the kraken microbiome count table (original, without filtering those with low counts) but not on the resistome table
setdiff(colnames(otu_table_ARG), sample_names(phyloseq.samples.bacteria)) #The resistome count table includes H21_0120, plus the mock communities
setdiff(sample_names(phyloseq.samples.bacteria), colnames(otu_table_ARG)) #


##PHYLOSEQ####
OTU_ARG <-phyloseq::otu_table(otu_table_ARG, 
                              taxa_are_rows = TRUE)
TAX_ARG <-phyloseq::tax_table(tax.table.ARG)
phyloseq_ARG <- phyloseq(OTU_ARG, TAX_ARG, sampledata_phyloseq) 
phyloseq_ARG ##8481 taxa and 208 (samples + mock communities)

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU_ARG), metadata$SampleID) #No

##PREPROCESSING ####
phyloseq_ARG #8481 taxa and 208 samples 
min(sample_sums(phyloseq_ARG)) # 156 (H21_0120)
max(sample_sums(phyloseq_ARG)) # 1,197,153  (ZymoMock1a_S142) 
mean(sample_sums(phyloseq_ARG)) #226362.71
median(sample_sums(phyloseq_ARG)) #11275.5
sort(sample_sums(phyloseq_ARG))

###Zymo#####
### pulling out samples from ZYMOs and EB, NTC and those samples with low OTUs
phyloseq_ARG.controls <- subset_samples(phyloseq_ARG, 
                                        grepl("NTC|EB|Zymo", sample_names(phyloseq_ARG)))
phyloseq_ARG.controls <- prune_taxa(taxa_sums(phyloseq_ARG.controls) > 0, phyloseq_ARG.controls) 
phyloseq_ARG.controls #1907 taxa and 3 samples (only zymos)

###Samples#####
##New phyloseq of just samples
phyloseq_ARG.samples <- subset_samples(phyloseq_ARG, 
                                       !grepl("NTC|EB|Zymo", sample_names(phyloseq_ARG)))
phyloseq_ARG.samples #8481 taxa and 205 samples 
#Taking out those with low counts:
sort(sample_sums(phyloseq_ARG.samples)) #H21_0120 has low counts, dropping it here as i did with kraken taxonomy
phyloseq_ARG.samples <- prune_samples(sample_sums(phyloseq_ARG.samples) > 1000, phyloseq_ARG.samples) 
phyloseq_ARG.samples <- prune_taxa(taxa_sums(phyloseq_ARG.samples) > 0, phyloseq_ARG.samples)
phyloseq_ARG.samples #7182 taxa and 204 samples (dropped H21_0120)
sort(sample_sums(phyloseq_ARG.samples)) #OK

#####H21####
phyloseq_ARG.samples_H21 <- subset_samples(phyloseq_ARG.samples, Enclosure == "H21")
phyloseq_ARG.samples_H21 <- prune_taxa(taxa_sums(phyloseq_ARG.samples_H21) > 0, 
                                             phyloseq_ARG.samples_H21)
phyloseq_ARG.samples_H21 #6228 taxa and 87 samples (From the naive systm, H21_0109 did not even have enough counts for alignment, and H21_0120 was filtered out as well for ow counts here)
range(phyloseq_ARG.samples_H21@sam_data$Collection_Date)#OK, "2023-10-09" through "2024-02-20"

#####P1####
phyloseq_ARG.samples_P1 <- subset_samples(phyloseq_ARG.samples, Enclosure == "P1")
phyloseq_ARG.samples_P1 <- prune_taxa(taxa_sums(phyloseq_ARG.samples_P1) > 0, 
                                            phyloseq_ARG.samples_P1)
phyloseq_ARG.samples_P1 #8269 taxa and 117 samples (P1_0308 did not have any alignments)
range(phyloseq_ARG.samples_P1@sam_data$Collection_Date)#OK, "2023-11-14" through "2024-04-02"

#Am I missing samples that have a kraken taxonomic classification?
setdiff(sample_names(phyloseq_ARG.samples), sample_names(phyloseq.samples.bacteria)) #Nope

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

###RA PLOT CLASS LEVEL METALS####
phyloseq_ARG.samples.ra.class #4 types and 216 samples

#Out of this overall communities object, select only metals 
phyloseq_ARG.samples.ra.type.metals <- subset_taxa(phyloseq_ARG.samples.ra.class, 
                                                    Type == "Metals") 

phyloseq_ARG.samples.ra.type.metals # 16 classes in the Metals type and 216 samples 

#Melt to plot 
phyloseq_ARG.samples.ra.type.metals.melt <- psmelt(phyloseq_ARG.samples.ra.type.metals)

#Color palette
palette_metals_gene_groups <- distinctColorPalette(length(unique(phyloseq_ARG.samples.ra.type.metals.melt$Class)))
metals_names <- unique(phyloseq_ARG.samples.ra.type.metals.melt$Class)# Create a named vector for the palette, where the names correspond to phlyum names
palette_metals_gene_groups <- setNames((palette_metals_gene_groups)[1:length(palette_metals_gene_groups)], metals_names)

#Plot  
RA_enclosures_ARG_metals.plot <- ggplot(phyloseq_ARG.samples.ra.type.metals.melt,
                                        aes(x=factor(Date_num), y= Abundance, fill = Class)) +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  scale_fill_manual(values = palette_metals_gene_groups)+
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
RA_enclosures_ARG_metals.plot


###RA PLOT CLASS LEVEL COPPER ####
phyloseq_ARG.samples.ra.class #58 classes and 216 samples

#Out of this overall communities object, select only copper 
phyloseq_ARG.samples.ra.class.copper <- subset_taxa(phyloseq_ARG.samples.ra.class, 
                                                          Class == "Copper resistance") 

phyloseq_ARG.samples.ra.class.copper # 1 class and 216 samples 

#Melt to plot 
phyloseq_ARG.samples.ra.class.copper.melt <- psmelt(phyloseq_ARG.samples.ra.class.copper)

####SUPPLEMENTARY TABLE 6 #####
stable6 <-  phyloseq_ARG.samples.ra.class.copper.melt %>%
  mutate(System = ifelse(grepl("H21", Enclosure), "Naive", "Established"))%>%
  group_by(System, Date_num_phase_abbrv) %>%
  summarise(
    mean_abun = round(mean(Abundance, na.rm = TRUE), 2),
    sd_abun = round(sd(Abundance, na.rm = TRUE),3), 
    min_abun = round(min(Abundance, na.rm = TRUE), 2), 
    max_abun = round(max(Abundance, na.rm = TRUE), 2),
    .groups = "drop_last") %>%
  arrange(System,  Date_num_phase_abbrv, desc(mean_abun))%>%
  mutate(`Mean Relative Abundance (%) ± SD` = 
           paste0(
             mean_abun,
             " ± ", 
             sd_abun
           ), 
         `Min - Max Relative Abundance (%)` = paste0(
           min_abun,
           " - ",
           max_abun
         ))%>%
  rename(Phase = Date_num_phase_abbrv)%>%
  mutate(Class = "Copper Resistance")%>%
  select(Class, System, Phase, 
         `Mean Relative Abundance (%) ± SD`, 
         `Min - Max Relative Abundance (%)`
         )
stable6
write_xlsx(stable6, 
           "SupplementaryTable6.xlsx")


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



##Create color palette
palette_copper_gene_groups <- distinctColorPalette(length(unique(phyloseq_ARG.samples.ra.group.copper_filt.melt$Group)))
copper_filt_names <- unique(phyloseq_ARG.samples.ra.group.copper_filt.melt$Group)# Create a named vector for the palette, where the names correspond to phlyum names
palette_copper_gene_groups <- setNames((palette_copper_gene_groups)[1:length(palette_copper_gene_groups)], copper_filt_names)
#order.filt.palette <- unname(alphabet2())
palette_copper_gene_groups$'Others <0.1% RA' <- "grey95"

# #Just top most abundant groups to include in legend
# top_taxa_copper_groups <- top_taxa_legend(phyloseq_ARG.samples.ra.group.copper.melt, taxlevel = "Group", n = 23)

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


###Poster VERO ADVISORY COUNCIL#######
alpha_div_wq_date_num_factor_other_metadata_VERO <- ggplot(alpha_div_wq_time_long%>%
                                                        filter(Index %in% c("Copper_level_mg_L",
                                                                            "Shannon",
                                                                            "Ammonia_mg_L"
                                                                            #"Chlorine_mg_L", 
                                                                            #"Alkalinity_mg_L",
                                                                            # "Temperature_F",
                                                                            # "pH_spu",
                                                                            # "Salinity_ppt"
                                                        )),
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
             # #switch = "y", 
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive",
                                      "Copper_level_mg_L"= "Copper\n(mg/L)",
                                      "Shannon" = "Shannon", 
                                      "Ammonia_mg_L" = "Ammonia\n(mg/L)"
                                      # "Temperature_F"= "Temperature (F)",
                                      # "Salinity_ppt" = "Salinity (ppt)", 
                                      # "pH_spu" = "pH (spu)"
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
    legend.position = c(0.95, 1.5),
    legend.justification = c(0, 1),
    legend.title.position = "left",
    legend.direction = "horizontal",
    legend.text = element_text(size = 15, angle = 45,
                               hjust = 0.5, vjust = 0.5),
    legend.title = element_text(size = 20, face = "bold"),
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
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks.x = element_line(colour = "black", linewidth = 1),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    plot.title = element_text(colour = "black", size = 48, face = "bold"))
alpha_div_wq_date_num_factor_other_metadata_VERO

alpha_div_wq_date_num_factor_other_metadata_VERO_1 <- alpha_div_wq_date_num_factor_other_metadata_VERO +
  theme(
    plot.title = element_blank(),
    legend.position = c(0.92, 1.35),  
    axis.text.x = element_blank(),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 20),
    strip.text.y = element_text(size = 30))

RA_enclosures_ARG_copper_genegroup.plot_8 <- RA_enclosures_ARG_copper_genegroup.plot +
  labs(fill = "Copper ARG Group")+
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme(
    legend.position = c(1.06, 0.5), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 21, face = "bold"),
    axis.title.y = element_text(size = 23),
    axis.text.y = element_text(size = 25),
    axis.title.x = element_blank(),
    axis.text.x = element_blank())

RA_family_enclosures_nit_plot_datenum_4 <- RA_family_enclosures_nit_plot_datenum + 
  theme(
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 22, face = "bold"),
    legend.key.size = unit(0.5, "cm"),
    legend.position = c(1.049, 0.55), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 23),
    axis.text.y = element_text(size = 25)) 

RA_enclosures_nitrifiers_only_species.plot_4 <- RA_enclosures_nitrifiers_only_species.plot +
  scale_fill_manual(
    values = palette_nitrifiers_only_species,
    breaks = top_nitrifying_only_species,
    labels = function(x) str_wrap(x, width = 20),
    drop = FALSE
  )+
  theme(
    legend.position = c(1.045, 0.46), 
    axis.title.x = element_text(size = 28),
    axis.text.x = element_text(size = 28),
    axis.title.y = element_text(size = 23),
    axis.text.y = element_text(size = 25), 
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18, face = "bold")
  )


figure1_poster_vero <-
  alpha_div_wq_date_num_factor_other_metadata_VERO_1 /
  RA_enclosures_ARG_copper_genegroup.plot_8 /
  RA_family_enclosures_nit_plot_datenum_4 /
  RA_enclosures_nitrifiers_only_species.plot_4 +
  plot_layout(heights = c(0.6, 0.5, 0.5, 0.5))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 30, face = "bold"))
figure1_poster_vero

ggsave("/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/Posters/Figure1_veroposter.png", 
       figure1_poster_vero, 
       device = "png", 
       dpi = 600, 
       height = 18, 
       width = 28)

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
phyloseq.samples.bacteria.family.ra #6755 families 

#Out of this overall communities object, select only rhodobacteraceae
phyloseq.samples.bacteria.family.ra.rhodobacteraceae <- subset_taxa(phyloseq.samples.bacteria.family.ra, 
                                                                    Family == "Rhodobacteraceae") 
phyloseq.samples.bacteria.family.ra.rhodobacteraceae <- subset_samples(phyloseq.samples.bacteria.family.ra.rhodobacteraceae, 
                                                                       sample_sums(phyloseq.samples.bacteria.family.ra.rhodobacteraceae) > 0)
phyloseq.samples.bacteria.family.ra.rhodobacteraceae #Rhodobacteraceae (1 taxa) in 216 samples 

#Melt to long format
phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt <- psmelt(phyloseq.samples.bacteria.family.ra.rhodobacteraceae)

#Join with copper resistance class RA
#Copper resistance class melted object (long format)
phyloseq_ARG.samples.ra.class.copper.melt_join <- phyloseq_ARG.samples.ra.class.copper.melt %>%
  rename(Abundance_copper = Abundance)
phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt <- phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt %>%
  rename(Abundance_rhodobacteraceae = Abundance)

#Join 
rhodobacteraceae_copper.ra.melt <- merge(phyloseq.samples.bacteria.family.ra.rhodobacteraceae.melt%>%
                                           select(SampleID, Abundance_rhodobacteraceae, Enclosure, Date_num_colors_systems), 
                                         phyloseq_ARG.samples.ra.class.copper.melt_join%>%
                                           select(SampleID, Abundance_copper, Enclosure), 
                                         by = c("SampleID", "Enclosure"))

#Plot - separate by phase
copperARG_rhodobacteraceae_correlation_plot<- ggplot(rhodobacteraceae_copper.ra.melt,
                                                     aes(y = Abundance_copper, 
                                                         x = Abundance_rhodobacteraceae,
                                                         fill = Enclosure, 
                                                         color = Enclosure)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Enclosure)) +
  facet_wrap(~ Enclosure, 
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
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.6,
    alpha = 0.3
  )+
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

####FIGURE 2 #######
figure2 <- copperARG_rhodobacteraceae_correlation_plot
ggsave("Figure2.png", 
       figure2, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 16)

#Plot - separate by phase
copperARG_rhodobacteraceae_correlation_phase_plot<- ggplot(rhodobacteraceae_copper.ra.melt,
                                                       aes(y = Abundance_copper, 
                                                           x = Abundance_rhodobacteraceae,
                                                          fill = Date_num_colors_systems, 
                                                          color = Date_num_colors_systems)) +
  geom_point(size = 7, shape = 18, 
             aes(color = Date_num_colors_systems)) +
  facet_wrap(Date_num_colors_systems~ Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive", 
                                      "L_Naive" = "L", 
                                      "T1_Naive" = "T1", 
                                      "E1_Naive" = "E1",
                                      "T2_Naive" = "T2",
                                      "T3_Naive" = "T3", 
                                      "E2_Naive" = "E2", 
                                      "P_Naive" = "P", 
                                      "L_Established" = "L", 
                                      "T1_Established" = "T1", 
                                      "E_Established" = "E", 
                                      "P_Established" = "P"
                                      )))+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.15)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.03)))+
  labs(x = expression(italic(Rhodobacteraceae) ~ "RA (%)"), 
       y = "Copper ARG Class RA (%)") +
  theme_bw() +
  geom_smooth(
    aes(group = Date_num_colors_systems,
        color = Date_num_colors_systems),
    method = "lm",
    se = TRUE,
    linewidth = 0.6,
    alpha = 0.3
  )+
  stat_cor(method = "spearman",
           label.x.npc = "center",
           label.y.npc = "bottom",
           color = "black",
           size = 8) +
  scale_color_manual(values = phases_naive_established_palette)+
  scale_fill_manual(values = phases_naive_established_palette)+
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
copperARG_rhodobacteraceae_correlation_phase_plot


####Correlation of Nitrosopumilaceae with RA of Copper resistance class RA#########
phyloseq.samples.bacteria.family.ra #6755 families 

#Out of this overall communities object, select only Nitrosopumilaceae
phyloseq.samples.bacteria.family.ra.nitrosopumilaceae <- subset_taxa(phyloseq.samples.bacteria.family.ra, 
                                                                     Family == "Nitrosopumilaceae") 
phyloseq.samples.bacteria.family.ra.nitrosopumilaceae <- subset_samples(phyloseq.samples.bacteria.family.ra.nitrosopumilaceae, 
                                                                        sample_sums(phyloseq.samples.bacteria.family.ra.nitrosopumilaceae) > 0)
phyloseq.samples.bacteria.family.ra.nitrosopumilaceae #Nitrosopumilaceae (1 taxa) in 216 samples 

#Melt to long format
phyloseq.samples.bacteria.family.ra.nitrosopumilaceae.melt <- psmelt(phyloseq.samples.bacteria.family.ra.nitrosopumilaceae)

#Join with copper resistance class RA
#Copper resistance class melted object (long format)
phyloseq_ARG.samples.ra.class.copper.melt_join <- phyloseq_ARG.samples.ra.class.copper.melt %>%
  rename(Abundance_copper = Abundance)
phyloseq.samples.bacteria.family.ra.nitrosopumilaceae.melt <- phyloseq.samples.bacteria.family.ra.nitrosopumilaceae.melt %>%
  rename(Abundance_nitrosopumilaceae = Abundance)

#Join 
nitrosopumilaceae_copper.ra.melt <- merge(phyloseq.samples.bacteria.family.ra.nitrosopumilaceae.melt%>%
                                           select(SampleID, Abundance_nitrosopumilaceae, Enclosure), 
                                         phyloseq_ARG.samples.ra.class.copper.melt_join%>%
                                           select(SampleID, Abundance_copper, Enclosure), 
                                         by = c("SampleID", "Enclosure"))


#Plot
copperARG_nitrosopumilaceae_correlation_plot<- ggplot(nitrosopumilaceae_copper.ra.melt,
                                                     aes(y = Abundance_copper, 
                                                         x = Abundance_nitrosopumilaceae,
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
  labs(x = expression(italic(Nitrosopumilaceae) ~ "RA (%)"), 
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
copperARG_nitrosopumilaceae_correlation_plot
ggsave("copperARG_rhodobacteraceae_correlation_plot.png", 
       copperARG_rhodobacteraceae_correlation_plot, 
       device = "png", 
       dpi = 600, 
       height = 10, 
       width = 16)


###RA PLOT GROUP FIEF, PERR, RUVB, TUPC####
#These were groups i found to be present mostly in late phase MAGs but absent in early ones, gonna look at how the reads look
phyloseq_ARG.samples.ra.group #1372 taxa and 216 samples

grep("tup", ignore.case = T, value = T, metals_df$Group)
#"TUPC" "TUPA" "TUPB"
grep("fie", ignore.case = T, value = T, metals_df$Group)
#"FIEF"
grep("per", ignore.case = T, value = T, metals_df$Group)
#"PERR" "PERO" "PER" 
grep("ruv", ignore.case = T, value = T, metals_df$Group)
#"RUVB" "RUBVM" 

#Out of this overall communities object, select only copper 
phyloseq_ARG.samples.ra.group.AOA.dif <- subset_taxa(
  phyloseq_ARG.samples.ra.group,
  Group %in% c("TUPC", "TUPA", "TUPB", "FIEF", "PERR", "PERO", "PER", "RUVB", "RUVM"))
phyloseq_ARG.samples.ra.group.AOA.dif # 7 gene groups and 216 samples 

# #Merge low abundance groups into one
# phyloseq_ARG.samples.ra.group.AOA.dif.filt <- merge_low_abundance_ARG_ra(phyloseq_ARG.samples.ra.group.copper,
#                                                                         "Enclosure",
#                                                                         level = "Group",
#                                                                         threshold = 0.1)
# phyloseq_ARG.samples.ra.group.AOA.dif.filt #19 gene groups with mean RA above 0.1%

# #Melt to plot 
# phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt <- psmelt(phyloseq_ARG.samples.ra.group.AOA.dif.filt)
# #For some reason the function is leaving the Others group as NA, fixing that...
# phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt <-
#   phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt %>%
#   mutate(Group = ifelse(is.na(Group), "Others <0.1% RA", Group))
# 
# #Factoring so "Others" group is last
# phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt <- phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt%>%
#   mutate(Group = factor(Group,
#                         levels = c(setdiff(Group,
#                                            unique(grep("Others", Group, value = TRUE))),
#                                    unique(grep("Others", Group, value = TRUE)))))##Factoring the Group column so that "Others.." is the last category

#Melt to long format
phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt <- psmelt(phyloseq_ARG.samples.ra.group.AOA.dif)

# ##Create color palette
# palette_copper_gene_groups <- distinctColorPalette(length(unique(phyloseq_ARG.samples.ra.group.copper_filt.melt$Group)))
# copper_filt_names <- unique(phyloseq_ARG.samples.ra.group.copper_filt.melt$Group)# Create a named vector for the palette, where the names correspond to phlyum names
# palette_copper_gene_groups <- setNames((palette_copper_gene_groups)[1:length(palette_copper_gene_groups)], copper_filt_names)
# #order.filt.palette <- unname(alphabet2())
# palette_copper_gene_groups$'Others <0.1% RA' <- "grey95"

# #Just top most abundant groups to include in legend
# top_taxa_copper_groups <- top_taxa_legend(phyloseq_ARG.samples.ra.group.copper.melt, taxlevel = "Group", n = 23)

#Plot
RA_enclosures_ARG_AOA_dif_genegroup.plot <- ggplot(phyloseq_ARG.samples.ra.group.AOA.dif.filt.melt,
                                                  aes(x=factor(Date_num), y= Abundance, fill = Group)) +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.0015,0)) +
  # scale_fill_manual(values = palette_copper_gene_groups,
  #                   labels = function(x) str_wrap(x, width = 15)) +
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
    legend.position = "right",
    # legend.position = c(1.08, 0.5),  # x, y inside plot
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
RA_enclosures_ARG_AOA_dif_genegroup.plot 
