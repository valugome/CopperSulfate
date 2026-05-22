#TAXONOMIC CLASSIFICATIONS OF COPPER-ARG-ALIGNING READS##########
#Importing data from kraken output nt_core - counts_copper_reads will be classified reads#### 
counts_copper_reads <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/Copper_ARG_reads/By_sample/kraken_analytic_matrix_samples_copper.conf_0.1.csv')
##Separating into taxonomy levels
counts_copper_reads_separated_tax <- counts_copper_reads %>%
  separate(taxa, 
           into = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|",#splits the strings by the "|" symbol.
           fill = "right") #fill = "right", missing components are added as "NA" to the right (last columns) instead of to the left 

##Extracting just taxonomy  (columns 1:8 are taxonomy, the rest are counts_copper_reads)
tax.table_copper_reads<- counts_copper_reads_separated_tax %>%
  dplyr::select(1:8)
tax.table_copper_reads


##Filling up actual NAs and string "NA"s in the taxonomy table
filled_taxonomy_copper_reads <- fill_taxonomy(tax.table_copper_reads) ##apply the function to the taxonomy table
anyNA(filled_taxonomy_copper_reads) ##OK, no NAs now 
grep("^NA$", filled_taxonomy_copper_reads, value = T) ##OK, no "NA" strings now

##Now, to add the row names as "OTU1, OTU2, etc..." for phyloseq_copper_reads later on
filled_taxonomy_2_copper_reads<- filled_taxonomy_copper_reads %>%
  mutate(OTU = paste0("OTU", 1:nrow(filled_taxonomy_copper_reads))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##convert into matrix for phyloseq_copper_reads
filled_taxonomy_2_copper_reads

#Make a csv file for the kraken taxonomy table
write.csv(filled_taxonomy_2_copper_reads,
          "kraken_taxonomy_copper_reads.csv",
          row.names = F)

###OTU table #####
otu_table_copper_reads <- counts_copper_reads[, -1]%>% #Excludes the first column (taxonomy)
  mutate(OTU = paste0("OTU", 1:nrow(counts_copper_reads))) %>% ##add "OTU#" column
  column_to_rownames(var= "OTU") %>% ##Make OTU column into row names
  as.matrix() ##make into matrix so it is compatible with otu_table function from phyloseq_copper_reads
otu_table_copper_reads

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

##Making into phyloseq_copper_reads-compatible object
sampledata_phyloseq_copper_reads <- metadata %>%
  mutate(rows = SampleID)%>%
  column_to_rownames(var= "rows") %>%##Make sampleID column into row names, so they match sample_names() with OTU and TAX
  sample_data(metadata) ##use phyloseq_copper_reads function sample_data() to make metadata into phyloseq_copper_reads sample data object

#phyloseq_copper_reads####
#Make phyloseq_copper_reads object
OTU_copper_reads <-phyloseq::otu_table(otu_table_copper_reads, taxa_are_rows = TRUE)
TAX_copper_reads <-phyloseq::tax_table(filled_taxonomy_2_copper_reads)
phyloseq_copper_reads <- phyloseq(OTU_copper_reads, 
                                  TAX_copper_reads, 
                                  sampledata_phyloseq_copper_reads)

#Am I missing metadata for any sampleIDs?
setdiff(sample_names(OTU_copper_reads), metadata$SampleID) #Yes, "H21_1021a" and "H21_1021b". These will not be included in the phyloseq_copper_reads object.

#Are there samples in metadata that don't have sequencing data?
setdiff(metadata$SampleID, sample_names(OTU_copper_reads)) #Yes, "P1_1126", "P1_1203", 
#"P1_1216", "P1_1225", "P1_1228", "P1_0104", "P1_0112", "P1_0115", 
#"P1_0205", "P1_0212", "P1_0218", "H21_1021", "H21_1122b", Zymo mocks and NTC

#COLOR PALETTES#####
enclosure.palette <- c("H21" = "#fc8d62",  
                       "P1"  = "#8da0cb" )

#PREPROCESSING ####
phyloseq_copper_reads #2102 taxa and 223 samples 
      
##Selecting only Bacteria/Archaea####
phyloseq_copper_reads.bacteria <- subset_taxa(phyloseq_copper_reads, Domain=="Archaea" | Domain=="Bacteria")
phyloseq_copper_reads.bacteria #1934 taxa and 223 samples

##Selecting only viruses######
phyloseq_copper_reads.viruses <- subset_taxa(phyloseq_copper_reads, Domain=="Viruses")
taxanames_viruses <- c("Kingdom", "Realm", "Phylum", "Class", "Order", "Family", "Genus", "Species") ##they have a different classification system, updating it here
colnames(phyloseq_copper_reads.viruses@tax_table) <- taxanames_viruses #replacing col names of the tax_table for new ones
colnames(phyloseq_copper_reads.viruses@tax_table) #OK taxonomy ranks
phyloseq_copper_reads.viruses #3 taxa and 223 samples

##Selecting only eukaryota #####
phyloseq_copper_reads.eukaryota <- subset_taxa(phyloseq_copper_reads, Domain=="Eukaryota")
colnames(phyloseq_copper_reads.eukaryota@tax_table) ##These are OK taxonomy ranks
phyloseq_copper_reads.eukaryota #165 taxa and 223 samples

#WORKING ON BACTERIA/ARCHAEA ONLY####
# some QC checks of the "classified" reads per samples
min(sample_sums(phyloseq_copper_reads.bacteria)) # 22 (P1_0211)
max(sample_sums(phyloseq_copper_reads.bacteria)) # 511  (H21_1110) 
mean(sample_sums(phyloseq_copper_reads.bacteria)) #165.9507
median(sample_sums(phyloseq_copper_reads.bacteria)) # 143
sort(sample_sums(phyloseq_copper_reads.bacteria))

###Dropping samples before copper dosage started####
#What's the range of dates 
range(phyloseq_copper_reads@sam_data$Collection_Date)#"2023-04-20" "2024-04-30"

#Actual Copper dosing starts from 10/09/2023 (sampling ends on 03/02/2024) for naive system
#Actual copper dosing starts from 11/14/2023 (sampling ends on 04/30/2024) for established system
phyloseq_copper_reads.dates <- subset_samples(phyloseq_copper_reads, Collection_Date > "2023-10-05")
phyloseq_copper_reads.dates #2102 taxa and 218 samples
phyloseq_copper_reads.dates <- prune_taxa(taxa_sums(phyloseq_copper_reads.dates) > 0, phyloseq_copper_reads.dates) 
phyloseq_copper_reads.dates #2049 taxa and 218 samples
setdiff(sample_names(phyloseq_copper_reads), sample_names(phyloseq_copper_reads.dates)) #Dropped "H21_0912" "H21_1005" "P1_0420"  "P1_0427"  "P1_0504" 

#Also, have H21_1202a and 1202b. These were taken on Dec 2nd, 2023. A is before backwash, B is afterbackwash. 
#Have more reliable metadata for H21_1202a (before backwash)
#Also, Sample H21_0120 was dropped from the overall analysis, dropping it here too
phyloseq_copper_reads.dates <- subset_samples(phyloseq_copper_reads.dates,
                                              !(SampleID %in% c("H21_1202b", "H21_0120")))
phyloseq_copper_reads.dates #2049 taxa and 216 samples (Dopped H21_1202b and H21_0120)
phyloseq_copper_reads.dates <- prune_taxa(taxa_sums(phyloseq_copper_reads.dates) > 0, phyloseq_copper_reads.dates) 
phyloseq_copper_reads.dates #2040 taxa and 216 samples
setdiff(sample_names(phyloseq_copper_reads), sample_names(phyloseq_copper_reads.dates)) 

###H21####
phyloseq_copper_reads.dates_H21 <- subset_samples(phyloseq_copper_reads.dates, Enclosure == "H21")
phyloseq_copper_reads.dates_H21 <- prune_taxa(taxa_sums(phyloseq_copper_reads.dates_H21) > 0, 
                                                  phyloseq_copper_reads.dates_H21)
phyloseq_copper_reads.dates_H21 #1565 taxa and 92 samples
range(phyloseq_copper_reads.dates_H21@sam_data$Collection_Date)#OK, now "2023-10-09" through "2024-03-02"

###P1####
phyloseq_copper_reads.dates_P1 <- subset_samples(phyloseq_copper_reads.dates, Enclosure == "P1")
phyloseq_copper_reads.dates_P1 <- prune_taxa(taxa_sums(phyloseq_copper_reads.dates_P1) > 0, 
                                                 phyloseq_copper_reads.dates_P1)
phyloseq_copper_reads.dates_P1 #920 taxa and 124 samples
range(phyloseq_copper_reads.dates_P1@sam_data$Collection_Date)#OK, now "2023-11-14" through "2024-04-30"


##NITRIFYING TAXA####
nitrifiers_all_copper <- subset_taxa(phyloseq_copper_reads.dates, Family == "Nitrosomonadaceae" | # AOB; some, plus a new one!
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
nitrifiers_all_copper #49 taxa and 216 samples
nitrifiers_copper <- subset_samples(nitrifiers_all_copper, sample_sums(nitrifiers_all_copper) > 0)
nitrifiers_copper #49 taxa and 142 samples (only in 142 samles were copper ARG reads classified within Nitrifying communities)

##QC checks again
min(sample_sums(phyloseq_copper_reads.dates)) #23 (P1_0211)
max(sample_sums(phyloseq_copper_reads.dates)) #524 (H21_1012) 
mean(sample_sums(phyloseq_copper_reads.dates)) #168.4213
median(sample_sums(phyloseq_copper_reads.dates)) #142.5
sort(sample_sums(phyloseq_copper_reads.dates)) 

#COMPARING CLASSIFIED READS BY KRAKEN#######
kraken_unclassified_reads <- readr::read_csv('/Users/valerialugo/Library/CloudStorage/OneDrive-TexasA&MUniversity/Documents/Projects/CuSo4/AMR_counts/Copper_ARG_reads/By_sample/unclassifieds_kraken_analytic_matrix_samples_copper.conf_0.1.csv')

#Metadata file phyloseq_copper_reads.dates
phyloseq_copper_reads.dates_metadata <- data.frame(phyloseq_copper_reads.dates@sam_data)

#Filtering just samples included in phyloseq_copper_reads.dates, adding metadata, calculating percentage classified
kraken_unclassified_copper_reads_samples_metadata <- kraken_unclassified_reads %>%
  dplyr::right_join(phyloseq_copper_reads.dates_metadata, by = "SampleID")%>%
  rename(Kraken2_Input_PairedEnd_Reads = Total, 
         Kraken2_Unclassified_PairedEnd_Reads = NumberUnclassified, 
         Kraken2_Unclassified_Percentage_Reads = PercentUnclassified)%>%
  mutate(Kraken2_Classified_Percentage_Reads = (100 - Kraken2_Unclassified_Percentage_Reads))
nrow(kraken_unclassified_copper_reads_samples_metadata) #Ok, 216 samples

###Kraken2 Classified Percentages Established vs Naive####
kraken2_classified_copper_read_percentages_P1vsH21<- ggplot(kraken_unclassified_copper_reads_samples_metadata, 
                               aes(x = Enclosure, 
                                   y= Kraken2_Classified_Percentage_Reads, 
                                   color = Enclosure, 
                                   fill = Enclosure)) +
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
kraken2_classified_copper_read_percentages_P1vsH21

###Kraken2 Classified Percentages Established and Naive over time#### 
kraken2_classified_copper_read_percentages_P1andH21_overtime <- 
  ggplot(kraken_unclassified_copper_reads_samples_metadata,
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
kraken2_classified_copper_read_percentages_P1andH21_overtime

#COMPARING SAMPLE SUMS (CLASSIFIED READS FROM KRAKEN)#######
##ALL TAXA#####
sample.sums_copper <- sample_sums(phyloseq_copper_reads.dates) #making a sample sums object
phyloseq_copper_reads.dates.samplessums.df <- cbind(phyloseq_copper_reads.dates@sam_data, 
                                                    sample.sums_copper) #combining sample sums with metaphyloseq_copper_reads
phyloseq_copper_reads.dates.samplessums.df
phyloseq_copper_reads.dates.samplessums.df$sampleID <- rownames(phyloseq_copper_reads.dates.samplessums.df) ##making a sampleID column


###Established vs Naive####
bacteria_archaea_samplesums_copper_reads_P1vsH21<- ggplot(phyloseq_copper_reads.dates.samplessums.df, 
                                  aes(x = Enclosure, y= sample.sums_copper, 
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
bacteria_archaea_samplesums_copper_reads_P1vsH21

##Stats 
wilcox_test(phyloseq_copper_reads.dates.samplessums.df, sample.sums_copper~Enclosure) #S. p = 1.98e-13

###Established and Naive over time####
bacteria_archaea_samplesums_P1andH21_overtime<- ggplot(phyloseq_copper_reads.dates.samplessums.df, 
                                    aes(x = factor(Date_num), 
                                        y= sample.sums_copper, 
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

##NITRIFIERS#####
sample.sums.copper.nit <- sample_sums(nitrifiers_copper) #making a sample sums object
nitrifiers_copper.samplesums.df <- cbind(nitrifiers_copper@sam_data, 
                                      sample.sums.copper.nit) #combining sample sums with metadata
nitrifiers_copper.samplesums.df
nitrifiers_copper.samplesums.df$SampleID <- rownames(nitrifiers_copper.samplesums.df) ##making a sampleID column

###Established vs Naive####
nitrifier_bacteria_archaea_samplesums_P1vsH21<- ggplot(nitrifiers_copper.samplesums.df, 
                                  aes(x = Enclosure, y= sample.sums.copper.nit, 
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
wilcox_test(nitrifiers_copper.samplesums.df, sample.sums.copper.nit~Enclosure) #S. p = 0.00228

#RELATIVE ABUNDANCE####
any(sample_sums(phyloseq_copper_reads.dates)== 0) ## no samples with 0 OTUs
phyloseq_copper_reads.dates.ra <- transform_sample_counts(phyloseq_copper_reads.dates, 
                                                        function(x) x/sum(x)*100) ##Relative abundance from normalized data

##CLASSIFICATION PERCENTAGES AT DIFFERENT TAXONOMIC LEVELS####
###PHYLUM######
phyloseq_copper_reads.dates_phylum.ra <- tax_glom(phyloseq_copper_reads.dates.ra, taxrank = "Phylum", NArm = F) 
phyloseq_copper_reads.dates_phylum.ra #79 phyla and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_phylum.ra)[, "Phylum"])) #79 phyla (so No duplicates)

Unknown_phylum_abundance <- phyloseq_copper_reads.dates_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Phylum, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_phylum_abundance ##0.187% abundance by Unknown Phyla

Unclassified_phylum_abundance <- phyloseq_copper_reads.dates_phylum.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Phylum, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_phylum_abundance ##31.2% abundance by Unclassified Phyla

Classified_phylum_abundance <- phyloseq_copper_reads.dates_phylum.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Phylum, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_phylum_abundance ##68.7% abundance by Classified Phyla

#How many unclassified?
phyloseq_copper_reads.dates_phylum.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_phylum.ra)[, "Phylum"]),
  phyloseq_copper_reads.dates_phylum.ra)
phyloseq_copper_reads.dates_phylum.unclassified.ra #8 unclassified Phyla

#How many unknown?
phyloseq_copper_reads.dates_phylum.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq_copper_reads.dates_phylum.ra)[, "Phylum"]),
  phyloseq_copper_reads.dates_phylum.ra)
phyloseq_copper_reads.dates_phylum.unknown.ra #16 "unknown" Phyla

#Keep just classified Phyla
phyloseq_copper_reads.dates_phylum.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_phylum.ra)[, "Phylum"]),
  phyloseq_copper_reads.dates_phylum.ra)
phyloseq_copper_reads.dates_phylum.classified.ra ##55 classified (not unknown or unclassified) Phyla

###CLASS#####
phyloseq_copper_reads.dates_class.ra <- tax_glom(phyloseq_copper_reads.dates.ra, taxrank = "Class", NArm = F) 
phyloseq_copper_reads.dates_class.ra #176 taxa and 216 samples

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_class.ra)[, "Class"])) #176 classes (so No duplicates)

Unknown_class_abundance <- phyloseq_copper_reads.dates_class.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Class, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_class_abundance #0.471% Abundance by Unknown classes

Unclassified_class_abundance <- phyloseq_copper_reads.dates_class.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Class, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_class_abundance ##43.4% Abundance by Unclassified Classes

Classified_class_abundance <- phyloseq_copper_reads.dates_class.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Class, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_class_abundance ##56.2% Abundance by Classified classes

#How many unclassified?
phyloseq_copper_reads.dates_class.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_class.ra)[, "Class"]),
  phyloseq_copper_reads.dates_class.ra)
phyloseq_copper_reads.dates_class.unclassified.ra #28 unclassified classes

#How many unknown?
phyloseq_copper_reads.dates_class.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq_copper_reads.dates_class.ra)[, "Class"]),
  phyloseq_copper_reads.dates_class.ra)
phyloseq_copper_reads.dates_class.unknown.ra #51 "unknown" classes

#Keep just classified Classes
phyloseq_copper_reads.dates_class.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_class.ra)[, "Class"]),
  phyloseq_copper_reads.dates_class.ra)
phyloseq_copper_reads.dates_class.classified.ra #97 classified classes

###ORDER######
phyloseq_copper_reads.dates_order.ra <- tax_glom(phyloseq_copper_reads.dates.ra, taxrank = "Order", NArm = F) 
phyloseq_copper_reads.dates_order.ra #334 orders

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_order.ra)[, "Order"])) #334 orders (no duplicates)

Unknown_order_abundance <- phyloseq_copper_reads.dates_order.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Order, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_order_abundance ##0.652% abundance by Unknown Orders

Unclassified_order_abundance <- phyloseq_copper_reads.dates_order.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Order, ignore.case = TRUE)) %>%  # Filter unclassified phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_order_abundance ##57.5% abundance by Unclassified Orders

Classified_order_abundance <- phyloseq_copper_reads.dates_order.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Order, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_order_abundance ##41.9% abundance by Classified orders

#How many unclassified?
phyloseq_copper_reads.dates_order.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_order.ra)[, "Order"]),
  phyloseq_copper_reads.dates_order.ra)
phyloseq_copper_reads.dates_order.unclassified.ra #47 unclassified orders

#How many unknown?
phyloseq_copper_reads.dates_order.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq_copper_reads.dates_order.ra)[, "Order"]),
  phyloseq_copper_reads.dates_order.ra)
phyloseq_copper_reads.dates_order.unknown.ra #73 "unknown" orders

#Keep just classified Orders
phyloseq_copper_reads.dates_order.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_order.ra)[, "Order"]),
  phyloseq_copper_reads.dates_order.ra)
phyloseq_copper_reads.dates_order.classified.ra #214 classified orders
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_order.classified.ra)[, "Order"])) ##214 classified orders (unique - without duplicates)

###FAMILY######
phyloseq_copper_reads.dates_family.ra <- tax_glom(phyloseq_copper_reads.dates.ra, taxrank = "Family", NArm = F) 
phyloseq_copper_reads.dates_family.ra #569 families

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_family.ra)[, "Family"])) #569 taxa (5 duplicates)

Unknown_family_abundance <- phyloseq_copper_reads.dates_family.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_family_abundance #0.949% abundance by Unknown Families

Unclassified_family_abundance <- phyloseq_copper_reads.dates_family.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Family, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_family_abundance ##66.3% abundance by Unclassified Families

Classified_family_abundance <- phyloseq_copper_reads.dates_family.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Family, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_family_abundance ##32.8% abundance by Classified Families

#How many unclassified?
phyloseq_copper_reads.dates_family.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_family.ra)[, "Family"]),
  phyloseq_copper_reads.dates_family.ra)
phyloseq_copper_reads.dates_family.unclassified.ra #93 unclassified families

#How many unknown?
phyloseq_copper_reads.dates_family.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq_copper_reads.dates_family.ra)[, "Family"]),
  phyloseq_copper_reads.dates_family.ra)
phyloseq_copper_reads.dates_family.unknown.ra #96 "unknown" families

#Keep just classified Families
phyloseq_copper_reads.dates_family.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_family.ra)[, "Family"]),
  phyloseq_copper_reads.dates_family.ra)
phyloseq_copper_reads.dates_family.classified.ra #380 classified families

###GENUS ######
phyloseq_copper_reads.dates_genus.ra <- tax_glom(phyloseq_copper_reads.dates.ra, taxrank = "Genus", NArm = F) 
phyloseq_copper_reads.dates_genus.ra #1105 genera
#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_genus.ra)[, "Genus"])) #1105 taxa (77 duplicates)

Unknown_genus_abundance <- phyloseq_copper_reads.dates_genus.ra %>%
  psmelt()%>%
  filter(grepl("unknown", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unknown_sum = sum(OTU_Abundance))   # Sum across OTUs
Unknown_genus_abundance ##1.32%  abundance by unknown genera

Unclassified_genus_abundance <- phyloseq_copper_reads.dates_genus.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Genus, ignore.case = TRUE)) %>%  # Filter unknown<tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance)) #Sum across OTUs
Unclassified_genus_abundance ##73.6% abundance by unclassified genera

Classified_genus_abundance <- phyloseq_copper_reads.dates_genus.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Genus, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_genus_abundance ##25.1% abundance by Classified Genera

#How many unclassified?
phyloseq_copper_reads.dates_genus.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_genus.ra)[, "Genus"]),
  phyloseq_copper_reads.dates_genus.ra)
phyloseq_copper_reads.dates_genus.unclassified.ra #189 unclassified genera

#How many unknown?
phyloseq_copper_reads.dates_genus.unknown.ra <- prune_taxa(
  grepl("unknown", phyloseq::tax_table(phyloseq_copper_reads.dates_genus.ra)[, "Genus"]),
  phyloseq_copper_reads.dates_genus.ra)
phyloseq_copper_reads.dates_genus.unknown.ra #128 "unknown" genera

#Keep just classified Genera
phyloseq_copper_reads.dates_genus.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_genus.ra)[, "Genus"]),
  phyloseq_copper_reads.dates_genus.ra)
phyloseq_copper_reads.dates_genus.classified.ra #788 classified genera
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_genus.classified.ra)[, "Genus"])) ##788 classified genera (unique - without duplicates)


###SPECIES######
phyloseq_copper_reads.dates.ra ##2040 Species- OTUs
phyloseq_copper_reads.dates_species.ra <- phyloseq_copper_reads.dates.ra
phyloseq_copper_reads.dates_species.ra #2040 Species

#Are there duplicates? 
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_species.ra)[, "Species"])) #2040 species (no duplicates)

Unclassified_species_abundance <- phyloseq_copper_reads.dates_species.ra %>%
  psmelt()%>%
  filter(grepl("unclassified", Species, ignore.case = TRUE)) %>%  # Filter unclassified <tax_rank> phyla
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance), .groups = "drop") %>%  # Mean abundance per OTU
  summarize(Unclassified_sum = sum(OTU_Abundance))   # Sum across OTUs
Unclassified_species_abundance ##82.1% abundance by unclassified species

Classified_species_abundance <- phyloseq_copper_reads.dates_species.ra %>%
  psmelt()%>%
  filter(!grepl("unclassified|unknown", Species, ignore.case = TRUE)) %>%  ##filter out unclassified and unknown
  group_by(OTU) %>%  #group by OTU
  summarize(OTU_Abundance = mean(Abundance, .groups = "drop")) %>%  # Mean abundance per OTU
  summarize(Classified_sum = sum(OTU_Abundance))   # Sum across OTUs
Classified_species_abundance ##17.9% abundance by Classified Species

#How many unclassified?
phyloseq_copper_reads.dates_species.unclassified.ra <- prune_taxa(
  grepl("unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_species.ra)[, "Species"]),
  phyloseq_copper_reads.dates_species.ra)
phyloseq_copper_reads.dates_species.unclassified.ra #399 unclassified species
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_species.unclassified.ra)[, "Species"])) ##399 unclassified species (unique - without duplicates)

#Keep just classified Genera
phyloseq_copper_reads.dates_species.classified.ra <- prune_taxa(
  !grepl("unknown|unclassified", phyloseq::tax_table(phyloseq_copper_reads.dates_species.ra)[, "Species"]),
  phyloseq_copper_reads.dates_species.ra)
phyloseq_copper_reads.dates_species.classified.ra #2641 classified species
length(unique(phyloseq::tax_table(phyloseq_copper_reads.dates_species.classified.ra)[, "Species"])) ##1641 classified species (unique - without duplicates)

#RELATIVE ABUNDANCE######
#ALL TAXA######
## ORDER #####
phyloseq_copper_reads.dates_order.ra #334 taxa and 216 samples 

#Grouping the low abundance orders into one category
phyloseq_copper_reads.dates.order.filt <- merge_low_abundance_grouped_ra(phyloseq_copper_reads.dates_order.ra, 
                                                                        "Enclosure", 
                                                                        level = "Order", 
                                                                        threshold = 0.5)
phyloseq_copper_reads.dates.order.filt #26 orders over 0.5% mean RA
phyloseq_copper_reads.dates.order.filt.melt <- psmelt(phyloseq_copper_reads.dates.order.filt)%>%
  mutate(Order = factor(Order, 
                         levels = c(setdiff(Order, 
                                            unique(grep("Others", Order, value = TRUE))), 
                                    unique(grep("Others", Order, value = TRUE)))))##Factoring the Order column so that "Others.." is the last category
levels(phyloseq_copper_reads.dates.order.filt.melt$Order) ##ok

##Create color palette - using the main order palette from the main order plot
#Taxa in the main plot (orders) and taxa just in copper reads
taxa_orders_main <- unique(phyloseq_copper_reads.dates.order.filt.melt$Order)
taxa_orders_copper <- unique(phyloseq_copper_reads.dates.order.filt.melt$Order)

#Find shared and unique taxa
shared_taxa_orders_plots <- intersect(taxa_orders_main, taxa_orders_copper)
new_taxa_orders_plots <- setdiff(taxa_orders_copper, taxa_orders_main)

#Make shared_taxa into character
shared_taxa_orders_plots <- as.character(shared_taxa_orders_plots)

#Start with colors for shared taxa (from the original palette)
order_copper_named_palette <- order_named_palette[shared_taxa_orders_plots]
order_copper_named_palette
#Generate new colors for taxa not in the original palette
new_colors_orders <- setNames(
  scales::hue_pal()(length(new_taxa_orders_plots)),
  new_taxa_orders_plots)

#Combine them
order_copper_named_palette <- c(order_copper_named_palette, new_colors_orders)
order_copper_named_palette <- order_copper_named_palette[taxa_orders_copper]

##Apply the function to obtain top orders (n=15)
top_orders_copper <- top_taxa_legend(phyloseq_copper_reads.dates.order.filt.melt, 
                              taxlevel = "Order", n = 15)
top_orders_copper

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_order_overall_taxa_copper_reads_plot <- ggplot(phyloseq_copper_reads.dates.order.filt.melt,
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
  scale_fill_manual(values = order_copper_named_palette,
                    breaks = top_orders_copper,
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
RA_order_overall_taxa_copper_reads_plot

## FAMILY #####
phyloseq_copper_reads.dates_family.ra #569 families and 216 samples 

phyloseq_copper_reads.dates.family.filt <- merge_low_abundance_grouped_ra(phyloseq_copper_reads.dates_family.ra, 
                                                                             "Enclosure", 
                                                                             level = "Family", 
                                                                             threshold = 0.5)
phyloseq_copper_reads.dates.family.filt #30 families over 0.5% mean RA
phyloseq_copper_reads.dates.family.filt.melt <- psmelt(phyloseq_copper_reads.dates.family.filt)%>%
  mutate(Family = factor(Family, 
                        levels = c(setdiff(Family, 
                                           unique(grep("Others", Family, value = TRUE))), 
                                   unique(grep("Others", Family, value = TRUE)))))##Factoring the Family column so that "Others.." is the last category
levels(phyloseq_copper_reads.dates.family.filt.melt$Family) ##ok

##Create color palette - based on families within the same order
palette_family_level_copper_df <- phyloseq_copper_reads.dates.family.filt.melt %>% 
  arrange(Order, Family) %>%   # ensure consistent shading order
  group_by(Order) %>%
  mutate(
    base_color = order_copper_named_palette[Order],
    shade = seq(0.02, 0.5, length.out = n()),  # avoid extremes
    color = darken(base_color, amount = shade)
  ) %>%
  ungroup()
palette_family_level_copper_df

#Set up final palette
family_named_copper_palette <- setNames(
  palette_family_level_copper_df$color,
  palette_family_level_copper_df$Family)
family_named_copper_palette
family_named_copper_palette$'Others <0.5% RA' <- "grey95"

##Apply the function to obtain top familys (n=15)
top_families_copper <- top_taxa_legend(phyloseq_copper_reads.dates.family.filt.melt, 
                              taxlevel = "Family", n = 15)
top_families_copper

### Plot RA at the order level with days since start (Date_num) as factor #####
RA_family_overall_taxa_copper_reads_plot <- ggplot(phyloseq_copper_reads.dates.family.filt.melt,
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
  scale_fill_manual(values = family_named_copper_palette,
                    breaks = top_families_copper,
                    labels = function(x) str_wrap(x, width = 20)) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.06, 0.5),  # x, y inside plot
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"), 
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
RA_family_overall_taxa_copper_reads_plot


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
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold"))
RA_enclosures_ARG_copper_genegroup.plot_7 <- RA_enclosures_ARG_copper_genegroup.plot +
  theme(
    legend.position = c(1.08, 0.5), 
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 17),
    axis.text.y = element_text(size = 17), 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank())

#Final plot
figure_alpha_overall_div_family_time_taxa_copper_reads <- 
  alpha_div_wq_date_num_factor_other_metadata_6 /
  RA_family_enclosures_overall_plot_datenum_6 /
  RA_enclosures_ARG_copper_genegroup.plot_7 /
  RA_family_overall_taxa_copper_reads_plot +
  plot_layout(heights = c(0.9, 0.6, 0.6, 0.6))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alpha_overall_div_family_time_taxa_copper_reads
ggsave("figure_alpha_overall_div_family_time_taxa_copper_reads.png", 
       figure_alpha_overall_div_family_time_taxa_copper_reads, 
       device = "png", 
       dpi = 600, 
       height = 23, 
       width = 26)


#NITRIFIERS WITHIN THE OVERALL COMMUNITY#####
## FAMILY #######
phyloseq_copper_reads.dates_family.ra #569 families

##Which families are nitrifiers_copper? 
nitrifiers_copper.melt <- psmelt(nitrifiers_copper)
unique(nitrifiers_copper.melt$Family) #"Nitrobacteraceae", "Nitrospinaceae" , "Nitrospiraceae", "Nitrosopumilaceae"     
#"Chromatiaceae", "Gallionellaceae", "Ectothiorhodospiraceae"

#Out of this overall communities object, select only nitrifiers_copper 
phyloseq_copper_reads.dates_family.ra.nitrifiers <- subset_taxa(phyloseq_copper_reads.dates_family.ra, 
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

#When I plot, need to have all the ticks for all samples, so I will not subset out the samples without reads classified as nitrifiers
# phyloseq_copper_reads.dates_family.ra.nitrifiers <- subset_samples(phyloseq_copper_reads.dates_family.ra.nitrifiers, 
#                                                                        sample_sums(phyloseq_copper_reads.dates_family.ra.nitrifiers) > 0)
# phyloseq_copper_reads.dates_family.ra.nitrifiers #7 nitrifying families in 142 samples 


#Melt to plot 
phyloseq_copper_reads.dates_family.ra.nitrifiers.melt <- psmelt(phyloseq_copper_reads.dates_family.ra.nitrifiers)


##Add a column for which type of  ammonia-nitrate group (AOA, AOB, NOB)
phyloseq_copper_reads.dates_family.ra.nitrifiers.melt <- phyloseq_copper_reads.dates_family.ra.nitrifiers.melt %>%
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
#Using the same one as for the overall nitrifyinf families plot
palette_nitrifiers_family

##Apply the function to obtain top orders (n=15)
top_nitrifying_families_copper <- top_taxa_legend(phyloseq_copper_reads.dates_family.ra.nitrifiers.melt, 
                                           n = 5)
top_nitrifying_families_copper <- c("Nitrosopumilaceae", #AOA
                             "Nitrobacteraceae",#NOB
                             "Nitrospinaceae", #NOB
                             "Nitrospiraceae",#NOB
                             "Chromatiaceae")#AOB
top_nitrifying_families_copper

#Plot
figure_alpha_overall_div_nit_time_taxa_copper_reads <- ggplot(phyloseq_copper_reads.dates_family.ra.nitrifiers.melt,
                                                aes(x=factor(Date_num), y= Abundance, fill = Family)) +
  theme_minimal() +
  labs(y= "Relative Abundance (%)", x = "Days") +
  facet_grid(~Enclosure, 
             scales = "free",
             labeller = as_labeller(c("P1" = "Established",
                                      "H21" = "Naive")))+
  geom_bar(stat = "summary", color = "black") +
  scale_y_continuous(expand = c(0.0015,0,0.03,0)) +
  # geom_vline(data = line_breaks_phases,
  #            aes(xintercept = Date_num),
  #            linetype = "dashed",
  #            color = "black",
  #            alpha = 0.8) +
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
                    breaks = top_nitrifying_families_copper
  ) +
  guides(fill=guide_legend(title.position="top", ncol = 1))+
  theme_bw()+
  theme(
    #legend.position = "right",
    legend.position = c(1.08, 0.5),  # x, y inside plot
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
    axis.title.x = element_text(colour = "black", size = 22),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text.y = element_text(colour = "black", size = 20),
    axis.ticks = element_line(colour = "black", linewidth = 0.8))
figure_alpha_overall_div_nit_time_taxa_copper_reads

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
    legend.position = c(1.08, 0.5), 
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 20)) 
RA_enclosures_ARG_copper_genegroup.plot_3 <- RA_enclosures_ARG_copper_genegroup.plot + 
  theme(
    legend.title = element_text(size = 16),
    axis.text.x = element_blank(),
    axis.title.x = element_blank()
    ) 

#Final plot
figure_alphadiv_nit_familyRA_ARG_taxa_copper_reads <-
  alpha_div_nit_wq_date_num_factor_other_metadata_3 /
  RA_family_enclosures_nit_plot_datenum_3 /
  RA_enclosures_ARG_copper_genegroup.plot_3 /
  figure_alpha_overall_div_nit_time_taxa_copper_reads+
  plot_layout(heights = c(1.1, 0.4, 0.4, 0.4))+
  plot_annotation(
    tag_levels = "A") &
  theme(plot.tag = element_text(size = 24, face = "bold"))
figure_alphadiv_nit_familyRA_ARG_taxa_copper_reads

#Saving figure
ggsave("figure_alphadiv_nit_familyRA_ARG_taxa_copper_reads.png", 
       figure_alphadiv_nit_familyRA_ARG_taxa_copper_reads, 
       device = "png", 
       dpi = 600, 
       height = 20, 
       width = 26)
