
########  CLEMENT DASTAN MLAY ####
#######MSC TAXONOMIC CLASSIFICATION ######
install.packages("remotes")
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
library(vegan)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ape)
library(plyr)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(scales)
library(data.table)
library(microbiome)
library(dplyr)
library(phyloseq)
library(DT) #for interactive tables
library(microbiomeutilities)

# Exploration of data from pavian (TSV files)
getwd()
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/taxonomic_classification")
kisite_data <- read.table("kisite_report.krk-report-220330.tsv", header = T)
kisite_krk <- read_table("kisite_report.krk")
colnames(kisite_krk)
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/kraken_reports")
kisite <- read_table("kisite_report.krk")
malindi <- read_table("malindi_report.krk")
mombasa <- read_table("mombasa_report.krk")
specie <- read_table("all_species.bracken")
s
########################## MSC TAXONOMIC CLASSIFICATION WORK ###########################

library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
#install.packages("ggtext")
library("ggtext")

# set working directory
setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/taxonomic_classification/kraken_reports")

# Import data to a phyloseq object
transcriptome_data <- import_biom("mydata.biom")
class(transcriptome_data)

#################### taxa data #######################
# reformat the names in each taxa (remain with the fourth character in each name)
transcriptome_data@tax_table@.Data <- substring(
  transcriptome_data@tax_table@.Data, 4)

# provide names in each taxonomic level
colnames(transcriptome_data@tax_table@.Data) <- c("Kingdom", 
                                                  "Phylum", "Class", "Order", "Family", 
                                                  "Genus", "Species")

# extract taxa data from phyloseq object
tax_data <- transcriptome_data@tax_table@.Data

y <- transcriptome_data@tax_table@.Data
x <- substring(transcriptome_data@tax_table@.Data, 4)

################## OTU data ##########################
# Obtain OTU of each site from phylosq object
otu_data <- transcriptome_data@otu_table@.Data
z <- t(otu_table(transcript_bact))
rarecurve(t(otu_table(transcript_bact)), step=50, label = FALSE)
rarecurve(otu_table(transcript_bact), step = 100, ylab = "Observed otu",
          xlab = "Number of reads", main = "Rarefaction curves", label = FALSE)

# subset data to obtain information of bacteria only 
transcript_bact <- subset_taxa(transcriptome_data, Kingdom == "Bacteria")
transcript_bact
sample_sums(transcript_bact)

# Raefy data without replacement
transcript.rarefied = rarefy_even_depth(transcript_bact, rngseed=1, sample.size=0.9*min(sample_sums(transcript_bact)), replace=F)
transcript.phylum = tax_glom(percentages, taxrank="Phylum", NArm=FALSE)
plot_bar(transcript.phylum, fill="Phylum")

# Determine alpha diversity
plot_richness(transcript.rarefied)
rich = estimate_richness(transcript_bact)
rich

plot_richness(physeq = transcript_bact, 
              measures = c("Simpson")) 

plot_richness(physeq = transcript_bact, 
              measures = c("Shannon"))

plot_richness(physeq = transcript_bact, 
              measures = c("InvSimpson")) 

plot_richness(physeq = transcript_bact, 
              measures = c("Fisher")) 

plot_richness(physeq = transcript_bact, 
              measures = c("Chao1")) 

# obtain number of reads from each site
deep <- data.frame(Samples = sample_names(transcript_bact),
                     Reads = sample_sums(transcript_bact))


# Plot a graph to show number of reads
ggplot(data = deep, mapping = aes(x = Samples,y = Reads)) +
geom_col()

# Determine unidentified reads 
summary(transcriptome_data@tax_table@.Data== "")

# Convert reads count into relative abundance (percentages)
head(transcript_bact@otu_table@.Data)
percentages <- transform_sample_counts(transcript_bact, 
                                       function(x) x*100 / sum(x))
head(percentages@otu_table@.Data)

# Determine beta diversity
distanceMethodList

# PCoA plot using the unweighted UniFrac as distance
wunifrac_dist = phyloseq::distance(transcript_bact, method="unifrac", weighted=F)

# generate an object where the distances between our samples will be allocated 
meta.ord <- ordinate(physeq = percentages, method = "NMDS", 
                     distance = "bray")

# Plot ordination to see beta diversity
plot_ordination(physeq = percentages, ordination = meta.ord)

# group all the OTUs that have the same taxonomy at a certain taxonomic rank
glom <- tax_glom(percentages, taxrank = 'Phylum')
glom1 <- tax_glom(percentages, taxrank = "Class")
View(glom@tax_table@.Data)
hig_glom <- tax_glom(hig_perc, taxrank = 'Phylum')
# melts phyloseq objects into a data frame so as manipulate 
# Relative data
percentages <- psmelt(glom)
percentages1 <- psmelt(glom1)
str(percentages)
hig_perc <- psmelt(hig_glom)
percentages
# Original Data
raw <- tax_glom(physeq = metatranscriptomics_data, taxrank = "Phylum")
raw.data <- psmelt(raw)
str(raw.data)

raw.plot <- ggplot(data=raw.data, aes(x=Sample, y=Abundance, fill=Phylum))+ 
geom_bar(aes(), stat="identity", position="stack")
raw.plot

class_plot <- ggplot(data=percentages1, aes(x=Sample, y=Abundance, fill=Class))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  geom_col()+
  labs(x=NULL,
       y= "Relative abundance %")+
  theme_classic()+
  theme(legend.text = element_text(face = "italic"),
        legend.key.size = unit(10, "pt"))
class_plot

rel.plot <- ggplot(data=percentages, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  geom_col()+
  labs(x=NULL,
       y= "Relative abundance %")+
  theme_classic()+
  theme(legend.text = element_text(face = "italic"),
        legend.key.size = unit(10, "pt"))
rel.plot
percentages

# Plit the most dominant taxa
percentage2 <- percentages
percentage2$Phylum[percentage2$Abundance[1:15]]
unique(percentage2$Phylum[percentage2$Abundance > 2.201315e-02]) 
#g.cyanos$Genus[g.cyanos$Abundance < 10] <- "Genera < 10.0 abund"

rel.plot <- ggplot(data=hig_perc, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  geom_col()+
  labs(x=NULL,
       y= "Mean relative abundance %")+
  theme_classic()+
  theme(legend.text = element_text(face = "italic"),
        legend.key.size = unit(10, "pt"))
rel.plot


#top5P.names = sort(tapply(taxa_sums(percentages), tax_table(percentages)[, "Phylum"], sum), TRUE)[1:5]
#Cut down the physeq.tree data to only the top 10 Phyla
#top5P = subset_taxa(physeq.tree, Phylum %in% names(top5P.names))
#Plot
#plot_bar(top5P, x="AgeGroup", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


#plot_bar(physeq.tree, x="AgeGroup", fill="Phylum") + 
 # geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#rel1.plot <- ggplot(data=percentages, aes(x=Sample, y=Abundance, fill=Phylum))+ 
# geom_bar(aes(), stat="identity", position="stack")
#rel1.plot

# DATA FOR MEGAN
raw_met <- tax_glom(physeq = metatranscriptomics_data, taxrank = "Species")
raw_met_data <- psmelt(raw_met)
str(raw_met_data)
write.table(raw_met_data, file = "all_sites_data", row.names = TRUE, col.names = TRUE)
table1 <- read.table("all_sites_data")  
  
library(vegan)
data(dune)
ord <- decorana(dune)
View(metatranscriptomics_data)

ord2 <- metatranscriptomics_data


tab <- matrix(c(7501735640,2606169198), ncol=2, byrow=TRUE)
colnames(tab) <- c('reads before','reads after')
rownames(tab) <- ('reads')
tab <- as.table(tab)
7501735640 - 2606169198

install.packages("plotrix")
library(plotrix)
slices <- c(70, 30 )
lbls <- c("filtered(70%)", "host free(30%)")
pie3D(slices,labels=lbls,explode=0.5)


# TAXONOMIC CLASSIFICATION OF INDIVIDUAL SAMPLES

setwd("/Users/cley/Desktop/Msc Project/Msc_data_analysis/taxonomic_classification")

DATA <- import_biom("all_samples.biom")
DATA@tax_table@.Data <- substring(DATA@tax_table@.Data, 4)
colnames(DATA@tax_table@.Data) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

DATA <- subset_taxa(DATA, Kingdom == "Bacteria")
DATA
plot_richness(physeq = DATA, 
              measures = c("Observed","Simpson", "Shannon")) 

READS <- data.frame(Samples = sample_names(DATA),
                            Reads = sample_sums(DATA))

READS

# Plot a graph to show number of reads
ggplot(data = READS, mapping = aes(x = Samples,y = Reads)) +
  geom_col()

READSperc <- transform_sample_counts(DATA, function(x) x*100 / sum(x))

DATA.ord <- ordinate(physeq = READSperc, method = "NMDS", 
                     distance = "bray")

DATA.glom <- tax_glom(READSperc, taxrank = 'Phylum')
READSperc1 <- psmelt(DATA.glom)

DATA.plot <- ggplot(data=READSperc1, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")+
  geom_col()+
  labs(x=NULL,
       y= "Relative abundance %")+
  theme_classic()+
  theme(legend.text = element_text(face = "italic"),
        legend.key.size = unit(10, "pt"))
DATA.plot
