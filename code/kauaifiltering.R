setwd("~/Documents/GitHub/Diet_Analysis/kauaidiet/kauaidiet/raw")


library(tidyverse)
install.packages("qiime2R")
library(qiime2R)
library(phyloseq)
library(ggplot2)

####Bringing in Data#####

#Read in OTU table
otu <- read.table(file = "kauaiotu_table.txt", header = TRUE, sep = '\t', fill = T, check.names = FALSE)
head(otu)
#read in taxonomy table
tax <- read.table(file = "kauaitaxonomy.tsv", sep = '\t', header = TRUE)
head(tax)

####Import into Phyloseq####

otu_table_in <- otu
rownames(otu_table_in) <- otu_table_in$OTU
otu_table_in$OTU <- NULL
otu_table_in <- as.matrix(otu_table_in)
head(otu_table_in)

rm(otu)

#fix taxonomy file and then export
tax <- as.matrix(read.delim("kauaitaxonomy.tsv", row.names = 1, fill=TRUE, stringsAsFactors = FALSE))
tax <- gsub("k__","", tax) #to remove unwanted "d_" in taxa names
tax <- gsub("p__","", tax)
tax <- gsub("c__","", tax)
tax <- gsub("o__","", tax)
tax <- gsub("f__","", tax)
tax <- gsub("g__","", tax)
tax <- gsub("s__","", tax)
write.table(tax, file="kauaitaxonomy_edited.tsv", sep = "\t", 
            col.names = TRUE, row.names = TRUE)
taxonomy <- read.table("kauaitaxonomy_edited.tsv", sep = "\t", header=T,  check.names = FALSE)
taxonomy <- as.matrix(taxonomy)
head(taxonomy)

#read in metadata
metadata <- read.table("dietmetamerged.txt", row.names = 1, header=T,
                       sep = "\t")
head(metadata)

# Import all as phyloseq objects
OTU <- otu_table(otu_table_in, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy)
META <- sample_data(metadata)

taxa_names(TAX)
taxa_names(OTU)


sample_names(OTU)
sample_names(META)
head(OTU)
head(META)

length(colnames(OTU))

#Put it all together in phyloseq object 
kps <- phyloseq(OTU, TAX, META) #this actually has all islands still
kps

# check meta data and tax table
head(sample_data(kps))
summary(sample_data(kps))
head(tax_table(kps))

# remove and clean up unnecessary objects
rm(OTU, TAX, META) 
rm(otu_table_in, taxonomy)

# return some metadata sampleID to a column for easier use later in adonis
sample_data(kps)$SampleID <- rownames(sample_data(kps))

####Check blank samples####

# calculate depth
metadata <- as.data.frame(sample_data(kps)) 
metadata$LibrarySize <- sample_sums(kps)
metadata <- metadata[order(metadata$LibrarySize),]
head(metadata)
metadata$Index <- seq(nrow(metadata))

meta.mat <- as.matrix(metadata) #so i can look at the lowest library sizes and see if i need to throw any out
#there are several samples below 1000 reads

# plot library sizes by Sample Type
ggplot() + 
  geom_point(data=metadata[metadata$Sample.or.Control == "Sample" ,], 
             aes(x=Index, y=LibrarySize), size=2.5, color= "lightskyblue") +
  geom_point(data=metadata[metadata$Sample.or.Control == "Control" ,], 
             aes(x=Index, y=LibrarySize), size=2.5, color= "darkorange3") +
  theme_bw() +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 10)), 
        axis.title.y = element_text(size=16, margin = margin(r = 10)),
        axis.text.y  = element_text(size=13, color= "black"), 
        axis.text.x  = element_text(size=13, color= "black")) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank()) +
  theme(panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 16),
        legend.key.size = unit(0.9, "cm"),
        legend.text  = element_text(size = 16)) 

#### DECONTAM ####

###First trying Orourke method on GitHub

## Filtering phyloseq object to remove any instances where ASV in just single sample
## function from here: https://github.com/joey711/phyloseq/issues/917

library(decontam)
#install.packages("maditr")
library(maditr)
psf <- filter_taxa(kps, function (x) {sum(x > 0) > 1}, prune=TRUE)

sample_data(psf)$is.neg <- sample_data(psf)$Sample.or.Control =="Control"

basic_contam.function <- function(threshold, label){
  tmp <- isContaminant(psf, method="prevalence", neg="is.neg", threshold=threshold)
  tmp <- tmp %>% mutate(ASVid=row.names(tmp))
  tmp %>% mutate(threshold=label) %>% mutate(batch="basic")
}

cp_basic01 <- basic_contam.function(0.1, "0.1")
cp_basic02 <- basic_contam.function(0.2, "0.2")
cp_basic03 <- basic_contam.function(0.3, "0.3")
cp_basic05 <- basic_contam.function(0.5, "0.5")

basic_contam.prev <- rbind(cp_basic01, cp_basic02, cp_basic03, cp_basic05)
rm(cp_basic01, cp_basic02, cp_basic03, cp_basic05)

basic_contam_table <- basic_contam.prev %>% 
  group_by(threshold, contaminant) %>%
  tally() %>% 
  dcast(., threshold~contaminant, value.var = "n")

ggplot(data = basic_contam.prev %>% filter(threshold=="0.1"), aes(p)) + 
  geom_histogram(bins=100, color='black', fill="gray50") +
  theme_bw() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=15)) +
  #geom_vline(xintercept = 0.2, color="red", linetype="dashed") +
  labs(x="decontam Score", y="Number ASVs")
