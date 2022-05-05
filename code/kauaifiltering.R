setwd("~/Documents/PhDData/Diet_Analysis/honeycreeperdiet/rawdata")


library(tidyverse)
install.packages("qiime2R")
library(qiime2R)
library(phyloseq)

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

test 