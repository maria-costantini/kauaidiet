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

##Having some trouble following along with this may come back, gonna use my method

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("decontam")
library(decontam); packageVersion("decontam")
library(ggplot2)

df <- as.data.frame(sample_data(kps))

df$LibrarySize <- sample_sums(kps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(df, aes(x=Index, y=LibrarySize, color=Sample.or.Control)) + geom_point()
sample_data(kps)$is.neg <- sample_data(kps)$Sample.or.Control =="Control"
contamdf.prev <- isContaminant(kps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant) #table of how many samples are contaminants or not
head(which(contamdf.prev$contaminant)) 
#2056
contamdf.prev05 <- isContaminant(kps, method="prevalence", neg="is.neg", threshold=0.5) #moreaggressive, means it will identify as contaminants all sequences that are more prevalent in negative controls than in positive samples
table(contamdf.prev05$contaminant)

kps.pa <- transform_sample_counts(kps, function(abund) 1*(abund>0))
kps.pa.neg <- prune_samples(sample_data(kps.pa)$Sample.or.Control == "Control", kps.pa)
kps.pa.pos <- prune_samples(sample_data(kps.pa)$Sample.or.Control == "Sample", kps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(kps.pa.pos), pa.neg=taxa_sums(kps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")



contams <- which(contamdf.prev05$contaminant) #remove these from data?

which(sample_sums(kps) < 500)
which(sample_sums(kps) < 1000)
kps.pa.neg
kps.pa.pos
prune_taxa(taxa_sums(kps.pa.neg) > 0, kps.pa.neg)
prune_taxa(taxa_sums(kps.pa.pos) > 0, kps.pa.pos)

#filter out contaminants 
kps.nocontam <- prune_taxa(!contamdf.prev05$contaminant, kps)
kps.nocontam

# then remove blank samples + OTUs with 0 remaining reads
kps.sample <- prune_samples(sample_data(kps.nocontam)$Sample.or.Control == "Sample", kps.nocontam)
kps.sample <- prune_taxa(taxa_sums(kps.sample) > 0, kps.sample)
kps.sample

rm(kps.pa, kps.pa.neg, kps.pa.pos)
rm(contamdf.prev, contamdf.prev05)

# Just checking which SAMPLES are likely to be removed later based on low nr of reads
table(sample_sums(kps.sample) < 1000)
which(sample_sums(kps.sample) < 1000)
table(sample_sums(kps.sample) < 500)
which(sample_sums(kps.sample) < 500)

head(tax_table(kps.sample))

kps.sample <- subset_taxa(kps.sample, Phylum!= "Chordata")

table(taxa_sums(kps.sample) < 1)

# prune on prevalence - using fast_melt function from above
install.packages("remotes")
remotes::install_github("alexpiper/seqateurs")
library(seqateurs)

mdt <- fast_melt(kps.sample)
prevdt <- mdt[, list(Prevalence = sum(count > 0), 
                     TotalCounts = sum(count)), by = taxaID]
head(prevdt, 10)

#prevalence equals how many samples the otu is present in
#count equals how many times that otu shows up in a given sample
#total count equals how many times an otu shows up in all samples
hist(prevdt$Prevalence, breaks=30)
hist(prevdt$Prevalence, xlim = c(0,10))

kps.sample
kps.sample.3 <- prune_taxa(prevdt[(Prevalence >= 3 ), taxaID], kps.sample)
kps.sample.5 <- prune_taxa(prevdt[(Prevalence >= 5 ), taxaID], kps.sample)
kps.sample.3
kps.sample.5
kps.sample.1 <- prune_taxa(prevdt[(Prevalence > 1), taxaID], kps.sample)
kps.sample.1
rm(kps.sample.1)
rm(kps.sample.3)
rm(kps.sample.5)
#filtering using above code seems very agrresive. I am going to go forward by just removing singletons, aka TotalCounts >1

##### prune / filter samples based on singletons / total counts ######
# total counts
library("data.table"); packageVersion("data.table")
tdt <- data.table(tax_table(kps.sample),  
                  TotalCounts = taxa_sums(kps.sample),  
                  OTU = taxa_names(kps.sample))
head(tdt)

library(ggplot2)

tdt[(TotalCounts <= 0), .N]
tdt[(TotalCounts <= 1), .N]
# How many doubletons?
tdt[(TotalCounts <= 5), .N]
#
tdt[(TotalCounts <= 2), .N]
# taxa cumulative sum
taxsum = tdt[, .N, by = TotalCounts]
setkey(taxsum, TotalCounts)
taxsum[, CumSum := cumsum(N)]
p <- ggplot(taxsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold") +
  theme_bw() 

p

# Zoom-in
p + xlim(0, 100) + geom_vline(xintercept = c(5,10))

p + xlim(0, 100) + geom_vline(xintercept = c(5,15))
kps.single <- prune_samples(sample_sums(kps.sample) >= 100, kps.sample)
kps.single
kps.sample <- prune_taxa(tdt$OTU[tdt$TotalCounts > 1], kps.sample)
kps.sample

table(sample_sums(kps.sample) < 1000)
which(sample_sums(kps.sample) < 1000)
which(sample_sums(kps.sample )< 500)
#remove samples with too low number of reads
sample_sums(kps.sample)[order(-sample_sums(kps.sample))]
#kps.clean1 <- prune_samples(sample_sums(kps.sample) >= 100, kps.sample)

kps.clean <- prune_samples(sample_sums(kps.sample) >= 1000, kps.sample)
kps.clean <- prune_taxa(taxa_sums(kps.clean) > 0, kps.clean)
kps.clean

hist(sample_sums(kps.sample)[order(-sample_sums(kps.sample))], breaks=20)
summary(sample_data(kps.clean))

## ----see-depths AFTER contaminant removal--------------------------------------------------------

sample_data(kps.clean)$LibrarySizeFilt <- sample_sums(kps.clean)

df <- as.data.frame(sample_data(kps.clean)) 
df <- df[order(df$LibrarySize),]
# or
df <- df[order(df$LibrarySizeFilt),]
head(df)
df$Index <- seq(nrow(df))
df$Date.Collected
ggplot(data=metadata, aes(x=Index, y=LibrarySize, color=Island)) + 
  geom_point(size=2.5) +
  theme_bw() +
  theme(axis.title.x= element_text(size=16, margin = margin(t = 10)), 
        axis.title.y= element_text(size=16, margin = margin(r = 10)),
        axis.text.y = element_text(size=13, color= "black"), 
        axis.text.x = element_text(size=13, color= "black")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 16),
        legend.key.size = unit(0.9, "cm"),
        legend.text  = element_text(size = 16)) 

# after filtering
ggplot(data=df, aes(x=Index, y=LibrarySizeFilt, color=Island)) + 
  geom_point(size=2.5) +
  theme_bw() +
  theme(axis.title.x= element_text(size=16, margin = margin(t = 10)), 
        axis.title.y= element_text(size=16, margin = margin(r = 10)),
        axis.text.y = element_text(size=13, color= "black"), 
        axis.text.x = element_text(size=13, color= "black")) +
  theme(panel.grid.minor=element_blank()) +
  theme(panel.border = element_rect(colour = "grey29")) +
  theme(legend.title = element_text(size = 16),
        legend.key.size = unit(0.9, "cm"),
        legend.text  = element_text(size = 16)) 

df = data.frame(nreads = sort(taxa_sums(kps.clean), TRUE), sorted = 1:ntaxa(kps.clean), 
                type = "OTUs")
df = rbind(df, data.frame(nreads = sort(sample_sums(kps.clean), TRUE), 
                          sorted = 1:nsamples(kps.clean), type = "Samples"))

ggplot(df, aes(x = sorted, y = nreads)) + 
  geom_bar(stat = "identity") +
  ggtitle("Total nr of reads") + 
  scale_y_log10() + 
  facet_wrap(~type, 1, scales = "free") + 
  theme_bw()

###Need to get only Kauai samples

kps <- prune_samples(sample_data(kps.clean)$Island == "Kauai", kps.clean)
#kps.single <- prune_samples(sample_data(kps.single)$Island == "Kauai", kps.single)

##make sure only samples are in data set

kps <- prune_samples(sample_data(kps)$Sample.or.Control == "Sample", kps)

# Transform to RELATIVE ABUNDANCES for weighted Unifrac and bray-curtis 
# since they are sensitive to differences in total counts. 
kps
kps.clean.rel <- transform_sample_counts(kps, function(x) x / sum(x))
kps.clean.rel

###############################################################
### RAREFY ### 

set.seed(386)
kps.clean.R = rarefy_even_depth(kps, sample.size = min(sample_sums(kps)))
kps.clean.R
get_taxa_unique(kps, "Phylum")
get_taxa_unique(kps, "Class")

# Top otus
# class 20
top20otus = names(sort(taxa_sums(kps), TRUE)[1:20])
kps.clean20 <- prune_taxa(top20otus, kps) 
tax_table(kps.clean20)



rm(basic_contam_table,basic_contam.prev,cp_basic01,cp_basic02,cp_basic03,cp_basic05)
rm(kps.clean1, kps.single,p)
