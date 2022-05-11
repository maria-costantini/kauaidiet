####Exploratory Analyses####


###############################################################
# NETWORK plot
# prune to just the top 100 most abundant OTUs across all samples.
set.seed(386)

# NOTE: this takes lots of memory/time so only use a few OTUs, e.g. top 100.
kps100 = prune_taxa(names(sort(taxa_sums(kps.clean.rel), TRUE))[1:100], kps.clean.rel)
kps20 = prune_taxa(names(sort(taxa_sums(kps.clean.rel), TRUE))[1:20], kps.clean.rel)

ig = make_network(kps100, type = "samples", distance = "jaccard", max.dist=0.8)
plot_network(ig, kps100, type = "samples", color = "PCR.Plate", line_weight = 0.4, label = NULL)
#doesnt look like theyre clustering by sequening run
plot_network(ig, kps100, type = "samples", color = "Species", line_weight = 0.4, label = NULL)

plot_network(ig, kps100, type = "samples", color = "Foraging.Guild", line_weight = 0.4, label = NULL)

ig2 = make_network(kps100, type = "samples", distance = "bray", max.dist=0.8)
plot_network(ig2, kps100, type = "samples", color = "Species", line_weight = 0.4, label = NULL)


#####ORDINATIONS####

library(plyr)
egen_col <- c("#F8766D",
              "#7fb9d7")

# choose one of the following distances measures and dataset
ord_psbray = ordinate(kps.clean.rel, "NMDS", "bray") #Still getting that stress is nearly zero
#ord_pswuni = ordinate(kps.clean.rel, "NMDS", "wunifrac")
ord_psjaccard = ordinate(kps.clean.R, "NMDS", "jaccard")

pcoabray = ordinate(kps.clean.rel, "PCoA", "bray")
pcoabrayR = ordinate(kps.clean.R, "PCoA", "bray")
pcoajac = ordinate(kps.clean.rel, "PCoA", "jaccard")

#pcoabrayR = ordinate(kps.clean.R, "PCoA", "bray")
#pcoajacR = ordinate(kps.clean.R, "PCoA", "jaccard")
#ord_ps = ordinate(coips.clean.rel, "PCoA", "wunifrac")
#ord_ps = ordinate(coips.clean.rel, "PCoA", "unifrac")
#ord_ps = ordinate(coips.clean.R)

# plot ordination to object
g <- plot_ordination(kps.clean.rel, pcoabrayR, color = "Species")
g
g <- plot_ordination(kps.clean.rel, pcoabray, color = "Foraging.Guild")
g
g1 <- plot_ordination(kps.clean.rel, pcoajac, color = "Species")
g1

g2 <- plot_ordination(kps.clean.rel, ord_psbray, color = "Species")
g2


# ALPHA DIVERSITY 
# (and RICHNESS)
library("ggplot2")
# You must use untrimmed, non-normalized count data for meaningful results, as many of these estimates 
# are highly dependent on the number of singletons. You can always trim the data later on if needed, 

# phyloseq default plot
alpha <- plot_richness(kps.sample, x = "Species", measures="Shannon") + geom_boxplot()
plot_richness(kps, x = "Species", measures="Shannon") + geom_boxplot()
alpha


####Doing Alpha diversity analyses in Vegan###

library(vegan)

##Trying to export of physeq so I can use in vegan##
coiasv = as(otu_table(kps), "matrix")
coiotu <- otu_table(kps, taxa_are_rows = FALSE)
# transpose if necessary
#if(taxa_are_rows(physeq1)){OTU1 <- t(OTU1)}
# Coerce to data.frame

coiasv <- t(coiasv)

coiasv = as.data.frame(coiasv)
coiasv <- cbind(rownames(coiasv), coiasv)
rownames(coiasv) <- NULL
colnames(coiasv) <- c(names(coiasv)) #to not write all the column names
colnames(coiasv)[1] <- "SampleID" 



#coiOTUtable=otu_table(coiphyseq)
#coiOTUtable <- as.data.frame(coiOTUtable)
#coiOTUtable <- t(coiOTUtable) #need to transpose 
#export your metadata table for other programs (like vegan)
#metadata_table=sample_data(coiphyseq)
#coimeta <- as.data.frame(metadata_table)
#export your taxonomy table for other programs (like vegan)
coitax=tax_table(kps)
coitax <- as.data.frame(coitax)

##Alpha in vegan##
library(tibble)
alpha <-join(coiasv, metadata, by = "SampleID", type = "left", match = "first") #join asvs and metadata
wide2 <-alpha[,c(1:4819)] #trim to get rid of metadata
head(colnames(wide2))



wide2 <- column_to_rownames(wide2, 'SampleID')

wide2 <- wide2 %>% 
rownames_to_column(var = "SampleID")
#wide2 <- as.data.frame(wide2)

#names(wide2) <- as.matrix(wide2[1, ])
#wide2 <- wide2[-1, ]
#wide2[] <- lapply(wide2, function(x) type.convert(as.character(x)))
#wide2
#wide2<- as_tibble(wide2)
#wide2 <-t(wide2)
#wide2 <- as.data.frame(wide2)
#names(wide2) <- as.matrix(wide2[1, ])
#wide2 <- wide2[,-1]
#rownames(wide2) <- wide2[,1]
#wide2[] <- lapply(wide2, function(x) type.convert(as.character(x)))
#wide2

metacoi$simpson <- diversity(wide2[,2:4819], MARGIN=1, index = "simpson") #does the same thing, but should be more direct
metacoi$shannon <- diversity(wide2[,2:4819], MARGIN=1, index = "shannon") 
metacoi$invsimpson <- diversity(wide2[,2:4819], MARGIN=1, index = "invsimpson")

####making some boxplots 
#cbbPalette <- c("#E69F00", "#0072B2") #for colorblind palette
#akik <- subset(metacoi, Species=="Akikiki")
ggplot(metacoi, aes(x = Species, y = shannon, fill = Species)) + 
  geom_boxplot() +
  labs(
    y = "Shannon Diversity"
  ) +
  theme_classic() +
  scale_fill_brewer(palette="Set1")
)
#cbbPalette2 <- c("#D55E00", "#999999")
#labels <- c("Breeding", "Non-breeding")
#plot <- ggplot(metatabcoi, aes(x = Species, y = shannon, fill=Season))+ 
  #geom_boxplot() +
 # labs(
  #  y = "Shannon Diversity"
 # ) +
 # theme_classic() +
#  scale_fill_manual(values=cbbPalette2, labels=labels) #scalefillmanual to change to cbbpalette
#)   

#plot #+ scale_x_discrete(labels=labels)  





#calculate sequencing depth 
long <- t(wide2)
long <- as.data.frame(long)
colnames(long) <- long[1,] #copies the first row to column names
long <- long[-1,] #deletes the first row so there isnt a duplicated row

#coidepth <- colSums(long2) #NOT WORKKING for some reason
#metacoi$depth <- coidepth

#rarefaction curve
#S <- specnumber(coiotu)
tcoiotu <- t(coiasv)
tcoiotu <- as.data.frame(tcoiotu)
colnames(tcoiotu) <- tcoiotu[1,] #copies the first row to column names
tcoiotu<- tcoiotu[-1,]
raremax <- min(rowSums(coiasv[,2:4819]))
#Srare <- rarefy(coiotu, raremax)
coiasv <-column_to_rownames(coiasv, 'SampleID')
rarecurve(coiasv, step = 20, sample = raremax,  col = "blue")

####glm for alpha
#having issues running for some reason
#alpha <- as.data.frame(metacoi)
#mod3 = glm(shannon ~ Species + Location + Species*Location, data=alpha)
#summary(mod3)
#Anova(mod3)

#### HEATMAPS ####
###making a heatmap###

families=tax_glom(kps.clean.rel, taxrank="Family",NArm=TRUE)

class= tax_glom(kps.clean.rel, taxrank="Class", NArm=TRUE)

genus= tax_glom(kps.clean.rel, taxrank="Genus", NArm=TRUE)

phylum=tax_glom(kps.clean.rel, taxrank="Phylum", NArm=TRUE)

order=tax_glom(kps.clean.rel, taxrank="Order",NArm=TRUE)


#heatmap by species for order
plotheat<-plot_heatmap(order, "NMDS", "jaccard", taxa.label="Order",low="#66CCFF", high="#000033", na.value="white")
q1 <- plotheat + facet_grid(~Foraging.Guild, scales= "free_x", switch = "x")
q2 <- q1 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
) 
plot(q2)

#TopNOTUs = names(sort(taxa_sums(families), TRUE)[1:35])
#top50= prune_taxa(TopNOTUs, families)


plotheat2 <-plot_heatmap(order, "NMDS", "jaccard", taxa.label="Order",low="#66CCFF", high="#000033", na.value="white")
q3 <- plotheat2 + facet_grid(~sample_Species, scales= "free_x", switch = "x")
q4 <- q3 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
) 
plot(q4)

#### Abundance Plots ####

######## ABUNDANCE PLOTS ################3

kps.clean <- subset_samples(kps.clean.rel, Island=="Kauai") #just kauai
abund <- transform_sample_counts(kps.clean.rel, function(x) x / sum(x) )

glom <- tax_glom(abund, taxrank = 'Order') 
glom # should list # taxa as # order
data <- psmelt(glom) # create dataframe from phyloseq object
data$Order <- as.character(data$Order) #convert to character
#simple way to rename phyla with < 1% abundance

data$Order[data$Abundance < 0.15] <- "Other"
medians <- ddply(data, ~Order, function(x) c(median=median(x$Abundance)))

remainder <- medians[medians$median <= 0.15,]$Order
data[data$Order %in% remainder,]$Order <- "Other"
data$Order[data$Abundance < 0.15] <- "Other"
p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Order))
p + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c( "gold" ,"darkorchid3","darkgreen","green3",
                                "seagreen1" ,"magenta" , "dodgerblue1","pink","brown1",
                                "firebrick" , "gray1","cyan1","deeppink", "darkorange1","lightpink", "blue", "snow3" ,"blue3","chocolate2", "lemonchiffon", "coral2", "cadetblue", "maroon4", "thistle1")) +
   #scale_fill_brewer(palette="Set3") +
   theme(legend.position="bottom", axis.text.x=element_blank(),#panel.border = element_blank(), 
        panel.spacing.x = unit(0,"line"), 
        axis.ticks.x=element_blank()) +  
  guides(fill=guide_legend(nrow=5)) + 
  facet_grid(~sample_Species, scales="free", space="free") 

