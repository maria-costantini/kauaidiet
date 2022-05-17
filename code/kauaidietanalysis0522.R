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
ord_psjaccard = ordinate(kauai.rel, "NMDS", "jaccard")

pcoabray = ordinate(kps.clean.rel, "PCoA", "bray")
pcoabrayR = ordinate(kps.clean.R, "PCoA", "bray")
pcoajac = ordinate(kauai.rel, "PCoA", "jaccard")

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
g1 <- plot_ordination(kauai, pcoajac, color = "Foraging.Guild")
g1

g2 <- plot_ordination(kps.clean.rel, ord_psbray, color = "Species")
g2

g3 <- plot_ordination(kauai.rel, ord_psjaccard, color = "Species")
g3

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

remotes::install_github("vmikk/metagMisc")

families=tax_glom(kauai.rel, taxrank="Family",NArm=TRUE)

class= tax_glom(kauai.rel, taxrank="Class", NArm=TRUE)

genus= tax_glom(kauai.rel, taxrank="Genus", NArm=TRUE)

phylum=tax_glom(kauai.rel, taxrank="Phylum", NArm=TRUE)

order=tax_glom(kauai.rel, taxrank="Order",NArm=TRUE) 
#order2 <- phyloseq_rm_na_tax(order) #this removes unused taxonomic ranks
#head(tax_table(order2))

#heatmap by species for order
plotheat<-plot_heatmap(class, "NMDS", "jaccard", taxa.label="Class",low="#66CCFF", high="#000033", na.value="white")
q1 <- plotheat + facet_grid(~sample_Species, scales= "free_x", switch = "x")
q2 <- q1 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
) 
plot(q2)

#TopNOTUs = names(sort(taxa_sums(families), TRUE)[1:35])
#top50= prune_taxa(TopNOTUs, families)


plotheat2 <-plot_heatmap(order, "PCoA", "jaccard", taxa.label="Order",low="#66CCFF", high="#000033", na.value="white")
q3 <- plotheat2 + facet_grid(~sample_Species, scales= "free_x", switch = "x")
q4 <- q3 + theme(
  axis.text.x = element_blank(),
  axis.ticks = element_blank()
) 
plot(q4)

#### Abundance Plots ####

######## ABUNDANCE PLOTS ################3

#kps.clean <- subset_samples(kps.clean.rel, Island=="Kauai") #just kauai
abund <- transform_sample_counts(kauai.rel, function(x) x / sum(x) )

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
  facet_grid(~Foraging.Guild, scales="free", space="free") 

##Phyloseq barplot
plot_bar(kauai, fill = "Class", x="sample_Species")

#Top 25 taxa
#Barplot of top 25 OTUs
Top25OTUs = names(sort(taxa_sums(kauai), TRUE)[1:25])
comparetop25 = prune_taxa(Top25OTUs, kauai)
relative  = transform_sample_counts(comparetop25, function(OTU) OTU / sum(OTU))
plot_bar(physeq = relative, fill = "Order", facet_grid =~Year, x = "sample_Species")

birdspecies<-kauai@sam_data$Species
dupes<-duplicated(birdspecies)
kauai@sam_data$duplicated<-dupes
multiple.spp.obs <- prune_samples(kauai@sam_data$duplicated==TRUE, kauai)
multiple.spp.obs <- prune_samples(sample_sums(multiple.spp.obs)>=1, multiple.spp.obs)

#Evaluate within-species variation

####not really needed to do
install.packages("remotes")
install.packages("adespatial")
remotes::install_github("umerijaz/microbiomeSeq")
library(microbiomeSeq)
library(ggplot2)
kauai2 <- normalise_data(order, norm.method = "proportion")
p <- plot_taxa(kauai2, grouping_column = "Species", method="hellinger")
relative<- transform_sample_counts(kauai2, function(OTU) OTU / sum(OTU))
p<- plot_taxa(relative, grouping_column = "Order", method = 'hellinger', number.taxa = 10)
pdf(file = "LCBD_multispp.pdf", width = 10, height = 5)
p
dev.off()

#Agglomerate by species
#Remove confidence from Tax table (messes up agglomeration)
tax_table(kauai) <- tax_table(kauai)[,1:5]
mergedspp<-merge_samples(x = kauai, group = "Species")
relative.merged  = transform_sample_counts(mergedspp, function(OTU) OTU / sum(OTU))

#Alpha diversity
Obs.MOTU<-plot_richness(physeq = kauai, color = "Foraging.Guild", measures = "Observed", x = "Species") + geom_boxplot(outlier.color = "black") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))
Obs.MOTU

#Check correlation between # of species replicate & # MOTUs, controlling for outliers
nobs<-as.data.frame(multispp)
Obs.MOTU<-estimate_richness(mergedspp, measures = c("Observed", "Shannon"))
Obs.MOTU$Var1<-rownames(Obs.MOTU)
samplesize_testdr<-merge.data.frame(nobs, Obs.MOTU, by = "Var1")

outliers <- boxplot(samplesize_testdr$Observed, plot=TRUE)$out
outliers

samplesize_testdr<-samplesize_testdr[-c(15, 17, 18), ]

#Observed MOTUs
cor.test(samplesize_testdr$Freq, samplesize_testdr$Observed, method = "kendall", exact = FALSE)
plot(samplesize_testdr$Freq, samplesize_testdr$Observed)
#theres definitely a correlation between sample size and richness

cor.test(samplesize_testdr$Freq, samplesize_testdr$Shannon, method = "kendall", exact = FALSE)
#less so with shannon but still pretty high

############# DATA TRANSFORMATION AND ORDINATION #################
#Center Log transform the data (don't rarefy)
require(microbiome)
physeq.trans<-microbiome::transform(kauai, "log10")


funx.ord <- ordinate(
  physeq = physeq.trans, 
  method = "NMDS", 
  distance = "bray"
)

funx.ord.jacP <- ordinate(
  physeq = physeq.trans, 
  method = "PCoA", 
  distance = "jaccard"
)

# Calculate distance matrix
study.bray <- phyloseq::distance(physeq.trans, method = "bray")
study.jaccard<- phyloseq::distance(physeq.trans, method = "jaccard")

noPU <- subset_samples(physeq.trans, Species!="Puaiohi ")
study.jaccardnoPU <- phyloseq::distance(noPU, method = "jaccard")


# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq.trans))


library(vegan)
dis <- vegdist(study.bray)
beta <- betadisper(study.jaccard, sampledf$Species, sqrt.dist = T)
permutest(beta)
(beta.HSD <- TukeyHSD(beta))
plot(beta.HSD)
plot(beta)
plot(beta, ellipse = TRUE, hull = FALSE)
par(cex.axis=0.4)
boxplot(beta)

sampledfnoPU<-subset(sampledf, sampledf$Species!="Puaiohi ")
beta <- betadisper(study.jaccardnoPU, sampledfnoPU$Species, sqrt.dist = T)
permutest(beta)
(beta.HSD <- TukeyHSD(beta))
plot(beta.HSD)
plot(beta)
plot(beta, ellipse = TRUE, hull = FALSE)
par(cex.axis=0.4)
boxplot(beta)

require(viridis)
pdf(file = "NMDS_Diet.pdf", width = 7, height = 5)
plot_ordination(
  physeq = physeq.trans,
  ordination = funx.ord.jacP,
  axes = c(1,2), 
  color = "Foraging.Guild", 
  title = "NMDS of Dietary MOTUs") +
 # scale_color_viridis_d()+
  geom_point(aes(color = Foraging.Guild), size = 2.5) +
  #scale_shape_manual(values = 0:7) +
  theme_bw() +
  stat_ellipse(fill="Foraging.Guild") + geom_jitter()
dev.off()


################## BIPARTITE NETWORK AND FOOD WEB #####################
######## use package bipartite
#install.packages('bipartite')
require(bipartite)

######## CONVERT PHYLOSEQ MOTU TABLE AND TAX TABLE INTO BIPARTITE INTERACTION MATRIX 
merged_birds_Order<-tax_glom(physeq = mergedspp, taxrank = "Order")
merged.motus<-as.data.frame(merged_birds_Order@otu_table)
merged.tax<-as.matrix((tax_table(merged_birds_Order)))
merged.tax<-data.frame(merged.tax[,1:4])
colnames(merged.motus)<-merged.tax[,4]

interactions<-as.matrix(merged.motus)
t.interactions <- t(interactions)
t.interactions <- as.data.frame(t.interactions)

#Visualize web
#pdf(file = 'FinalBipartite Network.pdf', width = 16, height = 13)
plotweb(sortweb(t.interactions, sort.order="inc"), text.rot = 90, labsize = 0.7)
#dev.off()
visweb(sortweb(interactions, sort.order="inc"), type ="nested")

#Null model comparison
#To test network metrics (e.g. H2 specialization index) against the Vazquez null
Iobs <- networklevel(web = t.interactions, index = "H2") #0.217
nulls <- nullmodel(web=t.interactions, N=1000, method='vaznull') # takes a while!
Inulls <- sapply(nulls, networklevel, index="H2")
plot(density(Inulls), xlim=c(0, 1), lwd=2, main="H'2")
abline(v=Iobs, col="red", lwd=2)

shapiro.test(Inulls) # p = =8,28x10-5, not signif diff from normal distribution
res <- t.test(Inulls, mu = 0.217) #One sample t-test comparing nulls to observed mean
#t = -1796.3, df = 99, p-value < 2.2e-16

t.inter.noPU <- t.interactions[,c(1:7,9)]

cmod <- computeModules(t.inter.noPU)
plotModuleWeb(cmod, labsize = 0.3)

##links per species

net.metrics.links <- networklevel(web = t.interactions, index = "links per species")
# Make null models 
nullslinks <- nullmodel(web=t.interactions, N=500, method='r2dtable') # takes a while!
Inullslinks <- sapply(nullslinks, networklevel, index="links per species")
plot(density(Inullslinks), xlim=c(0, 1), lwd=2, main="Links per species")
abline(v=net.metrics.links, col="red", lwd=2,)

shapiro.test(Inullslinks) # p = =8,28x10-5, not signif diff from normal distribution
res <- t.test(Inullslinks, mu = 2.85) 

obs.spec <- specieslevel(t.interactions, index = "ALLBUTD")

#####piankas overlap######3
install.packages("EcoSimR")
library(EcoSimR)
#merged_birds <-tax_glom(physeq = mergedspp, taxrank = "Order")
spp.asvs <-as.data.frame(mergedspp@otu_table)
spp.tax<-as.matrix((tax_table(mergedspp)))
#merged.tax<-data.frame(merged.tax[,1:4])
#colnames(spp.asvs)<-merged.tax
aks <- spp.asvs[1:2,]
ekwhite <- spp.asvs[c(1,9), ]
ikwhite <- spp.asvs[c(2,9),]
aniwhite <- spp.asvs[c(3,9),]
apwhite <- spp.asvs[c(4,9),]
iwwhite <- spp.asvs[c(5,9),]
kawhite <- spp.asvs[c(6,9),]
kaewhite <- spp.asvs[c(7,9),]
apiw <- spp.asvs[c(4,5),]
ekiw <- spp.asvs[c(1,5),]
ikani <- spp.asvs[c(2,3),]

pianka(m = spp.asvs) #0.1155
pianka(m=aks) #0.2059
pianka(ekwhite) #0.0068
pianka(m=ikwhite) #0.0646
pianka(aniwhite) #0.4126
pianka(apwhite) #0.0113
pianka(iwwhite)#0.2555
pianka(kawhite) #0.0016
pianka(kaewhite) #0.0088
pianka(apiw) #0.1099
pianka(ekiw) #0.097
pianka(ikani) #0.4297



pian <- niche.overlap(spp.asvs, method = "levins")

