#Gijsbert Werner, Balliol College, University of Oxford
#April 2019

#Script to explore new visualisations for the updated TRY manuscript 

#Load libraries
library(ape)
library(phytools)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(viridis)
library(diversitree)
library(corHMM)
library(Rphylopars)

#####Data reading

# Data reading ------------------------------------------------------------

#Read the tree
tree<-read.tree("./Data/v0.1/ALLMB.tre")
tree

#Read in data directly from the zipped file
dat<-as.data.table(
  read.csv(unz("./Data/ALLMB_Consolidated_FleshedGenera_TRY_Traits.csv.zip",
             "ALLMB_Consolidated_FleshedGenera_TRY_Traits.csv")))
head(dat)
nrow(dat)

dat_gf<-as.data.table(
  read.csv(unz("./Data/ALLMB_Consolidated_FleshedGenera_TRY_Traits_GIFT_GrowthForm.csv.zip",
               "ALLMB_Consolidated_FleshedGenera_TRY_Traits_GIFT_GrowthForm.csv")))
head(dat_gf)
nrow(dat_gf) #Data not present for all tree tips, about 150k fewer. All missing ones treat as NA?



# Data formatting ---------------------------------------------------------

#We'll use ALLMB_OrigSpecies for matching. 
dat$match_col<-gsub(pattern=" ","_",dat$ALLMB_OrigSpecies)
#Are the all present in the tree, and vice versa? 
length(which(dat$match_col %in% tree$tip.label))
length(which(tree$tip.label %in% dat$match_col))
dat[which(!dat$match_col %in% tree$tip.label),]

#Ok, so there about 142 non-overlapping species names. These are spelling/typo problems. We will drop these. 
dat<-dat %>% filter(dat$match_col %in% tree$tip.label)
tree<-drop.tip(tree,tree$tip.label[!tree$tip.label %in% dat$match_col])
tree

#Now consider the growth form dataset
table(dat_gf$GIFT_PlantGrowthForm,useNA = "ifany")
#What column to match them by? There appears to not be the ALLMB_OrigSpecies present
head(dat_gf)
names(dat_gf)
#We'll use TRY_ACC_SpeciesID to match
dat$GIFT_PlantGrowthForm<-
  dat_gf$GIFT_PlantGrowthForm[match(dat$TRY_AccSpeciesID,dat_gf$TRY_AccSpeciesID)]
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")
#Ok, this looks all good
rm(dat_gf)
gc()

#Some more data organisation
summary(dat$Number.of.Traits)
#NA's here represent absence, i.e. zero
dat$Number.of.Traits<-ifelse(is.na(dat$Number.of.Traits),0,dat$Number.of.Traits)
summary(dat$Number.of.Traits) #This looks convincing
#Log of numbers
dat$Log.Number.of.Traits<-ifelse(dat$Number.of.Traits>0,log(dat$Number.of.Traits),0) #I am not sure if it's appropriate to only take the log of the positive numbers
summary(dat$Log.Number.of.Traits)
#Break up in categories for where we want to use it. 
dat <- dat %>% mutate(trait_num_bins=cut(Number.of.Traits,
                                         breaks = c(-Inf,1,6,11,101,Inf),
                                         labels = c("Absence","# 1-5","# 6-10","# 11-100","# > 100"),right=F))
table(dat$trait_num_bins)

#Create the proper categories for GF
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")
dat$GIFT_PlantGrowthForm<-gsub(pattern="/","&",dat$GIFT_PlantGrowthForm) #We do this so that with the quantitative reconstruction they will be properly read as being uknown
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")
dat %>% filter(GIFT_PlantGrowthForm=="other") #Drop the 'others' to limit computation -> Discuss with Jens
dat <- dat %>% filter(GIFT_PlantGrowthForm!="other" | is.na(GIFT_PlantGrowthForm))
dat$GIFT_PlantGrowthForm<-ifelse(is.na(dat$GIFT_PlantGrowthForm),"herb&shrub&tree",dat$GIFT_PlantGrowthForm)
table(dat$GIFT_PlantGrowthForm,useNA = "ifany") #Discuss this strategy with Jens

# Visualisation - Smallest ------------------------------------------------

##Set up overall analysis for very small tree (1000 species, i.e. 0.2%)
#Let's make some small dataset for trial code.
set.seed(01865)
small_dat<-dat[sample(nrow(dat),size=1000),]
small_tree<-drop.tip(tree,
                     tree$tip.label[!tree$tip.label %in% small_dat$match_col])
plot.phylo(small_tree,type="f",cex = 0.15)
small_tree

##ASRs

#Two potential semi-quantitative ones: 1. log of number / or numbers, or 2. bins of numbers (ordinal) 
#A discrete one (gf)

#Quantiative reconstruction, log numbers
summary(small_dat$Log.Number.of.Traits)
#Create vector
small_trait_num_log<-small_dat$Log.Number.of.Traits
names(small_trait_num_log)<-small_dat$match_col
#For ease of plotting, order vector same order as in tree
small_trait_num_log<-small_trait_num_log[match(small_tree$tip.label,names(small_trait_num_log))]

system.time(
  small_tree_rec_num_log<-anc.recon(trait_data = small_trait_num_log,tree = small_tree)
)

#Plot the results
plot.phylo(small_tree,type="fan",cex=0.2,
           tip.color = viridis(100)[cut(small_trait_num_log,breaks=100)],
           edge.color = viridis(100)[cut(small_tree_rec_num_log[match(small_tree$edge[,1],names(small_tree_rec_num_log[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(small_tree_rec_num_log[,1],breaks=100)],pch=16)

