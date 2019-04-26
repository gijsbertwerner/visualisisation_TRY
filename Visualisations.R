#Gijsbert Werner, Balliol College, University of Oxford
#April 2019

#Script to explore new visualisations for the updated TRY manuscript 

#Load libraries
library(ape)
library(phytools)
library(dplyr)
library(data.table)

#####Data reading

#Read in data directly from the zipped file
dat<-as.data.table(
  read.csv(unz("./Data/ALLMB_Consolidated_FleshedGenera_TRY_Traits.csv.zip",
             "ALLMB_Consolidated_FleshedGenera_TRY_Traits.csv")))
head(dat)
dat2<-as.data.table(
  read.csv(unz("./Data/ALLMB_JK_AccSpec_6.csv.zip",
                   "ALLMB_JK_AccSpec_6.csv")))
head(dat2)
#Slightly different dataset

tree<-read.tree("./Data/v0.1/ALLMB.tre")
tree

######Data Cleaning and Organising

#Combine the two datases
dat_comb<-merge(dat2,dat,by="TRY_AccSpeciesID",all.x = T)
#Problem with this is the duplicated species. Don't do for now
length(which(duplicated(dat$TRY_AccSpeciesID)))
#There are even 140k of them!!
#Consider how to deal with. 
#An option is to match them by ALLMB_OrigSpecies
length(which(duplicated(dat$ALLMB_OrigSpecies)))
#Ok, that's good no duplicatinos here. We'll use this one for matching. 
dat$match_col<-gsub(pattern=" ","_",dat$ALLMB_OrigSpecies)
length(which(dat$match_col %in% tree$tip.label))
length(which(tree$tip.label %in% dat$match_col))
#Ok, so there about 142 non-overlapping species names -> why? 
dat[which(!dat$match_col %in% tree$tip.label),]


####Visualisation trying

#Let's make some small dataset for trial code.
set.seed(01865)
small_dat<-dat[sample(nrow(dat),size=75),]
small_tree<-drop.tip(tree,
                     tree$tip.label[!tree$tip.label %in% small_dat$match_col])
plot.phylo(small_tree,type="f",cex = 0.75)

#Explore some potential variables, first number of traits with data
summary(small_dat$Number.of.Traits)
#Prep vector
small_trait_num<-small_dat$Number.of.Traits
names(small_trait_num)<-small_dat$match_col
small_trait_num[is.na(small_trait_num)]<-0 #Because NA means 0 traits here
plotTree.wBars(small_tree,x = small_trait_num,border="white")
small_trait_num_log<-ifelse(small_trait_num>0,log(small_trait_num),0)
plotTree.wBars(small_tree,x = small_trait_num_log,border="white")

##Next step: colour gradient on branches and in coloured tip bars..
