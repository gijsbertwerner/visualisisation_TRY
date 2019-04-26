#Gijsbert Werner, Balliol College, University of Oxford
#April 2019

#Script to explore new visualisations for the updated TRY manuscript 

#Load libraries
library(ape)
library(dplyr)
library(data.table)

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

#Combine the two datases
dat_comb<-merge(dat2,dat,by="TRY_AccSpeciesID",all.x = T)
#Problem with this is the duplicated species. Don't do for now
length(which(duplicated(dat$TRY_AccSpeciesID)))
#There are even 140k of them!!

#Consider how to deal with. 
#An option is to match them by ALLMB_OrigSpecies
