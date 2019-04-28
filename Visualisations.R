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

####Binary yes/no
head(small_dat)
binary_traits<-small_dat %>% select(Number.of.Traits,SLA,All.six.traits..Diaz.et.al.2016.)
names(binary_traits)[1]<-"Presence"
names(binary_traits)[3]<-"All.six"
binary_traits$Presence<-ifelse(is.na(binary_traits$Presence),0,1)
binary_traits$All.six<-ifelse(is.na(binary_traits$All.six),0,1)
binary_traits$SLA<-ifelse(is.na(binary_traits$SLA),0,1)
rownames(binary_traits)<-small_dat$match_col
head(binary_traits)

trait.plot(tree = small_tree,dat = binary_traits,cols = list(Presence=c("gray90","red"),
                                                             SLA=c("gray90","green"),
                                                             All.six=c("gray90","blue")),
           cex.lab=0.01)

#Explore some potential variables, first number of traits with data
summary(small_dat$Number.of.Traits)
#Prep vector
small_trait_num<-small_dat$Number.of.Traits
names(small_trait_num)<-small_dat$match_col
small_trait_num[is.na(small_trait_num)]<-0 #Because NA means 0 traits here
#Visualise this in two ways
plotTree.wBars(small_tree,x = small_trait_num,border="white",type="fan")
plotTree.wBars(small_tree,x = small_trait_num,border="white",type="phylogram")
#Take natural log of values
small_trait_num_log<-ifelse(small_trait_num>0,log(small_trait_num),0)
plotTree.wBars(small_tree,x = small_trait_num_log,border="white",type="fan")


#Colour the bars based on the matching colour
#Create a function to generate a continuous color palette 
#inspired by https://stackoverflow.com/questions/9946630/colour-points-in-a-plot-differently-depending-on-a-vector-of-values and contMap which used rainbow

#Try with viridis gradients
plotTree.wBars(small_tree,x = small_trait_num_log,border="white",type="fan",
               col = viridis(100)[as.numeric(cut(small_trait_num_log[match(small_tree$tip.label,names(small_trait_num_log))],
                                              breaks = 100))])
add.color.bar(100,viridis(100),title = "Log of trait #",prompt = F,
              lims = c(min(small_trait_num_log),max(small_trait_num_log)),
              x=-185,y=-185)
plotTree.wBars(small_tree,x = small_trait_num,border="white",type="fan",
               col = viridis(100)[as.numeric(cut(small_trait_num[match(small_tree$tip.label,names(small_trait_num))],
                                                 breaks = 100))])
add.color.bar(100,viridis(100),title = "Trait #",prompt = F,
              lims = c(min(small_trait_num),max(small_trait_num)),
              x=-185,y=-185)

###Do some ancestral state reconstructing using contMap
small_tree_rec_num<-contMap(tree = small_tree,x = small_trait_num,outline = F,
                            type="fan",col="white",fsize=c(0,1))
plot(setMap(small_tree_rec_num,viridis(100)),fsize=c(0,1),type="fan") #Only with coloured branches. 
small_tree_rec_num_viridis<-setMap(small_tree_rec_num,viridis(100))
#And with a coloured branch and bars
plotTree.wBars(tree=small_tree_rec_num_viridis$tree,x = small_trait_num,
               method = "plotSimmap",
               colors=small_tree_rec_num_viridis$cols,
               fsize=c(0,1),type="fan",border="white",
               col = viridis(100)[as.numeric(cut(small_trait_num[match(small_tree$tip.label,names(small_trait_num))],
                                                                                        breaks = 100))])

small_tree_rec_num_log<-contMap(tree = small_tree,x = small_trait_num_log,outline = F,
                            type="fan",col="white",fsize=c(0,1))
plot(setMap(small_tree_rec_num_log,viridis(100)),fsize=c(0,1),type="fan") #Only with coloured branches. 
small_tree_rec_num_log_viridis<-setMap(small_tree_rec_num_log,viridis(100))
#And with a coloured branch and bars
plotTree.wBars(tree=small_tree_rec_num_log_viridis$tree,x = small_trait_num_log,
               method = "plotSimmap",
               colors=small_tree_rec_num_viridis$cols,
               fsize=c(0,1),type="fan",border="white",
               col = viridis(100)[as.numeric(cut(small_trait_num_log[match(small_tree$tip.label,names(small_trait_num_log))],
                                                 breaks = 100))])

#

