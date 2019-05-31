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
library(geiger)
library(castor)

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
dat$Log.Number.of.Traits<-log(dat$Number.of.Traits+1) #Do +1 so that zero traits present becomes log(0+1) = 0. 
summary(dat$Log.Number.of.Traits)
#Break up in categories for where we want to use it. 
dat <- dat %>% mutate(trait_num_bins=cut(Number.of.Traits,
                                         breaks = c(-Inf,1,6,11,101,Inf),
                                         labels = c("Absence","# 1-5","# 6-10","# 11-100","# > 100"),right=F))
table(dat$trait_num_bins)

#Create the proper categories for GF
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")
dat$GIFT_PlantGrowthForm<-gsub(pattern="/","&",dat$GIFT_PlantGrowthForm) #We do this so that with the quantitative reconstruction they will be properly read as being uknown, when using rayDISC in the corHMM package
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")
dat %>% filter(GIFT_PlantGrowthForm=="other") #Drop the 'others' to limit computation -> Discuss with Jens
dat <- dat %>% filter(GIFT_PlantGrowthForm!="other" | is.na(GIFT_PlantGrowthForm))
dat$GIFT_PlantGrowthForm<-ifelse(is.na(dat$GIFT_PlantGrowthForm),"herb&shrub&tree",dat$GIFT_PlantGrowthForm)
table(dat$GIFT_PlantGrowthForm,useNA = "ifany") #Discuss this strategy with Jens

#Code presence/absence properly
table(dat$SLA,useNA = "ifany")
dat$SLA<-ifelse(is.na(dat$SLA),0,1)
table(dat$SLA)

table(dat$Leaf.Area,useNA = "ifany")
dat$Leaf.Area<-ifelse(is.na(dat$Leaf.Area),0,1)
table(dat$Leaf.Area)

table(dat$Leaf.Nitrogen.Content.Per.Dry.Mass,useNA = "ifany")
dat$Leaf.Nitrogen.Content.Per.Dry.Mass<-ifelse(is.na(dat$Leaf.Nitrogen.Content.Per.Dry.Mass),0,1)
table(dat$Leaf.Nitrogen.Content.Per.Dry.Mass)

table(dat$Seed.Dry.Mass,useNA = "ifany")
dat$Seed.Dry.Mass<-ifelse(is.na(dat$Seed.Dry.Mass),0,1)
table(dat$Seed.Dry.Mass)

table(dat$Plant.Height,useNA = "ifany")
dat$Plant.Height<-ifelse(is.na(dat$Plant.Height),0,1)
table(dat$Plant.Height)

table(dat$Stem.Specific.Density..SSD.,useNA = "ifany")
dat$Stem.Specific.Density..SSD.<-ifelse(is.na(dat$Stem.Specific.Density..SSD.),0,1)
table(dat$Stem.Specific.Density..SSD.)

table(dat$All.six.traits..Diaz.et.al.2016.,useNA = "ifany")
dat$All.six.traits..Diaz.et.al.2016.<-ifelse(is.na(dat$All.six.traits..Diaz.et.al.2016.),0,1)
table(dat$All.six.traits..Diaz.et.al.2016.)

# Visualisation - Smallest ------------------------------------------------

##Set up overall analysis for very small tree (1000 species, i.e. 0.2%)
#Let's make some small dataset for trial code.
set.seed(01865)
small_dat<-dat[sample(nrow(dat),size=10),]
small_tree<-drop.tip(tree,
                     tree$tip.label[!tree$tip.label %in% small_dat$match_col])
plot.phylo(small_tree,type="f",cex = 0.15)
small_tree

##ASRs

#Two potential semi-quantitative ones: 1. log of number / or numbers, or 2. bins of numbers (ordinal) 
#A discrete one (gf)

####Quantiative reconstruction, log numbers
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
nodelabels(col=viridis(100)[cut(small_tree_rec_num_log[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 

#So this models num trait = 0, as the same as num-trait = 1. 

###Quantiative reconstruction, absolute numbers
summary(small_dat$Number.of.Traits)
#Create vector
small_trait_num<-small_dat$Number.of.Traits
names(small_trait_num)<-small_dat$match_col
#For ease of plotting, order vector same order as in tree
small_trait_num<-small_trait_num[match(small_tree$tip.label,names(small_trait_num))]

system.time(
  small_tree_rec_num<-anc.recon(trait_data = small_trait_num,tree = small_tree)
)

#Plot the results
plot.phylo(small_tree,type="fan",cex=0.2,
           tip.color = viridis(100)[cut(small_trait_num,breaks=100)],
           edge.color = viridis(100)[cut(small_tree_rec_num[match(small_tree$edge[,1],names(small_tree_rec_num[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(small_tree_rec_num[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 

#This retains the distinction between absence and number of traits = 1, but you can't actually see it really. 

#Bins reconstructed, ordinal
table(small_dat$trait_num_bins)
#Create vector
small_trait_num_bins<-small_dat$trait_num_bins
names(small_trait_num_bins)<-small_dat$match_col
#For ease of plotting, order vector same order as in tree
small_trait_num_bins<-small_trait_num_bins[match(small_tree$tip.label,names(small_trait_num_bins))]
table(small_trait_num_bins)

#turn tree into bifurcating one
small_tree_bifurc<-multi2di(small_tree)

####Model as an (ordered) character 
#ordered using meristic in fitDiscrete 
system.time(
  small_tree_trait_num_bins_meristic<-fitDiscrete(phy = small_tree,dat=small_trait_num_bins,
                                                  model = "meristic")
) #Only works with a bifurcating tree

#Do it using raydisc (Note if wanted to do could constraint rate matrix so that this becomes ordered)
small_dat_trait_num_bins <- small_dat %>% select(match_col,trait_num_bins)
table(small_dat_trait_num_bins$trait_num_bins)

system.time(
  small_tree_trait_num_bins_ER<-rayDISC(phy = small_tree_bifurc,data = small_dat_trait_num_bins,ntraits = 1,
                                        model="ER",node.states = "marginal",root.p="yang",
                                        verbose = T)
)

system.time(
  small_tree_trait_num_bins_SYM<-rayDISC(phy = small_tree_bifurc,data = small_dat_trait_num_bins,ntraits = 1,
                                         model="SYM",node.states = "marginal",root.p="yang",
                                         verbose = T)
)

system.time(
  small_tree_trait_num_bins_ARD<-rayDISC(phy = small_tree_bifurc,data = small_dat_trait_num_bins,ntraits = 1,
                                         model="ARD",node.states = "marginal",root.p="yang",
                                         verbose = T)
)

small_tree_trait_num_bins_ER$AICc
small_tree_trait_num_bins_SYM$AICc
small_tree_trait_num_bins_ARD$AICc

#Plot with pies
plot.phylo(small_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[small_trait_num_bins])
nodelabels(pie = small_tree_trait_num_bins_ARD$states,piecol = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90"),cex=0.25)

ASR_small_tree_trait_num_bins_ARD_vec<-apply(small_tree_trait_num_bins_ARD$states,1,which.max)
names(ASR_small_tree_trait_num_bins_ARD_vec)<-1:small_tree$Nnode+Ntip(small_tree)
#Plot with colours
plot.phylo(small_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[small_trait_num_bins],
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_small_tree_trait_num_bins_ARD_vec[match(small_tree$edge[,1],names(ASR_small_tree_trait_num_bins_ARD_vec))]])

#I will use only ARD

#Now growth form

table(small_dat$GIFT_PlantGrowthForm)
#Create vector
small_gf<-small_dat$GIFT_PlantGrowthForm
names(small_gf)<-small_dat$match_col
#For ease of plotting, order vector same order as in tree
small_gf<-small_gf[match(small_tree$tip.label,names(small_gf))]
table(small_gf)

table(small_dat$GIFT_PlantGrowthForm)
small_dat_gf <- small_dat %>% select(match_col,GIFT_PlantGrowthForm)
table(small_dat_gf$GIFT_PlantGrowthForm)

system.time(
  small_tree_gf_ARD<-rayDISC(phy = small_tree_bifurc,data = small_dat_gf,ntraits = 1,
                             model="ARD",node.states = "marginal",root.p="yang",
                             verbose = T)
)

#Plot with nodes
plot.phylo(small_tree,type="f",cex=1)
nodelabels(pie = small_tree_gf_ARD$states,piecol = c("lightgreen","darkgreen","brown"),cex=0.25)

ASR_small_tree_gf_ARD_vec<-apply(small_tree_gf_ARD$states,1,which.max)
names(ASR_small_tree_gf_ARD_vec)<-1:small_tree$Nnode+Ntip(small_tree)
#Plot with colours
plot.phylo(small_tree,type="f",cex=0.25,
           tip.color = c("lightgreen","darkgreen","brown")[small_gf],
           edge.color = c("lightgreen","darkgreen","brown")[ASR_small_tree_gf_ARD_vec[match(small_tree$edge[,1],names(ASR_small_tree_gf_ARD_vec))]])


####Combine everything (for potential combinations)

#Generate baseplot
names(small_dat)
small_dat_plotting_traits <- small_dat %>% select(trait_num_bins,
                                                  Leaf.Area,
                                                  SLA,
                                                  Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                  Seed.Dry.Mass,
                                                  Plant.Height,
                                                  Stem.Specific.Density..SSD.,
                                                  All.six.traits..Diaz.et.al.2016.)
head(small_dat_plotting_traits)
names(small_dat_plotting_traits)[1]<-"Presence"
names(small_dat_plotting_traits)[4]<-"Leaf.N"
names(small_dat_plotting_traits)[7]<-"SSD"
names(small_dat_plotting_traits)[8]<-"All_Diaz"
rownames(small_dat_plotting_traits)<-small_dat$match_col

#RColorbrewer 9-class Set3, seelction
system.time(
  small_base_plot<-
    trait.plot(tree = small_tree,dat = small_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                             Leaf.Area=c("gray90","#8dd3c7"),
                                                                             SLA=c("gray90","#bebada"),
                                                                             Leaf.N=c("gray90","#fb8072"),
                                                                             Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                             Plant.Height=c("gray90","#fdb462"),
                                                                             SSD=c("gray90","#b3de69"),
                                                                             All_Diaz=c("gray90","#fccde5")),
               legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5)
)

##Plot baseplot
Sys.time()
pdf(file="./small_1k_spec_base_plot.pdf",width = 8.2,height = 8.2)
trait.plot(tree = small_tree,dat = small_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5)

dev.off()
Sys.time()


#Baseplot with log ASR
Sys.time()
pdf(file="./small_1k_spec_trait_num_log_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = small_tree,dat = small_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(small_tree_rec_num_log[match(small_tree$edge[,1],names(small_tree_rec_num_log[,1])),1],breaks=100)])

dev.off()
Sys.time()


#Baseplot with absolute ASR
Sys.time()
pdf(file="./small_1k_spec_trait_asbolute_num_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = small_tree,dat = small_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(small_tree_rec_num[match(small_tree$edge[,1],names(small_tree_rec_num[,1])),1],breaks=100)])

dev.off()
Sys.time()

#Baseplot with categorical trait numbers
Sys.time()
pdf(file="./small_1k_spec_trait_num_bins_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = small_tree,dat = small_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_small_tree_trait_num_bins_ARD_vec[
             match(small_tree$edge[,1],names(ASR_small_tree_trait_num_bins_ARD_vec))]])

dev.off()
Sys.time()


#Baseplot with growth forms
Sys.time()
pdf(file="./small_1k_spec_gf_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = small_tree,dat = small_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_small_tree_gf_ARD_vec[match(small_tree$edge[,1],names(ASR_small_tree_gf_ARD_vec))]])

dev.off()
Sys.time()







# Visualisation - Intermediate ------------------------------------------------

##Set up overall analysis for very intermediate tree (35k species, i.e. 10%)
#Let's make some intermediate dataset for trial code.
set.seed(01865)
intermediate_dat<-dat[sample(nrow(dat),size=35000),]
intermediate_tree<-drop.tip(tree,
                            tree$tip.label[!tree$tip.label %in% intermediate_dat$match_col])
plot.phylo(intermediate_tree,type="f",cex = 0.15)
intermediate_tree

##ASRs

#Two potential semi-quantitative ones: 1. log of number / or numbers, or 2. bins of numbers (ordinal) 
#A discrete one (gf)

####Quantiative reconstruction, log numbers
summary(intermediate_dat$Log.Number.of.Traits)
#Create vector
intermediate_trait_num_log<-intermediate_dat$Log.Number.of.Traits
names(intermediate_trait_num_log)<-intermediate_dat$match_col
#For ease of plotting, order vector same order as in tree
intermediate_trait_num_log<-intermediate_trait_num_log[match(intermediate_tree$tip.label,names(intermediate_trait_num_log))]

system.time(
  intermediate_tree_rec_num_log<-anc.recon(trait_data = intermediate_trait_num_log,tree = intermediate_tree)
)

#Plot the results
plot.phylo(intermediate_tree,type="fan",cex=0.2,
           tip.color = viridis(100)[cut(intermediate_trait_num_log,breaks=100)],
           edge.color = viridis(100)[cut(intermediate_tree_rec_num_log[match(intermediate_tree$edge[,1],names(intermediate_tree_rec_num_log[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(intermediate_tree_rec_num_log[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 

#So this models num trait = 0, as the same as num-trait = 1. 

###Quantiative reconstruction, absolute numbers
summary(intermediate_dat$Number.of.Traits)
#Create vector
intermediate_trait_num<-intermediate_dat$Number.of.Traits
names(intermediate_trait_num)<-intermediate_dat$match_col
#For ease of plotting, order vector same order as in tree
intermediate_trait_num<-intermediate_trait_num[match(intermediate_tree$tip.label,names(intermediate_trait_num))]

system.time(
  intermediate_tree_rec_num<-anc.recon(trait_data = intermediate_trait_num,tree = intermediate_tree)
)

#Plot the results
plot.phylo(intermediate_tree,type="fan",cex=0.2,
           tip.color = viridis(100)[cut(intermediate_trait_num,breaks=100)],
           edge.color = viridis(100)[cut(intermediate_tree_rec_num[match(intermediate_tree$edge[,1],names(intermediate_tree_rec_num[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(intermediate_tree_rec_num[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 

#This retains the distinction between absence and number of traits = 1, but you can't actually see it really. 

#Bins reconstructed, ordinal
table(intermediate_dat$trait_num_bins)
#Create vector
intermediate_trait_num_bins<-intermediate_dat$trait_num_bins
names(intermediate_trait_num_bins)<-intermediate_dat$match_col
#For ease of plotting, order vector same order as in tree
intermediate_trait_num_bins<-intermediate_trait_num_bins[match(intermediate_tree$tip.label,names(intermediate_trait_num_bins))]
table(intermediate_trait_num_bins)

#turn tree into bifurcating one
intermediate_tree_bifurc<-multi2di(intermediate_tree)

####Model as an (ordered) character 
#ordered using meristic in fitDiscrete 
system.time(
  intermediate_tree_trait_num_bins_meristic<-fitDiscrete(phy = intermediate_tree,dat=intermediate_trait_num_bins,
                                                         model = "meristic")
) #Only works with a bifurcating tree

#Do it using raydisc (Note if wanted to do could constraint rate matrix so that this becomes ordered)
intermediate_dat_trait_num_bins <- intermediate_dat %>% select(match_col,trait_num_bins)
table(intermediate_dat_trait_num_bins$trait_num_bins)

system.time(
  intermediate_tree_trait_num_bins_ER<-rayDISC(phy = intermediate_tree_bifurc,data = intermediate_dat_trait_num_bins,ntraits = 1,
                                               model="ER",node.states = "marginal",root.p="yang",
                                               verbose = T)
)

system.time(
  intermediate_tree_trait_num_bins_SYM<-rayDISC(phy = intermediate_tree_bifurc,data = intermediate_dat_trait_num_bins,ntraits = 1,
                                                model="SYM",node.states = "marginal",root.p="yang",
                                                verbose = T)
)

system.time(
  intermediate_tree_trait_num_bins_ARD<-rayDISC(phy = intermediate_tree_bifurc,data = intermediate_dat_trait_num_bins,ntraits = 1,
                                                model="ARD",node.states = "marginal",root.p="yang",
                                                verbose = T)
)

intermediate_tree_trait_num_bins_ER$AICc
intermediate_tree_trait_num_bins_SYM$AICc
intermediate_tree_trait_num_bins_ARD$AICc

#Plot with pies
plot.phylo(intermediate_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[intermediate_trait_num_bins])
nodelabels(pie = intermediate_tree_trait_num_bins_ARD$states,piecol = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90"),cex=0.25)

ASR_intermediate_tree_trait_num_bins_ARD_vec<-apply(intermediate_tree_trait_num_bins_ARD$states,1,which.max)
names(ASR_intermediate_tree_trait_num_bins_ARD_vec)<-1:intermediate_tree$Nnode+Ntip(intermediate_tree)
#Plot with colours
plot.phylo(intermediate_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[intermediate_trait_num_bins],
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_intermediate_tree_trait_num_bins_ARD_vec[match(intermediate_tree$edge[,1],names(ASR_intermediate_tree_trait_num_bins_ARD_vec))]])

#I will use only ARD

#Now growth form

table(intermediate_dat$GIFT_PlantGrowthForm)
#Create vector
intermediate_gf<-intermediate_dat$GIFT_PlantGrowthForm
names(intermediate_gf)<-intermediate_dat$match_col
#For ease of plotting, order vector same order as in tree
intermediate_gf<-intermediate_gf[match(intermediate_tree$tip.label,names(intermediate_gf))]
table(intermediate_gf)

table(intermediate_dat$GIFT_PlantGrowthForm)
intermediate_dat_gf <- intermediate_dat %>% select(match_col,GIFT_PlantGrowthForm)
table(intermediate_dat_gf$GIFT_PlantGrowthForm)

system.time(
  intermediate_tree_gf_ARD<-rayDISC(phy = intermediate_tree_bifurc,data = intermediate_dat_gf,ntraits = 1,
                                    model="ARD",node.states = "marginal",root.p="yang",
                                    verbose = T)
)

#Plot with nodes
plot.phylo(intermediate_tree,type="f",cex=1)
nodelabels(pie = intermediate_tree_gf_ARD$states,piecol = c("lightgreen","darkgreen","brown"),cex=0.25)

ASR_intermediate_tree_gf_ARD_vec<-apply(intermediate_tree_gf_ARD$states,1,which.max)
names(ASR_intermediate_tree_gf_ARD_vec)<-1:intermediate_tree$Nnode+Ntip(intermediate_tree)
#Plot with colours
plot.phylo(intermediate_tree,type="f",cex=0.25,
           tip.color = c("lightgreen","darkgreen","brown")[intermediate_gf],
           edge.color = c("lightgreen","darkgreen","brown")[ASR_intermediate_tree_gf_ARD_vec[match(intermediate_tree$edge[,1],names(ASR_intermediate_tree_gf_ARD_vec))]])


####Combine everything (for potential combinations)

#Generate baseplot
names(intermediate_dat)
intermediate_dat_plotting_traits <- intermediate_dat %>% select(trait_num_bins,
                                                                Leaf.Area,
                                                                SLA,
                                                                Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                                Seed.Dry.Mass,
                                                                Plant.Height,
                                                                Stem.Specific.Density..SSD.,
                                                                All.six.traits..Diaz.et.al.2016.)
head(intermediate_dat_plotting_traits)
names(intermediate_dat_plotting_traits)[1]<-"Presence"
names(intermediate_dat_plotting_traits)[4]<-"Leaf.N"
names(intermediate_dat_plotting_traits)[7]<-"SSD"
names(intermediate_dat_plotting_traits)[8]<-"All_Diaz"
rownames(intermediate_dat_plotting_traits)<-intermediate_dat$match_col

#RColorbrewer 9-class Set3, seelction
system.time(
  intermediate_base_plot<-
    trait.plot(tree = intermediate_tree,dat = intermediate_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                                           Leaf.Area=c("gray90","#8dd3c7"),
                                                                                           SLA=c("gray90","#bebada"),
                                                                                           Leaf.N=c("gray90","#fb8072"),
                                                                                           Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                                           Plant.Height=c("gray90","#fdb462"),
                                                                                           SSD=c("gray90","#b3de69"),
                                                                                           All_Diaz=c("gray90","#fccde5")),
               legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5)
)

##Plot baseplot
Sys.time()
pdf(file="./intermediate_35k_spec_base_plot.pdf",width = 8.2,height = 8.2)
trait.plot(tree = intermediate_tree,dat = intermediate_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                                       SLA=c("gray90","#bebada"),
                                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                                       SSD=c("gray90","#b3de69"),
                                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5)

dev.off()
Sys.time()


#Baseplot with log ASR
Sys.time()
pdf(file="./intermediate_35k_spec_trait_num_log_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = intermediate_tree,dat = intermediate_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                                       SLA=c("gray90","#bebada"),
                                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                                       SSD=c("gray90","#b3de69"),
                                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(intermediate_tree_rec_num_log[match(intermediate_tree$edge[,1],names(intermediate_tree_rec_num_log[,1])),1],breaks=100)])

dev.off()
Sys.time()


#Baseplot with absolute ASR
Sys.time()
pdf(file="./intermediate_35k_spec_trait_asbolute_num_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = intermediate_tree,dat = intermediate_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                                       SLA=c("gray90","#bebada"),
                                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                                       SSD=c("gray90","#b3de69"),
                                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(intermediate_tree_rec_num[match(intermediate_tree$edge[,1],names(intermediate_tree_rec_num[,1])),1],breaks=100)])

dev.off()
Sys.time()

#Baseplot with categorical trait numbers
Sys.time()
pdf(file="./intermediate_35k_spec_trait_num_bins_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = intermediate_tree,dat = intermediate_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                                       SLA=c("gray90","#bebada"),
                                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                                       SSD=c("gray90","#b3de69"),
                                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_intermediate_tree_trait_num_bins_ARD_vec[
             match(intermediate_tree$edge[,1],names(ASR_intermediate_tree_trait_num_bins_ARD_vec))]])

dev.off()
Sys.time()


#Baseplot with growth forms
Sys.time()
pdf(file="./intermediate_35k_spec_gf_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = intermediate_tree,dat = intermediate_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                                       SLA=c("gray90","#bebada"),
                                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                                       SSD=c("gray90","#b3de69"),
                                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_intermediate_tree_gf_ARD_vec[match(intermediate_tree$edge[,1],names(ASR_intermediate_tree_gf_ARD_vec))]])

dev.off()
Sys.time()





# Visualisation - Full ------------------------------------------------

##Set up overall analysis for very full tree (1000 species, i.e. 0.2%)
#Let's make some full dataset for trial code.
set.seed(01865)
full_dat<-dat
full_tree<-tree
plot.phylo(full_tree,type="f",cex = 0.15)
full_tree

##ASRs

#Two potential semi-quantitative ones: 1. log of number / or numbers, or 2. bins of numbers (ordinal) 
#A discrete one (gf)

####Quantiative reconstruction, log numbers
summary(full_dat$Log.Number.of.Traits)
#Create vector
full_trait_num_log<-full_dat$Log.Number.of.Traits
names(full_trait_num_log)<-full_dat$match_col
#For ease of plotting, order vector same order as in tree
full_trait_num_log<-full_trait_num_log[match(full_tree$tip.label,names(full_trait_num_log))]

system.time(
  full_tree_rec_num_log<-anc.recon(trait_data = full_trait_num_log,tree = full_tree)
)

#Plot the results
plot.phylo(full_tree,type="fan",cex=0.2,
           tip.color = viridis(100)[cut(full_trait_num_log,breaks=100)],
           edge.color = viridis(100)[cut(full_tree_rec_num_log[match(full_tree$edge[,1],names(full_tree_rec_num_log[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(full_tree_rec_num_log[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 

#So this models num trait = 0, as the same as num-trait = 1. 

###Quantiative reconstruction, absolute numbers
summary(full_dat$Number.of.Traits)
#Create vector
full_trait_num<-full_dat$Number.of.Traits
names(full_trait_num)<-full_dat$match_col
#For ease of plotting, order vector same order as in tree
full_trait_num<-full_trait_num[match(full_tree$tip.label,names(full_trait_num))]

system.time(
  full_tree_rec_num<-anc.recon(trait_data = full_trait_num,tree = full_tree)
)

#Plot the results
plot.phylo(full_tree,type="fan",cex=0.2,
           tip.color = viridis(100)[cut(full_trait_num,breaks=100)],
           edge.color = viridis(100)[cut(full_tree_rec_num[match(full_tree$edge[,1],names(full_tree_rec_num[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(full_tree_rec_num[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 

#This retains the distinction between absence and number of traits = 1, but you can't actually see it really. 

#Bins reconstructed, ordinal
table(full_dat$trait_num_bins)
#Create vector
full_trait_num_bins<-full_dat$trait_num_bins
names(full_trait_num_bins)<-full_dat$match_col
#For ease of plotting, order vector same order as in tree
full_trait_num_bins<-full_trait_num_bins[match(full_tree$tip.label,names(full_trait_num_bins))]
table(full_trait_num_bins)

#turn tree into bifurcating one
full_tree_bifurc<-multi2di(full_tree)

####Model as an (ordered) character 
#ordered using meristic in fitDiscrete 
system.time(
  full_tree_trait_num_bins_meristic<-fitDiscrete(phy = full_tree,dat=full_trait_num_bins,
                                                 model = "meristic")
) #Only works with a bifurcating tree

#Do it using raydisc (Note if wanted to do could constraint rate matrix so that this becomes ordered)
full_dat_trait_num_bins <- full_dat %>% select(match_col,trait_num_bins)
table(full_dat_trait_num_bins$trait_num_bins)

system.time(
  full_tree_trait_num_bins_ER<-rayDISC(phy = full_tree_bifurc,data = full_dat_trait_num_bins,ntraits = 1,
                                       model="ER",node.states = "marginal",root.p="yang",
                                       verbose = T)
)

system.time(
  full_tree_trait_num_bins_SYM<-rayDISC(phy = full_tree_bifurc,data = full_dat_trait_num_bins,ntraits = 1,
                                        model="SYM",node.states = "marginal",root.p="yang",
                                        verbose = T)
)

system.time(
  full_tree_trait_num_bins_ARD<-rayDISC(phy = full_tree_bifurc,data = full_dat_trait_num_bins,ntraits = 1,
                                        model="ARD",node.states = "marginal",root.p="yang",
                                        verbose = T)
)

full_tree_trait_num_bins_ER$AICc
full_tree_trait_num_bins_SYM$AICc
full_tree_trait_num_bins_ARD$AICc

#Plot with pies
plot.phylo(full_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[full_trait_num_bins])
nodelabels(pie = full_tree_trait_num_bins_ARD$states,piecol = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90"),cex=0.25)

ASR_full_tree_trait_num_bins_ARD_vec<-apply(full_tree_trait_num_bins_ARD$states,1,which.max)
names(ASR_full_tree_trait_num_bins_ARD_vec)<-1:full_tree$Nnode+Ntip(full_tree)
#Plot with colours
plot.phylo(full_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[full_trait_num_bins],
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_full_tree_trait_num_bins_ARD_vec[match(full_tree$edge[,1],names(ASR_full_tree_trait_num_bins_ARD_vec))]])

#I will use only ARD

#Now growth form

table(full_dat$GIFT_PlantGrowthForm)
#Create vector
full_gf<-full_dat$GIFT_PlantGrowthForm
names(full_gf)<-full_dat$match_col
#For ease of plotting, order vector same order as in tree
full_gf<-full_gf[match(full_tree$tip.label,names(full_gf))]
table(full_gf)

table(full_dat$GIFT_PlantGrowthForm)
full_dat_gf <- full_dat %>% select(match_col,GIFT_PlantGrowthForm)
table(full_dat_gf$GIFT_PlantGrowthForm)

system.time(
  full_tree_gf_ARD<-rayDISC(phy = full_tree_bifurc,data = full_dat_gf,ntraits = 1,
                            model="ARD",node.states = "marginal",root.p="yang",
                            verbose = T)
)

#Plot with nodes
plot.phylo(full_tree,type="f",cex=1)
nodelabels(pie = full_tree_gf_ARD$states,piecol = c("lightgreen","darkgreen","brown"),cex=0.25)

ASR_full_tree_gf_ARD_vec<-apply(full_tree_gf_ARD$states,1,which.max)
names(ASR_full_tree_gf_ARD_vec)<-1:full_tree$Nnode+Ntip(full_tree)
#Plot with colours
plot.phylo(full_tree,type="f",cex=0.25,
           tip.color = c("lightgreen","darkgreen","brown")[full_gf],
           edge.color = c("lightgreen","darkgreen","brown")[ASR_full_tree_gf_ARD_vec[match(full_tree$edge[,1],names(ASR_full_tree_gf_ARD_vec))]])


####Combine everything (for potential combinations)

#Generate baseplot
names(full_dat)
full_dat_plotting_traits <- full_dat %>% select(trait_num_bins,
                                                Leaf.Area,
                                                SLA,
                                                Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                Seed.Dry.Mass,
                                                Plant.Height,
                                                Stem.Specific.Density..SSD.,
                                                All.six.traits..Diaz.et.al.2016.)
head(full_dat_plotting_traits)
names(full_dat_plotting_traits)[1]<-"Presence"
names(full_dat_plotting_traits)[4]<-"Leaf.N"
names(full_dat_plotting_traits)[7]<-"SSD"
names(full_dat_plotting_traits)[8]<-"All_Diaz"
rownames(full_dat_plotting_traits)<-full_dat$match_col

#RColorbrewer 9-class Set3, seelction
system.time(
  full_base_plot<-
    trait.plot(tree = full_tree,dat = full_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                           Leaf.Area=c("gray90","#8dd3c7"),
                                                                           SLA=c("gray90","#bebada"),
                                                                           Leaf.N=c("gray90","#fb8072"),
                                                                           Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                           Plant.Height=c("gray90","#fdb462"),
                                                                           SSD=c("gray90","#b3de69"),
                                                                           All_Diaz=c("gray90","#fccde5")),
               legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5)
)

##Plot baseplot
Sys.time()
pdf(file="./full_full_spec_base_plot.pdf",width = 8.2,height = 8.2)
trait.plot(tree = full_tree,dat = full_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                       SLA=c("gray90","#bebada"),
                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                       SSD=c("gray90","#b3de69"),
                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5)

dev.off()
Sys.time()


#Baseplot with log ASR
Sys.time()
pdf(file="./full_full_spec_trait_num_log_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = full_tree,dat = full_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                       SLA=c("gray90","#bebada"),
                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                       SSD=c("gray90","#b3de69"),
                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(full_tree_rec_num_log[match(full_tree$edge[,1],names(full_tree_rec_num_log[,1])),1],breaks=100)])

dev.off()
Sys.time()


#Baseplot with absolute ASR
Sys.time()
pdf(file="./full_35k_spec_trait_asbolute_num_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = full_tree,dat = full_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                       SLA=c("gray90","#bebada"),
                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                       SSD=c("gray90","#b3de69"),
                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(full_tree_rec_num[match(full_tree$edge[,1],names(full_tree_rec_num[,1])),1],breaks=100)])

dev.off()
Sys.time()

#Baseplot with categorical trait numbers
Sys.time()
pdf(file="./full_full_spec_trait_num_bins_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = full_tree,dat = full_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                       SLA=c("gray90","#bebada"),
                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                       SSD=c("gray90","#b3de69"),
                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_full_tree_trait_num_bins_ARD_vec[
             match(full_tree$edge[,1],names(ASR_full_tree_trait_num_bins_ARD_vec))]])

dev.off()
Sys.time()


#Baseplot with growth forms
Sys.time()
pdf(file="./full_full_spec_gf_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = full_tree,dat = full_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                       Leaf.Area=c("gray90","#8dd3c7"),
                                                                       SLA=c("gray90","#bebada"),
                                                                       Leaf.N=c("gray90","#fb8072"),
                                                                       Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                       Plant.Height=c("gray90","#fdb462"),
                                                                       SSD=c("gray90","#b3de69"),
                                                                       All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_full_tree_gf_ARD_vec[match(full_tree$edge[,1],names(ASR_full_tree_gf_ARD_vec))]])

dev.off()
Sys.time()




