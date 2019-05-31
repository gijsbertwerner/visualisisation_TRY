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
                                         labels = c("Absence","Trait# 1-5","Trait# 6-10","Trait# 11-100","Trait# > 100"),right=F))
table(dat$trait_num_bins)

#Create the proper categories for GF
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")
dat$GIFT_PlantGrowthForm<-gsub(pattern="/","&",dat$GIFT_PlantGrowthForm) #We do this so that with the quantitative reconstruction they will be properly read as being uknown, when using rayDISC in the corHMM package
table(dat$GIFT_PlantGrowthForm,useNA = "ifany")

dat %>% filter(GIFT_PlantGrowthForm=="other") #Drop the 'others' to limit computation -> Discuss with Jens
dat <- dat %>% filter(GIFT_PlantGrowthForm!="other" | is.na(GIFT_PlantGrowthForm))
tree<-drop.tip(tree,tree$tip.label[!tree$tip.label %in% dat$match_col])
tree

dat$GIFT_PlantGrowthForm<-ifelse(is.na(dat$GIFT_PlantGrowthForm),"herb&shrub&tree",dat$GIFT_PlantGrowthForm)
table(dat$GIFT_PlantGrowthForm,useNA = "ifany") #Discuss this strategy with Jens

#See help function of ?asr_mk_model
#Generate tip state pror based on these.
#Let's code herb as state 1, shrub as state 2 and tree as state 3
dat$state1_herb<-NA
dat$state2_shrub<-NA
dat$state3_tree<-NA

#For those were all possible, model as if all states are equally likely

#Model the herb onlies
for(i in 1:nrow(dat)){
  if(dat$GIFT_PlantGrowthForm[i]=="herb"){
    dat$state1_herb[i]<-1
    dat$state2_shrub[i]<-0
    dat$state3_tree[i]<-0
  }
  if(dat$GIFT_PlantGrowthForm[i]=="herb&shrub"){
    dat$state1_herb[i]<-1/2
    dat$state2_shrub[i]<-1/2
    dat$state3_tree[i]<-0
  }
  if(dat$GIFT_PlantGrowthForm[i]=="herb&shrub&tree"){
    dat$state1_herb[i]<-1/3
    dat$state2_shrub[i]<-1/3
    dat$state3_tree[i]<-1/3
  }
  if(dat$GIFT_PlantGrowthForm[i]=="herb&tree"){
    dat$state1_herb[i]<-1/2
    dat$state2_shrub[i]<-0
    dat$state3_tree[i]<-1/1
  }
  if(dat$GIFT_PlantGrowthForm[i]=="shrub"){
    dat$state1_herb[i]<-0
    dat$state2_shrub[i]<-1
    dat$state3_tree[i]<-0
  }
  if(dat$GIFT_PlantGrowthForm[i]=="shrub&tree"){
    dat$state1_herb[i]<-0
    dat$state2_shrub[i]<-1/2
    dat$state3_tree[i]<-1/2
  }
  if(dat$GIFT_PlantGrowthForm[i]=="tree"){
    dat$state1_herb[i]<-0
    dat$state2_shrub[i]<-0
    dat$state3_tree[i]<-1
  }
}

head(dat %>% select(GIFT_PlantGrowthForm,state1_herb,state2_shrub,state3_tree))
tail(dat %>% select(GIFT_PlantGrowthForm,state1_herb,state2_shrub,state3_tree))

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

###Write analysed tree and data to folder
write.csv(x = dat,file = "./Data/analysed_data.csv",row.names = F)
write.tree(phy = tree,file = "./Data/analysed_tree.tre")


# Visualisation - Trial  ------------------------------------------------

##Set up overall analysis for very small tree (1000 species, i.e. 0.2%)
#Let's make some small dataset for trial code.
set.seed(01865)
trial_dat<-dat[sample(nrow(dat),size=250),]
trial_tree<-drop.tip(tree,
                     tree$tip.label[!tree$tip.label %in% trial_dat$match_col])
plot.phylo(trial_tree,type="f",cex = 0.15)
trial_tree

##ASRs

#Two potential semi-quantitative ones: 1. log of number / or absolute numbers, or 2. bins of numbers (ordinal) 
#A discrete one (gf)

####Quantiative reconstruction, log numbers
summary(trial_dat$Log.Number.of.Traits)
#Create vector
trial_trait_num_log<-trial_dat$Log.Number.of.Traits
names(trial_trait_num_log)<-trial_dat$match_col
#For ease of plotting, order vector same order as in tree
trial_trait_num_log<-trial_trait_num_log[match(trial_tree$tip.label,names(trial_trait_num_log))]

system.time(
  trial_tree_rec_num_log<-anc.recon(trait_data = trial_trait_num_log,tree = trial_tree)
)

#Plot the results
plot.phylo(trial_tree,type="fan",cex=0.5,
           tip.color = viridis(100)[cut(trial_trait_num_log,breaks=100)],
           edge.color = viridis(100)[cut(trial_tree_rec_num_log[match(trial_tree$edge[,1],names(trial_tree_rec_num_log[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(trial_tree_rec_num_log[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 
gc()

###Quantiative reconstruction, absolute numbers
summary(trial_dat$Number.of.Traits)
#Create vector
trial_trait_num<-trial_dat$Number.of.Traits
names(trial_trait_num)<-trial_dat$match_col
#For ease of plotting, order vector same order as in tree
trial_trait_num<-trial_trait_num[match(trial_tree$tip.label,names(trial_trait_num))]

system.time(
  trial_tree_rec_num<-anc.recon(trait_data = trial_trait_num,tree = trial_tree)
)

#Plot the results
plot.phylo(trial_tree,type="fan",cex=0.5,
           tip.color = viridis(100)[cut(trial_trait_num,breaks=100)],
           edge.color = viridis(100)[cut(trial_tree_rec_num[match(trial_tree$edge[,1],names(trial_tree_rec_num[,1])),1],breaks=100)])
nodelabels(col=viridis(100)[cut(trial_tree_rec_num[,1],breaks=100)],pch=16) #Probably leave this out of the eventual one. 
gc()

#Bins reconstructed, ordinal
table(trial_dat$trait_num_bins)
#Create vector
trial_trait_num_bins<-trial_dat$trait_num_bins
names(trial_trait_num_bins)<-trial_dat$match_col
#For ease of plotting, order vector same order as in tree
trial_trait_num_bins<-trial_trait_num_bins[match(trial_tree$tip.label,names(trial_trait_num_bins))]
table(trial_trait_num_bins)

#turn tree into bifurcating one
trial_tree_bifurc<-multi2di(trial_tree)

####Model as an (ordered) character 

#ARD
mapped_trial_trait_num_bins<-map_to_state_space(raw_states = trial_trait_num_bins)
system.time(
  trial_tree_trait_num_bins_ARD<-
    asr_mk_model(tree = trial_tree,Nstates = 5,tip_states = mapped_trial_trait_num_bins$mapped_states,
             rate_model = "ARD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)

trial_tree_trait_num_bins_ARD

#Plot with pies
plot.phylo(trial_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[trial_trait_num_bins])
nodelabels(pie = trial_tree_trait_num_bins_ARD$ancestral_likelihoods,piecol = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),cex=0.25)

#Plot with colours
ASR_trial_tree_trait_num_bins_ARD_vec<-apply(trial_tree_trait_num_bins_ARD$ancestral_likelihoods,1,which.max)
names(ASR_trial_tree_trait_num_bins_ARD_vec)<-1:trial_tree$Nnode+Ntip(trial_tree)
plot.phylo(trial_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[trial_trait_num_bins],
           edge.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[ASR_trial_tree_trait_num_bins_ARD_vec[match(trial_tree$edge[,1],names(ASR_trial_tree_trait_num_bins_ARD_vec))]])

#SRD
system.time(
  trial_tree_trait_num_bins_SRD<-
    asr_mk_model(tree = trial_tree,Nstates = 5,tip_states = mapped_trial_trait_num_bins$mapped_states,
                 rate_model = "SRD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)

trial_tree_trait_num_bins_SRD

#Plot with pies
plot.phylo(trial_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[trial_trait_num_bins])
nodelabels(pie = trial_tree_trait_num_bins_SRD$ancestral_likelihoods,piecol = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),cex=0.25)

#Plot with colours
ASR_trial_tree_trait_num_bins_SRD_vec<-apply(trial_tree_trait_num_bins_SRD$ancestral_likelihoods,1,which.max)
names(ASR_trial_tree_trait_num_bins_SRD_vec)<-1:trial_tree$Nnode+Ntip(trial_tree)
plot.phylo(trial_tree,type="f",cex=0.25,
           tip.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[trial_trait_num_bins],
           edge.color = c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026")[ASR_trial_tree_trait_num_bins_SRD_vec[match(trial_tree$edge[,1],names(ASR_trial_tree_trait_num_bins_SRD_vec))]])


#####And growth form
table(trial_dat$GIFT_PlantGrowthForm)
#Create vector
trial_gf<-trial_dat$GIFT_PlantGrowthForm
names(trial_gf)<-trial_dat$match_col
#For ease of plotting, order vector same order as in tree
trial_gf<-trial_gf[match(trial_tree$tip.label,names(trial_gf))]
table(trial_gf)

table(trial_dat$GIFT_PlantGrowthForm)
trial_dat_gf <- trial_dat %>% select(match_col,GIFT_PlantGrowthForm)
table(trial_dat_gf$GIFT_PlantGrowthForm)

system.time(
  trial_tree_gf_ARD<-rayDISC(phy = trial_tree_bifurc,data = trial_dat_gf,ntraits = 1,
                             model="ARD",node.states = "marginal",root.p="yang",
                             verbose = T)
)

#Plot with nodes
plot.phylo(trial_tree,type="f",cex=1)
nodelabels(pie = trial_tree_gf_ARD$states,piecol = c("lightgreen","darkgreen","brown"),cex=0.25)

ASR_trial_tree_gf_ARD_vec<-apply(trial_tree_gf_ARD$states,1,which.max)
names(ASR_trial_tree_gf_ARD_vec)<-1:trial_tree$Nnode+Ntip(trial_tree)
#Plot with colours
plot.phylo(trial_tree,type="f",cex=0.25,
           tip.color = c("lightgreen","darkgreen","brown")[trial_gf],
           edge.color = c("lightgreen","darkgreen","brown")[ASR_trial_tree_gf_ARD_vec[match(trial_tree$edge[,1],names(ASR_trial_tree_gf_ARD_vec))]])


####Combine everything (for potential combinations)

#Generate baseplot
names(trial_dat)
trial_dat_plotting_traits <- trial_dat %>% select(trait_num_bins,
                                                  Leaf.Area,
                                                  SLA,
                                                  Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                  Seed.Dry.Mass,
                                                  Plant.Height,
                                                  Stem.Specific.Density..SSD.,
                                                  All.six.traits..Diaz.et.al.2016.)
head(trial_dat_plotting_traits)
names(trial_dat_plotting_traits)[1]<-"Presence"
names(trial_dat_plotting_traits)[4]<-"Leaf.N"
names(trial_dat_plotting_traits)[7]<-"SSD"
names(trial_dat_plotting_traits)[8]<-"All_Diaz"
rownames(trial_dat_plotting_traits)<-trial_dat$match_col

#RColorbrewer 9-class Set3, seelction
system.time(
  trial_base_plot<-
    trait.plot(tree = trial_tree,dat = trial_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
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
pdf(file="./trial_1k_spec_base_plot.pdf",width = 8.2,height = 8.2)
trait.plot(tree = trial_tree,dat = trial_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
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
pdf(file="./trial_1k_spec_trait_num_log_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = trial_tree,dat = trial_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(trial_tree_rec_num_log[match(trial_tree$edge[,1],names(trial_tree_rec_num_log[,1])),1],breaks=100)])

dev.off()
Sys.time()


#Baseplot with absolute ASR
Sys.time()
pdf(file="./trial_1k_spec_trait_asbolute_num_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = trial_tree,dat = trial_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = viridis(100)[cut(trial_tree_rec_num[match(trial_tree$edge[,1],names(trial_tree_rec_num[,1])),1],breaks=100)])

dev.off()
Sys.time()

#Baseplot with categorical trait numbers
Sys.time()
pdf(file="./trial_1k_spec_trait_num_bins_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = trial_tree,dat = trial_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_trial_tree_trait_num_bins_ARD_vec[
             match(trial_tree$edge[,1],names(ASR_trial_tree_trait_num_bins_ARD_vec))]])

dev.off()
Sys.time()


#Baseplot with growth forms
Sys.time()
pdf(file="./trial_1k_spec_gf_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = trial_tree,dat = trial_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.25,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_trial_tree_gf_ARD_vec[match(trial_tree$edge[,1],names(ASR_trial_tree_gf_ARD_vec))]])

dev.off()
Sys.time()
