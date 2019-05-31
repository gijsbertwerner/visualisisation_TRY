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
dat<-
  dat %>% mutate(
    state1_herb = case_when(
      GIFT_PlantGrowthForm=="herb" ~ 1,
      GIFT_PlantGrowthForm=="herb&shrub" ~ 1/2,
      GIFT_PlantGrowthForm=="herb&shrub&tree" ~ 1/3,
      GIFT_PlantGrowthForm=="herb&tree" ~ 1/2,
      GIFT_PlantGrowthForm=="shrub" ~ 0,
      GIFT_PlantGrowthForm=="shrub&tree" ~ 0,
      GIFT_PlantGrowthForm=="tree" ~ 0
    ),
    state2_shrub = case_when(
      GIFT_PlantGrowthForm=="herb" ~ 0,
      GIFT_PlantGrowthForm=="herb&shrub" ~ 1/2,
      GIFT_PlantGrowthForm=="herb&shrub&tree" ~ 1/3,
      GIFT_PlantGrowthForm=="herb&tree" ~ 0,
      GIFT_PlantGrowthForm=="shrub" ~ 1,
      GIFT_PlantGrowthForm=="shrub&tree" ~ 1/2,
      GIFT_PlantGrowthForm=="tree" ~ 0
    ),
    state3_tree = case_when(
      GIFT_PlantGrowthForm=="herb" ~ 0,
      GIFT_PlantGrowthForm=="herb&shrub" ~ 0,
      GIFT_PlantGrowthForm=="herb&shrub&tree" ~ 1/3,
      GIFT_PlantGrowthForm=="herb&tree" ~ 1/2,
      GIFT_PlantGrowthForm=="shrub" ~ 0,
      GIFT_PlantGrowthForm=="shrub&tree" ~ 1/2,
      GIFT_PlantGrowthForm=="tree" ~ 1
    )
  )
head(dat %>% select(GIFT_PlantGrowthForm,state1_herb,state2_shrub,state3_tree))
tail(dat %>% select(GIFT_PlantGrowthForm,state1_herb,state2_shrub,state3_tree))
dat$state1_herb+dat$state2_shrub+dat$state3_tree
#All looks good. 

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


# Analyses - run_87.5k  ------------------------------------------------

analysis_start<-Sys.time()

##Set up overall analysis for very small tree (1000 species, i.e. 0.2%)
#Let's make some small dataset for run_87.5k code.
set.seed(01865)
run_87.5k_dat<-dat[sample(nrow(dat),size=87500),]
run_87.5k_tree<-drop.tip(tree,
                     tree$tip.label[!tree$tip.label %in% run_87.5k_dat$match_col])
run_87.5k_tree

##ASRs

#Two potential semi-quantitative ones: 1. log of number / or absolute numbers, or 2. bins of numbers (ordinal) 
#A discrete one (gf)

####Quantiative reconstruction, log numbers
summary(run_87.5k_dat$Log.Number.of.Traits)
#Create vector
run_87.5k_trait_num_log<-run_87.5k_dat$Log.Number.of.Traits
names(run_87.5k_trait_num_log)<-run_87.5k_dat$match_col
#For ease of plotting, order vector same order as in tree
run_87.5k_trait_num_log<-run_87.5k_trait_num_log[match(run_87.5k_tree$tip.label,names(run_87.5k_trait_num_log))]

system.time(
  run_87.5k_tree_rec_num_log<-anc.recon(trait_data = run_87.5k_trait_num_log,tree = run_87.5k_tree)
)
save(run_87.5k_tree_rec_num_log,file = "./Models/run_87.5k_tree_rec_num_log")

gc()

###Quantiative reconstruction, absolute numbers
summary(run_87.5k_dat$Number.of.Traits)
#Create vector
run_87.5k_trait_num<-run_87.5k_dat$Number.of.Traits
names(run_87.5k_trait_num)<-run_87.5k_dat$match_col
#For ease of plotting, order vector same order as in tree
run_87.5k_trait_num<-run_87.5k_trait_num[match(run_87.5k_tree$tip.label,names(run_87.5k_trait_num))]

system.time(
  run_87.5k_tree_rec_num<-anc.recon(trait_data = run_87.5k_trait_num,tree = run_87.5k_tree)
)
save(run_87.5k_tree_rec_num,file = "./Models/run_87.5k_tree_rec_num")

gc()

#Bins reconstructed, ordinal
table(run_87.5k_dat$trait_num_bins)
#Create vector
run_87.5k_trait_num_bins<-run_87.5k_dat$trait_num_bins
names(run_87.5k_trait_num_bins)<-run_87.5k_dat$match_col
#For ease of plotting, order vector same order as in tree
run_87.5k_trait_num_bins<-run_87.5k_trait_num_bins[match(run_87.5k_tree$tip.label,names(run_87.5k_trait_num_bins))]
table(run_87.5k_trait_num_bins)

####Model as an (ordered) character 

#ARD
mapped_run_87.5k_trait_num_bins<-map_to_state_space(raw_states = run_87.5k_trait_num_bins)
system.time(
  run_87.5k_tree_trait_num_bins_ARD<-
    asr_mk_model(tree = run_87.5k_tree,Nstates = 5,tip_states = mapped_run_87.5k_trait_num_bins$mapped_states,
                 rate_model = "ARD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)
save(run_87.5k_tree_trait_num_bins_ARD,file="./Models/run_87.5k_tree_trait_num_bins_ARD")

ASR_run_87.5k_tree_trait_num_bins_ARD_vec<-apply(run_87.5k_tree_trait_num_bins_ARD$ancestral_likelihoods,1,which.max)
names(ASR_run_87.5k_tree_trait_num_bins_ARD_vec)<-1:run_87.5k_tree$Nnode+Ntip(run_87.5k_tree)

#SRD
system.time(
  run_87.5k_tree_trait_num_bins_SRD<-
    asr_mk_model(tree = run_87.5k_tree,Nstates = 5,tip_states = mapped_run_87.5k_trait_num_bins$mapped_states,
                 rate_model = "SRD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)

save(run_87.5k_tree_trait_num_bins_SRD,file="./Models/run_87.5k_tree_trait_num_bins_SRD")


ASR_run_87.5k_tree_trait_num_bins_SRD_vec<-apply(run_87.5k_tree_trait_num_bins_SRD$ancestral_likelihoods,1,which.max)
names(ASR_run_87.5k_tree_trait_num_bins_SRD_vec)<-1:run_87.5k_tree$Nnode+Ntip(run_87.5k_tree)


#####And growth form
table(run_87.5k_dat$GIFT_PlantGrowthForm)
#Create vector
run_87.5k_gf<-run_87.5k_dat$GIFT_PlantGrowthForm
names(run_87.5k_gf)<-run_87.5k_dat$match_col
#For ease of plotting, order vector same order as in tree
run_87.5k_gf<-run_87.5k_gf[match(run_87.5k_tree$tip.label,names(run_87.5k_gf))]
table(run_87.5k_gf)

##Select the data matrix containint the tip priors
run_87.5k_gf_priors<-run_87.5k_dat %>% select(match_col,state1_herb,state2_shrub,state3_tree)
rownames(run_87.5k_gf_priors)<-run_87.5k_gf_priors$match_col
run_87.5k_gf_priors<-run_87.5k_gf_priors[match(run_87.5k_tree$tip.label,run_87.5k_gf_priors$match_col),]
run_87.5k_gf_priors<-as.matrix(run_87.5k_gf_priors[,2:4])
head(run_87.5k_gf_priors)

#ARD
system.time(
  run_87.5k_tree_gf_ARD<-
    asr_mk_model(tree = run_87.5k_tree,Nstates = 3,tip_states = NULL,tip_priors = run_87.5k_gf_priors,
                 rate_model = "ARD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)
save(run_87.5k_tree_gf_ARD,file="./Models/run_87.5k_tree_gf_ARD")

ASR_run_87.5k_tree_gf_ARD_vec<-apply(run_87.5k_tree_gf_ARD$ancestral_likelihoods,1,which.max)
names(ASR_run_87.5k_tree_gf_ARD_vec)<-1:run_87.5k_tree$Nnode+Ntip(run_87.5k_tree)

#SRD
system.time(
  run_87.5k_tree_gf_SRD<-
    asr_mk_model(tree = run_87.5k_tree,Nstates = 3,tip_states = NULL,tip_priors = run_87.5k_gf_priors,
                 rate_model = "SRD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)
save(run_87.5k_tree_gf_SRD,file="./Models/run_87.5k_tree_gf_SRD")

ASR_run_87.5k_tree_gf_SRD_vec<-apply(run_87.5k_tree_gf_SRD$ancestral_likelihoods,1,which.max)
names(ASR_run_87.5k_tree_gf_SRD_vec)<-1:run_87.5k_tree$Nnode+Ntip(run_87.5k_tree)

#SRD - herb ancestor
herb_anc_vec<-c(1,0,0)
names(herb_anc_vec)<-colnames(run_87.5k_gf_priors)
herb_anc_vec
system.time(
  run_87.5k_tree_gf_SRD_herb_anc<-
    asr_mk_model(tree = run_87.5k_tree,Nstates = 3,tip_states = NULL,tip_priors = run_87.5k_gf_priors,root_prior = herb_anc_vec,
                 rate_model = "SRD",include_ancestral_likelihoods = T,reroot = T,Ntrials = 24,Nthreads = 6)
)
save(run_87.5k_tree_gf_SRD_herb_anc,file="./Models/run_87.5k_tree_gf_SRD_herb_anc")

ASR_run_87.5k_tree_gf_SRD_herb_anc_vec<-apply(run_87.5k_tree_gf_SRD_herb_anc$ancestral_likelihoods,1,which.max)
names(ASR_run_87.5k_tree_gf_SRD_herb_anc_vec)<-1:run_87.5k_tree$Nnode+Ntip(run_87.5k_tree)
gc()

analysis_end_plot_start<-Sys.time()

# Plotting ----------------------------------------------------------------

####Combine everything (four potential combinations)

#Generate baseplot
names(run_87.5k_dat)
run_87.5k_dat_plotting_traits <- run_87.5k_dat %>% select(trait_num_bins,
                                                  Leaf.Area,
                                                  SLA,
                                                  Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                  Seed.Dry.Mass,
                                                  Plant.Height,
                                                  Stem.Specific.Density..SSD.,
                                                  All.six.traits..Diaz.et.al.2016.)
head(run_87.5k_dat_plotting_traits)
names(run_87.5k_dat_plotting_traits)[1]<-"Presence"
names(run_87.5k_dat_plotting_traits)[4]<-"Leaf.N"
names(run_87.5k_dat_plotting_traits)[7]<-"SSD"
names(run_87.5k_dat_plotting_traits)[8]<-"All_Diaz"
rownames(run_87.5k_dat_plotting_traits)<-run_87.5k_dat$match_col

#RColorbrewer 9-class Set3, seelction

##Plot baseplot
Sys.time()
pdf(file="./Figures/run_87.5k_base_plot.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5)

dev.off()
Sys.time()
gc()

#Baseplot with log ASR
Sys.time()
pdf(file="./Figures/run_87.5k_trait_num_log_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = viridis(100)[cut(run_87.5k_tree_rec_num_log[match(run_87.5k_tree$edge[,1],names(run_87.5k_tree_rec_num_log[,1])),1],breaks=100)])
add.color.bar(100,viridis(100),title = "Log of trait #",prompt = F,
              lims = c(min(run_87.5k_trait_num_log),max(run_87.5k_trait_num_log)),fsize=0.5,
              x=-100,y=-50)
dev.off()
Sys.time()
gc()

#Baseplot with absolute ASR
Sys.time()
pdf(file="./Figures/run_87.5k_trait_asbolute_num_ASR.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = viridis(100)[cut(run_87.5k_tree_rec_num[match(run_87.5k_tree$edge[,1],names(run_87.5k_tree_rec_num[,1])),1],breaks=100)])
add.color.bar(100,viridis(100),title = "Trait #",prompt = F,
              lims = c(min(run_87.5k_trait_num),max(run_87.5k_trait_num)),fsize=0.5,
              x=-100,y=-50)
dev.off()
Sys.time()
gc()

#Baseplot with categorical trait numbers
Sys.time()
pdf(file="./Figures/run_87.5k_trait_num_bins_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_run_87.5k_tree_trait_num_bins_ARD_vec[
             match(run_87.5k_tree$edge[,1],names(ASR_run_87.5k_tree_trait_num_bins_ARD_vec))]])

dev.off()
Sys.time()
gc()

Sys.time()
pdf(file="./Figures/run_87.5k_trait_num_bins_ASR_SRD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = c("#bd0026","#f03b20","#fd8d3c","#fecc5c","gray90")[ASR_run_87.5k_tree_trait_num_bins_SRD_vec[
             match(run_87.5k_tree$edge[,1],names(ASR_run_87.5k_tree_trait_num_bins_SRD_vec))]])
dev.off()
Sys.time()
gc()

#Baseplot with growth forms
Sys.time()
pdf(file="./Figures/run_87.5k_gf_ASR_ARD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_run_87.5k_tree_gf_ARD_vec[match(run_87.5k_tree$edge[,1],names(ASR_run_87.5k_tree_gf_ARD_vec))]])

dev.off()
Sys.time()
gc()

#Baseplot with growth forms
Sys.time()
pdf(file="./Figures/run_87.5k_gf_ASR_SRD.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_run_87.5k_tree_gf_SRD_vec[match(run_87.5k_tree$edge[,1],names(ASR_run_87.5k_tree_gf_SRD_vec))]])

dev.off()
Sys.time()
gc()

#Baseplot with growth forms
Sys.time()
pdf(file="./Figures/run_87.5k_gf_ASR_SRD_herc_anc.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_87.5k_tree,dat = run_87.5k_dat_plotting_traits,cols = list(Presence=c("gray90","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
                                                                         Leaf.Area=c("gray90","#8dd3c7"),
                                                                         SLA=c("gray90","#bebada"),
                                                                         Leaf.N=c("gray90","#fb8072"),
                                                                         Seed.Dry.Mass=c("gray90","#80b1d3"),
                                                                         Plant.Height=c("gray90","#fdb462"),
                                                                         SSD=c("gray90","#b3de69"),
                                                                         All_Diaz=c("gray90","#fccde5")),
           legend=T,cex.lab=0.0001,edge.width=0.2,cex.legend = 0.5,
           edge.color = c("lightgreen","darkgreen","brown")[ASR_run_87.5k_tree_gf_SRD_herc_anc_vec[match(run_87.5k_tree$edge[,1],names(ASR_run_87.5k_tree_gf_SRD_herc_anc_vec))]])

dev.off()
Sys.time()
gc()

plotting_end<-Sys.time()

analysis_start
analysis_end_plot_start
plotting_end
save.image()


