#Gijsbert Werner, Balliol College, University of Oxford
#April-June 2019

#Script to explore new visualisations for the updated TRY manuscript 

#Load libraries
library(ape)
library(phytools)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(viridis)
library(diversitree)
library(Rphylopars)
library(castor)
library(ggplot2)
library(plotrix)

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
dat$Double.Log.Number.of.Traits<-log(dat$Log.Number.of.Traits+1) #Do +1 so that zero traits present becomes log(0+1) = 0. 
summary(dat$Log.Number.of.Traits)
summary(dat$Double.Log.Number.of.Traits)
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

# Analyses - run_thin5perc  ------------------------------------------------

analysis_start<-Sys.time()

#Let's make some small dataset for run_thin5perc code.
set.seed(01865)
run_thin5perc_dat<-dat[sample(nrow(dat),size=round(0.05*nrow(dat),0)),]
#run_thin5perc_dat<-dat[sample(nrow(dat),size=1000),] for trialling
run_thin5perc_tree<-drop.tip(tree,
                             tree$tip.label[!tree$tip.label %in% run_thin5perc_dat$match_col])
run_thin5perc_tree

##ASRs

####Quantiative reconstruction, double log numbers
summary(run_thin5perc_dat$Double.Log.Number.of.Traits)
#Create vector
run_thin5perc_trait_num_double_log<-run_thin5perc_dat$Double.Log.Number.of.Traits
names(run_thin5perc_trait_num_double_log)<-run_thin5perc_dat$match_col
#For ease of plotting, order vector same order as in tree
run_thin5perc_trait_num_double_log<-run_thin5perc_trait_num_double_log[match(run_thin5perc_tree$tip.label,names(run_thin5perc_trait_num_double_log))]

system.time(
  run_thin5perc_tree_rec_num_double_log<-anc.recon(trait_data = run_thin5perc_trait_num_double_log,tree = run_thin5perc_tree)
)
save(run_thin5perc_tree_rec_num_double_log,file = "./Models/run_thin5perc_tree_rec_num_double_log")
gc()


# Plotting ----------------------------------------------------------------

####Combine everything (four potential combinations)

#Generate baseplot
names(run_thin5perc_dat)
run_thin5perc_dat_plotting_traits <- run_thin5perc_dat %>% select(trait_num_bins,
                                                                  Leaf.Area,
                                                                  SLA,
                                                                  Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                                  Seed.Dry.Mass,
                                                                  Plant.Height,
                                                                  Stem.Specific.Density..SSD.)
head(run_thin5perc_dat_plotting_traits)
names(run_thin5perc_dat_plotting_traits)[1]<-"Presence"
names(run_thin5perc_dat_plotting_traits)[4]<-"Leaf.N"
names(run_thin5perc_dat_plotting_traits)[7]<-"SSD"
rownames(run_thin5perc_dat_plotting_traits)<-run_thin5perc_dat$match_col

#RColorbrewer 9-class Set3, seelction
#Suggestions Jens, June 19
viridis(100)[1]
table(run_thin5perc_dat$trait_num_bins)
length(which(run_thin5perc_dat$Number.of.Traits<1))

#Ok, do the median of each of these in double log terms
log(log(median(c(1,5))+1)+1)
log(log(median(c(6,10))+1)+1)
log(log(median(c(11,100))+1)+1)
log(log(median(c(101,max(run_thin5perc_dat$Number.of.Traits)))+1)+1)
#and the max is
log(log(max(run_thin5perc_dat$Number.of.Traits)+1)+1)

#So that means bin number:
round(100*(log(log(median(c(1,5))+1)+1)/log(log(max(run_thin5perc_dat$Number.of.Traits)+1)+1)))
round(100*(log(log(median(c(6,10))+1)+1)/log(log(max(run_thin5perc_dat$Number.of.Traits)+1)+1)))
round(100*(log(log(median(c(11,100))+1)+1)/log(log(max(run_thin5perc_dat$Number.of.Traits)+1)+1)))
round(100*(log(log(median(c(101,max(run_thin5perc_dat$Number.of.Traits)))+1)+1)/log(log(max(run_thin5perc_dat$Number.of.Traits)+1)+1)))

#And the corresponding colours
viridis(100)[44]
viridis(100)[59]
viridis(100)[82]
viridis(100)[96]

Sys.time()
pdf(file="./Figures/Main_text_fig_raw.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_thin5perc_tree,dat = run_thin5perc_dat_plotting_traits,cols = list(Presence=c("#440154FF","#26818EFF","#1FA287FF","#7FD34EFF","#DDE318FF"),
                                                                                         Leaf.Area=c("#440154FF","#8dd3c7"),
                                                                                         SLA=c("#440154FF","#bebada"),
                                                                                         Leaf.N=c("#440154FF","#fb8072"),
                                                                                         Seed.Dry.Mass=c("#440154FF","#80b1d3"),
                                                                                         Plant.Height=c("#440154FF","#fdb462"),
                                                                                         SSD=c("#440154FF","#b3de69")),
           legend=T,cex.lab=0.0001,edge.width=0.1,cex.legend = 0.5,plot=F,
           edge.color = viridis(100)[cut(run_thin5perc_tree_rec_num_double_log[match(run_thin5perc_tree$edge[,1],names(run_thin5perc_tree_rec_num_double_log[,1])),1],breaks=100)])
draw.circle(x=0,y=0,radius=325.0508,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-50,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-100,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-150,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-200,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-250,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-300,col="#80808025",border="#80808025")
par(new=T)
trait.plot(tree = run_thin5perc_tree,dat = run_thin5perc_dat_plotting_traits,cols = list(Presence=c("#440154FF","#26818EFF","#1FA287FF","#7FD34EFF","#DDE318FF"),
                                                                                         Leaf.Area=c("#440154FF","#8dd3c7"),
                                                                                         SLA=c("#440154FF","#bebada"),
                                                                                         Leaf.N=c("#440154FF","#fb8072"),
                                                                                         Seed.Dry.Mass=c("#440154FF","#80b1d3"),
                                                                                         Plant.Height=c("#440154FF","#fdb462"),
                                                                                         SSD=c("#440154FF","#b3de69")),
           legend=T,cex.lab=0.0001,edge.width=0.1,cex.legend = 0.5,plot=T, 
           edge.color = viridis(100)[cut(run_thin5perc_tree_rec_num_double_log[match(run_thin5perc_tree$edge[,1],names(run_thin5perc_tree_rec_num_double_log[,1])),1],breaks=100)])
add.color.bar(150,viridis(100),title = "To change manually",prompt = F,
              lims = c(min(run_thin5perc_tree_rec_num_double_log[,1]),max(run_thin5perc_tree_rec_num_double_log[,1])),fsize=0.5,
              x=-100,y=-50)
dev.off()
Sys.time()
gc()


pdf(file="./Figures/Main_text_fig_onlynodenames.pdf",width = 20,height = 20)
plot.phylo(x = run_thin5perc_tree,type = "f",show.tip.label = F,show.node.label = T,edge.width = 0.1,edge.color = "gray",cex=0.6)
dev.off()


pdf(file="./Figures/Main_text_fig_nodes_and_tips.pdf",width = 100,height = 100)
plot.phylo(x = run_thin5perc_tree,type = "f",show.tip.label = T,show.node.label = T,edge.width = 0.1,edge.color = "gray",cex = 0.1)
dev.off()


# Analyses - Supplementary  ------------------------------------------------

analysis_start<-Sys.time()

#Let's make some small dataset for run_fullk code.
set.seed(01865)
run_fullk_dat<-dat[sample(nrow(dat),size=nrow(dat)),]
#run_fullk_dat<-dat[sample(nrow(dat),size=1000),] for trialling
run_fullk_tree<-drop.tip(tree,
                         tree$tip.label[!tree$tip.label %in% run_fullk_dat$match_col])
run_fullk_tree

##ASRs

####Quantiative reconstruction, double log numbers
summary(run_fullk_dat$Double.Log.Number.of.Traits)
#Create vector
run_fullk_trait_num_double_log<-run_fullk_dat$Double.Log.Number.of.Traits
names(run_fullk_trait_num_double_log)<-run_fullk_dat$match_col
#For ease of plotting, order vector same order as in tree
run_fullk_trait_num_double_log<-run_fullk_trait_num_double_log[match(run_fullk_tree$tip.label,names(run_fullk_trait_num_double_log))]

system.time(
  run_fullk_tree_rec_num_double_log<-anc.recon(trait_data = run_fullk_trait_num_double_log,tree = run_fullk_tree)
)
save(run_fullk_tree_rec_num_double_log,file = "./Models/run_fullk_tree_rec_num_double_log")
gc()


# Plotting ----------------------------------------------------------------

#Generate baseplot
names(run_fullk_dat)
run_fullk_dat_plotting_traits <- run_fullk_dat %>% select(trait_num_bins,
                                                          Leaf.Area,
                                                          SLA,
                                                          Leaf.Nitrogen.Content.Per.Dry.Mass,
                                                          Seed.Dry.Mass,
                                                          Plant.Height,
                                                          Stem.Specific.Density..SSD.)
head(run_fullk_dat_plotting_traits)
names(run_fullk_dat_plotting_traits)[1]<-"Presence"
names(run_fullk_dat_plotting_traits)[4]<-"Leaf.N"
names(run_fullk_dat_plotting_traits)[7]<-"SSD"
rownames(run_fullk_dat_plotting_traits)<-run_fullk_dat$match_col

#RColorbrewer 9-class Set3, seelction
#Suggestions Jens, June 19
viridis(100)[1]
table(dat$trait_num_bins)
length(which(dat$Number.of.Traits<1))

#Ok, do the median of each of these in double log terms
log(log(median(c(1,5))+1)+1)
log(log(median(c(6,10))+1)+1)
log(log(median(c(11,100))+1)+1)
log(log(median(c(101,max(dat$Number.of.Traits)))+1)+1)
#and the max is
log(log(max(dat$Number.of.Traits)+1)+1)

#So that means bin number:
round(100*(log(log(median(c(1,5))+1)+1)/1.994106))
round(100*(log(log(median(c(6,10))+1)+1)/1.994106))
round(100*(log(log(median(c(11,100))+1)+1)/1.994106))
round(100*(log(log(median(c(101,max(dat$Number.of.Traits)))+1)+1)/1.994106))

#And the corresponding colours
viridis(100)[44]
viridis(100)[58]
viridis(100)[81]
viridis(100)[96]

Sys.time()
pdf(file="./Figures/FullSupFigure_FullColours_NoNames.pdf",width = 8.2,height = 8.2)
trait.plot(tree = run_fullk_tree,dat = run_fullk_dat_plotting_traits,cols = list(Presence=c("#440154FF","#26818EFF","#1FA287FF","#7FD34EFF","#E4E419FF"),
                                                                                 Leaf.Area=c("#440154FF","#8dd3c7"),
                                                                                 SLA=c("#440154FF","#bebada"),
                                                                                 Leaf.N=c("#440154FF","#fb8072"),
                                                                                 Seed.Dry.Mass=c("#440154FF","#80b1d3"),
                                                                                 Plant.Height=c("#440154FF","#fdb462"),
                                                                                 SSD=c("#440154FF","#b3de69")),
           legend=F,cex.lab=0.00001,edge.width=0.1,cex.legend = 0.5,plot=F,
           edge.color = viridis(100)[cut(run_fullk_tree_rec_num_double_log[match(run_fullk_tree$edge[,1],names(run_fullk_tree_rec_num_double_log[,1])),1],breaks=100)])
draw.circle(x=0,y=0,radius=325.0508,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-50,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-100,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-150,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-200,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-250,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-300,col="#80808025",border="#80808025")
par(new=T)
trait.plot(tree = run_fullk_tree,dat = run_fullk_dat_plotting_traits,cols = list(Presence=c("#440154FF","#26818EFF","#1FA287FF","#7FD34EFF","#E4E419FF"),
                                                                                 Leaf.Area=c("#440154FF","#8dd3c7"),
                                                                                 SLA=c("#440154FF","#bebada"),
                                                                                 Leaf.N=c("#440154FF","#fb8072"),
                                                                                 Seed.Dry.Mass=c("#440154FF","#80b1d3"),
                                                                                 Plant.Height=c("#440154FF","#fdb462"),
                                                                                 SSD=c("#440154FF","#b3de69")),
           legend=F,cex.lab=0.00001,edge.width=0.1,cex.legend = 0.5,plot=T,
           edge.color = viridis(100)[cut(run_fullk_tree_rec_num_double_log[match(run_fullk_tree$edge[,1],names(run_fullk_tree_rec_num_double_log[,1])),1],breaks=100)])
dev.off()
Sys.time()
gc()

Sys.time()
pdf(file="./Figures/FullSupFigure_Names.pdf",width = 200,height = 200)
plot.phylo(x = run_fullk_tree,show.node.label = F,show.tip.label = T,cex=0.1,plot=F,edge.width = 0.1)
draw.circle(x=0,y=0,radius=325.0508,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-50,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-100,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-150,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-200,col="#80808025",border="#80808025")
draw.circle(x=0,y=0,radius=325.0508-250,col="#FFFFFFFF",border="#FFFFFFFF")
draw.circle(x=0,y=0,radius=325.0508-300,col="#80808025",border="#80808025")
par(new=T)
plot.phylo(x = run_fullk_tree,show.node.label = F,show.tip.label = T,cex=0.1,plot=T,edge.width = 0.1)
dev.off()
Sys.time()
gc()

#merge

