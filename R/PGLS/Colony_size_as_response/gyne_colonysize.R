library(emmeans)
library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)
library(performance)

#PGLS for gyne data

#load tree
FLTree <- read.tree("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Phylogenetic_Data/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP(1).tre")
#read in species data
df<-readRDS("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Trait_Data/species_level_data_final.rds")
#here we can fiddle and see if we can substitute certain species for one another
#get species name in right format
df$Genus.Species<-gsub("\\.", "_",df$Genus.Species)

#not sure about matching leviscula
df$Genus.Species[df$Genus.Species == "Crematogaster_atkinsoni"]<-"Crematogaster_pilosa"
df$Genus.Species[df$Genus.Species == "Colobopsis_impressa"]<-"Colobopsis_saundersi"
df$Genus.Species[df$Genus.Species == "Dorymyrmex_bureni"]<-"Dorymyrmex_elegans"

df$Genus.Species[df$Genus.Species == "Lasius_americanus"]<-"Lasius_alienus"
df$Genus.Species[df$Genus.Species == "Neivamyrmex_carolinensis"]<-"Neivamyrmex_texanus"
df$Genus.Species[df$Genus.Species == "Proceratium_pergandei"]<-"Proceratium_avium"
df$Genus.Species[df$Genus.Species == "Pseudoneoponera_stigma"]<-"Pseudoponera_stigma"

df$Genus.Species[df$Genus.Species == "Solenopsis_pergandei"]<-"Solenopsis_carolinensis"
df$Genus.Species[df$Genus.Species == "Tapinoma_litorale"]<-"Tapinoma_melanocephalum"

morph<-df %>%
  dplyr::select(Genus.Species, Gyne, median_size) %>%
  drop_na()

#get tip labels
tiplabelsFLtree=FLTree$tip.label

#species in the data but not in the phylogeny
sort(setdiff(morph$Genus.Species, FLTree$tip.label))
#species in the phylogeny but not in the data
sort(setdiff(FLTree$tip.label, morph$Genus.Species))
#remove species that are not in the phylogeny from the main dataset
new_morph_df<-morph[which(morph$Genus.Species%in%
                            tiplabelsFLtree),]


#trim tree to species dataset
pruned_FLtree<-drop.tip(FLTree,setdiff(FLTree$tip.label,new_morph_df$Genus.Species ))
length(pruned_FLtree$tip.label)
length(unique(new_morph_df$Genus.Species))



new_morph_df2<-new_morph_df %>%
  group_by(Genus.Species, Gyne) %>%
  summarise(
    median_size = mean(median_size, na.rm =T)) %>%
  ungroup()

length(unique(new_morph_df2$Genus.Species))

#running pgls with caper
new_morph_df2$morphism<-as.factor(new_morph_df2$Gyne)
row.names(new_morph_df2)<-new_morph_df2$Genus.Species

orderc<-row.names(new_morph_df2)

bm3<-corBrownian(phy = pruned_FLtree, form = ~orderc)
bm3
gls1<-gls(log10(median_size) ~ morphism,
          data= new_morph_df2,
          correlation = bm3)

summary(gls1)

bm4<-corPagel(1, pruned_FLtree, form = ~orderc)
bm4
gls2<-gls(log10(median_size) ~ morphism,data= new_morph_df2,
          correlation = bm4)
bbmle::AIC(gls1, gls2)

summary(gls2)
rr2::R2(gls2)

pairs1<-emmeans(gls2, specs = c("morphism"), type = "response")
pairs1
pairs(pairs1)


get_model_vals_QWDg<-function(x, directory){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = "morphism", type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_colony_size_gyne.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_colony_size_gyne.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Colony_size_as_response/Colony~gyne_summary/"
get_model_vals_QWDg(gls2, directory = dir)



