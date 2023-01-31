library(emmeans)
library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)
library(performance)

#PGLS for colony, native

#load tree
FLTree <- read.tree("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Phylogenetic_Data/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP(1).tre")
#read in species data
df<-readRDS("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Trait_Data/species_level_data_final.rds")
df$Genus.Species<-gsub("\\.", "_",df$Genus.Species)

#here we can fiddle and see if we can substitute certain species for one another
df$Genus.Species[df$Genus.Species == "Dolichoderus_mariae"]<-"Dolichoderus_pustulatus"
df$Genus.Species[df$Genus.Species == "Pheidole_lamia"]<-"Pheidole_pilifera"
df$Genus.Species[df$Genus.Species == "Leptogenys_manni"]<-"Leptogenys_diminuta"
df$Genus.Species[df$Genus.Species == "Crematogaster_atkinsoni"]<-"Crematogaster_pilosa"
df$Genus.Species[df$Genus.Species == "Colobopsis_impressa"]<-"Colobopsis_saundersi"
df$Genus.Species[df$Genus.Species == "Dorymyrmex_bureni"]<-"Dorymyrmex_elegans"
df$Genus.Species[df$Genus.Species == "Eurhopalothrix_floridana"]<-"Eurhopalothrix_australis"
df$Genus.Species[df$Genus.Species == "Formica_archboldi"]<-"Formica_incerta"
df$Genus.Species[df$Genus.Species == "Lasius_americanus"]<-"Lasius_alienus"
df$Genus.Species[df$Genus.Species == "Neivamyrmex_carolinensis"]<-"Neivamyrmex_texanus"
df$Genus.Species[df$Genus.Species == "Proceratium_crassicorne"]<-"Proceratium_silaceum"
df$Genus.Species[df$Genus.Species == "Proceratium_pergandei"]<-"Proceratium_avium"
df$Genus.Species[df$Genus.Species == "Pseudoneoponera_stigma"]<-"Pseudoponera_stigma"
df$Genus.Species[df$Genus.Species == "Strumigenys_clypeata"]<-"Strumigenys_laevinasis"
df$Genus.Species[df$Genus.Species == "Temnothorax_bradleyi"]<-"Temnothorax_schaumii"
df$Genus.Species[df$Genus.Species == "Monomorium_trageri"]<-"Monomorium_floricola"
df$Genus.Species[df$Genus.Species == "Prionopelta_antillana"]<-"Prionopelta_amabilis"
df$Genus.Species[df$Genus.Species == "Tetramorium_lanuginosum"]<-"Tetramorium_bicarinatum"

morph<-df %>%
  dplyr::select(Genus.Species, HW,HL,ML,EL,WL, status) %>%
  drop_na() %>%
  rowwise() %>%
  mutate(HW = HW/WL,
         HL = HL/WL,
         ML = ML/WL,
         EL = EL/WL,
         WL = WL)

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
  group_by(Genus.Species, status) %>%
  summarise(
    HW = mean(HW, na.rm =T),
    HL = mean(HL, na.rm =T),
    ML = mean(ML, na.rm =T),
    EL = mean(EL, na.rm =T),
    WL = mean(WL, na.rm =T)) %>%
  ungroup()

length(unique(new_morph_df2$Genus.Species))

row.names(new_morph_df2) <-new_morph_df2$Genus.Species
ordercc = new_morph_df2$Genus.Species

hist(log(new_morph_df2$HW))
hist(log(new_morph_df2$HL))
hist(sqrt(new_morph_df2$ML))
hist(sqrt(new_morph_df2$EL))
hist(log(new_morph_df2$WL))

bm<-corBrownian(1, pruned_FLtree, form = ~ordercc)
bm

bm1<-corPagel(1, pruned_FLtree, form = ~ordercc)
bm1

gls1<-gls(log(HW) ~ status,data= new_morph_df2,
          correlation = bm)

gls2<-gls(log(HW) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_worker<-function(x, directory, spec){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = spec, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_hw_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_hw_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Worker_trait_as_response/HW~status_model_summary/"
get_model_vals_worker(gls2, directory = dir, spec ="status")



##########head length
gls1<-gls(log(HL) ~ status,data= new_morph_df2,
          correlation = bm)

gls2<-gls(log(HL) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_worker<-function(x, directory, spec){
  setwd(directory)
  pairs1<-emmeans(x, specs = spec, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_hl_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_hl_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Worker_trait_as_response/HL~status_model_summary/"
get_model_vals_worker(gls2, directory = dir, spec ="status")


#######mandible length
gls1<-gls(sqrt(ML) ~ status,data= new_morph_df2,
          correlation = bm)
gls2<-gls(sqrt(ML) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_worker<-function(x, directory, spec){
  setwd(directory)
  pairs1<-emmeans(x, specs = spec, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_ml_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_ml_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Worker_trait_as_response/ML~status_model_summary/"
get_model_vals_worker(gls2, directory = dir, spec ="status")



######eye length
gls1<-gls((EL) ~ status,data= new_morph_df2,
          correlation = bm)

gls2<-gls((EL) ~ status,data= new_morph_df2,
          correlation = bm1)


AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_worker<-function(x, directory, spec){
  setwd(directory)
  pairs1<-emmeans(x, specs = spec, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_el_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_el_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Worker_trait_as_response/EL~status_model_summary/"
get_model_vals_worker(gls2, directory = dir, spec ="status")


#### Weber length
gls1<-gls(log(WL) ~ status,data= new_morph_df2,
          correlation = bm)


gls2<-gls(log(WL) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_worker<-function(x, directory, spec){
  setwd(directory)
  pairs1<-emmeans(x, specs = spec, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_wl_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_wl_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Worker_trait_as_response/WL~status_model_summary/"
get_model_vals_worker(gls2, directory = dir, spec ="status")

