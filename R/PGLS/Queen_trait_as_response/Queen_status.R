library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)
library(performance)
library(emmeans)

#PGLS for queen traits


#load tree
FLTree <- read.tree("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Phylogenetic_Data/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP(1).tre")
#read in species data
df<-readRDS("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Trait_Data/species_level_data_final.rds")
df$Genus.Species<-gsub("\\.", "_",df$Genus.Species)


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
  dplyr::select(Genus.Species, hw,hl,ml,el,wl, status) %>%
  drop_na() %>%
  rowwise() %>%
  mutate(HW = hw/wl,
         HL = hl/wl,
         ML = ml/wl,
         EL = el/wl,
         WL = wl)

#get tip labels
tiplabelsFLtree=FLTree$tip.label


#species in the data but not in the phylogeny
sort(setdiff(morph$Genus.Species, FLTree$tip.label))
#species in the phylogeny but not in the data
sort(setdiff(FLTree$tip.label, morph$Genus.Species))
#remove species that are not in the phylogeny from the main dataset
new_morph_df<-morph[which(morph$Genus.Species %in%
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

#running pgls with caper
row.names(new_morph_df2) <-new_morph_df2$Genus.Species

#check distribution of response variable
hist((new_morph_df2$HW))
hist(new_morph_df2$HL)
hist((new_morph_df2$ML))
hist(sqrt(new_morph_df2$EL))
hist(log(new_morph_df2$WL))

ordercc<-new_morph_df2$Genus.Species

########head width pgls
bm<-corBrownian(phy =pruned_FLtree, form = ~ordercc)
bm
gls1<-gls((HW) ~ status,data= new_morph_df2,
          correlation = bm)

bm1<-corPagel(1, pruned_FLtree, form = ~ordercc)
bm1
gls2<-gls((HW) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenhw_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenhw_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QHW~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")

#################Weber length

gls1<-gls(log10(WL) ~ status,data= new_morph_df2,
          correlation = bm)


gls2<-gls(log10(WL) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenhw_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenhw_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QWL~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")

########head width pgls
bm<-corBrownian(phy = pruned_FLtree, form = ~ordercc)
bm
gls1<-gls((HW) ~ status,data= new_morph_df2,
          correlation = bm)

bm1<-corPagel(1, pruned_FLtree, form = ~ordercc)
bm1
gls2<-gls((HW) ~ status,data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenhw_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenhw_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QHW~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")



##############head length pgls

gls1<-gls((HL) ~ status,data= new_morph_df2,
          correlation = bm)

gls2<-gls((HL) ~ status,data= new_morph_df2,
          correlation = bm1)


AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenhl_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenhl_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}
dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QHL~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")

################Mandible length pgls

gls1<-gls((ML) ~ status,data= new_morph_df2,
          correlation = bm)

gls2<-gls((ML) ~ status,data= new_morph_df2,
          correlation = bm1)


AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenml_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenml_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}
dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QML~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")


################Eye length pgls

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


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenel_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenel_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}
dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QEL~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")

################Mandible length pgls

gls1<-gls((ML) ~ status,data= new_morph_df2,
          correlation = bm)

gls2<-gls((ML) ~ status,data= new_morph_df2,
          correlation = bm1)


AIC(gls1,gls2)

summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "status", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDg<-function(x, directory, specs){
  setwd(directory)
  pairs1<-emmeans(gls2, specs = specs, type = "response")

  write.csv(as.data.frame(pairs1), "Emmeans_output_queenml_status.csv")


  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_queenml_status.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}
dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Queen_trait_as_response/QML~status_model_summary/"
get_model_vals_QWDg(gls2, directory = dir, specs = "status")
setwd("/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/")





