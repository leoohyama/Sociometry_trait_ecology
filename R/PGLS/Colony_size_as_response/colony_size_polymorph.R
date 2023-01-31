library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)
library(emmeans)



#load tree
FLTree <- read.tree("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Phylogenetic_Data/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP(1).tre")
#read in species data
df<-readRDS("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Trait_Data/species_level_data_final.rds")

df$Genus.Species<-gsub("\\.", "_",df$Genus.Species)

#not sure about matching leviscula
df$Genus.Species[df$Genus.Species == "Crematogaster_atkinsoni"]<-"Crematogaster_pilosa"
df$Genus.Species[df$Genus.Species == "Colobopsis_impressa"]<-"Colobopsis_saundersi"
df$Genus.Species[df$Genus.Species == "Dorymyrmex_bureni"]<-"Dorymyrmex_elegans"
df$Genus.Species[df$Genus.Species == "Dorymyrmex_bossutus"]<-"Dorymyrmex_elegans"
df$Genus.Species[df$Genus.Species == "Dorymyrmex_reginicula"]<-"Dorymyrmex_elegans"
df$Genus.Species[df$Genus.Species == "Eurhopalothrix_floridana"]<-"Eurhopalothrix_australis"
df$Genus.Species[df$Genus.Species == "Formica_dolosa"]<-"Formica_incerta"
df$Genus.Species[df$Genus.Species == "Formica_pallidefulva"]<-"Formica_incerta"
df$Genus.Species[df$Genus.Species == "Formica_archboldi"]<-"Formica_incerta"
df$Genus.Species[df$Genus.Species == "Lasius_americanus"]<-"Lasius_alienus"
df$Genus.Species[df$Genus.Species == "Neivamyrmex_carolinensis"]<-"Neivamyrmex_opacithorax"
df$Genus.Species[df$Genus.Species == "Proceratium_crassicorne"]<-"Proceratium_silaceum"
df$Genus.Species[df$Genus.Species == "Proceratium_croceum"]<-"Proceratium_silaceum"
df$Genus.Species[df$Genus.Species == "Proceratium_pergandei"]<-"Proceratium_avium"
df$Genus.Species[df$Genus.Species == "Pseudoneoponera_stigma"]<-"Pseudoponera_stigma"

df$Genus.Species[df$Genus.Species == "Solenopsis_pergandei"]<-"Solenopsis_carolinensis"
df$Genus.Species[df$Genus.Species == "Strumigenys_clypeata"]<-"Strumigenys_laevinasis"
df$Genus.Species[df$Genus.Species == "Strumigenys_epinotalis"]<-"Strumigenys_laevinasis"
df$Genus.Species[df$Genus.Species == "Tapinoma_litorale"]<-"Tapinoma_melanocephalum"
df$Genus.Species[df$Genus.Species == "Temnothorax_bradleyi"]<-"Temnothorax_schaumii"
df$Genus.Species[df$Genus.Species == "Temnothorax_palustris"]<-"Temnothorax_schaumii"
df$Genus.Species[df$Genus.Species == "Temnothorax_texanus"]<-"Temnothorax_schaumii"


morph<-df %>%
  dplyr::select(Genus.Species, morphism, median_size) %>%
  drop_na(median_size)
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
  group_by(Genus.Species, morphism) %>%
  summarise(
            median_size = mean(median_size, na.rm =T)) %>%
  ungroup()

new_morph_df2$morphism<-as.factor(new_morph_df2$morphism)

row.names(new_morph_df2)<-new_morph_df2$Genus.Species

ordercc<-new_morph_df2$Genus.Species

new_morph_df2$morphism<-relevel(new_morph_df2$morphism, ref = "Polymorphic")

bm<-corBrownian(1, pruned_FLtree, form = ~ordercc)
bm
gls1<-gls(log(median_size) ~ morphism,
          data= new_morph_df2,
          correlation = bm)

bm1<-corPagel(1, pruned_FLtree, form = ~ordercc)
bm1
gls2<-gls(log(median_size) ~ morphism,
          data= new_morph_df2,
          correlation = bm1)

AIC(gls1,gls2)
summary(gls2)
rr2::R2(gls2)

pair1<-emmeans::emmeans(gls2, specs = "morphism", type = "response")
pair1
pairs(pair1)


get_model_vals_QWDpoly<-function(x, directory, spec){

  setwd(dir)

  pairs1<-emmeans(gls2, specs = spec, type = "response")

  write.csv(as.data.frame(pairs1),
            "Emmeans_output_colony_size_morph.csv")

  write.csv(as.data.frame(pairs(pairs1)),
            "Emmeans_pairs_output_colony_size_morph.csv")

  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/Colony_size_as_response/Colony~morphism/"
get_model_vals_QWDpoly(gls2, directory = dir, spec = "morphism")





