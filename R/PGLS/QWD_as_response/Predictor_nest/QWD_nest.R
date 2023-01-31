library(emmeans)
library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)

#PGLS for QWD vs. nst

#load tree
FLTree <- read.tree("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Phylogenetic_Data/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP(1).tre")
#read in species data
df<-readRDS("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Trait_Data/species_level_data_final.rds")
df$Genus.Species<-gsub("\\.", "_",df$Genus.Species)

#here we can fiddle and see if we can substitute certain species for one another
df$Genus.Species[df$Genus.Species == "Strumigenys_carolinensis"]<-"Strumigenys_creightoni"
df$Genus.Species[df$Genus.Species == "Crematogaster_cerasi"]<-"Crematogaster_laeviuscula"
df$Genus.Species[df$Genus.Species == "Crematogaster_atkinsoni"]<-"Crematogaster_pilosa"
df$Genus.Species[df$Genus.Species == "Colobopsis_impressa"]<-"Colobopsis_saundersi"
df$Genus.Species[df$Genus.Species == "Dorymyrmex_bureni"]<-"Dorymyrmex_elegans"
df$Genus.Species[df$Genus.Species == "Eurhopalothrix_floridana"]<-"Eurhopalothrix_australis"
df$Genus.Species[df$Genus.Species == "Formica_archboldi"]<-"Formica_incerta"
df$Genus.Species[df$Genus.Species == "Lasius_americanus"]<-"Lasius_alienus"
df$Genus.Species[df$Genus.Species == "Neivamyrmex_carolinensis"]<-"Neivamyrmex_texanus"
df$Genus.Species[df$Genus.Species == "Proceratium_pergandei"]<-"Proceratium_avium"
df$Genus.Species[df$Genus.Species == "Pseudoneoponera_stigma"]<-"Pseudoponera_stigma"
df$Genus.Species[df$Genus.Species == "Strumigenys_clypeata"]<-"Strumigenys_laevinasis"
df$Genus.Species[df$Genus.Species == "Temnothorax_bradleyi"]<-"Temnothorax_schaumii"
df$Genus.Species[df$Genus.Species == "Monomorium_trageri"]<-"Monomorium_floricola"
df$Genus.Species[df$Genus.Species == "Prionopelta_antillana"]<-"Prionopelta_amabilis"

morph<-df %>%
  dplyr::select(Genus.Species, QWD, nesting_strata) %>%
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


new_morph_df$nesting_strata[new_morph_df$Genus.Species == "Aphaenogaster_fulva"][2]<-"Ground"
new_morph_df2<-new_morph_df %>%
  group_by(Genus.Species, nesting_strata) %>%
  summarise(
    QWD = mean(QWD, na.rm =T)) %>%
  ungroup()

length(unique(new_morph_df2$Genus.Species))

#running pgls with caper
row.names(new_morph_df2) <-new_morph_df2$Genus.Species
ordercc<-new_morph_df2$Genus.Species
new_morph_df2$nesting_strata<-as.factor(new_morph_df2$nesting_strata)


bm<-corBrownian(1, pruned_FLtree, form = ~ordercc)
bm
bm1<-corPagel(1, pruned_FLtree, form = ~ordercc)
bm1
gls1<-gls(log(QWD) ~ nesting_strata,data= new_morph_df2,
          correlation = bm)

gls2<-gls(log(QWD) ~ nesting_strata,data= new_morph_df2,
          correlation = bm1)



AIC(gls1, gls2)
summary(gls2)
rr2::R2(gls2)


pairs1<-emmeans(gls2, specs = "nesting_strata", type = "response")
pairs1
pairs(pairs1)
contrast(regrid(pairs1))

get_model_vals_QWDr<-function(x, directory, spec){
  setwd(directory)
  pairs1<-emmeans(x, specs = spec, type = "response")
  write.csv(as.data.frame(pairs1), "Emmeans_output_QWD.csv")
  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output_QWD.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}
dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/QWD_as_response/Predictor_nest/QWD~nest_model_summary/"
get_model_vals_QWDr(gls2, directory = dir, spec = "nesting_strata")

setwd("/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/")

