library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)
library(emmeans)
#PGLS for QWD_EL vs. diet

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
df$Genus.Species[df$Genus.Species == "Tapinoma_litorale"]<-"Tapinoma_melanocephalum"
df$Genus.Species[df$Genus.Species == "Temnothorax_bradleyi"]<-"Temnothorax_schaumii"
df$Genus.Species[df$Genus.Species == "Monomorium_trageri"]<-"Monomorium_floricola"
df$Genus.Species[df$Genus.Species == "Prionopelta_antillana"]<-"Prionopelta_amabilis"
df$Genus.Species[df$Genus.Species == "Tetramorium_lanuginosum"]<-"Tetramorium_bicarinatum"

morph<-df %>%
  dplyr::select(Genus.Species, QWD_EL, diet) %>%
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
  group_by(Genus.Species, diet) %>%
  summarise(
    QWD_EL = mean(QWD_EL, na.rm =T)) %>%
  ungroup()

length(unique(new_morph_df2$Genus.Species))

#running pgls with caper
row.names(new_morph_df2)<-new_morph_df2$Genus.Species
new_morph_df2$diet<-as.factor(new_morph_df2$diet)
ordercc<-new_morph_df2$Genus.Species


hist((new_morph_df2$QWD_EL)^(1/3))
new_morph_df2$QWDEF_cube<-new_morph_df2$QWD_EL^(1/3)

which(is.na(new_morph_df2$QWDEF_cube))

new_morph_df2$QWD_EL[106]

new_morph_df2$QWDEF_cube[8]<-1.561243
new_morph_df2$QWDEF_cube[86]<-1.672151
new_morph_df2$QWDEF_cube[106]<-2.422315

bm<-corBrownian(1, phy = pruned_FLtree, form = ~ordercc)
bm

gls1<-gls((QWDEF_cube) ~ diet ,data= new_morph_df2,
          correlation = bm)

bm1<-corPagel(1, phy= pruned_FLtree, form = ~ordercc)

gls2<-gls((QWDEF_cube) ~ diet ,data= new_morph_df2,
          correlation = bm1)


AIC(gls1, gls2)
summary(gls2)
rr2::R2(gls2)


pairs1<-emmeans(gls2, specs = "diet", type = "response")
pairs1
pairs(pairs1)

get_model_vals_QWDr<-function(x, directory, spec){
  setwd(directory)

  pairs1<-emmeans(x, specs = spec, type = "response")
  write.csv(as.data.frame(pairs1), "Emmeans_output.csv")
  write.csv(as.data.frame(pairs(pairs1)), "Emmeans_pairs_output.csv.csv")
  r2<-rr2::R2(gls2)
  write.csv(as.data.frame(r2), "R2.csv")

}

dir<-"/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/R/PGLS/QWD_as_response/Predictor_diet/QWD_EL_diet_model_summary/"
get_model_vals_QWDr(gls2, directory = dir, spec = "diet")
setwd("/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/")

