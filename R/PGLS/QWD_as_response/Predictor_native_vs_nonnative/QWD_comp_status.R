
library(ape)
library(picante)
library(phytools)
library(caper)
library(tidyverse)


#PGLS for queen traits

#load tree
FLTree <- read.tree("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Phylogenetic_Data/RAxML_bestTree.result_treepl_CV_195tips_STRUMIGENYSADDED_INGROUP(1).tre")
#read in species data
df<-readRDS("~/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/Data/Trait_Data/species_level_data_final.rds")
df$Genus.Species<-gsub("\\.", "_",df$Genus.Species)

morph<-df %>%
  dplyr::select(Genus.Species, QWD, QWD_EL, QWD_HW,QWD_ML,
                QWD_HL, status) %>%
  drop_na()
#get tip labels
tiplabelsFLtree=FLTree$tip.label
sort(tiplabelsFLtree)
#get species name in right format
morph$Genus.Species<-gsub("\\.", "_",morph$Genus.Species)


#species in the data but not in the phylogeny
sort(setdiff(morph$Genus.Species, FLTree$tip.label))
#species in the phylogeny but not in the data
sort(setdiff(FLTree$tip.label, morph$Genus.Species))

#here we can fiddle and see if we can substitute certain species for one another
morph$Genus.Species[morph$Genus.Species == "Crematogaster_cerasi"]<-"Crematogaster_laeviuscula"

morph$Genus.Species[morph$Genus.Species == "Crematogaster_atkinsoni"]<-"Crematogaster_pilosa"
morph$Genus.Species[morph$Genus.Species == "Colobopsis_impressa"]<-"Colobopsis_saundersi"
morph$Genus.Species[morph$Genus.Species == "Dorymyrmex_bureni"]<-"Dorymyrmex_elegans"
morph$Genus.Species[morph$Genus.Species == "Eurhopalothrix_floridana"]<-"Eurhopalothrix_australis"
morph$Genus.Species[morph$Genus.Species == "Formica_archboldi"]<-"Formica_incerta"
morph$Genus.Species[morph$Genus.Species == "Lasius_americanus"]<-"Lasius_alienus"
morph$Genus.Species[morph$Genus.Species == "Neivamyrmex_carolinensis"]<-"Neivamyrmex_texanus"
morph$Genus.Species[morph$Genus.Species == "Proceratium_crassicorne"]<-"Proceratium_silaceum"
morph$Genus.Species[morph$Genus.Species == "Proceratium_pergandei"]<-"Proceratium_avium"
morph$Genus.Species[morph$Genus.Species == "Pseudoneoponera_stigma"]<-"Pseudoponera_stigma"
morph$Genus.Species[morph$Genus.Species == "Strumigenys_clypeata"]<-"Strumigenys_laevinasis"
morph$Genus.Species[morph$Genus.Species == "Temnothorax_bradleyi"]<-"Temnothorax_schaumii"
morph$Genus.Species[morph$Genus.Species == "Monomorium_trageri"]<-"Monomorium_floricola"
morph$Genus.Species[morph$Genus.Species == "Prionopelta_antillana"]<-"Prionopelta_amabilis"
morph$Genus.Species[morph$Genus.Species == "Tetramorium_lanuginosum"]<-"Tetramorium_bicarinatum"


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
    QWD = mean(QWD, na.rm =T),
    QWD_EL = mean(QWD_EL, na.rm =T),
    QWD_HW = mean(QWD_HW, na.rm =T),
    QWD_HL = mean(QWD_HL, na.rm =T),
    QWD_ML = mean(QWD_ML, na.rm =T)) %>%
  ungroup()

length(unique(new_morph_df2$Genus.Species))

#running pgls with caper
row.names(new_morph_df2) <-new_morph_df2$Genus.Species


hist(log(new_morph_df2$QWD))
hist((new_morph_df2$QWD_EL))
hist((new_morph_df2$QWD_ML))
hist((new_morph_df2$QWD_HW))
hist((new_morph_df2$QWD_HL))


bm<-corPagel(1, pruned_FLtree)
bm
gls1<-gls(log(QWD) ~ status,data= new_morph_df2,
          correlation = bm)

summary(gls1)
rr2::R2(gls1)

pair1<-emmeans::emmeans(gls1, specs = "status", type = "response")
pair1
pairs(pair1)


gls1<-gls((QWD_EL) ~ status,data= new_morph_df2,
          correlation = bm)

summary(gls1)
rr2::R2(gls1)

pair1<-emmeans::emmeans(gls1, specs = "status", type = "response")
pair1
pairs(pair1)


gls1<-gls((QWD_ML) ~ status,data= new_morph_df2,
          correlation = bm)

summary(gls1)
rr2::R2(gls1)

pair1<-emmeans::emmeans(gls1, specs = "status", type = "response")
pair1
pairs(pair1)


gls1<-gls((QWD_HW) ~ status,data= new_morph_df2,
          correlation = bm)

summary(gls1)
rr2::R2(gls1)

pair1<-emmeans::emmeans(gls1, specs = "status", type = "response")
pair1
pairs(pair1)

gls1<-gls((QWD_HL) ~ status,data= new_morph_df2,
          correlation = bm)

summary(gls1)
rr2::R2(gls1)

pair1<-emmeans::emmeans(gls1, specs = "status", type = "response")
pair1
pairs(pair1)




setwd("/Users/leoohyama/Google Drive/Colony_Size_FL_Project/FL_colony_size_proj/")

