library(tidyverse)

#make sure we are in project directory
getwd()

#Load colony data
colony<-read.csv("Data/Trait_Data/FL_colony_size_data.csv", encoding = "UTF-8",
                 na.strings = "")

#look at column types
colony<-as_tibble(colony)
colnames(colony)

#also load in gyne data
gyne<-read.csv("Data/Trait_Data/gyne_datav2.csv")

head(colony)

#join colony and gyne data by species
colony<-left_join(colony, gyne, by = c("Genus.Species" = "Genus.Species"))

#get a dataset where each row is a species with its median colony size
# this also includes species without colony data
colony_df<-colony %>%
  mutate(colony_size = as.numeric(colony_size)) %>%
  group_by(Genus.Species,diet, nesting_strata, morphism , status, Gyne) %>%
  summarise(median_size = (median(colony_size, na.rm = T)),
            ) %>%
  ungroup()

#add categories for size classes
colony_df<-colony_df %>% mutate(size_class = case_when(
  median_size < 100 ~ "0-100",
  median_size >= 100 & median_size < 1000 ~ "100-1,000",
  median_size >= 1000 & median_size <10000 ~ "1000-10,000",
  median_size >= 10000 ~ ">10,000"
))

#choose columns needed
colony_df<-colony_df %>%
  dplyr::select(Genus.Species,diet,nesting_strata,morphism, status,
         median_size,size_class, Gyne)


#read morphological trait data again
dft<-read.csv("Data/Trait_Data/worker_traits.csv")
dft$WL<-(dft$Weber_length)
dft$HW<-((dft$Head_width))
dft$HL<-((dft$Head_length))
dft$ML<-((dft$Mandible_length))
dft$EL<-((dft$Eye.length_1 + dft$Eye.length_2)/ 2)

#calculate number of specimens per sp. unit
intra<-dft %>%
  group_by(Genus.Species) %>%
  summarise(n = n())

#trait dataframe
dft1<-dft %>%
  dplyr::select(Genus.Species,Caste, HW,HL,ML,WL,EL ) %>%
  group_by(Genus.Species) %>%
  summarise(
    HW = mean(HW, na.rm =T),
    HL = mean(HL, na.rm =T),
    ML = mean(ML, na.rm =T),
    EL = mean(EL, na.rm =T),
    WL = mean(WL, na.rm =T)
  )%>%
  mutate_if(is.character, factor)


#read in queen traits
queen<-read.csv("Data/Trait_Data/queen_traits.csv", encoding = "UTF-8")


queen<-as_tibble(queen)

#get queen trait by sp units
queen_df<-queen %>% dplyr::select(Genus.Species,hw,
                                  hl, ml, el1, el2, wl) %>%
  rowwise() %>%
  mutate(el = mean(c(el1, el2), na.rm = T)) %>%
  dplyr::select(-c(el1, el2))
queen_df




#combine everything together
full_Set<-left_join(colony_df, dft1, by = "Genus.Species")
full_Set<-left_join(full_Set, queen_df, by = "Genus.Species")


#Now calculate QWD metrics here
full_Set1<-full_Set %>%
  mutate(QWD = (wl/WL)*100,
         QWD_n = 100 * ( (2 * (wl - WL)) / (wl + WL)),
         QWD_HL = 100 * ( (2 * (hl - HL)) / (hl + HL)),
         QWD_EL = 100 * ( (2 * (el - EL)) / (el + EL)),
         QWD_ML = 100 * ( (2 * (ml - ML)) / (ml + ML)),
         QWD_HW = 100 * ( (2 * (hw - HW)) / (hw + HW)),

         QWD_HL2 = ((hl/wl)/(HL/WL)),
         QWD_EL2 = ((el/wl)/(EL/WL)),
         QWD_ML2 = ((ml/wl)/(ML/WL)),
         QWD_HW2 = ((hw/wl)/(HW/WL)),

         )


#let's get climate data



#First read in point data from Booher_et_al_2023_Data
tot_ants<-read.csv("Data/Spatial_Data/AntsFloridaTotalDataPaper.csv")

#Clean up the data geographically by unique combinations of:
#'
#'Source,locations, time, date of collection (start or end), location by gps
#'and species


Ant_species_fl<-tot_ants %>%
  group_by(collection.from,All.Locality,LocatedAt,DateCollectedStart,DateCollectedEnd, Day, Month,
           LocLatitude, LocLongitude, Genus, Species) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  dplyr::select(Genus, Species, LocLatitude, LocLongitude, n)  %>%
  drop_na(LocLongitude,LocLatitude) %>%
  filter(!LocLatitude < 19) %>%
  filter(!LocLongitude > 0) %>%
  mutate(Genus_species = paste(Genus, Species, sep = ".")) %>%
  dplyr::select(Genus_species, LocLongitude, LocLatitude)

#Here we clean up the names
Ant_species_fl<-Ant_species_fl %>%
  filter(! Genus_species %in% c(".", "Aphaenogaster.sp.",
                                "Brachymyrmex.sp", "Brachymyrmex.Species",
                                "Brachymyrmex.sp.", "Brachymyrmex.steinheili?",
                                "Cardiocondyla.wroughtonii_cf",
                                "Cardiocondyla.sp.","Colobopsis.impressus?",
                                "Crematogaster.", "Crematogaster.pinicola?",
                                "Dolichoderus.", "Dorymyrmex.", "Dorymyrmex.olsoni?",
                                "Dorymyrmex.flavus?", "Dorymyrmex.sp." , "Dorymyrmex.sp.",
                                "Forelius.sp.nov.", "Forelius.sp.",
                                "Forelius.sp.nov.", "missing.", "Monomorium.",
                                "Myrmecina.Sp." ,"no ants." ,"no ants.no ants",
                                "Nylanderia.", "Nylanderia.parvula?", "Nylanderia.pubens?",
                                "Nylanderia.sp", "Nylanderia.sp.","Nylanderia.sp.2",
                                "Nylanderia.sp.4", "Pachycondyla.sti", "Pheidole.",
                                "Pheidole.?male", "Pheidole.bilimeki?", "Pheidole.crassicornis_nr" ,
                                "Pheidole.sp.", "Polyergus.", "Proceratium.silaceum_nr", "Solenopsis.",
                                "Solenopsis.?", "Solenopsis.bd", "Solenopsis.andi", "Solenopsis.carolinensis_nr",
                                "Solenopsis.sp.", "Solenopsis.sp.1", "Solenopsis.sp.2", "Solenopsis.sp.A",
                                "Solenopsis.sp.B", "Solenopsis.texana_nr", "Strumigenys.", "Strumigenys.?",
                                "Strumigenys.(headmissing)", "Strumigenys.boltoni?", "Strumigenys.sp.",
                                "Temnothorax.texana_nr","Aphaenogaster.miamaiana/carolinensis",
                                "Aphaenogaster.miamiana/carolinensis", "Camponotus.pavidus/decipiens",
                                "Lasius.sp", "Pheidole.floridana_nr", "Solenopsis.diplo",
                                "Solenopsis.carolinensis/abdita",
                                "Strumigenys.inopina?","Strumigenys.inopina?",
                                "Strumigenys.talpa/deyrupi" ))

#using sort unique as a counter of sp in data
sort(unique(Ant_species_fl$Genus_species))


Ant_species_fl$Genus_species[Ant_species_fl$Genus_species == "Aphaenogaster.floridanus"]<-"Aphaenogaster.floridana"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species == "Camponotus.abddominalis"]<-"Camponotus.abdominalis"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species == "Camponotus.abdominalis floridanus"]<-"Camponotus.abdominalis"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Camponotus.flordanus", "Camponotus.floricanus",
                                                                 "Camponotus.floridana")]<-"Camponotus.floridanus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species == "Camponotus.inequalis"]<-"Camponotus.inaequalis"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species == "Camponotus.pennsylvanica"]<-"Camponotus.pennsylvanicus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species == "Camponotus.snelllingi"]<-"Camponotus.snellingi"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Colobopsis.impressus",
                                                                 "Colobobsis.impressus")]<-"Colobopsis.impressa"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Crematogaster.missuriensis")]<-"Crematogaster.missouriensis"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Cyphomyrmex.minut")]<-"Cyphomyrmex.minutus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Dorymyrmex.bossuta")]<-"Dorymyrmex.bossutus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Dorymyrmex.flavopecta")]<-"Dorymyrmex.flavopectus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Eurhopalothrix.floridanus","Eurhopalothrix.foridana")]<-"Eurhopalothrix.floridana"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Forelius.pruniosus")]<-"Forelius.pruinosus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Formica.palledefulva")]<-"Formica.pallidefulva"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Myrmecina.punctiventris" )]<-"Myrmica.punctiventris"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Nyanderia.arenivaga" )]<-"Nylanderia.arenivaga"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Paratrechina.longicornus" )]<-"Paratrechina.longicornis"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Pheidole.denatata",
                                                                 "Pheidole.dent")]<-"Pheidole.dentata"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Pheidole.flavus")]<-"Pheidole.flavens"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Pheidole.floidana",
                                                                 "Pheidole.floriana",
                                                                 "Pheidole.floridanus")]<-"Pheidole.floridana"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Prinopelta.antillana")]<-"Prionopelta.antillana"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Pseudomyrmex.pallida")]<-"Pseudomyrmex.pallidus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c( "Solenopsis.einvicta" )]<-"Solenopsis.invicta"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c( "Strumigenys.egg")]<-"Strumigenys.eggersi"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c( "Strumigenys.emm")]<-"Strumigenys.emmae"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Strumigenys.laevineasis")]<-"Strumigenys.laevinasis"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Tapinoma.littoralis")]<-"Tapinoma.litorale"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Temnothorax.texana")]<-"Temnothorax.texanus"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Wasmania.auropunctata")]<-"Wasmannia.auropunctata"
Ant_species_fl$Genus_species[Ant_species_fl$Genus_species %in% c("Xenomyrmex.floridana")]<-"Xenomyrmex.floridanus"
sort(unique(Ant_species_fl$Genus_species))





#Now get climate data for the geographical point data
library(raster)
MAT = raster("Data/Climate_data/World_clim/wc2.1_30s_bio_1.tif")#get mean annual temp
ATR = raster("Data/Climate_data/World_clim/wc2.1_30s_bio_7.tif")#get mean annual temp range
MAP = raster("Data/Climate_data/World_clim/wc2.1_30s_bio_12.tif")#get precipitation

Ant_species_fl$MAT<-raster::extract(MAT, Ant_species_fl[,2:3])
Ant_species_fl$ATR<-raster::extract(ATR, Ant_species_fl[,2:3])
Ant_species_fl$MAP<-raster::extract(MAP, Ant_species_fl[,2:3])

#save cleaned geographic point data
saveRDS(Ant_species_fl, "Data/Point_Data/cleaned_ant_point_data.rds")

#save trait data
saveRDS(full_Set1, "Data/trait_clean.rds")

#create environmental species -based data
Ant_species_fl_means<-Ant_species_fl %>%
  group_by(Genus_species) %>%
  summarise(
    mean_mat = mean(MAT, na.rm =T),
    mean_atr = mean(ATR, na.rm =T),
    mean_map = mean(MAP, na.rm =T),
    mean_cv_evi = mean(cv_EVI, na.rm =T),
    range_mat = max(MAT, na.rm = T)-min(MAT, na.rm=T),
    range_map = max(MAP, na.rm = T)-min(MAP, na.rm=T),
    wind = mean(wind, na.rm = T)

  )



#final data for all species
full_Set2<-left_join(full_Set1, Ant_species_fl_means, by = c("Genus.Species"= "Genus_species"))

#save species-level data
saveRDS(full_Set2, "Data/Trait_Data/species_level_data_final.rds")

