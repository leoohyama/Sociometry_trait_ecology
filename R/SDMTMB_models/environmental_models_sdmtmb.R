#Models for environmental predictors

library(tidyverse)
library(sf)
library(dismo)
library(ggsflabel)
library(bbmle)
library(ape)
library(sdmTMB)
library(performance)
library(patchwork)
library(exactextractr)

#load Florida map

fl_map<-st_read("~/Downloads/Detailed_Florida_State_Boundary/Detailed_Florida_State_Boundary.shp")
fl<-fl_map #make a safe copy
fl <- fl %>% st_transform(crs = 2031)

#load in ant point data
ant_pts<-readRDS("Data/Point_Data/cleaned_ant_point_data.rds")

#load ant traits
full_Set1<-readRDS("Data/trait_clean.rds")

ant_pts<-fortify(ant_pts) %>% st_as_sf(coords = c("LocLongitude", "LocLatitude"), crs = 4326) %>%
  st_transform(crs = 2031)

####make the hexgrid

make_grid2 <- function(x, cell_diameter, cell_area, clip = FALSE) {
  if (missing(cell_diameter)) {
    if (missing(cell_area)) {
      stop("Must provide cell_diameter or cell_area")
    } else {
      cell_diameter <- sqrt(2 * cell_area / sqrt(3))
    }
  }
  ext <- as(raster::extent(x) + cell_diameter, "SpatialPolygons")
  raster::projection(ext) <- raster::projection(x)
  # generate array of hexagon centers
  g <- st_make_grid(ext, cellsize = cell_diameter, crs = st_crs(ext),
                    square = FALSE) %>%
    st_transform(2031)

  # clip to boundary of study area
  if (clip) {
    g <- st_intersection(g, st_make_valid(x))
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs

  return(g)
}

#make hexgrid at 2K m2 note that units are in 'meters' based on proj.
hex_grid<-make_grid2(fl, cell_area = 2e9, clip = TRUE)

hex_grid<-st_as_sf(hex_grid) %>%
  mutate(ID = 1:nrow(.))
hex_grid<-hex_grid[-1,]

#recast this as a multipolgyon rather than single polygon layer
hex_grid<-hex_grid %>% st_cast("MULTIPOLYGON")

#remove hexgrids that are less than 25% of area of normal hexbin
hex_grid$area<-as.numeric(st_area(hex_grid))

#now intersect the points with grid
ant_int2<-st_intersection(ant_pts, hex_grid)


MAT = raster("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_1.tif")#get mean annual temp
ATR = raster("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_7.tif")#get mean annual temp range
MAP = raster("~/Downloads/wc2.1_30s_bio/wc2.1_30s_bio_12.tif")#get precipitation
cv_EVI = raster("~/Downloads/drive-download-20210713T202837Z-001/cv_01_05_1km_uint16_resampled.tif")
wind = raster("~/Downloads/wc2.1_30s_wind/wc2.1_30s_wind_07.tif")

hex_grid$MAT<- exact_extract(MAT, hex_grid, c("mean"))
hex_grid$ATR<- exact_extract(ATR, hex_grid, c("mean"))
hex_grid$MAP<- exact_extract(MAP, hex_grid, c("mean"))
hex_grid$cv_EVI<- exact_extract(cv_EVI, hex_grid, c("mean"))
hex_grid$wind<- exact_extract(wind, hex_grid, c("mean"))


#get rid of singletons
single<-ant_int2 %>%
  group_by(Genus_species) %>%
  summarise(n = n()) %>%
  filter(n %in% c(1)) %>%
  pull(Genus_species)

#filter out duplicate species for each bin
ant_assemblage<-ant_int2 %>%
  filter(!Genus_species %in% single) %>%
  group_by(ID) %>%
  filter(!duplicated(Genus_species)) %>%
  dplyr::select(ID, Genus_species) %>%
  left_join(., full_Set1, by = c("Genus_species" = "Genus.Species")) %>%
  st_drop_geometry()



srs<-ant_assemblage %>%
  group_by(ID) %>%
  summarise(SR = n())

#let's do some subsampling for native ants!
#first find distribution of SR
ant_assemblage_n<-ant_assemblage %>%
  filter(status == "Native")


#remove any cells with less than 5 species observed
removecells<-ant_assemblage_n %>%
  filter(status == "Native") %>%
  group_by(ID) %>%
  summarise(SR = n()) %>%
  arrange(SR) %>%
  filter(SR<5) %>%
  pull(ID)

#species list for native ants
species_list<-ant_assemblage_n %>%
  filter(!Genus_species %in% single) %>%
  drop_na(QWD) %>%
  dplyr::select(ID, Genus_species) %>%
  st_drop_geometry() %>%
  ungroup() %>%
  filter(!ID %in% removecells)


sp_list<- unique(species_list$Genus_species)

#so to account for spatial pseudreplication from wide-ranging species
#we only let a species contribute to half of the range size (hex bins)



set.seed(1236)
ant_prev<-ant_assemblage_n %>%
  filter(!ID %in% removecells) %>%
  group_by(Genus_species) %>%
  summarise(n = n(),
    perc = n()/length(unique(ant_assemblage$ID))) %>%
  ungroup() %>%
  mutate(cutoff = .5 * perc)


get_one<-function(n){
  sps<-subset(species_list, Genus_species == sp_list[n])
  if(nrow(sps)<2){
    size<-ceiling(nrow(sps)/2)
    sample_ant<-sample(1:nrow(sps), size, replace = F)
    sample_ant = sps[sample_ant,]

    return(sample_ant)
    rm(sps)
    rm(size)
  }else{
    #establish cutoff
    cutoff = ant_prev$cutoff[ant_prev$Genus_species ==  sp_list[n]]
    size<-ceiling(cutoff * nrow(sps))
    sample_ant<-sample(1:nrow(sps), size, replace = F)
    sample_ant = sps[sample_ant,]
    return(sample_ant)
    rm(sps)
    rm(size)
  }

}



list_ants<-as.list(1:1000)
for (i in 1:1000) {
  interm<-bind_rows(lapply(1:113, function(x){get_one(x)}))
  interm$iter = i
  list_ants[[i]]<-interm
}
subsampled_data<-bind_rows(list_ants)


ss_geo<-left_join(subsampled_data, full_Set1, by = c("Genus_species"= "Genus.Species"))


#QWD + colony size Figure for manuscript
qwdn<-as_tibble(ss_geo) %>%
  filter(status == "Native") %>%
  group_by(ID) %>%
  summarise(SR = n(),
            QWD_v= (sd(QWD, na.rm=T)),
            QWD = median(QWD, na.rm = T),
            median_size = median(median_size, na.rm = T)) %>%
  left_join(hex_grid, ., by = "ID") %>%
  ggplot() +
  geom_sf(aes(fill = QWD, color = QWD)) +
  scale_fill_viridis_c(option = "G",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_viridis_c(option = "G") +
  coord_sf()+
  theme_minimal() +
  theme(legend.position = 'bottom')


csnc<-as_tibble(ss_geo) %>%
  filter(status == "Native") %>%
  group_by(ID) %>%
  summarise(SR = n(),
            QWD_v= (sd(QWD, na.rm=T)),
            QWD = median(QWD, na.rm = T),
            median_size = median(median_size, na.rm = T)) %>%
  left_join(hex_grid, ., by = "ID") %>%
  ggplot() +
  geom_sf(aes(fill = log10(median_size),
              color =log10(median_size))) +
  scale_fill_viridis_c(option = "G",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_viridis_c(option = "G") +
  coord_sf()+
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(fill = "Median Colony Size (log)",
       color = "Median Colony Size (log)")
csnc + qwdn


#now do this with non natives
#first find distribution of SR
ant_assemblage_nn<-ant_assemblage %>%
  filter(!status == "Native")
#remove any cells with less than ten species observed
removecellsnn<-ant_assemblage_nn %>%
  group_by(ID) %>%
  summarise(SR = n()) %>%
  arrange(SR) %>%
  filter(SR<5) %>%
  pull(ID)



species_listnn<-ant_assemblage_nn %>%
  filter(!Genus_species %in% single) %>%
  drop_na(QWD) %>%
  dplyr::select(ID, Genus_species) %>%
  st_drop_geometry() %>%
  ungroup() %>%
  filter(!ID %in% removecellsnn)


sp_listnn<- unique(species_listnn$Genus_species)


ant_prevnn<-ant_assemblage_nn %>%
  filter(!ID %in% removecellsnn) %>%
  group_by(Genus_species) %>%
  summarise(n = n(),
            perc = n()/length(unique(ant_assemblage_nn$ID))) %>%
  ungroup() %>%
  mutate(cutoff = .5 * perc)

get_onenn<-function(n){

  sps<-subset(species_listnn, Genus_species == sp_listnn[n])
  if(nrow(sps)<2){
    size<-ceiling(nrow(sps)/2)
    sample_ant<-sample(1:nrow(sps), size, replace = F)
    sample_ant = sps[sample_ant,]

    return(sample_ant)
    rm(sps)
    rm(size)
  }else{
    #establish cutoff
    cutoff = ant_prev$cutoff[ant_prev$Genus_species ==  sp_list[n]]
    size<-ceiling(cutoff * nrow(sps))
    sample_ant<-sample(1:nrow(sps), size, replace = F)
    sample_ant = sps[sample_ant,]
    return(sample_ant)
    rm(sps)
    rm(size)
  }

}


list_antsnn<-as.list(1:1000)
for (i in 1:1000) {
  interm<-bind_rows(lapply(1:49, function(x){get_onenn(x)}))
  interm$iter = i
  list_antsnn[[i]]<-interm
}
subsampled_datann<-bind_rows(list_antsnn)


ss_geon_nonative<-left_join(subsampled_datann, full_Set1, by = c("Genus_species"= "Genus.Species"))


#QWD  + colony size for non natives
qwdnn<-as_tibble(ss_geon_nonative) %>%
  group_by(ID) %>%
  summarise(SR = n(),
            QWD_v= (sd(QWD, na.rm=T)),
            QWD = median(QWD, na.rm = T),
            median_size = median(median_size, na.rm = T)) %>%
  left_join(hex_grid, ., by = "ID") %>%
  ggplot() +
  geom_sf(aes(fill = QWD, color = QWD)) +
  scale_fill_viridis_c(option = "G",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_viridis_c(option = "G") +
  coord_sf()+
  theme_minimal()  +
  theme(legend.position = 'bottom')

csncnn<-as_tibble(ss_geon_nonative) %>%
  group_by(ID) %>%
  summarise(SR = n(),
            QWD_v= (sd(QWD, na.rm=T)),
            QWD = median(QWD, na.rm = T),
            median_size = median(median_size, na.rm = T)) %>%
  left_join(hex_grid, ., by = "ID") %>%
  ggplot() +
  geom_sf(aes(fill = log10(median_size),
              color =log10(median_size))) +
  scale_fill_viridis_c(option = "G",
                       guide = guide_colorbar(frame.colour = "black",
                                              ticks.colour = "black")) +
  scale_color_viridis_c(option = "G") +
  coord_sf()+
  theme_minimal() +
  theme(legend.position = 'bottom') +
  labs(fill = "Median Colony Size (log)",
       color = "Median Colony Size (log)")

#figure 4
csncnn + qwdnn


#######################SDMTMB MODELS
ss_geo_plot<-as_tibble(ss_geo) %>%
  group_by(iter,ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(median_size, na.rm =T)) %>%
  group_by(ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD_v, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(size, na.rm =T))

#this will allow us to grab coordinates for the models
gps_coords<-left_join(hex_grid,  ss_geo_plot, by = "ID")


bothstat<-as_tibble(ss_geo) %>%
  drop_na(status) %>%
  group_by(iter,ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(median_size, na.rm =T)) %>%
  group_by(ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD_v, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(size, na.rm =T))

bothstat1<-left_join(hex_grid,  bothstat, by = "ID")


xycoords<-st_coordinates(st_centroid(gps_coords))
gps_coords<-cbind(gps_coords,xycoords)


#here we divide the coordinates by 1000 to change the units from
#meters to kilometers (works better when fitting the mesh from the sdmtmb models)
gps_coords2<-gps_coords %>% st_drop_geometry() %>%
  mutate(X = X/1000, Y= Y/1000) %>%
  drop_na(QWD) %>%
  mutate(size = round(size))


#check for spatial autocorrelation with data before models
#using the global moran's I
library(ape)
inv_dists <- as.matrix(dist(gps_coords2[,c("X","Y")]))
diag(inv_dists) <- 0

#Significant effects of space found
Moran.I(gps_coords2$QWD, inv_dists)



#make mesh for sdmtmb QWD
mesh <- make_mesh(gps_coords2, xy_cols = c("X", "Y"), n_knots = 75, type = "cutoff_search")
plot(mesh)


m1<-sdmTMB::sdmTMB(data = gps_coords2, mesh = mesh, (QWD) ~ scale(ATR),
                   spatial = "on",family = Gamma(link = 'log'))
saveRDS(m1, "R/SDMTMB_models/QWD_native/QWD_ATR_sdmtmb_native.rds")

m2<-sdmTMB::sdmTMB(data = gps_coords2, mesh = mesh, (QWD) ~ scale(MAT) ,
                   spatial = "on",family = Gamma(link = 'log'))
saveRDS(m2, "R/SDMTMB_models/QWD_native/QWD_MAT_sdmtmb_native.rds")

m3<-sdmTMB::sdmTMB(data = gps_coords2, mesh = mesh, (QWD) ~ scale(MAP),
                   spatial = "on",family = Gamma(link = 'log'))
saveRDS(m3, "R/SDMTMB_models/QWD_native/QWD_MAP_sdmtmb_native.rds")

nullmod<-sdmTMB::sdmTMB(data = gps_coords2, mesh = mesh, (QWD) ~ 1,
                   spatial = "off",family = Gamma(link = 'log'))
saveRDS(nullmod, "R/SDMTMB_models/QWD_native/null_native.rds")



AICtab(m1,m2,m3,nullmod,
       weights = T, delta = T, base = TRUE, sort = TRUE)
AICQWD_native<-AICtab(m1,m2,m3,nullmod,
       weights = T, delta = T, base = TRUE, sort = TRUE)

#write.csv(AICQWD_native, "AIC_tables/AICQWD_native.csv")

summary(m2)
tidy(m2, conf.int = T)
s_gamma <- simulate(m2, nsim = 500)
pred_fixed <- m2$family$linkinv(predict(m2)$est_non_rf)

#using Dharma to examine model residuals
r_pois <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = gps_coords2$QWD,
  fittedPredictedResponse = pred_fixed
)

#plot residuals
plot(r_pois)

#examining p values and model estimates
summary(m2$sd_report, select = "fixed", p.value = TRUE)


#check residuals for spatial auto
#check for spatial autocorrelation with matern model
inv_dists <- as.matrix(dist(gps_coords2[,c("X","Y")]))
gps_coords2$residuals<-residuals(m2)
diag(inv_dists) <- 0
Moran.I(gps_coords2$residuals, inv_dists)



#Non native sdmtmb
bothstat<-as_tibble(ss_geon_nonative) %>%
  drop_na(status) %>%
  group_by(iter,ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(median_size, na.rm =T)) %>%
  group_by(ID) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD_v, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(size, na.rm =T))

bothstat1nn<-left_join(hex_grid,  bothstat, by = "ID")


#colony size models
xycoords<-st_coordinates(st_centroid(gps_coords))
gps_coords3<-cbind(gps_coords,xycoords)
gps_coords3<-gps_coords3 %>% st_drop_geometry() %>%
  mutate(X = X/1000, Y= Y/1000) %>%
  mutate(size = round(size)) %>%
  drop_na(size)


#check for spatial autocorrelation with matern model
inv_dists <- as.matrix(dist(gps_coords3[,c("X","Y")]))
diag(inv_dists) <- 0
Moran.I(gps_coords3$size, inv_dists)

#make mesh for sdmtmb colony size
mesh <- make_mesh(gps_coords3, xy_cols = c("X", "Y"), n_knots = 75, type = "cutoff_search")
plot(mesh)


cs1<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, size ~ scale(ATR),
                   spatial = "on", family = poisson(link = 'log'))
saveRDS(cs1, "R/SDMTMB_models/Size_native/CS_poisson_ATR_sdmtmb_native.rds")

cs1nb<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, (size) ~ scale(ATR),
                      spatial = "on",family = nbinom2(link = 'log'))
saveRDS(cs1nb, "R/SDMTMB_models/Size_native/CS_nb_ATR_sdmtmb_native.rds")

cs2<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, (size) ~ scale(MAT),
                   spatial = "on",family = poisson(link = 'log'))
saveRDS(cs2, "R/SDMTMB_models/Size_native/CS_poisson_MAT_sdmtmb_native.rds")

cs2nb<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, (size) ~ scale(MAT),
                    spatial = "on",family = nbinom2(link = 'log'))
saveRDS(cs2nb, "R/SDMTMB_models/Size_native/CS_nb_MAT_sdmtmb_native.rds")

cs3<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, (size) ~ scale(MAP),
                   spatial = "on",family = poisson(link = 'log'))
saveRDS(cs3, "R/SDMTMB_models/Size_native/CS_poisson_MAP_sdmtmb_native.rds")

cs3nb<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, (size) ~ scale(MAP),
                      spatial = "on",family = nbinom2())
saveRDS(cs3nb, "R/SDMTMB_models/Size_native/CS_nb_MAP_sdmtmb_native.rds")

nullmod<-sdmTMB::sdmTMB(data = gps_coords3, mesh = mesh, (size) ~ 1,
                        spatial = "off",family = poisson(link = 'log'))
saveRDS(nullmod, "R/SDMTMB_models/Size_native/CS_nullmod_sdmtmb_native.rds")


AICtab(cs1,cs1nb, cs2,cs2nb,cs3,cs3nb,nullmod,
       weights = T, delta = T, base = TRUE, sort = TRUE)
AICCS_native<-AICtab(cs1,cs1nb, cs2,cs2nb,cs3,cs3nb,nullmod,
       weights = T, delta = T, base = TRUE, sort = TRUE)


write.csv(AICCS_native, "AIC_tables/AIC_CS_native.csv")


summary(cs2nb)
summary(cs2nb$sd_report, select = "fixed", p.value = TRUE)



tidy(cs2nb, conf.int = T,exponentiate = T)

visreg::visreg(cs2nb, xvar = "MAT", scale = "response")


s_gamma <- simulate(cs2nb, nsim = 500)
pred_fixed <- cs2nb$family$linkinv(predict(cs2nb)$est_non_rf)
r_pois <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = gps_coords3$size,
  fittedPredictedResponse = pred_fixed
)
plot(r_pois)

#Nonnatives sdmtmb
ss_geo_plot2<-as_tibble(ss_geon_nonative) %>%
  group_by(iter,ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(median_size, na.rm =T)) %>%
  group_by(ID, status) %>%
  summarise(SR = n(),
            QWD_v = sd(QWD_v, na.rm=T),
            QWD = median(QWD, na.rm = T),
            WL = median(WL, na.rm =T),
            wl = median(wl, na.rm = T),
            size = median(size, na.rm =T))


gps_coordsnn<-left_join(hex_grid,  ss_geo_plot2, by = "ID")


#SDMTMB MODELS
xycoords<-st_coordinates(st_centroid(gps_coordsnn))
gps_coordsnn2<-cbind(gps_coordsnn,xycoords)
st_crs(gps_coordsnn2)
gps_coordsnn2<-gps_coordsnn2 %>% st_drop_geometry() %>%
  mutate(X = X/1000, Y= Y/1000) %>%
  drop_na(QWD) %>%
  mutate(size = round(size))
#make mesh for sdmtmb QWD
mesh <- make_mesh(gps_coordsnn2, xy_cols = c("X", "Y"), n_knots = 75, type = "cutoff_search")
plot(mesh)




m1<-sdmTMB::sdmTMB(data = gps_coordsnn2, mesh = mesh, (QWD) ~ scale(ATR),
                   spatial = "on",family = Gamma(link = 'log'))
saveRDS(m1, "QWD_ATR_sdmtmb_nonnative.rds")

m2<-sdmTMB::sdmTMB(data = gps_coordsnn2, mesh = mesh, (QWD) ~ scale(MAT) ,
                   spatial = "on",family = Gamma(link = 'log'))
saveRDS(m2, "QWD_MAT_sdmtmb_nonnative.rds")

m3<-sdmTMB::sdmTMB(data = gps_coordsnn2, mesh = mesh, (QWD) ~ scale(MAP),
                   spatial = "on",family = Gamma(link = 'log'))
saveRDS(m3, "QWD_MAP_sdmtmb_nonnative.rds")

nullmod<-sdmTMB::sdmTMB(data = gps_coordsnn2, mesh = mesh, (QWD) ~ 1,
                        spatial = "off",family = Gamma(link = 'log'))
saveRDS(nullmod, "QWD_nullmod_sdmtmb_nonnative.rds")


AICtab(m1,m2,m3,nullmod,
         weights = T, delta = T, base = TRUE, sort = TRUE)
AICQWD_nonnative<-AICtab(m1,m2,m3,nullmod,
                      weights = T, delta = T, base = TRUE, sort = TRUE)
write.csv(AICQWD_nonnative, "AICQWD_nonnative.csv")

#colony size models

xycoords<-st_coordinates(st_centroid(gps_coordsnn))
gps_coords2nn<-cbind(gps_coordsnn,xycoords)
st_crs(gps_coords2nn)
gps_coords2nn<-gps_coords2nn %>% st_drop_geometry() %>%
  mutate(X = X/1000, Y= Y/1000) %>%
  mutate(size = round(size)) %>%
  drop_na(size)


#make mesh for sdmtmb colony size
mesh <- make_mesh(gps_coords2nn, xy_cols = c("X", "Y"), n_knots = 75, type = "cutoff_search")
plot(mesh)

gps_coords2nn$obs = c(1:nrow(gps_coords2nn))

cs1<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, size ~ scale(ATR),
                    spatial = "on", family = poisson(link = 'log'))

cs1nb<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, (size) ~ scale(ATR),
                      spatial = "on",family = nbinom2())

cs2<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, (size) ~ scale(MAT),
                    spatial = "on",family = poisson(link = 'log'))

cs2nb<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, (size) ~ scale(MAT),
                      spatial = "on",family = nbinom2())

cs3<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, (size) ~ scale(MAP),
                    spatial = "on",family = poisson(link = 'log'))

cs3nb<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, (size) ~ scale(MAP),
                      spatial = "on",family = nbinom2())


nullmod<-sdmTMB::sdmTMB(data = gps_coords2nn, mesh = mesh, (size) ~ 1,
                        spatial = "off",family = poisson(link = 'log'))


AICtab(cs1,cs1nb, cs2,cs2nb,cs3,cs3nb,nullmod, weights = T,
       delta = T, base = T)

summary(cs2nb)
summary(cs3nb)
summary(cs1nb)


summary(cs2nb)
summary(cs1nb$sd_report, select = "fixed", p.value = TRUE)



tidy(cs2nb, conf.int = T,exponentiate = F)
tidy(cs1nb, conf.int = T,exponentiate = F)

visreg::visreg(cs2nb, xvar = "MAT", scale = "response")


s_gamma <- simulate(cs2nb, nsim = 500)
pred_fixed <- cs2nb$family$linkinv(predict(cs2nb)$est_non_rf)
r_pois <- DHARMa::createDHARMa(
  simulatedResponse = s_gamma,
  observedResponse = wtf3nn$size,
  fittedPredictedResponse = pred_fixed
)
plot(r_pois)




