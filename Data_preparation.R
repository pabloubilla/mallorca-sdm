### ---- Libraries ----###
library(tidyverse)
library(sf)
library(terra)

setwd("~/Documents/mallorca")

### --- Load data --- ###

#### Biodibal
Biodibal_data=read_sf('/Users/ivancortesfernandez/Library/CloudStorage/OneDrive-UniversitatdelesIllesBalears/Estancias/Montpellier/Marzo 2026/biodibal_terrestrial_25831.gpkg')
Biodibal_data$coordinateUncertaintyInMeters=1 ## revisar
# biodibal NO contiene Floras
names(Biodibal_data)
str(Biodibal_data)
Biodibal_data$eventDate <- as.POSIXct(
  Biodibal_data$eventDate,
  format = "%b %d, %Y @ %H:%M:%OS",
  tz = "UTC"
)

### GBIF


GBIF_data=read_sf('/Users/ivancortesfernandez/Library/CloudStorage/OneDrive-UniversitatdelesIllesBalears/Estancias/Montpellier/Marzo 2026/GBIF_Mallorca.gpkg')
# 281299
GBIF_data$taxonRank
GBIF_data=GBIF_data%>%
  filter(taxonRank %in% c("SPECIES","SUBSPECIES","VARIETY"))
#278000
GBIF_data_clean=GBIF_data%>%
  filter(coordinateUncertaintyInMeters<1000 & !is.na(coordinateUncertaintyInMeters))## Ojo a este paso, es un hiperparámetro
GBIF_data_def=GBIF_data_clean%>%
  filter(coordinateUncertaintyInMeters!=707)%>%
  select(eventDate,scientificName,coordinateUncertaintyInMeters)

GBIF_data_def$eventDate=as.POSIXct(
  GBIF_data_def$eventDate,
  format = "%b %d, %Y @ %H:%M:%OS",
  tz = "UTC"
)

library(dplyr)

library(sf)
library(dplyr)

GBIF_data_filtrado <- GBIF_data_def %>%
  mutate(geom_wkt = st_as_text(geom)) %>%
  st_drop_geometry()%>%
  anti_join(
    Biodibal_data %>% mutate(geom_wkt = st_as_text(geom))%>%
      st_drop_geometry(),
    by = c(
      setdiff(names(GBIF_data_def), "geom"),
      "geom_wkt"
    )
  ) %>%
  select(-geom_wkt)### sin coincidencias, asi que podemos usar GBIF_data_def
  
### FLORAs

potential_flora_data=GBIF_data_clean%>%
  filter(coordinateUncertaintyInMeters==707)

potential_flora_data$eventDate=as.POSIXct(
  potential_flora_data$eventDate,
  format = "%b %d, %Y @ %H:%M:%OS",
  tz = "UTC"
)

str(potential_flora_data)

UTM_net=read_sf('/Users/ivancortesfernandez/Library/CloudStorage/OneDrive-UniversitatdelesIllesBalears/Estancias/Montpellier/Marzo 2026/flora_net.gpkg')
Buffer_centroids=read_sf('/Users/ivancortesfernandez/Library/CloudStorage/OneDrive-UniversitatdelesIllesBalears/Estancias/Montpellier/Marzo 2026/buffer_centroids.gpkg')
#Buffers centroids ensures that not other CRS points (ED50) are considered as PAU
Flora_net=UTM_net%>%
  filter(Flora==1)
crs(potential_flora_data)

Flora_data=st_join(potential_flora_data,Flora_net)%>%
  filter(Flora==1)%>%
  select(eventDate,scientificName,geom)
Not_finished_flora_data_raw=st_join(potential_flora_data,Flora_net)%>%
  filter(is.na(Flora))%>%
  select(eventDate,scientificName,geom)

# if we stop here, some ED50 data and other points can be included as PAU, 


#Not_finished_flora_data=Not_finished_flora_data_raw%>%
#  st_join(Buffer_centroids)%>%
#  filter(!is.na(UTMCODE1X1))%>%
#  select(eventDate,scientificName,geom)

# It could be explored if than we could weight PAU by number of species
# as PAU with 20 species is not equally sampled than PAU with 200 species

#Flora_POV_data=Not_finished_flora_data_raw%>%
#  st_join(Buffer_centroids)%>%
#  filter(is.na(UTMCODE1X1))%>%
#  select(eventDate,scientificName,geom)

# we turn Flora POV into GBIF
#str(GBIF_data_def)
#str(Flora_POV_data)
#GBIF_data_def=rbind(GBIF_data_def,Flora_POV_data)

### final datasets
#Inventory_data Future apporximation
Flora_data #55888 PAComplete
str(Flora_data)
Not_finished_flora_data_raw #47,537 PAUncomplete
str(Not_finished_flora_data)
Biodibal_data #90931 POValidated
str(Biodibal_data)
GBIF_data_def #127363 POUnvalidated
str(GBIF_data_def)

### coincidences 

normalize_sf <- function(x, digits = 5) {
  coords <- st_coordinates(x)
  x %>%
    mutate(
      lon = round(coords[,1], digits),
      lat = round(coords[,2], digits)
    ) %>%
    st_drop_geometry()
}

#Inventory_n  <- normalize_sf(Inventory_data)
Flora_n      <- normalize_sf(Flora_data)
NotFlora_n   <- normalize_sf(Not_finished_flora_data_raw)
Biodibal_n   <- normalize_sf(Biodibal_data)
GBIF_n       <- normalize_sf(GBIF_data_def)




matches <- list(

  Flora_Biodibal       = inner_join(Flora_n, Biodibal_n,
                                    by = c("scientificName","lon","lat")),
  Flora_GBIF           = inner_join(Flora_n, GBIF_n,
                                    by = c("scientificName","lon","lat")),
  NotFlora_Biodibal    = inner_join(NotFlora_n, Biodibal_n,
                                    by = c("scientificName","lon","lat")),
  NotFlora_GBIF        = inner_join(NotFlora_n, GBIF_n,
                                    by = c("scientificName","lon","lat")),
  Biodibal_GBIF        = inner_join(Biodibal_n, GBIF_n,
                                    by = c("scientificName","lon","lat"))
)

sapply(matches, nrow)

### --- Merge data --- ###

Flora_data$class="PAC"
Flora_data$coordinateUncertaintyInMeters=707 # ojo a las 100x100 que pueda haber
names(Flora_data)

Not_finished_flora_data_raw$class="PAU"
Not_finished_flora_data_raw$coordinateUncertaintyInMeters=707
names(Not_finished_flora_data_raw)

Biodibal_data$class="POV"
names(Biodibal_data)

GBIF_data_def$class="POU"
names(GBIF_data_def)

full_data=rbind(Flora_data,Not_finished_flora_data_raw,Biodibal_data,GBIF_data_def)%>%
  arrange(scientificName)



### --- Unify scientific names --- ###

# Sacamos el taxonkey de cada scientificname, y de ese el acceptedname
library(rgbif)
taxa_levels=levels(as.factor(full_data$scientificName))

## 1) Obtener el taxonKey para cada scientificName
taxonkey_list <- data.frame(
  initial_name = character(),
  accepted_key = numeric(),
  accepted_rank=character(),
  accepted_name = character(),
  accepted_order=character(),
  accepted_family=character(),
  stringsAsFactors = FALSE
)

#i=taxa_levels[185]

for (i in taxa_levels) {
  
  print(i)
  
  ## 1) Match inicial
  res <- name_backbone(
    name = i,
    rank = "species"
  )
  
  ## Si no hay coincidencia válida, saltar
  if (is.null(res$usageKey) || is.na(res$usageKey) || res$rank=="GENUS" || res$rank=="FAMILY" || res$rank=="CLASS") {
    print("falla")
    taxon_sel <- data.frame(
      initial_name = i,
      accepted_key = NA,
      accepted_rank="OTHER",
      accepted_name = NA,
      accepted_order=NA,
      accepted_family=NA,
      stringsAsFactors = FALSE
    )
    taxonkey_list <- bind_rows(taxonkey_list, taxon_sel)
    next
  }
  
  ## 2) Determinar la clave aceptada
  accepted_key <- if (!is.null(res$acceptedUsageKey) && !is.na(res$acceptedUsageKey)) {
    res$acceptedUsageKey
  } else {
    res$usageKey
  }
  
  ## 3) Recuperar el nombre aceptado a partir de la clave
  acc_res <- name_usage(key = accepted_key)
  
  accepted_name <- acc_res$data$scientificName
  
  ## 4) Guardar resultado
  taxon_sel <- data.frame(
    initial_name = i,
    accepted_key = accepted_key,
    accepted_rank=res$rank,
    accepted_name = accepted_name,
    accepted_order=res$order,
    accepted_family=res$family,
    stringsAsFactors = FALSE
  )
  
  taxonkey_list <- bind_rows(taxonkey_list, taxon_sel)
}


taxonkey_list_clean=taxonkey_list%>%
  filter(!accepted_rank %in% c("CLASS","OTHER","GENUS","FAMILY"))%>%
  filter(!accepted_order %in% c("Bryales","Bryopsidales","Ceramiales",
                                "Cladophorales",
                                "Corallinales","Dasycladales",
                                "Dicranales","Fossombroniales","Funariales",
                                "Hypnales","Jungermanniales","Marchantiales",
                                "Nemaliales","Pelliales","Peyssonneliales",
                                "Porellales","Pottiales","Scouleriales",
                                "Siphonocladales"))%>%
  filter(!accepted_family %in% c("Posidoniaceae","Potamogetonaceae",
                                 "Ulvaceae","Phyllophoraceae",
                                 "Ceratophyllaceae","Anadyomenaceae",
                                 "Zosteraceae","Ruppiaceae","Bonnemaisoniaceae",
                                 "Cymodoceaceae","Hydrocharitaceae"))%>%
  rename(scientificName=initial_name)%>%
  droplevels()

length(levels(as.factor(taxonkey_list_clean$scientificName)))
levels(as.factor(taxonkey_list_clean$accepted_family))
length(levels(as.factor(taxonkey_list_clean$accepted_name)))#2306 especies en total
levels(as.factor(taxonkey_list_clean$accepted_name))
full_data_clean=full_data%>% #321715 ocurrences
  filter(scientificName %in% taxonkey_list_clean$scientificName)%>%# 320784
  left_join(taxonkey_list_clean%>%
              select(scientificName,accepted_name),by="scientificName")%>%
  select(-scientificName,eventDate)%>%
  rename(scientificName=accepted_name)


### ----- Final Save data ---- ###

st_write(full_data_clean,"full_data_clean.gpkg",append=FALSE)
