### ---- Libraries ----###
library(tidyverse)
library(sf)
library(terra)

# Unifica eventDate (GBIF ISO, "NA" como texto, formato GBIF legacy) antes de rbind
normalize_event_date <- function(x) {
  if (inherits(x, "POSIXct")) {
    return(x)
  }
  x_chr <- as.character(x)
  x_chr[is.na(x_chr) | x_chr %in% c("NA", "NULL", "")] <- NA_character_

  out <- suppressWarnings(
    as.POSIXct(x_chr, format = "%b %d, %Y @ %H:%M:%OS", tz = "UTC")
  )
  iso <- is.na(out) & !is.na(x_chr)
  if (any(iso)) {
    x_iso <- sub("Z$", "", sub("T", " ", x_chr[iso]))
    out[iso] <- suppressWarnings(
      as.POSIXct(x_iso, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
    )
  }
  out
}

### --- Load data --- ###

#### Biodibal
Biodibal_data=read_sf('./data/raw/biodibal_terrestrial_25831.gpkg')
Biodibal_data$coordinateUncertaintyInMeters=1 ## revisar
# biodibal NO contiene Floras
names(Biodibal_data)
str(Biodibal_data)
Biodibal_data$eventDate <- normalize_event_date(Biodibal_data$eventDate)

### Saez

Saez_data=read.csv('./data/raw/Puig_Major.csv',
sep = ";",
check.names = FALSE)%>%
gather(quad,presence,c("31SDE8206","31SDE8207"))%>%
mutate(presence=case_when(presence==1~1,TRUE~0),
XCENTROIDE=case_when(quad=="31SDE8206"~482507,quad=="31SDE8207"~483512),
YCENTROIDE=case_when(quad=="31SDE8206"~4406457,quad=="31SDE8207"~4407435))%>%
filter(presence==1)

Saez_data%>%
group_by(XCENTROIDE)%>%
summarise(n=n())

Saez_data_sf=st_as_sf(Saez_data,coords=c("XCENTROIDE","YCENTROIDE"),crs=25831)%>%
select(SP,presence)%>%
rename(scientificName=SP)

#st_write(Saez_data_sf,"Saez_data_sf.gpkg",append=FALSE)

str(Saez_data_sf)
dim(Saez_data_sf) #902

### GBIF


GBIF_data=read_sf('./data/raw/GBIF_Mallorca.gpkg')
# 281299
GBIF_data$taxonRank
GBIF_data=GBIF_data%>%
  filter(taxonRank %in% c("SPECIES","SUBSPECIES","VARIETY"))
#278000

GBIF_data_clean=GBIF_data%>%
  filter(coordinateUncertaintyInMeters<1000 & !is.na(coordinateUncertaintyInMeters))## Ojo a este paso, es un hiperparámetro


GBIF_data_clean%>%
  filter(collectionCode=="FV-MALLORCA")%>%
  ggplot(aes(x=coordinateUncertaintyInMeters))+
  geom_histogram()+
  theme_minimal()

nrow(GBIF_data_clean)

GBIF_data_def=GBIF_data_clean%>%
  filter(coordinateUncertaintyInMeters!=707)%>%
  filter(collectionCode!="FV-MALLORCA")%>%
  select(eventDate,scientificName,coordinateUncertaintyInMeters)

nrow(GBIF_data_def)

### En GBIF hay datos de POV, PAC y PAU

GBIF_data_def$eventDate <- normalize_event_date(GBIF_data_def$eventDate)

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
dim(GBIF_data_def)  

### FLORAs

potential_flora_data=GBIF_data_clean%>%
  filter(coordinateUncertaintyInMeters==707)## Las floras que vienen de nuestro dataset de GBIF

potential_flora_data$eventDate <- normalize_event_date(potential_flora_data$eventDate)

str(potential_flora_data)

UTM_net=read_sf('./data/raw/Flora_net.gpkg')
UTM_net
## creeate centroids for UTM_net
UTM_net_centroids=st_centroid(UTM_net)
## create a buffer of 100 meters around the centroids
UTM_net_centroids_buffer=st_buffer(UTM_net_centroids,dist=100)
#Buffers centroids ensures that not other CRS points (ED50) are considered as PAU
Flora_net=UTM_net%>%
filter(Flora==1)

dim(potential_flora_data)[1]
dim(Saez_data_sf)[1]

# Estandarizamos esquema y geometría antes del bind_rows
potential_flora_data_bind <- potential_flora_data %>%
  mutate(source = "GBIF") %>%
  select(eventDate, scientificName, geom, source)

Saez_data_sf_bind <- Saez_data_sf %>%
  st_set_geometry("geom") %>%
  mutate(eventDate = as.POSIXct(NA), source = "Saez") %>%
  select(eventDate, scientificName, geom, source)

full_flora_data <- bind_rows(potential_flora_data_bind, Saez_data_sf_bind)
str(full_flora_data)
dim(full_flora_data)[1] == dim(potential_flora_data_bind)[1] + dim(Saez_data_sf_bind)[1]
table(full_flora_data$source, useNA = "ifany")

#write_sf(full_flora_data,"full_flora_data.gpkg",append=FALSE)

Flora_data <- st_join(full_flora_data, Flora_net) %>%
  filter(Flora == 1) %>%
  select(eventDate, scientificName, geom)

#st_write(Flora_data,"Flora_data_raw.gpkg",append=FALSE)

## check sum

str(Flora_data)


str(Flora_data)
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
str(Not_finished_flora_data_raw)
Biodibal_data #90931 POValidated
str(Biodibal_data)
GBIF_data_def #46386 POUnvalidated
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
# Esquema común antes de combinar (Flora/PAU solo traían eventDate, scientificName, geom)
MERGE_COLS <- c(
  "eventDate", "scientificName", "geom",
  "class", "coordinateUncertaintyInMeters"
)

Flora_data <- Flora_data %>%
  select(eventDate, scientificName, geom) %>%
  mutate(
    class = "PAC",
    coordinateUncertaintyInMeters = 707
  ) %>%
  select(all_of(MERGE_COLS))

Not_finished_flora_data_raw <- Not_finished_flora_data_raw %>%
  select(eventDate, scientificName, geom) %>%
  mutate(
    class = "PAU",
    coordinateUncertaintyInMeters = 707
  ) %>%
  select(all_of(MERGE_COLS))

Biodibal_data <- Biodibal_data %>%
  mutate(class = "POV") %>%
  select(all_of(MERGE_COLS))

GBIF_data_def <- GBIF_data_def %>%
  mutate(class = "POU") %>%
  select(all_of(MERGE_COLS))

# Mismo tipo en eventDate; si no, rbind convierte "NA"/ISO a POSIXct y falla
Flora_data$eventDate <- normalize_event_date(Flora_data$eventDate)
Not_finished_flora_data_raw$eventDate <- normalize_event_date(
  Not_finished_flora_data_raw$eventDate
)
Biodibal_data$eventDate <- normalize_event_date(Biodibal_data$eventDate)
GBIF_data_def$eventDate <- normalize_event_date(GBIF_data_def$eventDate)

# rbind exige las mismas columnas (a diferencia de bind_rows, que rellena con NA)
stopifnot(
  identical(names(Flora_data), MERGE_COLS),
  identical(names(Not_finished_flora_data_raw), MERGE_COLS),
  identical(names(Biodibal_data), MERGE_COLS),
  identical(names(GBIF_data_def), MERGE_COLS)
)

full_data <- rbind(
  Flora_data,
  Not_finished_flora_data_raw,
  Biodibal_data,
  GBIF_data_def
) %>%
  arrange(scientificName)

  full_data%>%
  count(class)



### --- Unify scientific names --- ###

# Sacamos el taxonkey de cada scientificname, y de ese el acceptedname
library(rgbif)
normalize_scientific_name <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

full_data <- full_data %>%
  mutate(
    scientificName_original = scientificName,
    scientificName = normalize_scientific_name(scientificName)
  )

taxa_levels=sort(unique(full_data$scientificName))

## 1) Obtener el taxonKey para cada scientificName
taxonkey_list <- data.frame(
  initial_name = character(),
  accepted_key = numeric(),
  accepted_rank=character(),
  accepted_name = character(),
  accepted_species_key = numeric(),
  accepted_species_name = character(),
  accepted_order=character(),
  accepted_family=character(),
  stringsAsFactors = FALSE
)


download_taxonkey <- function() {
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
      accepted_species_key = NA,
      accepted_species_name = NA,
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
  acc_res <- name_usage(key = accepted_key, data = "all")
  acc_data <- acc_res$data
  
  accepted_name <- acc_data$scientificName
  
  ## 3.1) Forzar especie aceptada (sin subsp./var./form)
  if (!is.null(acc_data$rank) && acc_data$rank == "SPECIES") {
    accepted_species_key <- accepted_key
    accepted_species_name <- accepted_name
  } else if (!is.null(acc_data$speciesKey) && !is.na(acc_data$speciesKey)) {
    species_res <- name_usage(key = acc_data$speciesKey, data = "all")
    accepted_species_key <- species_res$data$key
    accepted_species_name <- species_res$data$scientificName
  } else {
    accepted_species_key <- NA
    accepted_species_name <- NA
  }
  
  ## 4) Guardar resultado
  taxon_sel <- data.frame(
    initial_name = i,
    accepted_key = accepted_key,
    accepted_rank=res$rank,
    accepted_name = accepted_name,
    accepted_species_key = accepted_species_key,
    accepted_species_name = accepted_species_name,
    accepted_order=res$order,
    accepted_family=res$family,
    stringsAsFactors = FALSE
  )
  
  taxonkey_list <- bind_rows(taxonkey_list, taxon_sel)
}
write_xlsx(taxonkey_list,"taxonkey_list.xlsx")
}

#donwload_taxonkey()




## In this phase we manually revise both failed and erroneous assignments

library(readxl)
taxonkey_list=read_xlsx("data/raw/taxonkey_list_review.xlsx")
str(taxonkey_list)
levels(as.factor(taxonkey_list$accepted_species_name))


taxonkey_list_clean <- taxonkey_list %>%
  filter(!accepted_rank %in% c("CLASS", "OTHER", "GENUS", "FAMILY")) %>%
  filter(!is.na(accepted_species_key)) %>%
  filter(!accepted_order %in% c(
    "Bryales", "Bryopsidales", "Ceramiales", "Cladophorales",
    "Corallinales", "Dasycladales", "Dicranales", "Fossombroniales", "Funariales",
    "Hypnales", "Jungermanniales", "Marchantiales", "Nemaliales", "Pelliales",
    "Peyssonneliales", "Porellales", "Pottiales", "Scouleriales", "Siphonocladales"
  )) %>%
  filter(!accepted_family %in% c(
    "Posidoniaceae", "Potamogetonaceae", "Ulvaceae", "Phyllophoraceae",
    "Ceratophyllaceae", "Anadyomenaceae", "Zosteraceae", "Ruppiaceae",
    "Bonnemaisoniaceae", "Cymodoceaceae", "Hydrocharitaceae"
  )) %>%
  filter(!is.na(accepted_key))

# Tabla de unión: cada variante de entrada -> taxón aceptado en GBIF
name_lookup <- taxonkey_list_clean %>%
  select(initial_name, accepted_key, accepted_name, accepted_species_key, accepted_species_name)

# Una fila por taxón aceptado (clave canónica para contar especies)
taxa_accepted <- taxonkey_list_clean %>%
  group_by(accepted_species_key) %>%
  summarise(
    accepted_name = first(accepted_species_name),
    accepted_family = first(accepted_family),
    n_variants = n_distinct(initial_name),
    .groups = "drop"
  )

# Variantes de entrada distintas (p. ej. con/sin autor en paréntesis)
length(unique(name_lookup$initial_name))
# Taxones aceptados únicos en GBIF (esto es el conteo taxonómico real)
length(unique(name_lookup$accepted_species_key))
levels(as.factor(taxa_accepted$accepted_family))

# Ejemplos: varias entradas que apuntan al mismo accepted_key
name_lookup %>%
  group_by(accepted_species_key, accepted_species_name) %>%
  filter(n() > 1) %>%
  summarise(variants = paste(sort(unique(initial_name)), collapse = " | "), .groups = "drop")

full_data_clean <- full_data %>%
  left_join(name_lookup, by = c("scientificName" = "initial_name")) %>%
  filter(!is.na(accepted_species_key)) %>%
  left_join(
    taxa_accepted %>% select(accepted_species_key, accepted_name_canonical = accepted_name),
    by = "accepted_species_key"
  ) %>%
  select(-scientificName, -accepted_name, -accepted_species_name) %>%
  rename(scientificName = accepted_name_canonical) %>%
  select(-accepted_key, -accepted_species_key, everything())

# Especies únicas en el dataset final (ya unificadas por GBIF)
length(unique(full_data_clean$scientificName))




### ----- Final Save data ---- ###
library(writexl)
st_write(full_data_clean,"data/raw/full_data_clean_saez.gpkg",append=FALSE)






