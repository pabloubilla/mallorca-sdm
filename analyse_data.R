library(sf)
library(ggplot2)
library(dplyr)

# load data
st_layers("full_data_clean.gpkg")
data <- st_read("full_data_clean.gpkg", layer = "full_data_clean")

## unique species and frequency
data$scientificName %>% unique() %>% length()
species_count <- data %>%
  group_by(scientificName) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


## plot data
ggplot(data) +
  geom_sf(size = 1, alpha = 0.6) +
  theme_minimal() 


# by different class (sources)
ggplot(data) +
  geom_sf(aes(color = class), size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Observations in Mallorca by class",
    color = "Class"
  )



### separate one dataset per class, names are PAC, PAU, POV, POU
df_pac <- data %>% filter(class == "PAC")
df_pau <- data %>% filter(class == "PAU")
df_pov <- data %>% filter(class == "POV")
df_pou <- data %>% filter(class == "POU")



classes <- c("PAC", "PAU", "POV", "POU")

counts <- data %>%
  st_drop_geometry() %>%
  filter(class %in% classes) %>%
  count(class)

df_plot <- data %>%
  filter(class %in% classes) %>%
  left_join(counts, by = "class") %>%
  mutate(
    class = factor(class, levels = classes),
    class_label = paste0(class, " (n = ", n, ")")
  )

ggplot(df_plot) +
  geom_sf(aes(color = class), size = 1.2, alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ class_label, ncol = 2) +
  scale_color_manual(values = c(
    PAC = "blue",
    PAU = "orange",
    POV = "green",
    POU = "red"
  )) +
  theme_minimal() +
  labs(title = "Observations in Mallorca by class")


