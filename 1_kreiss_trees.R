library(terra)
library(dplyr)
library(ggplot2)
library(cowplot)
library(sf)
rm(list = ls())

setwd("C:/Users/madbro/Universität Zürich UZH/GEO_410/data/")


# read data
trees <- vect("data_trees/trees.shp")
kreis <- vect("data_kreiss/data/kreis.shp")




# join layers 
trees_kreis <- terra::intersect(trees, kreis)
trees_df <- as.data.frame(trees_kreis)



# calculate kreis area
kreis_area <- data.frame(
  kname    = kreis$kname,
  area_km2 = expanse(kreis, unit = "km")
) %>%
  group_by(kname) %>%
  summarise(area_km2 = sum(area_km2))


# tree per kreis
counts_kreis <- trees_df %>%
  count(kname, name = "n_trees")


#species per kreiss
counts_species <- trees_df %>%
  count(kname, baumnamela, name = "n_trees") %>%
  arrange(kname, desc(n_trees))


#full summary table 
kreis_summary <- trees_df %>%
  filter(!is.na(baumnamela)) %>%
  count(kname, baumnamela, name = "n") %>%
  group_by(kname) %>%
  mutate(p = n / sum(n)) %>%
  summarise(
    n_trees   = sum(n),                   # total trees 
    n_species = n(),                      # species number
    shannon   = -sum(p * log(p)),         # Shannon H'- more sensitive to rare
    simpson   = 1 - sum(p^2)              # Gini-Simpson -more sensitive to comoon 
  ) %>%
  left_join(kreis_area, by = "kname") %>%
  mutate(trees_per_km2 = n_trees / area_km2) %>%
  arrange(kname)

kreis_summary


# join summary to kreis polygons
kreis_map <- merge(kreis, kreis_summary, by = "kname")
kreis_sf <- st_as_sf(kreis_map)


# kreis colours
kreis_cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")

# order kreise by number and give each one a colour
kreis_order <- kreis_summary$kname[order(as.numeric(gsub("\\D", "", kreis_summary$kname)))]
names(kreis_cols) <- kreis_order

kreis_sf$kname <- factor(kreis_sf$kname, levels = kreis_order)
kreis_summary$kname <- factor(kreis_summary$kname, levels = kreis_order)


# plot theme, black outline and no grid lines
theme_kreis <- theme_cowplot() +
  theme(
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid      = element_blank(),
    legend.position = "none"
  )


# map of kreise with lat and lon coordinates, no axis titles
p_map <- ggplot(kreis_sf) +
  geom_sf(aes(fill = kname), colour = "black", linewidth = 0.5) +
  geom_sf_text(aes(label = gsub("\\D", "", kname)), size = 4, fontface = "bold") +
  scale_fill_manual(values = kreis_cols) +
  coord_sf(crs = 4326, datum = 4326) +                       # show axes in lat / lon
  scale_x_continuous(breaks = seq(8.45, 8.65, by = 0.05),
                     labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(breaks = seq(47.32, 47.44, by = 0.04),
                     labels = function(y) paste0(y, "°N")) +
  theme_kreis +
  theme(
    axis.title  = element_blank(),                           # no "Longitude" / "Latitude" titles
    axis.text.x = element_text(size = 8)
  )


# number of trees
p_trees <- ggplot(kreis_summary, aes(x = kname, y = n_trees, fill = kname)) +
  geom_col(colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = kreis_cols) +
  scale_x_discrete(labels = function(x) gsub("\\D", "", x)) +
  labs(x = "Kreis", y = "Number of trees") +
  theme_kreis


# trees per km2
p_density <- ggplot(kreis_summary, aes(x = kname, y = trees_per_km2, fill = kname)) +
  geom_col(colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = kreis_cols) +
  scale_x_discrete(labels = function(x) gsub("\\D", "", x)) +
  labs(x = "Kreis", y = "Trees per km²") +
  theme_kreis


# number of species
p_species <- ggplot(kreis_summary, aes(x = kname, y = n_species, fill = kname)) +
  geom_col(colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = kreis_cols) +
  scale_x_discrete(labels = function(x) gsub("\\D", "", x)) +
  labs(x = "Kreis", y = "Number of species") +
  theme_kreis


# shannon index
p_shannon <- ggplot(kreis_summary, aes(x = kname, y = shannon, fill = kname)) +
  geom_col(colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = kreis_cols) +
  scale_x_discrete(labels = function(x) gsub("\\D", "", x)) +
  labs(x = "Kreis", y = "Shannon index") +
  theme_kreis


# simpson index
p_simpson <- ggplot(kreis_summary, aes(x = kname, y = simpson, fill = kname)) +
  geom_col(colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = kreis_cols) +
  scale_x_discrete(labels = function(x) gsub("\\D", "", x)) +
  labs(x = "Kreis", y = "Simpson index") +
  theme_kreis


# combine bar plots into a 3 x 2 grid
p_bars <- plot_grid(p_trees, p_density, p_species, p_shannon, p_simpson,
                    ncol = 3, labels = c("B", "C", "D", "E", "F"))

# combine map and bar plots
p_final <- plot_grid(p_map, p_bars,
                     ncol = 2, rel_widths = c(1, 2), labels = c("A", ""))

p_final


# save
ggsave("kreis_panels.png", p_final, width = 16, height = 7, dpi = 300)

write.csv(kreis_summary,  "kreis_summary.csv",           row.names = FALSE)
write.csv(counts_species, "species_counts_by_kreis.csv", row.names = FALSE)