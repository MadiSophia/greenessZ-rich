library(terra)
library(dplyr)
library(ggplot2)
library(cowplot)
library(sf)
rm(list = ls())

setwd("C:/Users/madbro/Universität Zürich UZH/GEO_410/data/")


# read dat
trees <- vect("data_trees/trees.shp")
kreis <- vect("data_kreiss/data/kreis.shp")
sat   <- rast("sentinel/NDVI_EVI_NDWI_BAP_Zurich_July_2025.tif")

names(sat)                                   # check layer order first!

# one polygon per Kreis, everything in the CRS of the satellite image
kreis <- aggregate(kreis, by = "kname")
kreis <- project(kreis, crs(sat))
trees <- project(trees, crs(sat))


# thresholds derived- add citations 
ndwi_thr <- 0.25   # NDWI > thr => water (masked out)
ndvi_thr <- 0.25   # NDVI > thr => green
evi_thr  <- 0.25   # EVI  > thr => green


# ---- water mask ----
water <- sat[["NDWI"]] > ndwi_thr
ndvi  <- mask(sat[["NDVI"]], water, maskvalues = TRUE)
evi   <- mask(sat[["EVI"]],  water, maskvalues = TRUE)


# ---- remove trees that fall on water pixels ----
n_before   <- nrow(trees)
tree_water <- terra::extract(water, trees)[[2]]      # 1 = water, 0 = land, NA = outside raster
trees      <- trees[!(tree_water %in% TRUE), ]       # keep land + trees outside raster
cat("Trees removed (on water):", n_before - nrow(trees), "of", n_before, "\n")


# ---- join layers ----
trees_kreis <- terra::intersect(trees, kreis)
trees_df    <- as.data.frame(trees_kreis)


# ---- kreis area (total and land only, water excluded) ----
kreis_area <- data.frame(
  kname    = kreis$kname,
  area_km2 = expanse(kreis, unit = "km")
)

land_km2 <- cellSize(water, unit = "km") * (1 - water)   # area of land cells only
kreis_area$land_km2 <- terra::extract(land_km2, kreis, fun = sum, na.rm = TRUE)[[2]]


# ---- trees per kreis ----
counts_kreis <- trees_df %>%
  count(kname, name = "n_trees")


# ---- species per kreis ----
counts_species <- trees_df %>%
  count(kname, baumnamela, name = "n_trees") %>%
  arrange(kname, desc(n_trees))


# ---- greenness per kreis (water already masked) ----
stk <- c(water * 1,
         (ndvi > ndvi_thr) * 1,                 # NA over water => land only
         (evi  > evi_thr)  * 1)
names(stk) <- c("water", "green_ndvi", "green_evi")

green <- terra::extract(stk, kreis, fun = mean, na.rm = TRUE)
green <- data.frame(
  kname          = kreis$kname,
  water_pct      = green$water      * 100,
  green_ndvi_pct = green$green_ndvi * 100,
  green_evi_pct  = green$green_evi  * 100
)


# ---- full summary table ----
kreis_summary <- trees_df %>%
  filter(!is.na(baumnamela)) %>%
  count(kname, baumnamela, name = "n") %>%
  group_by(kname) %>%
  mutate(p = n / sum(n)) %>%
  summarise(
    n_trees   = sum(n),                   # total trees
    n_species = n(),                      # species number
    shannon   = -sum(p * log(p)),         # Shannon H' - more sensitive to rare
    simpson   = 1 - sum(p^2)              # Gini-Simpson - more sensitive to common
  ) %>%
  left_join(kreis_area, by = "kname") %>%
  left_join(green,      by = "kname") %>%
  mutate(trees_per_km2 = n_trees / land_km2) %>%   # density on land area only
  arrange(kname)

kreis_summary


# ---- join summary to kreis polygons ----
kreis_map <- merge(kreis, kreis_summary, by = "kname")
kreis_sf  <- st_as_sf(kreis_map)


# ---- kreis colours ----
kreis_cols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
                "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
                "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")

# order kreise by number and give each one a colour
kreis_order <- kreis_summary$kname[order(as.numeric(gsub("\\D", "", kreis_summary$kname)))]
names(kreis_cols) <- kreis_order

kreis_sf$kname      <- factor(kreis_sf$kname,      levels = kreis_order)
kreis_summary$kname <- factor(kreis_summary$kname, levels = kreis_order)


# ---- plot theme: black outline, no grid lines ----
theme_kreis <- theme_cowplot() +
  theme(
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid      = element_blank(),
    legend.position = "none"
  )


# ---- map of kreise with lat / lon coordinates, no axis titles ----
p_map <- ggplot(kreis_sf) +
  geom_sf(aes(fill = kname), colour = "black", linewidth = 0.5) +
  geom_sf_text(aes(label = gsub("\\D", "", kname)), size = 4, fontface = "bold") +
  scale_fill_manual(values = kreis_cols) +
  coord_sf(crs = 4326, datum = 4326) +
  scale_x_continuous(breaks = seq(8.45, 8.65, by = 0.05),
                     labels = function(x) paste0(x, "°E")) +
  scale_y_continuous(breaks = seq(47.32, 47.44, by = 0.04),
                     labels = function(y) paste0(y, "°N")) +
  theme_kreis +
  theme(
    axis.title = element_blank(),
    axis.text  = element_text(size = 8)
  )


# ---- helper for the bar plots ----
bar_kreis <- function(var, ylab) {
  ggplot(kreis_summary, aes(x = kname, y = .data[[var]], fill = kname)) +
    geom_col(colour = "black", linewidth = 0.5) +
    scale_fill_manual(values = kreis_cols) +
    scale_x_discrete(labels = function(x) gsub("\\D", "", x)) +
    labs(x = "Kreis", y = ylab) +
    theme_kreis
}

p_trees   <- bar_kreis("n_trees",        "Number of trees")
p_density <- bar_kreis("trees_per_km2",  "Trees per km² (land)")
p_species <- bar_kreis("n_species",      "Number of species")
p_shannon <- bar_kreis("shannon",        "Shannon index")
p_simpson <- bar_kreis("simpson",        "Simpson index")
p_ndvi    <- bar_kreis("green_ndvi_pct", "Green area, NDVI (% of land)")
p_evi     <- bar_kreis("green_evi_pct",  "Green area, EVI (% of land)")
p_water   <- bar_kreis("water_pct",      "Water (% of area)")


# ---- combine bar plots into a 4 x 2 grid ----
p_bars <- plot_grid(p_trees, p_density, p_species, p_shannon,
                    p_simpson, p_ndvi, p_evi, p_water,
                    ncol = 4, labels = c("B", "C", "D", "E", "F", "G", "H", "I"))

# combine map and bar plots
p_final <- plot_grid(p_map, p_bars,
                     ncol = 2, rel_widths = c(1, 2.6), labels = c("A", ""))

p_final


# ---- save (before any base-graphics checks, so a plotting error can't block outputs) ----
ggsave("kreis_panels.png", p_final, width = 20, height = 7, dpi = 300)

write.csv(kreis_summary,  "kreis_summary.csv",           row.names = FALSE)
write.csv(counts_species, "species_counts_by_kreis.csv", row.names = FALSE)


# ---- optional check: masked rasters ----
# Saved to file so it never depends on the size of the RStudio Plots pane
png("check_masked_rasters.png", width = 2400, height = 1200, res = 200)
plot(c(ndvi, evi), main = c("NDVI (water masked)", "EVI (water masked)"))
dev.off()

# To view on screen instead, open a separate window (Windows):
# graphics.off()
# windows(width = 12, height = 6)
# plot(c(ndvi, evi), main = c("NDVI (water masked)", "EVI (water masked)"))