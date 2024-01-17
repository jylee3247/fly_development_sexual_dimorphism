##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####


# Draw Figure 6D ----------------------------------------------------------

library(ComplexHeatmap)
library(tidyverse)
library(data.table)

# setwd("wd")

# Load sequence identity matrix from Clustal Omega
mat <- fread("sequence_identity_matrix.mat", skip = 1, sep = " ") %>% column_to_rownames("V1")

# Load sequence metadata file
md <- fread("sequence_information.txt")

# Set the color of each components
type_cols <- c(`Bacterial FALDH` = "black",
               `AO FALDH` = "red", 
               `Five Bacterial FALDH (function identified)` = "#0000FF", 
               `D. melanogaster ALDH` = "darkgreen")

# Heatmap annotation
al_preset = list(border = "black", title_position = "topcenter", 
                 title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 8))

column_ha <- HeatmapAnnotation(
  Type = md$type,
  col = list(Type = type_cols),
  annotation_legend_param = al_preset)

row_ha <- rowAnnotation(
  Type = md$type,
  col = list(Type = type_cols),
  annotation_legend_param = al_preset, show_legend = FALSE)

# Draw heatmap with the UPGMA clustering algorithm.
ht <- Heatmap(as.matrix(mat), bottom_annotation = column_ha, right_annotation = row_ha, name = "Identity(%)",
              show_row_names = FALSE, show_column_names = FALSE, 
              clustering_method_columns = "average",  clustering_method_rows = "average", 
              row_split = 2, column_split = 2, 
              column_names_gp = gpar(fontsize = 2), row_names_gp = gpar(fontsize = 2))
ht_draw <- draw(ht)

# Save heatmap as png file
png("Heatmap.png", width = 8, height = 6, units = "in", res = 300)
print(ht_draw)
decorate_column_dend("Identity(%)", {grid.rect(gp = gpar(fill = "#00FF0040", col = NA))}, slice = 1)
decorate_column_dend("Identity(%)", {grid.rect(gp = gpar(fill = "#FF000040", col = NA))}, slice = 2)
decorate_row_dend("Identity(%)", {grid.rect(gp = gpar(fill = "#00FF0040", col = NA))}, slice = 1)
decorate_row_dend("Identity(%)", {grid.rect(gp = gpar(fill = "#FF000040", col = NA))}, slice = 2)
dev.off()