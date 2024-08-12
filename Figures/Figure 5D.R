##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####


# Draw Figure 5D ----------------------------------------------------------

library(treeio)
library(ggtreeExtra)
library(ggtree)
library(ggpubr)
library(cowplot)
library(ape)
library(tidyverse)

# setwd("wd")

# Load phylogenetic tree of ALDH sequences

tree <- read.tree("tree.nwk")

# Load sequence metadata file
md <- fread("sequence_information.txt")

# Set the color of each components
type_cols <- c(`Bacterial FALDH` = "black",
               `AO FALDH` = "red", 
               `Five Bacterial FALDH (function identified)` = "#0000FF", 
               `D. melanogaster ALDH` = "darkgreen")

# Combine tree and metadata
tree <- treeio::full_join(as.treedata(tree), md, by = join_by(label == id))

# Draw phylogenetic tree
p <- ggtree(tree, layout = "circular") +
  geom_tippoint(aes(color = type), size = 1.5, show.legend = TRUE) +
  scale_color_manual(name = "Type", values = type_cols) +
  theme(
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", face = "bold")
  )

p.lgd <- as_ggplot(get_legend(p))
p <- p + theme(legend.position = "none")

ggsave("phylogenetic_tree.png", p,  width = 6, height = 15, dpi = 300)
ggsave("phylogenetic_tree_legend.png", p.lgd,  width = 5, height = 5, dpi = 300)
