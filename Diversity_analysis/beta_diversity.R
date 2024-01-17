##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####
# Related to Figure 7A and B

# Calcaulte beta diversity ----------------------------------------------------------

library(phyloseq)
library(vegan)
library(ape)
library(tidyverse)

# setwd("wd")

###### Functions #####

rarefy.table <- function(input.ps, min.depth = NULL, seed.num = 123){
  if ( is.null(min.depth) ) {
    min.depth <- floor(min(sample_sums(input.ps)))
    cat("Rarefy depth is not specified. Use minimum sample depth", min.depth, "as min.depth.\n")
  }
  
  res.ps <-  rarefy_even_depth(input.ps, min.depth, replace = FALSE, rngseed = seed.num)
  return(res.ps)
}


calc.distance <- function(t, d = "bray"){ 
  return(vegdist(data.frame(t(as.matrix(t))), method = d)) 
}


calc.mds <- function(d, nk = NULL){
  if (is.null(nk)){
    warning("Number of 'k' is not defined. Use [sample number]-1 as 'k' value.")
    nk <- nrow(as.matrix(d))-1
  }
  return(stats::cmdscale(d, k = nk, eig = TRUE))
}

##### Diversity analysis #####

# Load data
table <- otu_table(read_delim("table.tsv"), taxa_are_rows = T)
tax.df <- tax_table(read_delim("taxonomy.tsv"))
tree <- phy_tree(ape::read.tree("tree.nwk"))
md <- sample_data(read_delim("metadata.tsv"))
ps <- phyloseq(table, tree, md, tax.df)

# Rarefying table
ps <- rarefy.table(ps)

# Get table object
t <- data.frame(ps@otu_table)

# Calcaulte Bray-Curtis dissimilarity matrix
d <- calc.distance(t, "bray")

# Calcaulte PCoA
mds <- calc.mds(d)


##### Visualization of the PCoA was conducted using GraphPad Prism v8.0 #####