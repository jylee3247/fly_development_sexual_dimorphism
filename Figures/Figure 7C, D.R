##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####


# Draw Figure 7C, D ----------------------------------------------------------

library(tidyverse)

# setwd("wd")

# Load Maaslin2 result.
# Use result from taxonomic analysis or functional (Pfam) analysis
maaslin2.res <- read.table("maaslin2_result.tsv", sep="\t") 

# Set colors
color <- c(Male = "blue", Female = "red")

# Use results with FDR corrected p-value (qval) below 0.1.
maaslin2.res <- maaslin2.res %>% filter(qval < 0.1)

# Coefficient >0 indicate enrichment in females, and coefficient <0 indicate enrichment in males.
maaslin2.res$Sex = ifelse(coef>0, "Female", "Male")

# Draw ggplot
p <- ggplot(maaslin2.res, aes(x = fct_reorder(feature, coef), 
                                 y = coef, fill = Sex)) + 
  geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), 
                width = 0.3, color = "black", lwd = 0.5) + 
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_hline(yintercept = 0, lwd = 1) +
  labs(x = NULL, y = "Coefficient") + 
  theme_bw() +
  scale_fill_manual(name = "Sex", values = color) +
  coord_flip() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length = unit(2, "mm"),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(size = 1.5)
  ) 

ggsave("Maaslin2_plot.png", p, width = 15, height = 10, dpi = 300)

# Resulting plot is further refined using other image editing programs.