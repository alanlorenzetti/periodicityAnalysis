# alorenzetti 202008

# loading packs ########
if(!require(pacman)){
  install.packages("pacman")
}
library(pacman)

# # the following package is only available at GitHub
# install_github("LabTranslationalArchitectomics/riboWaltz",
#                dependencies = TRUE, 
#                build_opts = c("--no-resave-data", "--no-manual"))

# packs
packs = c("BiocManager",
          "tidyverse",
          "Rsamtools",
          "rtracklayer",
          "GenomicRanges",
          "ggthemes",
          "devtools",
          "riboWaltz",
          "data.table",
          "ggpubr",
          "BSgenome",
          "ggtext",
          "svglite")

p_load(char = packs)

# setting ggplot theme
theme_set(theme_bw())

# getting tab10 colors
tab10 = ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`
