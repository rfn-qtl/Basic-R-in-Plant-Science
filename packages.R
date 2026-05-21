# Before the course, we need to install all the packages 

install.packages(c(
  "tidyverse",
  "tidyplots",
  "janitor",
  "skimr",
  "patchwork",
  "GGally",
  "desplot",
  "lme4",
  "car",
  "bestNormalize",
  "ScottKnott",
  "drc",
  "jtools",
  "reshape2",
  "ggplot2",
  "ggpubr",
  "dplyr",
  "ggstatsplot",
  "patchwork"
))

# Mixed model package
install.packages("devtools")
library(devtools)
devtools::install_github("famuvie/breedR")
