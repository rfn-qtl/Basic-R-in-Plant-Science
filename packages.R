# let's install all the packages we will need for this training

cran_packages <- c(
  "devtools", "BiocManager", "lme4", "car", "bestNormalize",
  "ggplot2", "ggstatsplot", "ScottKnott", "jtools", "drc", "reshape2"
)

install.packages(cran_packages)

devtools::install_github('famuvie/breedR')
BiocManager::install("impute")