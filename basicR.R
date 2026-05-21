# =========================================
# BASIC R IN PLANT SCIENCE
# NCSU
# Author: Roberto Fritsche-Neto
# email: roberto.neto@ncsu.edu
# Last update: May 21, 2026
# =========================================

################################################################################

# 1. INSTALLING AND LOADING PACKAGES

################################################################################

# FROM CRAN
install.packages("devtools")
install.packages("BiocManager")

# FROM github
library(devtools)
devtools::install_github('famuvie/breedR')
BiocManager::install("impute")

# Load packages ---------------------------------------------------------------
library(breedR)
library(tidyverse)
library(tidyplots)
library(janitor)
library(skimr)
library(patchwork)
library(GGally)
library(desplot)
library(lme4)
library(car)
library(bestNormalize)
library(ScottKnott)
library(drc)
library(jtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggstatsplot)
library(patchwork)

# you can also use 
require(ggplot2) # silent way to call the package
# or call just a specific function using ::
#Usefull when two different packages have the same functions names
#snpReady::G.matrix()

################################################################################

# 2. IMPORTING AND ORGANIZING DATA

################################################################################

# Import TXT file -------------------------------------------------------------

# header = TRUE means the first row contains column names.

# na.strings defines missing values.

data <- read.table(
  "data.txt",
  header = TRUE,
  na.strings = c("NA", ".")
)

# Alternative formats ---------------------------------------------------------

# data <- read.csv("data.csv")

# data <- readRDS("data.rds")

# Basic inspection ------------------------------------------------------------

head(data)      # first rows
head(data, 10)  # first 10 rows

tail(data)      # last rows

dim(data)       # number of rows and columns

names(data)     # variable names

str(data)       # variable classes

# Clean variable names --------------------------------------------------------

# Makes names easier to use in analyses.

data <- clean_names(data)

names(data)

# Convert variables to factors ------------------------------------------------

# Experimental design variables are usually factors.

# adjusting to have only factors and numeric
data$plot <- as.factor(data$plot)
str(data)
# now, doing the same for many cols at the same time
data[ ,1:9] <- lapply(data[ ,1:9], as.factor)

str(data)

################################################################################

# 3. DESCRIPTIVE ANALYSIS

################################################################################

# General summary statistics --------------------------------------------------

summary(data)

# More detailed summary -------------------------------------------------------

skim(data)

# Mean of a trait -------------------------------------------------------------

mean(data$sdm, na.rm = TRUE)

# Standard deviation ----------------------------------------------------------

sd(data$sdm, na.rm = TRUE)

# Coefficient of variation ----------------------------------------------------

cv_sdm <- sd(data$sdm, na.rm = TRUE) /
  mean(data$sdm, na.rm = TRUE) * 100

cv_sdm

################################################################################

# 4. CORRELATIONS AMONG TRAITS

################################################################################

# Select only numeric variables ----------------------------------------------

num_data <- data[ ,10:12]
head(num_data)

# more efficient
str(data)
sapply(data, is.numeric)
num_data <- data[ ,sapply(data, is.numeric)]
head(num_data)
  
# Correlation matrix ----------------------------------------------------------

cor_matrix <- cor(num_data, use = "pairwise.complete.obs")

round(cor_matrix, 2)

# Correlation heatmap using tidyplots ----------------------------------------

cor_df <- as.data.frame(as.table(cor_matrix))
head(cor_df)
colnames(cor_df)[3] <- "Correlation"
  
cor_plot <- cor_df %>%
  tidyplot(x = Var1, y = Var2, color = Correlation) %>%
  add_heatmap() %>%
  add_title("Trait Correlation Among Traits")

cor_plot

################################################################################

# 5. CHECKING THE FIELD DESIGN

################################################################################

# Visualizing spatial distribution in the field -------------------------------

# Useful to detect gradients and experimental problems.

d1 <- desplot(
  data,
  sdm ~ row * col,
  out1 = rep,
  out2 = n,
  out2.gpar = list(col = "darkgreen", lwd = 1)
)

print(d1)

################################################################################

# 6. QUALITY CONTROL

################################################################################

# 6.1 Outlier detection -------------------------------------------------------

# Simple visualization

data %>%
  tidyplot(x = type, y = sdm) %>%
  add_boxplot()

# Linear model for residual inspection ---------------------------------------

fit <- lm(sdm ~ type + row + col + n + gid, data = data)

summary(fit)

# Statistical outlier detection ----------------------------------------------

outlier_test <- outlierTest(fit)
outlier_test

# Removing outliers -----------------------------------------------------------

# Here we replace outlier observations with NA.

outlier_names <- names(outlier_test$rstudent)

if(length(outlier_names) > 0){
  data[outlier_names, "sdm"] <- NA
}

################################################################################

# 6.2 NORMALITY ASSESSMENT

################################################################################

# Histograms help identify asymmetry.

data %>%
  tidyplot(x = sdm) %>%
  add_histogram(bins = 20)

# Shapiro-Wilk test -----------------------------------------------------------

# H0 = data are normally distributed.

shapiro.test(na.omit(data$sdm))

# Transformation --------------------------------------------------------------

# bestNormalize automatically searches for better transformations.

sdm_adj <- bestNormalize(
  data$sdm,
  standardize = FALSE,
  allow_orderNorm = TRUE,
  out_of_sample = FALSE
)

# Best transformation selected ------------------------------------------------

sdm_adj$chosen_transform

# Store transformed data ------------------------------------------------------

data$sdm_adj <- sdm_adj$x.t

# Test transformed data -------------------------------------------------------

shapiro.test(na.omit(data$sdm_adj))

################################################################################

# 6.3 RESIDUAL DIAGNOSTICS

################################################################################

# Compare original and transformed models ------------------------------------

fit1 <- lm(sdm ~ type + row + col + n + gid, data = data)
fit2 <- lm(sdm_adj ~ type + row + col + n + gid, data = data)

# QQ-plots -------------------------------------------------------------------

par(mfrow = c(1,2))

qqnorm(resid(fit1))
qqline(resid(fit1), col = "red")

qqnorm(resid(fit2))
qqline(resid(fit2), col = "blue")

par(mfrow = c(1,1))

################################################################################

# 7. EXPORTING CLEANED DATA

################################################################################

# Save processed datasets -----------------------------------------------------

write.table(data, "data_clean.txt", row.names = FALSE)
write.csv(data, "data_clean.csv", row.names = FALSE)

################################################################################

# 8. GRAPHICS WITH TIDYPLOTS

################################################################################

# Histogram ------------------------------------------------------------------

p1 <- data %>%
  tidyplot(x = sdm) %>%
  add_histogram(bins = 20) %>%
  add_title("Distribution of SDM")

p1

# Boxplot --------------------------------------------------------------------

p2 <- data %>%
  tidyplot(x = type, y = sdm, color = type) %>%
  add_boxplot() %>%
  add_title("SDM across genotypic types")

p2

# Scatter plot with regression line ------------------------------------------

p3 <- data |>
  tidyplot(x = n, y = sdm, color = gid) |>
  add_mean_line() |>
  add_mean_dot() |>
  add_sem_ribbon()

p3

# Bar plot -------------------------------------------------------------------

p4 <- data %>% 
  tidyplot(x = type, y = sdm, color = type) %>% 
  add_mean_bar(alpha = 0.4) %>% 
  add_sem_errorbar() %>% 
  add_data_points_beeswarm()

p4

# barplot with a t-test
data |> 
  tidyplot(x = n, y = sdm, color = n) |> 
  add_boxplot() |> 
  add_test_pvalue(ref.group = 1)


# Combine figures ------------------------------------------------------------

final_plot <- p1 + p2 + p3 + p4

final_plot

# Save figure ----------------------------------------------------------------

ggsave(
  filename = "Figure_1.pdf",
  plot = final_plot,
  width = 12,
  height = 10
)


# NON-PARAMETRIC TESTS AND GRAPHS
# ----------------------------
# Contrast ideal and low N or type 
kruskal.test(sdm ~ n, data = data)
kruskal.test(sdm ~ type, data = data)

# Create boxplot with Kruskal-Wallis test result shown
ggstatsplot::ggbetweenstats(
  data = data,
  x = n,
  y = sdm,
  type = "nonparametric",  # Automatically chooses Kruskal-Wallis for >2 groups
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",  # show only significant comparisons
  p.adjust.method = "BH",            # adjust p-values
  messages = FALSE
)


################################################################################

# 9. ANALYSIS OF VARIANCE (ANOVA)

################################################################################

# ANOVA evaluates whether treatments differ significantly.

# Simple model ---------------------------------------------------------------

anova_1 <- aov(sdm ~ gid, data = data)

summary(anova_1)

plot(anova_1)

# More realistic model -------------------------------------------------------

# Includes replication and interaction effects.

anova_2 <- aov(sdm ~ rep + gid + n + gid:n, data = data)

summary(anova_2)

plot(anova_2)

################################################################################

# 9.1 COEFFICIENT OF VARIATION (CV%)

################################################################################

# CV measures experimental precision.

ve <- summary(anova_2)[[1]][5,3]

cv <- sqrt(ve) /
  mean(na.omit(data$sdm)) * 100

cv

################################################################################

# 9.2 BROAD-SENSE HERITABILITY

################################################################################

# Heritability estimates the proportion of phenotypic variation

# explained by genetic effects.

n_reps <- length(unique(data$rep))

vg <- (summary(anova_2)[[1]][2,3] - ve) / n_reps

h2g <- vg / (vg + ve / n_reps)

h2g

################################################################################

# 9.3 TUKEY TEST

################################################################################

# Multiple comparison among treatment means.

# Nitrogen levels ------------------------------------------------------------

TukeyHSD(anova_2, "n")

# Genotypes ------------------------------------------------------------------

TukeyHSD(anova_2, "gid")

################################################################################

# 9.4 SCOTT-KNOTT GROUPING

################################################################################

# Widely used in plant breeding.

sk_gid <- SK(anova_2, which = "gid")
summary(sk_gid)

sk_n <- SK(anova_2, which = "n")
summary(sk_n)

################################################################################

# 10. REGRESSION ANALYSIS

################################################################################

# Regression models evaluate relationships between variables.

################################################################################

# 10.1 SIMPLE LINEAR REGRESSION

################################################################################

lm_model <- lm(sdm ~ sra, data = data)

summary(lm_model)

# R-squared ------------------------------------------------------------------

summary(lm_model)$r.squared

# Scatter plot with regression line ------------------------------------------

reg_plot <- data %>%
  ggplot(aes(x = sra, y = sdm)) +
  
  geom_point(alpha = 0.7, size = 2) +
  
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  
  stat_regline_equation(
    label.x.npc = "left",
    label.y.npc = 0.95
  ) +
  
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = 0.88
  ) +
  
  labs(
    title = "Relationship between SRA and SDM",
    x = "SRA",
    y = "SDM"
  ) +
  
  theme_classic(base_size = 14)

reg_plot

################################################################################

# 10.2 MULTIPLE REGRESSION

################################################################################

multi_model <- lm(sdm ~ sra + nae, data = data)

summary(multi_model)

# Added-variable plots -------------------------------------------------------

avPlots(multi_model)

################################################################################

# 10.3 QUADRATIC REGRESSION

################################################################################

# Useful when responses are nonlinear.

quad_model <- lm(sdm ~ sra + I(sra^2), data = data)

summary(quad_model)

# Generate prediction values -------------------------------------------------

sra_values <- seq(
  min(data$sra, na.rm = TRUE),
  max(data$sra, na.rm = TRUE),
  length.out = 200
)

pred_data <- data.frame(sra = sra_values)

pred_data$pred <- predict(quad_model, newdata = pred_data)

# Plot quadratic curve -------------------------------------------------------

quad_plot <- ggplot(data, aes(x = sra, y = sdm)) +
  geom_point(alpha = 0.7, size = 2) +
  
  geom_line(
    data = pred_data,
    aes(x = sra, y = pred),
    color = "blue",
    linewidth = 1
  ) +
  
  labs(
    title = "Observed vs Predicted Relationship",
    x = "SRA",
    y = "SDM"
  ) +
  
  theme_minimal(base_size = 14)

quad_plot


################################################################################

# 10.4 LOGISTIC REGRESSION

################################################################################

# Logistic regression is used for binary responses.

# Create binary variable -----------------------------------------------------

data$n_bi <- ifelse(data$n == "ideal", 1, 0)

# Fit model ------------------------------------------------------------------

mylogit <- glm(
  n_bi ~ sdm + sra + nae,
  family = "binomial",
  data = data
)

summary(mylogit)

anova(mylogit)

summ(mylogit)

################################################################################

# 10.5 SIGMOIDAL MODEL

################################################################################

# Sigmoidal models are common in biology.

sig_model <- drm(sdm ~ nae, data = data, fct = G.3())

summary(sig_model)

plot(sig_model)

################################################################################

# 11. MIXED MODELS

################################################################################

# Mixed models allow fixed and random effects simultaneously.

# Very important in plant breeding.

# Create interaction term ----------------------------------------------------

data$gxn <- paste0(data$gid, "_", data$n)

# Mixed model ----------------------------------------------------------------

fit_mixed <- remlf90(
  fixed = sdm ~ n,
  random = ~ gid + rep + gxn,
  data = data
)

# Variance components --------------------------------------------------------

fit_mixed$var

################################################################################

# 11.1 HERITABILITY FROM MIXED MODELS

################################################################################

n_reps <- length(unique(data$rep))
n_levels <- length(unique(data$n))

h2_mixed <- fit_mixed$var[1,1] /
  (
    fit_mixed$var[1,1] +
      fit_mixed$var[3,1] / n_levels +
      fit_mixed$var[4,1] / (n_reps * n_levels)
  )

h2_mixed

################################################################################

# 11.2 BLUPs

################################################################################

# Best Linear Unbiased Predictions.

blups <- fit_mixed$ranef$gid[[1]]

head(blups)

################################################################################

# 11.3 BLUP VISUALIZATION WITH TIDYPLOTS

################################################################################
# reorganizing the data
blup_df <- data.frame(
  gid = rownames(blups),
  value = blups$value,
  se = blups$s.e.
)

# estimating the confidence interval 
blup_df$dms <- blup_df$se * 1.96

# creating a plot
blup_plot <- blup_df %>%
  mutate(gid = reorder(gid, value)) %>%
  ggplot(aes(x = gid, y = value)) +
  
  # confidence interval (thicker, light color)
  geom_errorbar(
    aes(ymin = value - dms,
        ymax = value + dms),
    width = 0.15,
    linewidth = 0.6,
    alpha = 0.7
  ) +
  
  # BLUP points (more visible, centered)
  geom_point(
    color = "firebrick",
    size = 2.2
  ) +
  
  coord_flip() +
  
  labs(
    x = "Genotype (GID)",
    y = "BLUP",
    title = "BLUP estimates with uncertainty intervals"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold")
  )


blup_plot

################################################################################

# 11.4 MODEL COMPARISON

################################################################################

# Reduced model --------------------------------------------------------------

fit_reduced <- remlf90(
  fixed = sdm ~ n,
  random = ~ gid + rep,
  data = data
)

# AIC comparison -------------------------------------------------------------

fit_mixed$fit$AIC
fit_reduced$fit$AIC

# Likelihood ratio test ------------------------------------------------------

lrt <- abs(
  fit_mixed$fit$`-2logL` -
    fit_reduced$fit$`-2logL`
)

pchisq(lrt, 1, lower.tail = FALSE)

################################################################################

# 12. CREATING FUNCTIONS

################################################################################

# Custom functions simplify repetitive tasks.

# Mean function --------------------------------------------------------------

my_mean <- function(x){
  
  result <- sum(x, na.rm = TRUE) /
    length(na.omit(x))
  
  return(result)
}

apply(data[,10:12], 2, my_mean)

# Positive-value function ----------------------------------------------------

be_positive <- function(x){
  
  if(x >= 0){
    return(x)
  } else {
    return(0)
  }
}

# Example --------------------------------------------------------------------

aux <- rnorm(10)

sapply(aux, be_positive)



### CHALLENGE TIME
# make a function to identify if the median of a vector x
x <- rnorm(10)
mean(x)

# median?


################################################################################

# 13. LOOPS

################################################################################

# Loops automate repetitive analyses.

traits <- c("sdm", "sra", "nae")

# Empty list to store outputs ------------------------------------------------

output <- list()

# Loop across traits ---------------------------------------------------------

for(i in traits){
  
  formula_i <- as.formula(
    paste(i, "~ gid + rep")
  )
  
  model_i <- aov(formula_i, data = data)
  
  output[[i]] <- summary(model_i)
}

# Display results ------------------------------------------------------------

output

################################################################################

# 13.1 AUTOMATED HISTOGRAMS

################################################################################

# Generate one histogram per trait.

plot_list <- list()

for(i in traits){
  
  p <- data %>%
    tidyplot(x = .data[[i]]) %>%
    add_histogram(bins = 20) %>%
    add_title(paste("Histogram of", i))
  
  plot_list[[i]] <- p
}

# wrap all plots 
final_plot <- wrap_plots(plot_list, nrow = 1)
final_plot


### CHALLENGE TIME
# now, combine these loops to deploy all results together

################################################################################

# 14. BEST PROGRAMMING PRACTICES

################################################################################

# 1. Use meaningful variable names.

# 2. Avoid spaces in object names.

# 3. Prefer snake_case naming.

# 4. Save intermediate results.

# 5. Comment all important analyses.

# 6. Keep raw data unchanged.

# 7. Build reproducible workflows.

# 8. Always inspect your data before analysis.

# 9. Graph your data before modeling.

# 10. Interpret biological meaning, not only p-values.

################################################################################

# 15. USEFUL RESOURCES

################################################################################

# Data visualization

# [https://grafify.shenoylab.com/](https://grafify.shenoylab.com/)

# [https://rkabacoff.github.io/datavis/](https://rkabacoff.github.io/datavis/)

# [http://www.cookbook-r.com/Graphs/](http://www.cookbook-r.com/Graphs/)

# R search engine

# [https://rseek.org/](https://rseek.org/)

# Tidyverse documentation

# [https://www.tidyverse.org/](https://www.tidyverse.org/)

# Plant breeding statistics

# [https://cran.r-project.org/](https://cran.r-project.org/)

################################################################################

# END OF SCRIPT

################################################################################