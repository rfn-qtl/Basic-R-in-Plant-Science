# =========================================
# BASIC R IN PLANT SCIENCE
# NCSU
# Author: Roberto Fritsche-Neto
# email: roberto.neto@ncsu.edu
# Last update: June 26, 2025
# =========================================

################################### 1st part ###################################

# 1. INSTALL AND LOAD PACKAGES 
# --------------------
# FROM CRAN
install.packages("devtools")
install.packages("BiocManager")

# FROM github
library(devtools)
devtools::install_github('famuvie/breedR')
BiocManager::install("impute")

# loading them
library(breedR)

# 2. DATA IMPORTATION 
# --------------------
# Load dataset from txt
data <- read.table("data.txt", header = TRUE, na.strings = NA)
# the missing data usually is NA, but can be any format, more than one
# data <- read.table("data.txt", header = TRUE, na.strings = c(NA, "."))
# it will no change your original file

# # Load dataset from CSV
# data <- read.csv("data.csv", header = TRUE, sep = ",")

# # Load dataset from R file
# data <- readRDS("data")

head(data)  # view first few rows
tail(data) # view the last few rows
dim(data) # dimension of the dataset

# types of data and preparing your file
str(data)   # structure of the dataset

# adjusting to have only factors and numeric
data$plot <- as.factor(data$plot)
str(data)
# now, doing the same for many cols at the same time
data[ ,1:9] <- lapply(data[ ,1:9], as.factor)
str(data)

# 2. DESCRIPTIVE ANALYSIS
# -----------------------
summary(data)      # general summary
mean(data$SDM)    # mean of a variable
sd(data$SDM)      # standard deviation
boxplot(data$SDM) # boxplot
cor(data$SDM, data$SRA) # correlation between two traits
round(cor(data[ ,10:12], use = "pairwise"),2) # correlation among many traits
heatmap(round(cor(data[ ,10:12], use = "pairwise"),2))

# doing the same for many traits at the same time
apply(data[ ,10:12], 2, mean) # 2 apply in col, 1 in rows

# check the experimental design and spatial distribution
library(desplot)
d1 <- desplot(data, SDM ~ row*col, 
              out1 = rep, 
              out2 = N,
              out2.gpar=list(col = "green", lwd = 1, lty = 1))
print(d1)

# 3. QUALITY CONTROL
# -----------------------
# 3.1 identifying outliers
boxplot(data$SDM, col = "red")
#install.packages("lme4")
library(lme4)
# outlier detection and elimination
fit <- lm(SDM ~ type + row + col + N + gid, data = data)
#install.packages("car")
library(car)
(outlier <- names(outlierTest(fit)$p))
data[outlier, "SDM"] <- NA


# 3.2 testing for normality
# First lets check using patterns
length(data$SDM)
rnorm(220)
plot(rnorm(220))
plot(density(rnorm(220)))
shapiro.test(rnorm(220)) # normal distribution
shapiro.test(runif(220)) # uniform distribution
# then, 
shapiro.test(data$SDM)

#install.packages("bestNormalize")
require(bestNormalize)
SDMadj <- bestNormalize(data$SDM, standardize = FALSE, allow_orderNorm = TRUE, out_of_sample = FALSE)
SDMadj$chosen_transform
shapiro.test(SDMadj$x.t)
data$SDMadj <- SDMadj$x.t
head(data)

# What about the residuals?
# Quartile‐Quartile (Q‐Q) normality plot for residuals
fit <- lm(SDM ~ type + row + col + N + gid, data = data)
fit2 <- lm(SDMadj ~ type + row + col + N + gid, data = data)

par(mfrow = c(2,2)) # organize the plot window in 1 row and 2 col
qqnorm(resid(fit))
qqline(resid(fit), col = "red")
qqnorm(resid(fit2))
qqline(resid(fit2), col = "blue")
hist(data$SDM, col = "red", main = "SDM", xlab = "SDM")
hist(data$SDMadj, col = "blue", main = "Adjusted SDM", xlab = "Adjusted SDM")
dev.off()

# 3.3 saving the newest data files
str(data)
head(data)
write.table(data, "data_clean.txt")
write.csv(data, "data_clean.csv")

# 4. IMPORTANT TYPES OF GRAPHS
# ----------------------------
# install.packages("ggplot2")
library(ggplot2)

colnames(data)

# Histogram
p1 <- ggplot(data, aes(x = SDM)) + 
  geom_histogram(bins = 20, fill = "skyblue") +
  theme_minimal()
p1

# Boxplot
p2 <- ggplot(data, aes(x = type, y = SDM, fill = type)) + 
  geom_boxplot() +
  theme_minimal()
p2

# Scatter Plot
p3 <- ggplot(data, aes(x = SRA, y = SDM)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
p3 
  
# Bar Plot
p4 <- ggplot(data, aes(x = type, fill = type)) + 
  geom_bar() +
  theme_minimal()
p4

# Line Plot
ggplot(data, aes(y = SDM, x = N, group = gid, color = gid)) + 
  geom_line() +
  theme_minimal()

# combining and saving graphs
plot_list <- list(p1, p2, p3, p4)

p_final <- ggstatsplot::combine_plots(
  plotlist = plot_list,
  plotgrid.args = list(nrow = 2))

p_final

ggsave(filename = './Fig1.pdf',
       plot = p_final,
       device = 'pdf',
       width = 300,
       height = 300,
       units = 'mm',
       dpi = 300)

################################### 2nd part ###################################

# 5. ANOVA (Analysis of Variance) and comparison tests
# -----------------------------------------------------
# 5.1. ANOVA
# Useful to understand the significance and importance of factors in your response variable

# first model
anova_1 <- aov(SDM ~ gid, data = data)
summary(anova_1)
plot(anova_1)

# a more realistic model
anova_2 <- aov(SDM ~ rep + gid + N + gid*N, data = data)
summary(anova_2)
plot(anova_2)

# estimate CV% and heritability
summary(anova_2)

# CV% - it gives us an idea the precision of the whole experiment
# CV = sqrt(Ve) / mean *100
Ve <- summary(anova_2)[[1]][5,3]
CV <- sqrt(Ve) / mean(na.omit(data$SDM)) *100
CV

# heritability (fixed model)
# h2g = Vg / (Vg + Ve/rep) = Vg / Vp
#it gives us an idea the accuracy of the whole experiment
# the proportion explained by genonoities, 
# the correlation between phenotypes and genotypes
n.reps <- length(unique(data$rep))
n.reps
Ve <- summary(anova_2)[[1]][5,3]
Vg <- (summary(anova_2)[[1]][2,3] - Ve)/n.reps
h2g <- Vg / (Vg + Ve / n.reps) 
h2g

# 5.2. TUKEY HSD test for factor in the model
# N levels
TukeyHSD(anova_2, "N", ordered = TRUE, conf.level=.95)
plot(TukeyHSD(anova_2, "N", ordered = TRUE,conf.level=.95))
# genotypes
TukeyHSD(anova_2, "gid", ordered = TRUE, conf.level=.95)
plot(TukeyHSD(anova_2, "gid", ordered = TRUE,conf.level=.95))

# 5.3 SCOTT-KNOTT test
# --------------------
#install.packages("ScottKnott")
library(ScottKnott)
# Genotypes
sk1 <- SK(x = anova_2, which = "gid")
summary(sk1)
# N levels
sk2 <- SK(x = anova_2, which = "N")
summary(sk2)


# 6. REGRESSIONS
# --------------------
# useful to understand the significance and importance of numeric variable to explain your response variable
colnames(data)

# 6.1 Linear = y = a + bX + e
lm_model <- lm(SDM ~ SRA, data = data)
summary(lm_model)
model_summary <- summary(lm_model)

# Access R-squared
model_summary$r.squared

# Access coefficients table
model_summary$coefficients

# creating a simple graph
plot(data$SRA, data$SDM)
abline(lm_model, col = "blue")


# 6.2 Multiple = y = a + b1X1 + b2X2 + ... + e
# -------------------------
m_model <- lm(SDM ~ SRA + NAE, data = data)
summary(m_model)
model_summary <- summary(m_model)
# Access R-squared
model_summary$r.squared
# Access coefficients table
model_summary$coefficients
# creating a simple graph
library(car)
#produce added variable plots
avPlots(m_model)

# TIPS 
# using all variables without typing them
#m_model <- lm(SDM ~ ., data = data)

# 6.3 Quadratic model = y = a + bX + b^2X + e
data$SRA2 <- data$SRA^2
qd_model <- lm(SDM ~ SRA + SRA2, data = data)
summary(qd_model)

# Plotting
#create sequence of possible SRA 
range(data$SRA)
SRA_Values <- seq(min(data$SRA), max(data$SRA), 0.005)
#create list of predicted values levels using quadratic model
SDM_Predict <- predict(qd_model, list(SRA = SRA_Values, SRA2=SRA_Values^2))
#create scatterplot of original data values
plot(data$SRA, data$SDM, pch=16)
#add predicted lines based on quadratic regression model
lines(SRA_Values, SDM_Predict, col='blue')

# 6.4 Logistic - the response variable is TRUE or FALSE (binomial)
# let's create one in our dataset for levels of N
data$N
data$N == "ideal"
as.numeric(data$N == "ideal")
data$N_bi <- as.numeric(data$N == "ideal")
head(data)

mylogit <- glm(N_bi ~ SDM + SRA + NAE , family = "binomial", data = data)
summary(mylogit)
anova(mylogit)
mylogit$coefficients
library(jtools)
summ(mylogit)

# 6.4 Sigmoidal = y ~ a/(1 + exp(-b * (x-c)) ) + d
# -----------------------
# tends to represent the "biology" of things
library(drc)
fm <- drm(SDM ~ NAE, data = data, fct = G.3())
plot(fm)
summary(fm)

################################### 3rd part ###################################

# 7. MIXED MODELS
# ----------------------------
# modeling better, considering random and fixed effects
# when it is important?

library(breedR)
# for this package we need to "create" a col for the interaction or nested effects 
data$GxN <- paste0(data$gid, data$N)
head(data)

# using only the classical experimental design
fit1 <- remlf90(fixed = SDM ~ N,
             random = ~gid + rep + GxN,
             data = data)

# components of variance
fit1$var
n.reps <- length(unique(data$rep))
n.levels <- length(unique(data$N))

#CV%
sqrt(fit1$var[4,1]) / mean(na.omit(data$SDM)) *100

# Broad-sense heritability
h2g <- fit1$var[1,1] / (fit1$var[1,1] + fit1$var[3,1]/n.levels + fit1$var[4,1]/(n.reps*n.levels))
h2g

# Best Linear Unbiased predictions - BLUPS
# for genotypes
BLUPs <- fit1$ranef$gid[[1]]
BLUPs

# for the interaction G x N
GxN <- fit1$ranef$GxN[[1]]
GxN

# reliability (r2) = 1 - PEV / (Vg + Vg*Fii)
# near to heritability
(r <- mean(1 - BLUPs$s.e.^2 / fit1$var[1,1]))

# confidence interval for BLUPS
DMS <- BLUPs$s.e.*1.96
blups2 <- data.frame(gid = rownames(BLUPs), BLUPs, "DMS" = DMS)
head(blups2)

library(ggplot2)
limits <- aes(ymax = blups2$value + blups2$DMS,
              ymin = blups2$value - blups2$DMS)
p <- ggplot(data = blups2, aes(x = reorder(factor(gid), -value), y = value))
p + geom_jitter(stat = "identity", colour = "red") +
  geom_errorbar(limits, position = position_dodge(0.5),
                width = 0.10) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "gid", y = "BLUP")

# Comparing models - LRT and AIC
# to test model factor or compare model, we need to run other models eliminating / including factors 

# So, let's remove the GxN factor and test if it is significant
fit2 <- remlf90(fixed = SDM ~ N,
                random = ~gid + rep,
                data = data)

# how the model is good (AIC) - the smaller, the better the model
fit1$fit$AIC
fit2$fit$AIC

# LRT test,
LRT <- abs(fit1$fit$`-2logL` - fit2$fit$`-2logL`)
LRT
pchisq(abs(LRT), 1, lower.tail = F) # chisq test


# 8. CREATING MY FUNCTIONS
# ----------------------------

# first, an example of how to estimate the mean
my_mean <- function(x){
  med <- sum(x) / length(x)
  return(med)
}

apply(data[,10:12],2, my_mean)

# let's create a function to identify if a number is positive, 
# if not, replace it with 0

be_positive <- function(x){
  if(x >= 0){ return(x)} else(0)
}

aux <- as.matrix(rnorm(10))
aux
apply(aux, 1, be_positive)

### CHALLENGE TIME
# make a function to identify if a number is prime
x <- 1:10

prime <- function(x){

  
  
  
  }

# let's test it
apply(matrix(x), 1, prime)
x[apply(matrix(x), 1, prime)]

################################### 4th part ###################################

# 9. USING LOOPS
# ----------------------------
# first, a basic example
x <- rnorm(10)
x

for (i in 1:length(x)){
  
  print(be_positive(x[i]))
  
}

# now, we need to reorganize the file
colnames(data)
traits <- colnames(data)[10:12]
traits

library(reshape2)
# reorganize the data
data.melted <- reshape2::melt(data, measure.vars = traits)
head(data.melted)

# create an grid
grid <- expand.grid(traits) 
grid
# more than two can be added to create a grid, for instance:
grid2 <- expand.grid(trait = traits, N = unique(data$N))
grid2

# create a object to storage all results
output <- list()

# doing my first loop from scratch - ANOVA for all traits within N levels
head(data.melted)
i = 1

### CHALLENGE TIME
# now, adjust this loop to add a histogram for each one of the traits and N levels

output <- list()
output_figures <- list()




# 10. TIPS AND RESOURCES
# ----------------------------

# 9.1 tips for programming
# avoid use the same name of already present functions
# avoid mix capital and small letters
# try to define a pattern for complex names data.field or data_field
# always add the last update in the script
# always have a copy of the main script
# take care to overwrite files, objects, and/or clean the environment
# when you clean your environment, you are not changing the script

# 9.2 graphs
# https://grafify.shenoylab.com/#intro
# https://rkabacoff.github.io/datavis/
# http://www.cookbook-r.com/Graphs/

# 9.3 R google
# https://rseek.org/

# 9.4 using chatGPT
# https://chatgpt.com/
# make a function if a number if cousin
# improve a graph

# ================================
# THE END 
# ================================