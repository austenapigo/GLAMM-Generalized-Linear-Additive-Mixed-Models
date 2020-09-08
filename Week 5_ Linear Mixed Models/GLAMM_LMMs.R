rm(list = ls()) # clears workspace
# loading packages
Packages = c("ggplot2", "ggthemes", "tidyverse", "car", "lattice", "MASS", "lmerTest", "emmeans", "MuMIn")
lapply(Packages, library, character.only = TRUE)

# standard error function
se <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

# function to visualize the boxplot, histogram, QQ-plot, and kernel density estimate plot for residuals
normality.plots = function(x) {
  par(mfrow=c(2,2))
  hist(residuals(x), main = "Histogram", xlab = "Values")
  boxplot(residuals(x), main = "Boxplot", ylab = "Values")
  qqnorm(residuals(x))
  qqline(residuals(x))
  plot(density(residuals(x)), main = "Kernel Density Estimate")
}

# function to visaulize the fitted vs. residuals of model and the absolute residuals vs. fitted
var.plots = function(x) {
  par(mfrow=c(1,2))
  plot(x = fitted(x), y = residuals(x), xlab="Fitted", ylab="Residuals")
  abline(h=0)
  title("Residuals vs. Fitted")
  plot(x = fitted(x), y = abs(residuals(x)), xlab="Fitted", ylab="Absolute Residuals")
  abline(h=0)
  title("Absolute Residuals vs. Fitted")
}

##############################################################################################################################
##### OPENING DF
##############################################################################################################################

# opening Mortality csv file
mortality = read.csv(file="GLAMM_data_LMMs.csv", encoding = 'UTF-8', na.strings=c("", "NA"))

##############################################################################################################################
##### EXPLORING DF
##############################################################################################################################

# note that I already set up the df in another script, so there's not much to explore here... 
glimpse(mortality)

##############################################################################################################################
##### DATA ANALYSIS: MIXED EFFECTS MODEL FOR MORTALITY OVER TIME AFTER SNAIL ADDITION
##############################################################################################################################

# mixed effects model to test effects of Nutrient, Sediment, and Snail treatments on Mortality over time
# Nutrient, Sediment, and Snail are fixed effects
# Shell length (mm) as fixed effect - covariate included in the model
# Concrete Block is a random effect
# Colony is a random effect
# Date is a random effect (repeated measures)

# Treatments:
# Nutrient (Control, Enriched) - fixed effect
# Sediment (Control, Addition) - fixed effect
# Snail (Control, Present) - fixed effect
# Shell length (covariate) - fixed effect
# need to include Block as a random effect: + (1 | Block)
# need to include Colony as a random effect: + (1 | Colony)
# need to include Date as a random effect: + (1 | Date) # this is the repeated measures bit

# first looking at the distribution of the response variable
par(mfrow=c(1,2))
hist(mortality$PercMortality, main = "Histogram")
qqnorm(mortality$PercMortality)
qqline(mortality$PercMortality)
par(mfrow=c(1,1))
# data is not normal, looks 0 inflated

# spread of data across treatments
boxplot(mortality$PercMortality ~ mortality$Treatment)

# looking at spread of data across random effects
par(mfrow=c(1,1))
boxplot(mortality$PercMortality ~ mortality$Date) # spread of data across survey dates
boxplot(mortality$PercMortality ~ mortality$Colony) # spread of data across parent colony
boxplot(mortality$PercMortality ~ mortality$Block) # spread of data across concrete blocks
# data seems pretty evenly spread across survey dates
# looks like there's an effect of colony
# looks like there's an effect of block

# the distribution of the data aren't normally distributed
# looking to see if other distributions fit the data better
# need to first determine the best distribution for the data
# this code determines the probability distribution of my data 
# source: https://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
par(mfrow=c(2,2))
qqp(mortality$PercMortality, "norm", main = "Gaussian") # normal distribution
qqp(mortality$PercMortality, "lnorm", main = "Log Normal") # log normal distribution
nbinom <- fitdistr(mortality$PercMortality, "Negative Binomial") 
qqp(mortality$PercMortality, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]], main = "Negative Binomial") # negative binomial distribution
# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter.
poisson <- fitdistr(mortality$PercMortality, "Poisson")
qqp(mortality$PercMortality, "pois", lambda = poisson$estimate, main = "Poisson") # poisson distribution
# looks like log normal is the best fitting distribution

# Gaussian is clearly not a contender, so going to start with a LMM with log transformed data
# because there are zero values adding 1 to the data
mortality$log.mortality = log10(mortality$PercMortality + 1)
sort(mortality$log.mortality)

# distribution of the response variable
par(mfrow=c(1,2))
hist(mortality$log.mortality)
qqnorm(mortality$log.mortality)
qqline(mortality$log.mortality)

par(mfrow=c(1,1))
boxplot(mortality$log.mortality ~ mortality$Treatment)
# distribution still doesn't look great

# the first step in building a LMM is determining the random effects structure of the model
# NOTE ON RANDOM EFFECTS
# Date must be included in all versions of the model because it is the repeated measures component
# creating a full model with Colony and Block 
# creating a separate model with just Colony and just Block as the other random effect
# then I will test the fit of these different random effects structures using AICc

# NOTE ON FIXED EFFECTS
# because this experiment has a factorial design, all the fixed effects must be crossed
# and fixed effects may not be eliminated because it's dishonest to the original experimental design and a priori hypotheses
# thus the fixed effects structure for this model is set and should not change
# shell length is included as a covariate because snail size varied across treatment and I wanted to account for this variation in the model

# full model with all random effects (Colony, Block, and Date)
log.full.model = lmer(log.mortality ~ Snail * Nutrient * Sediment + Shell.length + (1|Colony) + (1|Block) + (1|Date), data = mortality, REML = TRUE)
# model with just Colony and Date as random effects
log.colony.model = lmer(log.mortality ~ Snail * Nutrient * Sediment + Shell.length + (1|Colony) + (1|Date), data = mortality, REML = TRUE)
# model with just Block and Date as random effects
log.block.model = lmer(log.mortality ~ Snail * Nutrient * Sediment + Shell.length + (1|Block) + (1|Date), data = mortality, REML = TRUE)

# need to check assumptions for random effects:
# 1) residuals are independently distributed
# 2) residuals are from a normal distribution
par(mfrow=c(1,3))
qqnorm(ranef(log.full.model)$Colony[,1], main = "Colony Residuals")
qqline(ranef(log.full.model)$Colony[,1])
qqnorm(ranef(log.full.model)$Block[,1], main = "Block Residuals")
qqline(ranef(log.full.model)$Block[,1])
qqnorm(ranef(log.full.model)$Date[,1], main = "Date Residuals")
qqline(ranef(log.full.model)$Date[,1])
# the residuals look relatively normal

# looking at the LRT for the random effects
# this tests the significance of the random effects
ranova(log.full.model) 
# colony is not significant (p = 0.52469)
# block is significant (p < 0.0001)
# date is significant (p = 0.01501)

# dotplot() function looks at the variation of the random effects
# note that it outputs a single plot for each random effect, so if you have multiple you have to scroll through them
dotplot(ranef(log.full.model, condVar = TRUE))
# this is a good visual to see how well the LRT test "matches" how variable the random effects are
# note that the variance for Colony overlaps with 0, suggesting it explains less variance
# the variance for Block is far from 0 for several of the groups, suggesting it explains more of the variance

# the next step is comparing the random effects structmure using AIC
# also using AICc because this experiment has small sample sizes (n = 9-11)
AIC(log.full.model, log.colony.model, log.block.model)
AICc(log.full.model, log.colony.model, log.block.model)
# there's a similar pattern between AIC and AICc
# the Colony model is the worst fitting, suggesting that Block is an important random effect
# more evidence for this comes from the above LRT test where Block is significant and Colony is not
# the full model and the Block model have similar AICc values
# a delta AICc value of <2 suggests that the models are the same, thus the fit of the full and Block model are similar
# I have chosen to drop Colony, forego using the full model, and use the Block model as my final model because
# it is the simpler of the two models and will save a degree of freedom

# checking residuals for normality
normality.plots(log.block.model) # residuals look normally distributed
# checking residuals for constant variance
var.plots(log.block.model) # residuals look evenly spread across the 0 line despite that strange band

# f-test for fixed effects
anova(log.block.model, type=3, ddf = "Kenward-Roger")
summary(log.block.model)

# post-hoc tests

# multiple comparisons
emmeans(log.block.model, list(pairwise ~ Snail * Nutrient), adjust = "tukey")
emmeans(log.block.model, list(pairwise ~ Sediment * Nutrient), adjust = "tukey")
emmeans(log.block.model, list(pairwise ~ Sediment * Snail), adjust = "tukey")

# within group comparisons
emmeans(log.block.model, pairwise ~ Snail | Nutrient)
emmeans(log.block.model, pairwise ~ Snail | Sediment)
emmeans(log.block.model, pairwise ~ Nutrient | Sediment)
