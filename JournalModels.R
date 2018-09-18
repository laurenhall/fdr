###############################################################################
## This file contains the models fit in the main text and the supplemental.  ##
## This code does not require any other source code to be run; simply        ##
## make sure the required data sets are available.                           ##
###############################################################################
# Load the required library
library(nlme)

######################
# Main text models  ##
######################

# Load in the data - the main text models use the nzns data set.
# This is the data set with both zero estimates and journals that switch OA status removed.

# Load in the data set
lmdatnzns = read.csv("data_final_nzns.csv", header=TRUE)

# Create data subsets by subject
oncdatnzns = subset(lmdatnzns, lmdatnzns$SUBJ=="ONC")
gendatnzns = subset(lmdatnzns, lmdatnzns$SUBJ=="GEN")

# Unstratified model selection - Full Model
lma1nzns = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ, random=~1|JRNL, data=lmdatnzns)
summary(lma1nzns)
# Remove OA*JIF*SUBJ
lma2nzns = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdatnzns)
summary(lma2nzns)
# Remove OA*SUBJ
lma3nzns = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ, random=~1|JRNL, data=lmdatnzns)
summary(lma3nzns)
# Remove OA*JIF
lma4nzns = lme(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdatnzns)
summary(lma4nzns)
# All remaining interactions are significant.

# Test the joint hypothesis that SUBJ and JIF*SUBJ significantly contribute to the model
# This gets us a better idea of whether SUBJ significantly affects EST.
# We use a likelihood ratio test for nested models.
# Fit the reduced model, without SUBJ and JIF*SUBJ
lma5nzns = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=lmdatnzns)
summary(lma5nzns)
# Compare this model to Model #4, the larger model.
# DF = 2 coefficients difference
cs = -2*(lma5nzns$logLik - lma4nzns$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified Model - General Medicine, full model
lmg1nzns = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=gendatnzns)
summary(lmg1nzns)
# General Medicine, final model (OA*JIF removed)
lmg2nzns = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=gendatnzns)
summary(lmg2nzns)

# Stratified Model - Oncology, full model
lmo1nzns = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=oncdatnzns)
summary(lmo1nzns)
# Oncology, final model (OA*JIF removed)
lmo2nzns = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=oncdatnzns)
summary(lmo2nzns)


###############################
# Supplementary text models  ##
###############################
# Here we fit the models in the supplemental material.
# Each of these uses a different data set.
# We repeat the analysis above for each of these sets.

###########################################################
# First supplemental data set - Full Data Set, none removed

lmdat = read.csv("data_final.csv", header=TRUE) # Read in data with JIF added
summary(lmdat)
# Data subsets
oncdat = subset(lmdat, lmdat$SUBJ=="ONC")
gendat = subset(lmdat, lmdat$SUBJ=="GEN")

# Ustratified model
lma1 = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ, random=~1|JRNL, data=lmdat)
summary(lma1)
lma2 = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdat)
summary(lma2)
lma3 = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ, random=~1|JRNL, data=lmdat)
summary(lma3)
lma4 = lme(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdat)
summary(lma4)
# Reduced model for likelihood ratio test
lma5 = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=lmdat)
summary(lma5)
# LRT
cs = -2*(lma5$logLik - lma4$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified - General Medicine
lmg1 = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=gendat)
summary(lmg1)
lmg2 = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=gendat)
summary(lmg2)

# Stratified - Oncology
lmo1 = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=oncdat)
summary(lmo1)
lmo2 = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=oncdat)
summary(lmo2)

########################################################
# Second supplemental data set - Zero estimates removed

lmdatnz = read.csv("data_final_nz.csv", header=TRUE)
# Data subsets
oncdatnz = subset(lmdatnz, lmdatnz$SUBJ=="ONC")
gendatnz = subset(lmdatnz, lmdatnz$SUBJ=="GEN")

# Unstratified model
lma1nz = lme(EST ~ YEAR +  OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ, random=~1|JRNL, data=lmdatnz)
summary(lma1nz)
lma2nz = lme(EST ~ YEAR +  OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdatnz)
summary(lma2nz)
lma3nz = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ, random=~1|JRNL, data=lmdatnz)
summary(lma3nz)
lma4nz = lme(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdatnz)
summary(lma4nz)
# Reduced model for likelihood ratio test
lma5nz = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=lmdatnz)
summary(lma5nz)
#LRT
cs = -2*(lma5nz$logLik - lma4nz$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified - General Medicine
lmg1nz = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=gendatnz)
summary(lmg1nz)
lmg2nz = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=gendatnz)
summary(lmg2nz)
# Stratified - Oncology
lmo1nz = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=oncdatnz)
summary(lmo1nz)
lmo2nz = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=oncdatnz)
summary(lmo2nz)

########################################################################
# Third supplemental data set - No journals that were both OA and not OA

lmdatns = read.csv("data_final_ns.csv", header=TRUE)
# Data subsets
oncdatns = subset(lmdatns, lmdatns$SUBJ=="ONC")
gendatns = subset(lmdatns, lmdatns$SUBJ=="GEN")
# Unstratified model
lma1ns = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ, random=~1|JRNL, data=lmdatns)
summary(lma1ns)
lma2ns = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdatns)
summary(lma2ns)
lma3ns = lme(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ, random=~1|JRNL, data=lmdatns)
summary(lma3ns)
lma4ns = lme(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ, random=~1|JRNL, data=lmdatns)
summary(lma4ns)
# Reduced model for likelihood ratio test
lma5ns = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=lmdatns)
summary(lma5ns)
# LRT
cs = -2*(lma5ns$logLik - lma4ns$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified - General Medicine
lmg1ns = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=gendatns)
summary(lmg1ns)
lmg2ns = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=gendatns)
summary(lmg2ns)

# Stratified - Oncology
lmo1ns = lme(EST ~ YEAR + OA + JIF + OA*JIF, random=~1|JRNL, data=oncdatns)
summary(lmo1ns)
lmo2ns = lme(EST ~ YEAR + OA + JIF, random=~1|JRNL, data=oncdatns)
summary(lmo2ns)

##
