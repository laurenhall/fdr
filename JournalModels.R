###############################################################################
## This file contains the models fit in the main text and the supplemental.  ##
## This code does not require any other source code to be run; simply        ##
## make sure the required data sets are available.                           ##
###############################################################################
# Load the required library
library(lme4)
library(lmerTest)
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
lma1nzns = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ + (1 | JRNL), data=lmdatnzns)
summary(lma1nzns)
# Remove OA*JIF*SUBJ
lma2nzns = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + (1|JRNL), data=lmdatnzns)
summary(lma2nzns)
# Remove OA*SUBJ
lma3nzns = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ + (1|JRNL), data=lmdatnzns)
summary(lma3nzns)
# Remove OA*JIF
lma4nzns = lmer(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ + (1|JRNL), data=lmdatnzns)
summary(lma4nzns)


# All remaining interactions are significant.

# Test the joint hypothesis that SUBJ and JIF*SUBJ significantly contribute to the model
# This gets us a better idea of whether SUBJ significantly affects EST.
# We use a likelihood ratio test for nested models.
# Fit the reduced model, without SUBJ and JIF*SUBJ
lma5nzns = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=lmdatnzns)
summary(lma5nzns)
# Compare this model to Model #4, the larger model.
anova(update(lma4nzns, . ~ ., REML = F),update(lma5nzns, . ~ ., REML = F))

# Stratified Model - General Medicine, full model
lmg1nzns = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=gendatnzns)
summary(lmg1nzns)
# General Medicine, final model (OA*JIF removed)
lmg2nzns = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=gendatnzns)
summary(lmg2nzns)

# Stratified Model - Oncology, full model
lmo1nzns = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=oncdatnzns)
summary(lmo1nzns)
# Oncology, final model (OA*JIF removed)
lmo2nzns = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=oncdatnzns)
summary(lmo2nzns)

############################################################
# Bootstrapping confidence bands/Creating interaction plot #
############################################################
# Here we create the main text interaction plot of Subject and JIF
# Using confidence intervals created with bootMer()

# Select the unstratified model
model = lma4nzns

# Set the number of bootstrap replicates
nsim = 1e4

# Create the new data to get predictions for - Covering all levels
JIF = seq(from=0.1, to=60, by=0.1)
OA = 0.5; subj=c("GEN","ONC"); year=2013
newdata = data.frame(expand.grid(JIF,OA,subj,year)); colnames(newdata) = c("JIF","OA","SUBJ","YEAR")

# Perform the bootstrapping
newdata$pred = predict(model, newdata=newdata, re.form=NA)
boot.out = bootMer(model, nsim=nsim, FUN=function(x)predict(x, newdata=newdata, re.form=NA), .progress="txt")
pred.int = apply(boot.out$t, 2, function(x) quantile(x,c(0.025, 0.975)))
newdata$low = pred.int[1,]
newdata$high = pred.int[2,]

newdata_extrapolated = subset(newdata, newdata$JIF > 26.5 & newdata$SUBJ == "ONC")
newdata_NE = subset(newdata, newdata$JIF < 26.5 | newdata$SUBJ=="GEN")

library(ggplot2)

interaction_plot = ggplot(newdata,aes(x=JIF,y=pred)) +
  geom_line(data=newdata_NE, aes(col = factor(SUBJ, labels=c("Medicine","Oncology"))), size=1, linetype="solid") +
  geom_line(data=newdata_extrapolated, col="#00BFC4", size=1, linetype="dashed") +
  geom_ribbon(aes(ymin=low, ymax=high, fill=factor(SUBJ, labels=c("Medicine","Oncology"))),alpha=0.3) +
  scale_y_continuous("Estimated FDR", breaks=scales::pretty_breaks(n=10))+
  scale_x_continuous("Journal Impact Factor", breaks=scales::pretty_breaks(n=10)) +
  coord_cartesian(ylim=c(0,.35)) +
  labs(fill="Subject", col="Subject")
  #+ ggtitle("Interaction of JIF and Subject, All Journals")

tiff(file="Fig1.tiff", height = 1250, width = 2000,
     compression = "lzw", res = 300)
interaction_plot
dev.off()

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
lma1 = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ + (1|JRNL), data=lmdat)
summary(lma1)
lma2 = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + (1|JRNL), data=lmdat)
summary(lma2)
lma3 = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ + (1|JRNL), data=lmdat)
summary(lma3)
lma4 = lmer(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ + (1|JRNL), data=lmdat)
summary(lma4)
# Reduced model for likelihood ratio test
lma5 = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=lmdat)
summary(lma5)
# LRT
cs = -2*(lma5$logLik - lma4$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified - General Medicine
lmg1 = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=gendat)
summary(lmg1)
lmg2 = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=gendat)
summary(lmg2)

# Stratified - Oncology
lmo1 = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=oncdat)
summary(lmo1)
lmo2 = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=oncdat)
summary(lmo2)

#####################################
# Create Supplementary Figure S3-S4 #
#####################################
# Since figures S3-S4 are based on this data set, we generate them here. 
library(ggplot2)

fdr_jif_subject = ggplot(lmdat,aes(x=JIF,y=EST)) +
  geom_point(aes(col = factor(SUBJ, labels=c("Medicine","Oncology")))) +
  scale_y_continuous("Estimated FDR", breaks=scales::pretty_breaks(n=10))+
  scale_x_continuous("Journal Impact Factor", breaks=scales::pretty_breaks(n=10)) +
  labs(col="Subject") 
  #+ ggtitle("Disribution of FDR Estimates by JIF and Subject, All Journals")
fdr_jif_subject

tiff(file="FigS2A.tiff", height = 1250, width = 2000,
     compression = "lzw", res = 300)
fdr_jif_subject
dev.off()

fdr_jif_access = ggplot(lmdat,aes(x=JIF,y=EST)) +
  geom_point(aes(col = factor(OA, labels = c("Closed","Open")))) +
  scale_y_continuous("Estimated FDR", breaks=scales::pretty_breaks(n=10))+
  scale_x_continuous("Journal Impact Factor", breaks=scales::pretty_breaks(n=10)) +
  labs(col="Access") 
  #+ ggtitle("Disribution of FDR Estimates by JIF and Access, All Journals")
fdr_jif_access

tiff(file="FigS2B.tiff", height = 1250, width = 2000,
     compression = "lzw", res = 300)
fdr_jif_access
dev.off()

########################################################
# Second supplemental data set - Zero estimates removed

lmdatnz = read.csv("data_final_nz.csv", header=TRUE)
# Data subsets
oncdatnz = subset(lmdatnz, lmdatnz$SUBJ=="ONC")
gendatnz = subset(lmdatnz, lmdatnz$SUBJ=="GEN")

# Unstratified model
lma1nz = lmer(EST ~ YEAR +  OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ + (1|JRNL), data=lmdatnz)
summary(lma1nz)
lma2nz = lmer(EST ~ YEAR +  OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + (1|JRNL), data=lmdatnz)
summary(lma2nz)
lma3nz = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ + (1|JRNL), data=lmdatnz)
summary(lma3nz)
lma4nz = lmer(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ + (1|JRNL), data=lmdatnz)
summary(lma4nz)
# Reduced model for likelihood ratio test
lma5nz = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=lmdatnz)
summary(lma5nz)
#LRT
cs = -2*(lma5nz$logLik - lma4nz$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified - General Medicine
lmg1nz = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=gendatnz)
summary(lmg1nz)
lmg2nz = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=gendatnz)
summary(lmg2nz)
# Stratified - Oncology
lmo1nz = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=oncdatnz)
summary(lmo1nz)
lmo2nz = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=oncdatnz)
summary(lmo2nz)

########################################################################
# Third supplemental data set - No journals that were both OA and not OA

lmdatns = read.csv("data_final_ns.csv", header=TRUE)
# Data subsets
oncdatns = subset(lmdatns, lmdatns$SUBJ=="ONC")
gendatns = subset(lmdatns, lmdatns$SUBJ=="GEN")
# Unstratified model
lma1ns = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + OA*JIF*SUBJ + (1|JRNL), data=lmdatns)
summary(lma1ns)
lma2ns = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + OA*SUBJ + JIF*SUBJ + (1|JRNL), data=lmdatns)
summary(lma2ns)
lma3ns = lmer(EST ~ YEAR + OA + JIF + SUBJ + OA*JIF + JIF*SUBJ + (1|JRNL), data=lmdatns)
summary(lma3ns)
lma4ns = lmer(EST ~ YEAR + OA + JIF + SUBJ + JIF*SUBJ + (1|JRNL), data=lmdatns)
summary(lma4ns)
# Reduced model for likelihood ratio test
lma5ns = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=lmdatns)
summary(lma5ns)
# LRT
cs = -2*(lma5ns$logLik - lma4ns$logLik)
pchisq(cs, 2, lower.tail = FALSE)

# Stratified - General Medicine
lmg1ns = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=gendatns)
summary(lmg1ns)
lmg2ns = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=gendatns)
summary(lmg2ns)

# Stratified - Oncology
lmo1ns = lmer(EST ~ YEAR + OA + JIF + OA*JIF + (1|JRNL), data=oncdatns)
summary(lmo1ns)
lmo2ns = lmer(EST ~ YEAR + OA + JIF + (1|JRNL), data=oncdatns)
summary(lmo2ns)
