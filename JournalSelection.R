###############################################################
## Code for Journal selection and P-value collection         ##
##                                                           ##
## Journal names were taken from InCites by Thompson-Reuters ##
###############################################################
# NOTE: Before running this code, ensure that you have the following libraries
# installed: RCurl, XML, genefilter, tm, stats4, and stringr.
# RCurl, XML, and genefilter may require more extensive work to install.
# See the documentation on these packages for more information.


# Load in the required function and libraries
source("Functions.R")

# The following functions were written under the assumption that a multicore computer
# is available. If one is not available, or this code is being run on a Windows machine,
# set cores = 1.
# Otherwise, set cores = the # of cores that you would like to use.

cores=16
# Read in the CSV with journal candidates
jdat = read.csv("journals.csv", header=TRUE)

#
names(jdat); length(jdat$Abbreviation)
journals=jdat$Abbreviation

# Years to assess
years = 2011:2015

searchyj = expand.grid(years, journals)

# Collect abstracts from journal candidates
absout = mclapply(seq_len(nrow(search)), function(i){
  at = c()
   tmp = getAbstractsPmids(searchyj[i,2], search[i,1])
   abt = cbind(rep(searchyj[i,2], length(tmp$abs)), rep(search[i,1], length(tmp$abs)),matrix(unlist(tmp), ncol=3))
   at = rbind(at, abt)
  return(at)}, mc.cores=cores)
absdat = matrix(unlist(absout), ncol=5, byrow=TRUE) # Convert to a matrix
colnames(absdat) = c("journal","year","abstract","pmid","title")

# Extract summary
npaper = matrix(NA,nrow=length(journals),ncol=5); colnames(npaper) = years; rownames(npaper) = journals
for(i in 1:length(journals)){
for(j in 1:length(years)){
tmp = subset(absdat, absdat[,1] == journals[i] & absdat[,2] == years[j])
npaper[i,j] = dim(tmp)[1]
}}

#save(absdat,npaper, file = "absdat.rda") # This is a very large file!

# Remove journals with no data for one or more years
nep = which(apply(npaper, 1, min) ==0)
jsv = journals[-nep]

# Create p-value data objects
pvaldat = matrix(NA, ncol=8)
colnames(pvaldat) = c("journal","year","pvalue","truncated","pmid","abstract","title","string")
nna = nmr = npval = np05 = npwp
ind = 1:length(jsv)
# Collect p-values from journals with enough papers

searchyj2 = expand.grid(years, jsv)
pdatout = mclapply(seq_len(nrow(searchyj2)), function(x){
pt = numeric(0)
tmp = subset(absdat, absdat[,1] == searchyj2[x,2] & absdat[,2] == searchyj[x,1])
ptmp = getPvalues(tmp[,3],tmp[,4])
# Some journals in this sample reported p = 0.0000. Here, we re-record these as p < 0.0001.
ptmp$trunc[which(ptmp$pvalues == 0)] = 1
ptmp$pvalues[which(ptmp$pvalues == 0)] = 0.001
npv = length(ptmp$pvalues)
match = match(ptmp$ids, tmp[,4])
tmat = cbind(rep(searchyj2[x,2], npv), rep(searchyj[x,1],npv), ptmp$pvalues, ptmp$trunc, ptmp$ids, tmp[match,3], tmp[match,5], ptmp$strs)
pt = rbind(pt, tmat)
return(pt)}, mc.cores=cores) 
pvaldat = matrix(unlist(pdatout), ncol=8, by.row=TRUE) 
colnames(pvaldat) = c("journal","year","pvalue","truncated","pmid","abstract","title","string") 

# Extract summaries from p-value data. These are useful to try to identify journals that are being misread.
# Legend: npwp = Number of papers with p-values
#         nna = Number of p-values recorded as NA
#         nmr = Number obviously misread (as < 0 or > 1)
#         npval = Total number of p-values
#         np05 = number of p-values under 0.05
npwp = matrix(NA, nrow=length(jsv), ncol=length(years)); rownames(npwp) = jsv; colnames(npwp) = years
nna = nmr = npval = np05 = npwp

for(i in 1:length(jsv)){
for(j in 1:length(years)){
tmp = subset(pvaldat, pvaldat[,1] == jsv[i] & pvaldat[,2] == years[j])
npwp[i,j] = length(unique(tmp[,5]))
nna[i,j] = sum(is.na(tmp[,3]))
nmr[i,j] = sum(as.numeric(tmp[,3]) < 0 | as.numeric(tmp[,3]) > 1)
npval[i,j] = length(tmp[!is.na(tmp[,3]),3])
np05[i,j] = sum(as.numeric(tmp[!is.na(tmp[,3]),3]) < 0.05)
}}

#save(pvaldat, npwp, nna, nmr, npval, np05, file="pvaldat.rda")

# Getting ready to estimate the SWFDR for the selected journals.
pvaldat = pvaldat[!is.na(pvaldat[,3]),] # Remove NA vaues
# Extract the p-value data as a numeric matrix
pdata = matrix(as.numeric(pvaldat[,c(2:5)]),ncol=4); rownames(pdata) = pvaldat[,1]; colnames(pdata) = colnames(pvaldat)[2:5]

pi0out = mclapply(seq_len(nrow(searchyj2)), function(x){
# Some of these will fail to produce estimates. In these cases, we want to return NA.
 tc = tryCatch({
 tmp = subset(pdata, pdata[,1]==searchyj2[x,2] & rownames(pdata) == searchyj2[x,1] & pdata[,2] < 0.05 & pdata[,2] > 0)
 trnc = tmp[,3]
 rnd = rep(0,length(trnc))
 rnd[trnc == 0] = (tmp[trnc==0,2] == round(tmp[trnc==0,2],2))
 pv = tmp[,2]
 est = calculateSwfdr(pValues = pv, truncated = trnc, rounded = rnd, numEmIterations=100)
 }, error=function(e) e)
 if(!inherits(tc, "error")){
 return(est$pi0)} else {return(NA)}
}, mc.cores=cores)
pi0JY = matrix(unlist(pi0out), ncol=length(years), byrow=TRUE)
colnames(pi0JY) = years; rownames(pi0JY) = jsv

# Remove journals that produced a NA estimate for one or more years
noest = which(apply(pi0JY, 1, function(x){sum(is.na(x))}) > 0)
jsv2 = jsv[-noest]

# Collect data for export - This writes the raw data to a CSV.
# The format is: Journal, Year, Estimate
est = as.vector(t(pi0JY[-noest,]))
jrnl = rep(jsv2, each=5)
yr = rep(2011:2015, times=length(jsv2))
expdat = data.frame("JRNL" = jrnl, "YEAR" = yr, "EST" = est)
write.csv(expdat, file="modeldat_base.csv")
# The data set above was prepared for analysis by the following:
# JIF, Open Access status, and Subject were added to the CSV manually.
# Journals with at least one NA estimate, or journals that did not have JIF data for one or more years, were removed.
# The dataset used to fit the text models can be found in "data_final.csv".
#
# This dataset was reduced to form the following data sets:
# data_final_nz.csv: The final data set, with journals that produced zero estimates removed.
# data_final_ns.csv: The final data set, with journals that switched Open Access status removed.
# data_final_nzns.csv: the final data set with zero estimates and switching journals removed.
