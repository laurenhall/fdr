##############################################################################################
## This document contains all the necessary R packages and functions to replicate analysis. ##
##############################################################################################

# Load required libraries for analysis

library(RCurl)
library(XML)
library(tm)
library(stats4)
library(genefilter)
library(stringr)

#################################################################################################
## Function are based on versions provided by Jager and Leek (2013)                            ##
## See main text for citation.                                                                 ##
## These functions were written on a Linux machine and have not been tested on Mac or Windows. ##
#################################################################################################

# A function to download abstracts from PubMed
# This is an unmodified version of the function provided by Jager and Leek

getAbstractsPmids = function(journaltitle,year){
url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
q = paste("db=pubmed&term=",gsub(" ","+",journaltitle),"[ta]+AND+",year,"[dp]&usehistory=y",sep="")
esearch <- xmlTreeParse(getURL(paste(url, q, sep="")), useInternal = T)
webenv  <- xmlValue(getNodeSet(esearch, "//WebEnv")[[1]])
key     <- xmlValue(getNodeSet(esearch, "//QueryKey")[[1]])

url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
q   <- "db=pubmed&retmode=xml&rettype=abstract"
efetch <- xmlTreeParse(getURL(paste(url, q, "&WebEnv=", webenv, "&query_key=", key, sep="")), useInternal = T)
r = xmlRoot(efetch)
n = xmlSize(r)
abstracts = pmid = titles = rep(NA,n)
for(i in 1:n){abstracts[i] =  xmlValue(r[[i]][[1]][["Article"]][["Abstract"]]); pmid[i] = xmlValue(r[[i]][[1]][["PMID"]]); titles[i] = xmlValue(r[[i]][[1]][["Article"]][["ArticleTitle"]]) }
return(list(abs=abstracts,pmid=pmid,title=titles))
}

# A function to standardize notation, punctuation, and encoding in an abstract so it is easier to read.
# Since different articles in different journals may have different notation or encoding characters on PubMed, this standardizes several common differences to make them easier to use with regular expressions.
standardize = function(abstract){
# Replace unusually coded spaces with ordinary spaces
abstract = gsub(" "," ",abstract)
abstract = gsub(" "," ",abstract)
abstract = gsub(" "," ", abstract)
# Replace >= or <= notation, which improperly flags as =
abstract=gsub(">=|> =|=>|= >",">", abstract)
abstract=gsub("<=|< =|=<|= <","<", abstract)
# Replace brackets on numbers with parentheses
abstract = gsub("\\[","(",abstract)
abstract = gsub("\\]",")",abstract)
# Replace 'small dot' with periods
abstract = gsub("·",".",abstract)
# Replace unusual < notation
abstract = gsub("<;","<",abstract)
# Replace a common string with the appropriate character
abstract = gsub("\302\2400", "< ", abstract)
return(abstract)}

# A function to convert a read string into a number (for collecting p-values)
convert = function(string){
# Trim string
while(str_count(string, "[a-zA-Z]") > 1 || str_count(string, "[.,]") > 1 || length(grep("[0-9]",strsplit(string,"")[[1]][nchar(string)])) == 0){
    if(string == substr(string,1,(nchar(string)-1))){break;}
    string = substr(string,1,(nchar(string)-1))
  }
# Convert to a number
 string = gsub(" ","",string)
 if(is.na(as.numeric(string)) == TRUE){
    if(length(grep("[×x(]",string))>0){
    tmp = strsplit(string,"[\\*×x(]")[[1]]
    return(as.numeric(tmp[1])*10^(as.numeric(tmp[length(tmp)])))
  } else {
    return(as.numeric(string))}
} else {
  return(as.numeric(string))}
}


# A function to extract p-values from abstracts
# This function also saves the read string of characters.
# This function is based on the getPvalues function provided by Jager and Leek.
pget = function(abs,pmid){
  pvalues = trunc = ids = strs = numeric(0)
  abs = abs[!is.na(abs)]
  pmid = pmid[!is.na(abs)]
  abs = standardize(abs)
  # Collect truncated p-values
  ind = grep("[Pp][[:space:]]?[<≤]",abs)
    if(length(ind) > 0){
    for(i in 1:length(ind)){
    tmp = strsplit(abs[ind[i]],"[[:space:](][Pp][[:space:]]?[<≤]")[[1]][-1]
    n = length(tmp)
    for(j in 1:n){
      if(length(grep("[·.0123456789]",substr(tmp[j],1,3))) > 0){
        if(length(grep("[A-Z]",substr(tmp[j],1,1)))>0){next;}
        strng = substr(tmp[j],1,15)
        tmp2 = convert(tmp[j])
        pvalues = c(pvalues,tmp2)
        trunc = c(trunc,1)
        ids = c(ids,pmid[ind[i]])
        strs = c(strs,strng)
          }
        }
      }
    }
  # Collect exact/rounded p-values
  ind = grep("[Pp][[:space:]]?=",abs)
  if(length(ind) > 0){
  for(i in 1:length(ind)){
    tmp = strsplit(abs[ind[i]],"[[:space:](][Pp][[:space:]]?=")[[1]][-1]
    n = length(tmp)
    for(j in 1:n){
      if(length(grep("[·.0123456789]",substr(tmp[j],1,3))) > 0){
        if(length(grep("[A-Z]",substr(tmp[j],1,1)))>0){next;}
        strng = substr(tmp[j],1,15)
        tmp2 = convert(tmp[j])
        pvalues = c(pvalues,tmp2)
        trunc = c(trunc,0)
        ids = c(ids,pmid[ind[i]])
        strs = c(strs, strng)
          }
        }
      }
    }
  return(list(pvalues=pvalues,ids=ids,trunc=trunc, strs=strs))
}


# The following function is exactly reproduced as used by Jager and Leek
# The p.values.lmer function is reproduced from the blog:
# http://blog.lib.umn.edu/moor0554/canoemoore/2010/09/lmer_p-values_lrt.html
# and is due to Chrisotpher Moore (moor0554@umn.edu)
p.values.lmer <- function(x) {
  summary.model <- summary(x)
  data.lmer <- data.frame(model.matrix(x))
  names(data.lmer) <- names(fixef(x))
  names(data.lmer) <- gsub(pattern=":", x=names(data.lmer), replacement=".", fixed=T)
  names(data.lmer) <- ifelse(names(data.lmer)=="(Intercept)", "Intercept", names(data.lmer))
  string.call <- strsplit(x=as.character(x@call), split=" + (", fixed=T)
  var.dep <- unlist(strsplit(x=unlist(string.call)[2], " ~ ", fixed=T))[1]
  vars.fixef <- names(data.lmer)
  formula.ranef <- paste("+ (", string.call[[2]][-1], sep="")
  formula.ranef <- paste(formula.ranef, collapse=" ")
  formula.full <- as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef, collapse=" + "),
                  formula.ranef))
  data.ranef <- data.frame(x@frame[,
                which(names(x@frame) %in% names(ranef(x)))])
  names(data.ranef) <- names(ranef(x))
  data.lmer <- data.frame(x@frame[, 1], data.lmer, data.ranef)
  names(data.lmer)[1] <- var.dep
  out.full <- lmer(formula.full, data=data.lmer, REML=F)
  p.value.LRT <- vector(length=length(vars.fixef))
  for(i in 1:length(vars.fixef)) {
    formula.reduced <- as.formula(paste(var.dep, "~ -1 +", paste(vars.fixef[-i],
                       collapse=" + "), formula.ranef))
    out.reduced <- lmer(formula.reduced, data=data.lmer, REML=F)
    print(paste("Reduced by:", vars.fixef[i]))
    print(out.LRT <- data.frame(anova(out.full, out.reduced)))
    p.value.LRT[i] <- round(out.LRT[2, 7], 3)
  }
  summary.model@coefs <- cbind(summary.model@coefs, p.value.LRT)
  summary.model@methTitle <- c("\n", summary.model@methTitle,
                           "\n(p-values from comparing nested models fit by maximum likelihood)")
  print(summary.model)
}



# Function for calculating the science-wise FDR using the EM algorithm.
# This function is unmodified from the version provided by Jager and Leek.

calculateSwfdr = function(pValues,truncated,rounded,pi0 = 0.5,alpha=1,beta=50,numEmIterations=100){
  pp = pValues
  tt = truncated
  rr = rounded

  ll = function(a,b){
    tmp1 = rep(0,length(pp))
    tmp1[tt==0 & rr==0] = log(dbeta(pp[tt==0 & rr==0],a,b)/pbeta(0.05,a,b))
    tmp1[tt > 0 & rr==0] = log(pbeta(pp[tt > 0 & rr==0],a,b)/pbeta(0.05,a,b))
    tmp1 = -sum((1-z)*tmp1)

    probvec = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),a,b) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),a,b))/pbeta(0.05,a,b)
    probvec = rev(probvec)
    tmp2 = sum(-n1*log(probvec))


    return(tmp1 + tmp2)
  }

  n = table(cut(pp[rr > 0],c(-0.01,0.005,0.015,0.025,0.035,0.045,0.051)))


  for(i in 1:numEmIterations){

    ## E-step

     probvec1 = (pbeta(c(0.05,0.045,0.035,0.025,0.015,0.005),alpha,beta) - pbeta(c(0.045,0.035,0.025,0.015,0.005,0),alpha,beta))/pbeta(0.05,alpha,beta)
     probvec1 = rev(probvec1)
     probvec0 = c(0.005,0.01,0.01,0.01,0.01,0.005)*20

     pij0 = pi0*probvec0/(probvec0*pi0 + probvec1*(1-pi0))
     n0 = n*pij0
     n1 = n - n0

     z = rep(0,length(pp))
     z[tt == 0 & rr ==0] <- pi0*20/(pi0*20 + (1-pi0)*dbeta(pp[tt==0 & rr == 0],alpha,beta)/pbeta(0.05,alpha,beta))
     z[tt > 0 & rr ==0] <- pi0*20*pp[tt > 0 & rr ==0]/(pi0*20*pp[tt > 0 & rr==0] + (1-pi0)*pbeta(pp[tt > 0 & rr==0],alpha,beta)/pbeta(0.05,alpha,beta))

     ## M-step

     pi0 = (sum(n0) + sum(z))/(sum(n) + sum(rr == 0))
     tmp = mle(ll,start=list(a=0.05,b=100),lower=c(0.001,1),upper=c(1,500),method="L-BFGS-B")
     alpha =coef(tmp)[1]
     beta = coef(tmp)[2]
   }
  return(list(pi0 = pi0, alpha=alpha, beta = beta, z=z,n0=n0,n=n))
}
