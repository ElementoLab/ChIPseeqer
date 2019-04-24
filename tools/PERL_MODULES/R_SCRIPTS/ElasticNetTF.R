library(glmnet)

source("/Users/olivier/PROGRAMS/PLASMODIUM/EXPRESSION/IDC/lasso.R");

	
getActivityCoefs <- function(suf, alpha=0.05, shuffled=F, sortby=1) {  

  ef <- sprintf("%s.expcount.exp.weighted", suf)
  if (shuffled == T) {
    ef <- sprintf("%s.shuffled", ef)
  }
  mf <- sprintf("%s.expcount.mot.weighted", suf)
  e <<- data.frame(read.csv(ef, sep="\t", row.names=1, header=T))
  m <<- as.matrix (read.csv(mf, sep="\t", row.names=1, header=T))     

  coefs <- PBM.glmnet.getCoefs(e, m, alpha=alpha)
  if (sortby != 0) {
    coefs <- as.matrix(coefs[order(coefs[,sortby]),])
  }
  return(coefs)
}