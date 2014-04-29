#########################
###### FA/PCA ###########
#########################
# Author: Kory Kantenga
# Date  : 25-04-2014
# Macroeconometrics IV

setwd("~/Dropbox/Work/2nd_Year/Macroeconometrics/PS7/Q2")
library('ifa')
library('xtable')

# DGP

set.seed(5485)

z   = rnorm(200)
e_1 = rnorm(200)
e_2 = rnorm(200)
e_3 = rnorm(200, sd=25)

X <- cbind(z+e_1,z+e_2,z+e_3)

# PCA and EM algorithm

pr.comp  <- princomp(X,center = FALSE)
estimate <- pr.comp$scores[,1]
EM1      <- ifa.em(X,ni=1)

X[,1]<-(X[,1]-mean(X[,1]))/var(X[,1]) # Standardize X1
X[,2]<-(X[,2]-mean(X[,2]))/var(X[,2]) # Standardize X2
X[,3]<-(X[,3]-mean(X[,3]))/var(X[,3]) # Standardize X3
pr.compST <- princomp(X,center = TRUE)
estimateST<- pr.compST$scores[,1]

CORR    <- cbind(cor(estimate,z),cor(estimateST,z))
EM      <- ifa.em(X,ni=1)

# Results
correlation <- xtable(CORR)
loadingsEM  <- xtable(EM$H)
loadingsEM1  <- xtable(EM1$H)

lm(X[,1]~estimate-1)
lm(X[,2]~estimate-1)
lm(X[,3]~estimate-1)

lm(X[,1]~estimateST-1)
lm(X[,2]~estimateST-1)
lm(X[,3]~estimateST-1)
