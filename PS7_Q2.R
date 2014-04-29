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

z   = rnorm(200)
e_1 = rnorm(200)
e_2 = rnorm(200)
e_3 = rnorm(200, sd=25)

X <- cbind(z+e_1,z+e_2,z+e_3)

# PCA and EM algorithm

pr.comp  <- princomp(X,center = FALSE)
estimate <- pr.compS$scores[,1]
llambda  <- lm(estimate~X-1)

pr.compST <- princomp(X,center = TRUE)
estimateST<- pr.compST$scores[,1]
llambdaST <- lm(estimateST~X-1)

LLAMBDA <- cbind(llambda$coefficients, llambdaST$coefficients)
CORR    <- cbind(cor(estimate,z),cor(estimateST,z))
EM      <- ifa.em(X,ni=1)

# Results
l_lambda    <- xtable(LLAMBDA)
correlation <- xtable(CORR)
loadingsEM  <- xtable(EM$H)
