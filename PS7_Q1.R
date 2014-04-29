#########################
###### SHRINKAGE ########
#########################
# Author: Kory Kantenga
# Date  : 25-04-2014
# Macroeconometrics IV

setwd("~/Dropbox/Work/2nd_Year/Macroeconometrics/PS7/Q1")

# Control parameters
k    = 30
T    = 100
xdim = k*T
edim = T
nSim = 500

# Designs
beta_A = t(t(c(runif(10,min=0.5,max=1),runif(20,min=0,max=0.3))))
beta_B = t(t(runif(30,min=0.5,max=1)))
beta_C = t(t(c(rep(0,20),runif(10,min=0.5,max=1))))

# Penalty Terms
llambda<-seq(0,30, by = 0.5) # Ridge
lllambda<-c(1.25,1.034132318,0.856926272,0.710085762,0.588407435,0.487579568,0.404029285,
            0.334795947,0.277426242,0.229887251,0.190494409,0.157851815,0.130802766,0.108388767
            ,0.089815568,0.074425021,0.061671756,0.051103855,0.042346840,0.035090403,0.029077409
            ,0.024094785,0.019965969,0.016544656,0.013709610,0.011360369,0.009413688,0.007800584
            ,0.006463898,0.005356262,0.004438428,0.003677871,0.003047641,0.002525406,0.002092659
            ,0.001734067,0.001436922,0.001190695,0.001)

# Ridge Estimator
ridge_p <- function(y,x,l,xnew){
  S <- svd(x)
  D <- diag(S$d)
  w <- matrix(data = NA, nrow = p, ncol = length(llambda))
  for (o in 1:length(l)) {
    w[,o] = xnew%*%(S$v)%*%solve(D%*%D+l[o]*diag(k))%*%D%*%t(S$u)%*%y
  }
  return(w)
}

SPE_RidA <- matrix(data = NA, nrow = nSim, ncol = length(llambda))
SPE_RidB <- matrix(data = NA, nrow = nSim, ncol = length(llambda))
SPE_RidC <- matrix(data = NA, nrow = nSim, ncol = length(llambda))

LSPE_A <- matrix(data = NA, nrow = nSim, ncol = length(lllambda))
LSPE_B <- matrix(data = NA, nrow = nSim, ncol = length(lllambda))
LSPE_C <- matrix(data = NA, nrow = nSim, ncol = length(lllambda))

PCSPE_A <- matrix(data = NA, nrow = nSim, ncol = k)
PCSPE_B <- matrix(data = NA, nrow = nSim, ncol = k)
PCSPE_C <- matrix(data = NA, nrow = nSim, ncol = k)

for (j in 1:nSim) {
  x  <-rnorm(xdim)
  e  <-rnorm(edim)
  x  <-matrix(data=x,nrow=T,ncol=k)
  e  <-matrix(data=e,nrow=T,ncol=1)
  
  train <- sample(1:nrow(x), size = nrow(x)/2)
  test <- (-train)
  
  y_A<- x%*%beta_A + e
  y_B<- x%*%beta_B + e
  y_C<- x%*%beta_C + e
  
  #   #RIDGE
  y_A_p       <- matrix(rep(y_A[test],length(llambda)), ncol = length(llambda))
  SPE_RidA[j,] <- colMeans((y_A_p - ridge_p(y_A[train],x[train,],llambda,x[test,]))^2)
  
  y_B_p       <- matrix(rep(y_B[test],length(llambda)), ncol = length(llambda))
  SPE_RidB[j,] <- colMeans((y_B_p - ridge_p(y_B[train],x[train,],llambda,x[test,]))^2)
  
  y_C_p       <- matrix(rep(y_C[test],length(llambda)), ncol = length(llambda))
  SPE_RidC[j,] <- colMeans((y_C_p - ridge_p(y_C[train],x[train,],llambda,x[test,]))^2)
  
  #LASSO
  lasso.modA <- glmnet(x[train,], y_A[train], alpha = 1,lambda=lllambda)
  lasso.modB <- glmnet(x[train,], y_B[train], alpha = 1,lambda=lllambda)
  lasso.modC <- glmnet(x[train,], y_C[train], alpha = 1,lambda=lllambda)
  
  lasso.predA <- predict(lasso.modA,newx = x[test,])
  lasso.predB <- predict(lasso.modB,newx = x[test,])
  lasso.predC <- predict(lasso.modC,newx = x[test,])
  
  yl_A_p      <- matrix(rep(y_A[test],length(lllambda)), ncol = length(lllambda))
  yl_B_p      <- matrix(rep(y_B[test],length(lllambda)), ncol = length(lllambda))
  yl_C_p      <- matrix(rep(y_C[test],length(lllambda)), ncol = length(lllambda))
  
  LSPE_A[j,] <- colMeans((yl_A_p - lasso.predA)^2)
  LSPE_B[j,] <- colMeans((yl_B_p - lasso.predB)^2)
  LSPE_C[j,] <- colMeans((yl_C_p - lasso.predC)^2)
  
  # Principle Components
  
  pcr.A <- pcr(y_A[train]~x[train,])
  pcr.B <- pcr(y_B[train]~x[train,])
  pcr.C <- pcr(y_C[train]~x[train,])
  
  pcr.predA <- predict(pcr.A,newx = x[test,])
  pcr.predB <- predict(pcr.B,newx = x[test,])
  pcr.predC <- predict(pcr.C,newx = x[test,])
  
  pcr.predA <- pcr.predA[,1,]
  pcr.predB <- pcr.predB[,1,]
  pcr.predC <- pcr.predC[,1,]
  
  y_A_PC      <- matrix(rep(y_A[test],k), ncol = k)
  PCSPE_A[j,] <- colMeans((y_A_PC - pcr.predA)^2)
  y_B_PC      <- matrix(rep(y_B[test],k), ncol = k)
  PCSPE_B[j,] <- colMeans((y_B_PC - pcr.predB)^2)
  y_A_PC      <- matrix(rep(y_A[test],k), ncol = k)
  PCSPE_C[j,] <- colMeans((y_A_PC - pcr.predC)^2)
  
}

MSPE_RidA <- colMeans(SPE_RidA)
MSPE_RidB <- colMeans(SPE_RidB)
MSPE_RidC <- colMeans(SPE_RidC)

setEPS()
postscript("Ridge.eps")
par(mfrow=c(3,1))
plot(llambda,MSPE_RidA,type="b",main="Design A",ylab="Ridge MSPE",xlab=expression(lambda))
abline(h=MSPE_RidA[1],col = "blue")
plot(llambda,MSPE_RidB,type="b",main="Design B",ylab="Ridge MSPE",xlab=expression(lambda))
abline(h=MSPE_RidB[1],col = "blue")
plot(llambda,MSPE_RidC,type="b",main="Design C",ylab="Ridge MSPE",xlab=expression(lambda))
abline(h=MSPE_RidC[1],col = "blue")
dev.off()

MSPE_LassA <- colMeans(LSPE_A)
MSPE_LassB <- colMeans(LSPE_B)
MSPE_LassC <- colMeans(LSPE_C)

setEPS()
postscript("Lasso.eps")
par(mfrow=c(3,1))
plot(lllambda,MSPE_LassA,type="b",main="Design A",ylab="Lasso MSPE",xlab=expression(lambda))
abline(h=MSPE_RidA[1],col = "blue")
plot(lllambda,MSPE_LassB,type="b",main="Design B",ylab="Lasso MSPE",xlab=expression(lambda))
abline(h=MSPE_RidB[1],col = "blue")
plot(lllambda,MSPE_LassC,type="b",main="Design C",ylab="Lasso MSPE",xlab=expression(lambda))
abline(h=MSPE_RidC[1],col = "blue")
dev.off()

MSPE.PCR.A <- colMeans(PCSPE_A)
MSPE.PCR.B <- colMeans(PCSPE_B)
MSPE.PCR.C <- colMeans(PCSPE_C)

setEPS()
postscript("PCR.eps")
par(mfrow=c(3,1))
plot(1:1:30,MSPE.PCR.A,type="b",main="Design A",ylab="PCR MSPE",xlab="Principle Components")
abline(h=MSPE_RidA[1],col = "blue")
plot(1:1:30,MSPE.PCR.B,type="b",main="Design B",ylab="PCR MSPE",xlab="Principle Components")
abline(h=MSPE_RidB[1],col = "blue")
plot(1:1:30,MSPE.PCR.C,type="b",main="Design C",ylab="PCR MSPE",xlab="Principle Components")
abline(h=MSPE_RidC[1],col = "blue")
dev.off()