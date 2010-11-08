source("logistic.R")
source("fr.R")

# 13 time points
tt <- seq(0,6,by=0.5)
# two groups of n each
n <- 200
grp <- c(rep(0,n),rep(1,n))

# natural spline basis with 6+1 df
psi6 <- ns(tt,df=6)

# parameter values modeled after Ma et. al. paper
# correlation is chosen to be consistent with Matern covariance
# function to be considered later

rr2 <- matern(0.5,1,2)
rr1 <- matern(0.5,1,1)
rr0 <- matern(0.5,1,0.5)

beta21 <- matrix(nrow=2,ncol=5)
beta21[1,] <- c(1,9,1,rr2,1/100)
beta21[2,] <- c(19/20,17/2,1,rr2,1/100)

beta11 <- beta21
beta11[,4] <- c(rr1,rr1)

beta01 <- beta21
beta01[,4] <- c(rr0,rr0)

# null parameter value is average of alternative
beta20 <- apply(beta21,2,mean)
beta20 <- rbind(beta20,beta20)

beta10 <- apply(beta11,2,mean)
beta10 <- rbind(beta10,beta10)

beta00 <- apply(beta01,2,mean)
beta00 <- rbind(beta00,beta00)

################################################

nsim <- 1000

#########################
out000 <- compareSimAR(beta00,nsim,tt,psi6,grp,df=0)
out004 <- compareSimAR(beta00,nsim,tt,psi6,grp,df=4)
out00m <- compareSimMat(beta00,nsim,tt,psi6,grp,phi=1,kappa=2)

out010 <- compareSimAR(beta01,nsim,tt,psi6,grp,df=0)
out014 <- compareSimAR(beta01,nsim,tt,psi6,grp,df=4)
out01m <- compareSimMat(beta01,nsim,tt,psi6,grp,phi=1,kappa=2)

#########################
out100 <- compareSimAR(beta10,nsim,tt,psi6,grp,df=0)
out104 <- compareSimAR(beta10,nsim,tt,psi6,grp,df=4)
out10m <- compareSimMat(beta10,nsim,tt,psi6,grp,phi=1,kappa=2)

out110 <- compareSimAR(beta11,nsim,tt,psi6,grp,df=0)
out114 <- compareSimAR(beta11,nsim,tt,psi6,grp,df=4)
out11m <- compareSimMat(beta11,nsim,tt,psi6,grp,phi=1,kappa=2)

#########################
out200 <- compareSimAR(beta20,nsim,tt,psi6,grp,df=0)
out204 <- compareSimAR(beta20,nsim,tt,psi6,grp,df=4)
out20m <- compareSimMat(beta20,nsim,tt,psi6,grp,phi=1,kappa=2)

out210 <- compareSimAR(beta21,nsim,tt,psi6,grp,df=0)
out214 <- compareSimAR(beta21,nsim,tt,psi6,grp,df=4)
out21m <- compareSimMat(beta21,nsim,tt,psi6,grp,phi=1,kappa=2)


# tests for the null distribution
ksTest(out000)
ksTest(out004)
ksTest(out00m)

ksTest(out100)
ksTest(out104)
ksTest(out10m)

ksTest(out200)
ksTest(out204)
ksTest(out20m)

# power in the alternative case
powerCalc(out000,out010)
powerCalc(out100,out110)
powerCalc(out200,out210)

powerCalc(out004,out014)
powerCalc(out104,out114)
powerCalc(out204,out214)

powerCalc(out00m,out01m)
powerCalc(out10m,out11m)
powerCalc(out20m,out21m)

