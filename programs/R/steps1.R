source("R/fr.R")
library(qtl)

########################
prob <- read.csv("data/N1_CTvsProb_BW6min.csv")
prob <- as.matrix(prob[,-1])
prob <- prob[,-(1:2)]
########################

# time, number of time points, and number of individuals
tt <- 1:220
hr <- tt/10
n <- nrow(prob)
m <- ncol(prob)

# cross data; read, re-estimate map, delete X chromosome
cr <- read.cross(file="data/N1_GenoTypes.csv",format="csv")
mp <- est.map(cr)
cr <- replace.map(cr,mp)
cr <- subset(cr,ind=1:89)
cr <- subset(cr,chr=1:19)

# calculate genotype probabilities
cr <- calc.genoprob(cr,step=5)

acttot <- apply(prob,1,mean)
actday <- apply(prob[,c(1:50,181:220)],1,mean)
acteve <- apply(prob[,51:71],1,mean)
actnit <- apply(prob[,71:160],1,mean)
actmor <- apply(prob[,161:180],1,mean)

out0 <- scanone(cr,pheno.col=acttot,method="hk")
out1 <- scanone(cr,pheno.col=actday,method="hk")
out2 <- scanone(cr,pheno.col=acteve,method="hk")
out3 <- scanone(cr,pheno.col=actnit,method="hk")
out4 <- scanone(cr,pheno.col=actmor,method="hk")

out0.perm <- scanone(cr,pheno.col=acttot,method="hk",n.perm=1000)
out1.perm <- scanone(cr,pheno.col=actday,method="hk",n.perm=1000)
out2.perm <- scanone(cr,pheno.col=acteve,method="hk",n.perm=1000)
out3.perm <- scanone(cr,pheno.col=actnit,method="hk",n.perm=1000)
out4.perm <- scanone(cr,pheno.col=actmor,method="hk",n.perm=1000)

plot(out0,out2,out3)
add.threshold(out0,perms=out0.perm)
plot(out1,out4,ylim=c(0,4))
add.threshold(out4,perms=out4.perm)

phi16 <- bs(hr,df=16,intercept=FALSE)
out16 <- funcScanone(prob,cr,phi16,crit="ss")
out16.perm <- funcScanonePerm(prob,cr,phi16,crit="ss",nperm=1000)

phi20 <- ns(hr,df=20,intercept=FALSE)
out20 <- funcScanone(prob,cr,phi20,crit="ss")
out20.perm <- funcScanonePerm(prob,cr,phi20,crit="ss",nperm=1000)

phi40 <- ns(hr,df=40,intercept=FALSE)
out40 <- funcScanone(prob,cr,phi40,crit="ss")
out40w <- funcScanone(prob,cr,phi40,crit="ss",weightPhi=1/diag(var(prob)))



put01 <- funcScanone(prob,cr,NULL,crit="ss")
put02 <- funcScanone(prob,cr,NULL,crit="qf")

psi5 <- bs(hr,df=5,degree=1,intercept=FALSE)
put51 <- funcScanone(prob,cr,psi5,crit="ss")
put52 <- funcScanone(prob,cr,psi5,crit="qf")

psi6 <- bs(hr,df=6,degree=1,intercept=FALSE)
put61 <- funcScanone(prob,cr,psi6,crit="ss")
put62 <- funcScanone(prob,cr,psi6,crit="qf")

psi7 <- bs(hr,df=7,degree=1,intercept=FALSE)
put71 <- funcScanone(prob,cr,psi7,crit="ss")
put72 <- funcScanone(prob,cr,psi7,crit="qf")

psi8 <- bs(hr,df=8,degree=1,intercept=FALSE)
put81 <- funcScanone(prob,cr,psi8,crit="ss")
put82 <- funcScanone(prob,cr,psi8,crit="qf")

psi10 <- bs(hr,df=10,degree=1,intercept=FALSE)
put101 <- funcScanone(prob,cr,psi10,crit="ss")
put102 <- funcScanone(prob,cr,psi10,crit="qf")

cv3 <- bsCV(prob,z=NULL,hr,df=3:150,nfold=10)
plot(3:52,cv3[1:50,1],ylim=c(40,80))
lines(cv3[1:50,1]-cv3[1:50,2])
lines(cv3[1:50,1]+cv3[1:50,2])

cv1 <- bsCV(prob,z=NULL,hr,df=3:150,nfold=10,degree=1)
plot(3:52,cv1[1:50,1],ylim=c(40,80))
lines(cv1[1:50,1]-cv1[1:50,2])
lines(cv1[1:50,1]+cv1[1:50,2])

which.min(cv[,1])
cv[115,1]+cv[115,2]
cv[,1]>cv[115,1]+cv[115,2]
  

yh16chr1 <- yPred(prob,cr$geno[[1]]$prob[,31,1],phi16,
                  matrix(c(0,1),ncol=1),phi16)
yh16chr4 <- yPred(prob,cr$geno[[4]]$prob[,32,1],phi16,
                  matrix(c(0,1),ncol=1),phi16)
yh16chr9 <- yPred(prob,cr$geno[[9]]$prob[,7,1],phi16,
                  matrix(c(0,1),ncol=1),phi16)
yh16chr10 <- yPred(prob,cr$geno[[10]]$prob[,14,1],phi16,
                  matrix(c(0,1),ncol=1),phi16)


plot(hr,yh16chr1[1,],ylim=c(0,1),type="n")
lines(hr,yh16chr1[1,])
lines(hr,yh16chr1[1,],col="red")
lines(hr,yh16chr1[2,],col="blue")

plot(hr,yh16chr4[1,],ylim=c(0,1),type="n")
lines(hr,yh16chr4[1,])
lines(hr,yh16chr4[1,],col="red")
lines(hr,yh16chr4[2,],col="blue")

plot(hr,yh16chr9[1,],ylim=c(0,1),type="n")
lines(hr,yh16chr9[1,])
lines(hr,yh16chr9[1,],col="red")
lines(hr,yh16chr9[2,],col="blue")

plot(hr,yh16chr10[1,],ylim=c(0,1),type="n")
lines(hr,yh16chr10[1,])
lines(hr,yh16chr10[1,],col="red")
lines(hr,yh16chr10[2,],col="blue")


qtl4 <- cbind(cr$geno[[1]]$prob[,31,1],cr$geno[[4]]$prob[,32,1],
              cr$geno[[9]]$prob[,7,1],cr$geno[[10]]$prob[,14,1])

yh20qtl4 <- devSS(prob,qtl4,phi20)
yh20chr1 <- yPred(prob,qtl4[,1],phi20,
                  matrix(c(0,1),ncol=1),phi20)
yh20chr4 <- yPred(prob,qtl4[,2],phi20,
                  matrix(c(0,1),ncol=1),phi20)
yh20chr9 <- yPred(prob,qtl4[,3],phi20,
                  matrix(c(0,1),ncol=1),phi20)

qtl42 <- cbind(cr$geno[[4]]$prob[,9,1],cr$geno[[4]]$prob[,29,1])
yh20qtl42 <- devSS(prob,qtl42,phi20)

plot(hr,yh20chr1[1,],ylim=c(0,1),type="n")
lines(hr,yh20chr1[1,])
lines(hr,yh20chr1[1,],col="red")
lines(hr,yh20chr1[2,],col="blue")

plot(hr,yh20chr4[1,],ylim=c(0,1),type="n")
lines(hr,yh20chr4[1,])
lines(hr,yh20chr4[1,],col="red")
lines(hr,yh20chr4[2,],col="blue")

plot(hr,yh20chr9[1,],ylim=c(0,1),type="n")
lines(hr,yh20chr9[1,])
lines(hr,yh20chr9[1,],col="red")
lines(hr,yh20chr9[2,],col="blue")

