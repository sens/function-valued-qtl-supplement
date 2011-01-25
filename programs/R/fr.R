#########################################
# Routines for functional regression
#########################################

# y = data matrix
# z = covariate matrix
# phi = smoothing basis

library(qtl)
library(splines)
library(Matrix)
library(corpcor)

varEst <- function(y)
  {
    var(t(y))
  }

bHat <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL)
  {
    ty <- t(y)
    if(is.null(phi))
      out1 <- lm(ty~1,weights=weightPhi)
    else
      {
        if(addPhiIntercept)
          out1 <- lm(ty~phi,weights=weightPhi)
        else
          out1 <- lm(ty~phi-1,weights=weightPhi)
      }
    tb <- coef(out1)
    ## print(dim(t(tb)))
    ## print(tb)
    if(is.null(z))
      out2 <- lm(t(tb)~1,weights=weightZ)
    else
      {
        if(addZIntercept)
          out2 <- lm(t(tb)~z,weights=weightZ)
        else
          out2 <- lm(t(tb)~z-1,weights=weightZ)
      }
    b <- coef(out2)
    b
  }

varBHat <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                    shrink=TRUE)
  {
    # estimate yhat
    yhat <- yHat(y,z,phi,addPhiIntercept=addPhiIntercept,
                 addZIntercept=addZIntercept)

    if(is.null(z))
      z <- rep(1,length=nrow(y))
    else
      {
        if(addZIntercept)
          z <- model.matrix(~z)
        else
          z <- model.matrix(~z-1)
      }

    # degrees of freedom
    zdf <- ncol(z)
    # sample size
    n <- nrow(y)
    
    # estimate within sample variance matrix, sigma
    if(shrink)
      {
        sigma <- cov.shrink(y-yhat,verb=FALSE)
      }
    else
      {
        sigma <- var(y-yhat)*(n-1)/(n-zdf)
      }
    # project sigma on phi
    pp <- bHat(sigma,phi,phi,addZIntercept=addPhiIntercept,
               addPhiIntercept=addPhiIntercept)

    # calculate (z'z)^(-1)
    zzinv <- solve(t(z)%*%z)
    # take kron prod
    kronecker(zzinv,pp)
  }

quadForm <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                     notIntercept=TRUE,weightPhi=NULL,shrink=TRUE)
  {
    # estimate yhat
    yhat <- yHat(y,z,phi,addPhiIntercept=addPhiIntercept,
                 addZIntercept=addZIntercept,weightPhi=weightPhi)
    # estimate within sample variance matrix, sigma
    # sigma <- var(y-yhat)
    # sigma <- cov.shrink(y,verb=FALSE)
    n <- nrow(y)
    if(shrink)
      {
        sigma <- cov.shrink(y-yhat,verb=FALSE)
      }
    else
      {
        if(addZIntercept)
          zdf <- ncol(model.matrix(~z))
        else
          zdf <- ncol(model.matrix(~z-1))
        sigma <- var(y-yhat)*(n-1)/(n-zdf)
      }
    
    # sigma <- var(y)
    # calculate bhat
    b <- bHat(y,z,phi,addPhiIntercept=addPhiIntercept,
                 addZIntercept=addZIntercept,weightPhi=weightPhi)
    # project sigma on phi
    pp <- bHat(sigma,phi,phi,addZIntercept=addPhiIntercept,
               addPhiIntercept=addPhiIntercept,weightPhi=weightPhi,
               weightZ=weightPhi)
    if(is.null(z))
      z <- rep(1,length=nrow(y))
    else
      {
        if(addZIntercept)
          z <- model.matrix(~z)
        else
          z <- model.matrix(~z-1)
      }
    # calculate (z'z)^(-1)
    zzinv <- solve(t(z)%*%z)
    # take kron prod
    vbhat <- kronecker(zzinv,pp)
    if(notIntercept)
      {
        bb <- c(t(b[-1,]))
        idx <- 1:length(b[1,])
        qq <- solve(vbhat[-idx,-idx],bb)
        qq <- sum(bb*qq)
      }
    else
      {
        bb <- c(t(b))
        qq <- solve(vbhat,bb)
        qq <- sum(bb*qq)
      }
    list(quadform=qq,vbhat=vbhat,b=bb)

  }

yHat <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL)
  {
    ty <- t(y)
    if(is.null(phi))
      {
        phi <- rep(1,length=ncol(y))
        out1 <- lm(ty~phi,weights=weightPhi)
      }
    else
      {
        if(addPhiIntercept)
          out1 <- lm(ty~phi,weights=weightPhi)
        else
          out1 <- lm(ty~phi-1,weights=weightPhi)
      }
    tyhat <- fitted(out1)
      
    if(is.null(z))
      {
        z <- rep(1,length=nrow(y))
        out2 <- lm(t(tyhat)~z)
      }
    else
      {
        if(addZIntercept)
          out2 <- lm(t(tyhat)~z)
        else
          out2 <- lm(t(tyhat)~z-1)
        }
    yhat <- fitted(out2)
    yhat
  }

devSS <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                  weightPhi=NULL)
{
  yh <- yHat(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
             weightPhi=weightPhi)
  if(!is.null(weightPhi))
    sum((y-yh)^2%*%diag(weightPhi))
  else
    sum((y-yh)^2)
}


funcFit <- function(y,z,phi,type="ss",calcNull=FALSE,addPhiIntercept=TRUE,
                    addZIntercept=TRUE)
  {
    yhat1 <- yHat(y,z,phi)
    ss1 <- sum((y-yhat1)^2)
    if(calcNull)
      {
        yhat0 <- yHat(y,z=NULL,phi)
        ss0 <- sum((y-yhat0)^2)
      }
    else
      {
        ss0 <- NULL
      }

    if(type=="all")
      {
        b <- bHat(y,z,phi)
        Vbhat <- varBHat(y,z,phi)
        out <- list(b=b,yhat=yhat,Vbhat=Vbhat,ss1=ss1,ss0=ss0)
      }
    else
      {
        out <- list(ss1=ss1,ss0=ss0)
      }
    out
  }

## funcScanone <- function(y,geno,phi)
##   {
##     fit0 <- funcFit(y,geno[,1],phi,calcNull=TRUE)
##     out <- rep(fit0$ss0,ncol(geno))
##     for( i in 1:ncol(geno))
##       {
##         # print(i)
##         fit1 <- funcFit(y,geno[,i],phi)
##         out[i] <- log(fit0$ss0) - log(fit1$ss1)
##       }
##     out
##   }

 funcScanone <- function(y,cr,phi,method="hk",crit="qf",weightPhi=NULL,
                         shrink=TRUE)
  {
    if(method=="hk")
      out <- funcScanone.hk(y,cr,phi,crit=crit,weightPhi=weightPhi,
                            shrink=shrink)
    else
      if(method=="imp")
        out <- funcScanone.imp(y,cr,phi,crit=crit,weightPhi=weightPhi)
      else
        error("Unknown method.")
    out
  }

 funcScanone.hk <- function(y,cr,phi,crit="qf",weightPhi=NULL,
                            shrink=TRUE)
   {
     if(!match("prob",names(cr$geno[[1]])))
       {
         warning("First running calc.genoprob.")
         cr <- calc.genoprob(cr)
       }
     # get the array ready
     out <- scanone(cr,pheno.col=y[,1])
     npseudo <- 0
     nr <- nrow(out)
     gg <- array(dim=c(nrow(y),nr,dim(cr$geno[[1]]$prob)[3]))
     
     for( i in 1:length(cr$geno) )
       {
         # all the genotype probabilities on the i-th chromosome
         ggChr <- cr$geno[[i]]$prob
         npseudoChr <- dim(ggChr)[2]
         # print(dim(ggChr))
         # print(c(npseudo,npseudoChr))
         gg[,(npseudo+1):(npseudo+npseudoChr),] <- ggChr
         npseudo <- npseudo+npseudoChr
       }

     if( crit=="qf" )
       {
         for( j in 1:nr )
           {
             qq <- quadForm(y,gg[,j,-1],phi,weightPhi=weightPhi,
                            shrink=shrink)
             out[j,3] <- qq$quadform/2/log(10)
           }
       }

     if( crit=="ss" )
       {
         ss0 <- devSS(y,NULL,phi,weightPhi=weightPhi)
         for( j in 1:nr )
           {
            ss  <- devSS(y,gg[,j,-1],phi,weightPhi=weightPhi)
            out[j,3] <- -log(ss/ss0)
           }
       }
     out
   }

 funcScanone.imp <- function(y,cr,phi)
   {
     if(!match("draws",names(cr$geno[[1]])))
       {
         warning("First running sim.geno.")
         cr <- sim.geno(cr)
       }
     # get the array ready
     out <- scanone(cr,pheno.col=y[,1],method="imp")
     
     npseudo <- 0
     for( i in 1:length(cr$geno) )
       {
         # all the imputations on the i-th chromosome
         gg <- cr$geno[[i]]$draws
         npseudoChr <- dim(gg)[2]
         for( j in 1:npseudoChr )
           {
             qq <- quadForm(y,gg[,j],phi)
             out[npseudo+j,3] <- qq$quadform/2/log(10)
           }
         npseudo <- npseudo+npseudoChr
       }
     out
   }


funcScanonePerm <- function(y,cr,phi,nperm,method="hk",crit="qf",
                            weightPhi=NULL,shrink=TRUE)
  {
    # if genotype probabilities have not been calculated, calculate them
     if(!match("prob",names(cr$geno[[1]])))
       {
         warning("First running calc.genoprob.")
         cr <- calc.genoprob(cr)
       }
     # get the array ready
     out <- scanone(cr,pheno.col=y[,1])
     out[,3:(2+nperm)] <- matrix(rep(0,nperm*nrow(out)),ncol=nperm)

     # number of individuals
     n <- nrow(y)
     m <- ncol(y)
     
     # ready genotype matrix
     npseudo <- 0
     # number pseudomarkers
     nr <- nrow(out)
     gg <- array(dim=c(nrow(y),nr,dim(cr$geno[[1]]$prob)[3]))
     
     for( i in 1:length(cr$geno) )
       {
         # all the genotype probabilities on the i-th chromosome
         ggChr <- cr$geno[[i]]$prob
         npseudoChr <- dim(ggChr)[2]
         # print(dim(ggChr))
         # print(c(npseudo,npseudoChr))
         gg[,(npseudo+1):(npseudo+npseudoChr),] <- ggChr
         npseudo <- npseudo+npseudoChr
       }


     if( crit=="ss" )
       {
         ss0 <- devSS(y,NULL,phi,weightPhi=weightPhi)
       }
     
     # perform permutations
     for( k in 1:nperm )
       {
         print(k)
         npseudo <- 0
         permidx <- sample(n)
         if( crit=="qf" )
           {
             for( j in 1:nr )
               {
                 qq <- quadForm(y[permidx,],gg[,j,-1],phi,weightPhi=weightPhi,
                                shrink=shrink)
                 out[j,2+k] <- qq$quadform/2/log(10)
           }
           }
         
         if( crit=="ss" )
           {
             for( j in 1:nr )
               {
                 ss  <- devSS(y[permidx,],gg[,j,-1],phi,weightPhi=weightPhi)
                 out[j,2+k] <- -log(ss/ss0)
               }
           }
       }
     # to conform to R/qtl scanone.perm
     out <- apply(out[,-(1:2)],2,max)
     out <- matrix(out,ncol=1)
     class(out) <- c("scanoneperm", "matrix")
     out
  }


###########################################################
# routine for cross validated integrated sequared error
###########################################################

cvSS <- function(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL,nfolds=10)
  {
    # divide the data into nfolds
    n <- nrow(y)
    idx <- sample(n)
    folds <- ((0:(n-1))%%nfolds) + 1
    
    if(addPhiIntercept)
      phi1 <- model.matrix(~phi)
    if(addZIntercept)
      {
        if(is.null(z))
          {
            z1 <- model.matrix(~rep(1,n)-1)
          }
        else
          {
            z1 <- model.matrix(~z)
          }
      }

    ss <- vector(length=nfolds)
    z1 <- as.matrix(z1)
    
    for( i in 1:nfolds )
      {
        # select the i-the fold
        fidx <- idx[folds==i]
        bb <- bHat(y[-fidx,],z1[-fidx,],phi=phi,addPhiIntercept=addPhiIntercept,
                  addZIntercept=FALSE,weightPhi=weightPhi,
                  weightZ=weightZ[-fidx])

        yh <- z1[fidx,]%*%bb%*%t(phi1)

        if(!is.null(weightPhi))
          ss[i] <- sum((y[fidx,]-yh)^2%*%diag(weightPhi))
        else
          ss[i] <- sum((y[fidx,]-yh)^2)
      }
    ss
  }


bsCV <- function(y,z,x,df=3:20,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL,nfolds=10,degree=3,
                 intercept=FALSE)
{
  cv <- matrix(nrow=length(df),ncol=2)
  for( i in 1:nrow(cv) )
    {
      phi <- bs(x,df=df[i],degree=degree,intercept=intercept)
      out <- cvSS(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL,nfolds=10)
      cv[i,1] <- mean(out)
      cv[i,2] <- sd(out)
    }
  cv
}


yPred <- function(y,z,phi,znew,phinew,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL)
  {
    bb <- bHat(y,z,phi,addPhiIntercept=TRUE,addZIntercept=TRUE,
                 weightPhi=NULL,weightZ=NULL)
    if(addZIntercept)
      znew <- model.matrix(~znew)
    else
      znew <- model.matrix(~znew-1)

    if(addPhiIntercept)
      phinew <- model.matrix(~phinew)
    else
      phinew <- model.matrix(~phinew-1)

    print(dim(bb))
    
    yp <- znew %*% bb %*% t(phinew)

    yp
  }
  
