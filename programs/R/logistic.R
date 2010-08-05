library(geoR)

maternCov <- function(u,phi,kappa)
  {
    toeplitz(matern(u,phi,kappa))
  }

rnormMatern <- function(n,u,phi,kappa)
  {
    ee <- rnorm(n*length(u))
    A <- Matrix::chol(maternCov(u,phi,kappa))
    ee <- matrix(ee,nrow=n)%*%A
    ee
  }

logisticFun <- function(beta,tt)
  {
    # a/(1+b*exp(-c*t))
    # initial value = a/(1+b)
    # asymptotic limit = a
    # rate of growth = c
    beta[1] / (1+ beta[2]*exp(-beta[3]*tt))
  }

logisticSim1 <- function(beta,tt,n,df=0)
  {
    # autocorrelation
    rr <- beta[4]
    # residual variance
    ss <- beta[5]

    m <- length(tt)

    ww <- rr/sqrt(1-rr^2)
    
    mm <- logisticFun(beta[1:3],tt)
    mm <-  kronecker( matrix(mm,nrow=1), matrix(rep(1,n),ncol=1) )
    if(df==0)
      {
        ee <- matrix(rnorm(n*m),nrow=n,ncol=m)
      }
    else
      {
        scaleFac <- n/(n-2)
        ee <- matrix(rnorm(n*m),nrow=n,ncol=m)/scaleFac
      }
    for( i in 2:m )
      {
        ee[,i] <- (ww*ee[,i-1] + ee[,i])/sqrt(1+ww^2)
      }
   mm + sqrt(ss)*ee
  }


logisticSim1Mat <- function(beta,tt,n,phi,kappa)
  {
    # residual variance
    ss <- beta[5]

    m <- length(tt)

    mm <- logisticFun(beta[1:3],tt)
    mm <-  kronecker( matrix(mm,nrow=1), matrix(rep(1,n),ncol=1) )


    ee <- rnormMatern(n,tt,phi,kappa)
    
    mm + sqrt(ss)*ee
  }


logisticSim2 <- function(beta,tt,grp,df=0)
  {
    ngrp <- length(table(grp))
    grp <- as.numeric(as.factor(grp))-1

    y <- matrix(nrow=length(grp),ncol=length(tt))
    
    if( nrow(beta)!= ngrp )
      stop("Number of groups does not match the dimension of parameters.")
    else
      {
        for( i in 1:ngrp )
          {
            idx <- which(grp==(i-1))
            y[idx,] <- logisticSim1(beta[i,],tt,length(idx),df=df)
          }

      }
    y
  }

logisticSim2Mat <- function(beta,tt,grp,phi,kappa)
  {
    ngrp <- length(table(grp))
    grp <- as.numeric(as.factor(grp))-1

    y <- matrix(nrow=length(grp),ncol=length(tt))
    
    if( nrow(beta)!= ngrp )
      stop("Number of groups does not match the dimension of parameters.")
    else
      {
        for( i in 1:ngrp )
          {
            idx <- which(grp==(i-1))
            y[idx,] <- logisticSim1Mat(beta[i,],tt,length(idx),
                                       phi=phi,kappa=kappa)
          }

      }
    y
  }




logisticLik <- function(beta,y,tt)
  {
    rr <- beta[4]
    ss <- beta[5]

    n <- nrow(y)
    m <- ncol(y)

    # mean
    mm <- logisticFun(beta[1:3],tt)

    # residuals
    ee <- y - kronecker( matrix(mm,nrow=1), matrix(rep(1,n),ncol=1) )

    exponent <- sum((ee[,-m]-rr*ee[,-1])^2)/(1-rr^2) + sum(ee[,m]^2)
    exponent <- (-1/2)*exponent/ss
 
    loglik <- exponent - m*n*(1/2)*log(2*pi*ss) - n*(m-1)*(1/2)*log(1-rr^2)

    loglik
  }


logisticLik1 <- function(beta,y,tt)
  {
    beta1 <- c(exp(beta[1:3]),tanh(beta[4]),exp(beta[5]))
    -logisticLik(beta1,y,tt)
  }

logisticMLE1 <- function(y,tt,beta0=c(0,0,0,0,0),loglik=FALSE,...)
  {
    out <- optim(beta0,logisticLik1,y=y,tt=tt,control=list(fnscale=-1),
                 method="BFGS")
    betahat <- c(exp(out$par[1:3]),tanh(out$par[4]),exp(out$par[5]))
    if(loglik)
      {
        loglik <- out$value
        val <- list(betahat,loglik)
      }
    else
      val <- betahat

    val
  }



logisticLik2 <- function(beta,y,tt,grp)
  {
    grp <- as.numeric(as.factor(grp))-1
    ngrp <- max(grp)+1

    if(length(beta)!=ngrp*3+2)
      stop("Wrong parameter vector length.")
    else
      {
        beta1 <- matrix(nrow=ngrp,ncol=5)
        beta1[,1:3] <- matrix(beta[1:(ngrp*3)],nrow=ngrp,byrow=T)
        beta1[,4:5] <- matrix(rep(beta[(ngrp*3+1):(ngrp*3+2)],ngrp),
                             nrow=ngrp,byrow=T)
        beta1[,c(1,2,3,5)] <- exp(beta1[,c(1,2,3,5)])
        beta1[,4] <- tanh(beta1[,4])
      }
    # print(beta1)
    beta <- beta1
    loglik <- 0
    for( i in 1:ngrp )
      {
        idx <- which(grp==(i-1))
        loglik <- loglik + logisticLik(beta[i,],y[idx,],tt)
        # print(logisticLik(beta[i,],y[idx,],tt))
      }
    -loglik
  }


logisticMLE2 <- function(y,tt,grp,beta0=NULL,
                         loglik=FALSE,...)
  {
    ngrp <- length(unique(grp))

    if(is.null(beta0))
      {
        aa <- max(y)
        if(min(y)<=0)
          bb <- mean(y[,1])
        else
          bb <- (aa-min(y))/min(y)
        cc <- 1
        rr <- 0
        ss <- mean(diag(var(y)))
        beta0 <- c(log(aa),log(bb),log(cc),atanh(rr),log(ss))
        beta0 <- c(rep(beta0[1:3],ngrp),beta0[4:5])
      }

    # out <- optim(beta0,logisticLik2,y=y,tt=tt,grp=grp,...)
    # out <- nlm(logisticLik2,beta0,y=y,tt=tt,grp=grp,...)
    out <- nlminb(beta0,logisticLik2,y=y,tt=tt,grp=grp,...)
    # out$par <- out$estimate
    # out$convergence <-  as.numeric(out$code>2)
    # out$value <- out$minimum
    out$value <-  out$objective

    betahat <- matrix(nrow=ngrp,ncol=5)
    betahat[,1:3] <- matrix(out$par[1:(ngrp*3)],ncol=3,byrow=T)
    betahat[,4:5] <- matrix(rep(out$par[(ngrp*3+1):(ngrp*3+2)],ngrp),
                            ncol=2,byrow=T)
    betahat[,c(1,2,3,5)] <- exp(betahat[,c(1,2,3,5)])
    betahat[,4] <- tanh(betahat[,4])
    if(loglik)
      {
        loglik <- out$value
        val <- list(betahat=betahat,loglik=loglik,convergence=out$convergence)
      }
    else
      val <- list(betahat=betahat,convergence=out$convergence)

    val
  }

logisticTest2 <- function(y,tt,grp,beta0=NULL,...)
  {
    out1 <- logisticMLE2(y,tt,grp,loglik=TRUE,beta0=beta0,...)

    out0 <- logisticMLE2(y,tt,grp=rep(0,length(grp)),loglik=TRUE,
                         beta0=apply(out1$betahat,2,mean),...)
 
    if((out0$convergence!=0)|(out1$convergence!=0))
      stop("Did not converge.")
    else
      {
        out <- list(beta0=out0$betahat,
                    beta1=out1$betahat,
                    loglik= out0$loglik - out1$loglik)
      }
    out
  }


# simulate data with autorgressive errors
compareSimAR <- function(beta,nsim,tt,psi,grp,df=0)
  {
    
    pval <- matrix(nrow=nsim,ncol=3)
    ngrp <- length(unique(grp))
    psidf <- ncol(model.matrix(~psi))*(ngrp-1)
          
    for( i in 1:nsim )
      {
        print(i)
        y <- logisticSim2(beta,tt,grp,df=df)
        out1 <- logisticTest2(y,tt,grp,
                              beta0=c(beta[1,1:3],beta[2,1:3],beta[1,4:5]))
        out2 <- quadForm(y,grp,psi,shrink=TRUE)
        out3 <- quadForm(y,grp,psi,shrink=FALSE)

        pval[i,1] <- pchisq(2*out1$loglik,3,lower=FALSE)
        pval[i,2] <- pchisq(out2$quadform,psidf,lower=FALSE)
        pval[i,3] <- pchisq(out3$quadform,psidf,lower=FALSE)        
      }
    pval
}

# simulate data with Matern errors
compareSimMat <- function(beta,nsim,tt,psi,grp,phi,kappa)
  {
    
    pval <- matrix(nrow=nsim,ncol=3)
    ngrp <- length(unique(grp))
    psidf <- ncol(model.matrix(~psi))*(ngrp-1)
    
    for( i in 1:nsim )
      {
        print(i)
        y <- logisticSim2Mat(beta,tt,grp,phi=phi,kappa=kappa)
        out1 <- logisticTest2(y,tt,grp,
                              beta0=c(beta[1,1:3],beta[2,1:3],beta[1,4:5]))
        out2 <- quadForm(y,grp,psi,shrink=TRUE)        
        out3 <- quadForm(y,grp,psi,shrink=FALSE)
        
        pval[i,1] <- pchisq(2*out1$loglik,3,lower=FALSE)
        pval[i,2] <- pchisq(out2$quadform,psidf,lower=FALSE)
        pval[i,3] <- pchisq(out3$quadform,psidf,lower=FALSE)        
      }
    pval
}

# kolmogorov-smirnov test applied to each column of p-values
ksTest <- function(x)
{
  z <- vector(length=3)
  for( i in 1:3 )
    z[i] <- ks.test(x[,i],"punif")$p.value
  z
}

powerCalc <- function(x0,x1,alpha=0.05)
{
  z <- vector(length=3)
  for( i in 1:3 )
    z[i] <- mean(x1[,i]<quantile(x0[,i],alpha))
  z
}
