# functions for implementing warping ideas

# x = point where function is evaluated
# x0 = x coordinates of knots
# y0 = y coordinates of warped function
# ends = interval of warping function

warpFun <- function(x,x0,y0,ends=c(0,1))
  {
    x0 <- c(ends[1],x0,ends[2])
    y0 <- c(ends[1],y0,ends[2])
    predict(interpSpline(x0,y0),x)
  }

warpFun1 <- function(y0,x,x0,ends=c(0,1))
  {
    x0 <- c(ends[1],x0,ends[2])
    y0 <- c(ends[1],y0,ends[2])
    predict(interpSpline(x0,y0),x)$y
  }


# x = point where function is evaluated
# x0 = x coordinates of knots
# y0 = y coordinates of amplitude function

amplFun <- function(x,x0,y0)
  {
    predict(interpSpline(x0,y0),x)$y
  }


# y = trait
# z = covariate (categorical)
# x = time points
# x0 = warp knots

fitWarp <- function(y,z,x,x0)
  {

  }

devWarp0 <- function(y)
  {
    n <- nrow(y)
    m <- ncol(y)
    y0 <- apply(y,2,mean)
    y0 <- matrix(rep(y0,n),byrow=T,nrow=n)
    sum((y-y0)^2)
  }

# function for vertical warping
# y = data matrix
# x = time co-ordinates
# z = covariates
# v1 = y values corresponding to knots v0
# beta0 = base function

devWarp1 <- function(v1,y,z,x,v0,beta0)
  {
    n <- nrow(y)
    
    ampl <- amplFun(x,v0,v1[-1])
    
    # z <- matrix(z,ncol=1)
    
    y0 <- predict(beta0,x=x)$y
    y0 <- outer(rep(v1[1],n),y0)
    
    y1 <- outer(z,ampl)
    
    # print(dim(y1))
    
    yhat <- y0+y1
    ss1 <- sum((y-yhat)^2)

    ss1
  }


# horizontal warping
# y = data matrix
# z = covariate
# x = time points
# w0 = time interior warp knots
# w1 =
# beta0 = base function

devWarp2 <- function(w1,y,z,x,w0,ends=c(0,1),beta0)
  {
    n <- nrow(y)

    # warp deviations at knots
    ww1 <- outer(w1,z)
    # warp base
    ww0 <- outer(w0,rep(1,n))

    # warped times at knots
    ww2 <- ww1+ww0
    
    # all warped times
    warp <- apply(ww2,2,warpFun1,x=x,x0=w0,ends=ends)
    
    # predict values at warped times
    yhat <- predict(beta0,x=c(warp))$y
    yhat <- matrix(yhat,nrow=n,byrow=T)
    
    # yhat <- y0+y1
    ss1 <- sum((y-yhat)^2)

    #list(ss0,ss1)
    ss1
  }

###########################################################################

# horizontal and vertical warping
# y = data matrix
# z = covariate
# x = time points
# w0 = time interior warp knots
# w1 =
# beta0 = base function

devWarp3 <- function(w1,y,z,x,w0,ends=c(0,1),beta0)
  {
    # number of rows in data matrix
    n <- nrow(y)

    # warp deviations at knots
    ww1 <- outer(w1,z)
    # warp base
    ww0 <- outer(w0,rep(1,n))

    # warped times at knots
    ww2 <- ww1+ww0
    
    # all warped times
    warp <- apply(ww2,2,warpFun1,x=x,x0=w0,ends=ends)
    
    # predict values at warped times
    yhat <- predict(beta0,x=c(warp))$y
    yhat <- matrix(yhat,nrow=n,byrow=T)
    
    # yhat <- y0+y1
    ss1 <- sum((y-yhat)^2)

    #list(ss0,ss1)
    ss1
  }



# date()
# optim(rep(0,20),devWarp1,y=prob,z=g,x=hr,v0=seq(0,22,length=20),
#      method="BFGS")
# date()
## g1 <- cr$geno[[1]]$data[,16]
## g4 <- cr$geno[[4]]$data[,15]
## g9 <- cr$geno[[9]]$data[,5]
 
