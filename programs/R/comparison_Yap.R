library(qtl)
source('logistic.R')

beta.coef <- matrix(c(30, 5.0, 0.5, 28.5, 5.0, 0.5, 27.5, 5., 0.5), byrow=T, nrow=3)
tt <- 1:10                              # PLACEHOLDER
len.tt <- 10                                 # time series length

structured.cov <- matrix( c(0.72, 0.39, 0.45, 0.48, 0.50, 0.53, 0.60, 0.64, 0.68, 0.68, 0.39, 1.06, 1.61, 1.60, 1.50, 1.48, 1.55, 1.47, 1.35, 1.29, 0.45, 1.61, 3.29, 3.29, 3.17, 3.09, 3.19, 3.04, 2.78, 2.53, 0.48, 1.60, 3.29, 3.98, 4.07, 4.01, 4.17, 4.18, 4.00, 3.69, 0.50, 1.50, 3.17, 4.07, 4.70, 4.68, 4.66, 4.78, 4.70, 4.36, 0.53, 1.48, 3.09, 4.07, 4.68, 5.56, 6.23, 6.87, 7.11, 6.92, 0.60, 1.55, 3.19, 4.17, 4.66, 6.23, 8.59, 10.16, 10.80, 10.70, 0.64, 1.47, 3.04, 4.18, 4.78, 6.87, 10.16, 12.74, 13.80, 13.80, 0.68, 1.35, 2.78, 4.00, 4.70, 7.11, 10.80, 13.80, 15.33, 15.35, 0.68, 1.29, 2.53, 3.69, 4.36, 6.92, 10.70, 13.80, 15.35, 15.77), byrow=T, nrow=10)

## generate data
gen.data <- function(sample.size, cov.fcn){
  popu.size <- 10000  # population size
  
  ## simulate genotypes
  mp <- sim.map(100, n.mar=6, include.x=F, eq.spacing=T) # simulate map
  md <- c(1,32,0,0)                       # one QTL at 32cM on chrom. 1
  cross <- sim.cross(map=mp, model=md, type='f2', n.ind=popu.size, keep.qtlgeno=T) # simulate 10,000 individuals

  ## retrieve sample genotypes
  ind <- sample(popu.size, sample.size)
  samples <- subset.cross(cross, ind=ind)
  ## subsetting qtlgeno by hand since subset.cross does not do it.
  qtlgeno <- samples$qtlgeno[ind]
  samples$qtlgeno <- qtlgeno

  ## simulate phenotypes
  ## 1. get the means
  mean.vals <- matrix(0., nrow=sample.size, ncol=len.tt)
  for (i in 1:3){
    mean.vals[i,] <- logisticFun(beta.coef[i,],tt)
  }
  sample.means <- mean.vals[samples$qtlgeno,]
  ## 2. get the noises
  if (cov.fcn == 'autocorr'){           # case (1)
    ee <- rnormAutocor(sample.size, tt, 0.6, 3.0) # rho = 0.6, sigma^2 = 3.
  }
  else if (cov.fcn == 'equicorr'){      # case (2)
    cov <- matrix(0.5, nrow=len.tt, ncol=len.tt)  # rho = 0.5
    diag(cov) <- 1.0
    cov <- 3.0 * cov                    # sigma^2 = 3.0
    ee <- rnormMulti(sample.size, cov)
  }
  else if (cov.fcn == 'structured'){    # case(3)
    ee <- rnormMulti(sample.size, structured.cov)    
  }
  else{
    stop('Unknown covariance function.')
  }
  Y <- sample.means + ee

  samples$pheno <- Y
  return(samples)
}



