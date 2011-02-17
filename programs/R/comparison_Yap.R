library(qtl)
library(ggplot2)
library(multicore)
source('logistic.R')
source('fr.R')

beta.coef <- matrix(c(30, 5.0, 0.5, 28.5, 5.0, 0.5, 27.5, 5., 0.5), byrow=T, nrow=3)
tt <- 0.:9
len.tt <- 10                                 # time series length

structured.cov <- matrix(
c(0.72, 0.39, 0.45, 0.48, 0.50, 0.53, 0.60, 0.64, 0.68, 0.68,
  0.39, 1.06, 1.61, 1.60, 1.50, 1.48, 1.55, 1.47, 1.35, 1.29,
  0.45, 1.61, 3.29, 3.29, 3.17, 3.09, 3.19, 3.04, 2.78, 2.53,
  0.48, 1.60, 3.29, 3.98, 4.07, 4.01, 4.17, 4.18, 4.00, 3.69,
  0.50, 1.50, 3.17, 4.07, 4.70, 4.68, 4.66, 4.78, 4.70, 4.36,
  0.53, 1.48, 3.09, 4.07, 4.68, 5.56, 6.23, 6.87, 7.11, 6.92,
  0.60, 1.55, 3.19, 4.17, 4.66, 6.23, 8.59, 10.16,10.80,10.70,
  0.64, 1.47, 3.04, 4.18, 4.78, 6.87, 10.16,12.74,13.80,13.80,
  0.68, 1.35, 2.78, 4.00, 4.70, 7.11, 10.80,13.80,15.33,15.35,
  0.68, 1.29, 2.53, 3.69, 4.36, 6.92, 10.70,13.80,15.35,15.77), byrow=T, nrow=10)

## generate both genotypes and phenotypes
gen.data <- function(sample.size, cov.fcn){
  ##popu.size <- 10000  # population size
  
  ## simulate genotypes
  mp <- sim.map(100, n.mar=6, include.x=F, eq.spacing=T) # simulate map
  md <- c(1,32,0,0)                       # one QTL at 32cM on chrom. 1
  samples <- sim.cross(map=mp, model=md, type='f2', n.ind=sample.size, keep.qtlgeno=T)

  ## retrieve sample genotypes
  ## ind <- sample(popu.size, sample.size)
  ## samples <- subset.cross(cross, ind=ind)
  ## ## subsetting qtlgeno by hand since subset.cross does not do it.
  ## qtlgeno <- samples$qtlgeno[ind]
  ## samples$qtlgeno <- qtlgeno

  ## simulate phenotypes
  ## 1. get the means
  mean.vals <- matrix(0., nrow=3, ncol=len.tt) # only 3 genotypes means only 3 mean curves
  for (i in 1:3){
    mean.vals[i,] <- logisticFun(beta.coef[i,],tt)
  }
  sample.means <- mean.vals[samples$qtlgeno,]
  ## 2. get the noises
  if (cov.fcn == 'autocorr'){           # case (1)
    ee <- rnormAutocor(sample.size, tt, 0.6, 3.0) # rho = 0.6, sigma^2 = 3
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

## plot a set of time series as rows of a matrix.
## time is the sequence of time points
## ID is the sample IDs of the time series;
## color.code is the code that distinguish groups among samples
plot.ts <- function(ts, time, ID, color.code,...){
  stopifnot(nrow(ts) == length(ID), length(ID) == length(color.code))
  stopifnot(ncol(ts) == length(time))

  header <- paste(ID, color.code, sep=' ')
  tmp <- data.frame(cbind(time,t(ts)), row.names=time)
  colnames(tmp) <- c('time', header)
  temp <- melt(tmp, id.vars='time')
  temp <- cbind(temp, colsplit(temp$variable, ' ', c('ID', 'color.code')))

  return(geom_line(data=temp, aes(time, value,group=ID, color=factor(color.code)),...))
}

## plot a covariance matrix as a heatmap
## cov.mat is the covariance matrix to be plotted
## time is the time points of which we assume the covariance matrix comes from
## returns a heatmap layer
plot.hm <- function(cov.mat, time){
  tmp <- data.frame(cbind(1:10, cov.mat))
  colnames(tmp) <- c('time', paste(time))
  tmp <- melt(tmp, id.vars='time')
  hmp <- geom_tile(data=tmp, mapping=aes(time, variable, fill=value)) 
  return(hmp)
}
## generate a series of diagonistic plots
diagnostics <- function(sample.size=200, file.name='diagnostics.pdf'){
  diag.one <- function(cov.fcn){
    samples <- gen.data(sample.size, cov.fcn)
    Y <- samples$pheno
    X <- samples$qtlgeno

    ## compare mean curve against true values
    lrg <- lm(Y~factor(X)-1)
    tmp <- fitted(lrg)
    sam.ave <- rbind(tmp[which(X==1)[1],], tmp[which(X==2)[1],], tmp[which(X==3)[1],])
    fig.sam <- plot.ts(sam.ave, tt, 1:3, 1:3, linetype=1)
    mean.vals <- matrix(0., nrow=3, ncol=len.tt) # only 3 genotypes means only 3 mean curves
    for (i in 1:3){
      mean.vals[i,] <- logisticFun(beta.coef[i,],tt)
    }
    fig.true <- plot.ts(mean.vals, tt, 1:3, 1:3, linetype=2)
    mean.curves <- ggplot()+fig.sam+fig.true+xlab('time')+ylab('phenotype')+opts(title=paste('mean curves with', cov.fcn))

    ## plot all samples along with their mean curves
    all.curves <- ggplot() + fig.sam + plot.ts(Y, tt, 1:sample.size, X, alpha=I(1/3.)) + xlab('time')+ylab('phenotype')+opts(title=paste('all curves with', cov.fcn))

    ## look at the heatmap of residual covariance
    hm <- ggplot() + plot.hm(cov(resid(lrg)), 1:10) + xlab('time') + ylab('time') +opts(title=cov.fcn)
    
    return(list(mean=mean.curves, all=all.curves, cov=hm))
  }
  
  pdf(file=file.name)
  autocorr <- diag.one('autocorr')
  equicorr <- diag.one('equicorr')
  structured <- diag.one('structured')
  tmp <- list(autocorr, equicorr, structured)
  ## plot the figures
  ## 1. mean curves
  lapply(tmp, function(x) print(x$mean));
  ## 2. all curves
  lapply(tmp, function(x) print(x$all));
  ## 3. covariance heatmap
  lapply(tmp, function(x) print(x$cov))
  ## print(autocorr$cov)
  ## print(equicorr$cov)
  ## print(structured$cov)
  dev.off()
}

## One simulation run
one.sim <- function(sample.size=100, cov.fcn='autocorr', basis.fcn, meth){
  ## generate data and calculate the genotype probability
  samples <- gen.data(sample.size, cov.fcn)
  Y <- samples$pheno
  tmp <- calc.genoprob(samples,  step=4)
  
  ## functional regression
  
  res <- funcScanone(Y, tmp, basis.fcn, crit=meth)
  return(res$pos[which.max(res$lod)])
}

## simulations
simulate.multi <- function(nrun=1000, sample.size=100, cov.fcn='structured', crit='ss'){

  res <- rep(0., nrun)
  basis <- bs(tt, df=4, intercept=FALSE)#TRUE)

  for (i in 1:nrun){
    res[i] <- one.sim(sample.size, cov.fcn, basis, crit)
  }

  return(res)
}

## compute various metrics
calc.metrics <- function(vec){
  return( list(mean(vec), sd(vec), sqrt((mean(vec)-32.)^2+sd(vec)^2)) )
}

param.df <- data.frame( list(sample.size=rep(c(100, 400), 3), cov.fcn=c('autocorr', 'autocorr', 'equicorr', 'equicorr', 'structured', 'structured')) )
#res.ss <- mclapply(1:6, function(i, params) simulate.multi(10000, params[i,1], params[i,2], 'ss'), param.df, mc.cores=4)
#res.qf <- mclapply(1:6, function(i, params) simulate.multi(10000, params[i,1], params[i,2], 'qf'), param.df, mc.cores=4)
