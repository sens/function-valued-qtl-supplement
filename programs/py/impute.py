#!/usr/bin/env python

from glob import iglob
import csv
import numpy as np
from scipy.stats import tmean, tvar, scoreatpercentile

import rpy2.robjects as robj
from rpy2.robjects import r as R
import rpy2.robjects.numpy2ri

from FDA.bspline import Bspline
from FDA.fqtl import compute_qtl, add_ones
from FDA.raw_regress_qr import select, Regression

def impute(geno_fl, n=16):
    """
    Given a genotype file impute n times missing genotypes.

    The output file names are new_geno_i.csv, where i is the count.
    """
    R("""library('qtl')
    library('abind')
    geno <- read.cross(file='%s', format='csv', estimate.map=TRUE) # allow for estimating genetic map
    g2 <- sim.geno(geno, n.draws=%i, step=5)
    map <- c()
    for (x in g2$geno) {
       if (exists('pulled')) {
          pulled <- abind(pulled, x$'draws', along=2)
       }
       else {
          pulled <- x$'draws'
       }
       map <- c(map, attr(x$draws, 'map'))
    }
    stopifnot(length(map) == dim(pulled)[2])
    IDs <- geno$pheno$id
    """ % (geno_fl, n))
    R("""
    for (i in 1:%i){
       tmp  <- pulled[,,i]
       stopifnot(names(map) == colnames(tmp))
       write.csv(tmp, sep=',', file=paste('new_geno_', i,'.csv', sep=''), row.names=IDs, col.names=colnames(tmp))
       }""" % n)

def compute_one(flnm, gen_model, pheno, Psi, time_points, selMat):
    """
    Compute statistics for one genotype file.

    Parameter:
    flnm : the name of genotype file.
    gen_model : a function that given a sequence of genotypes will generate design matrix.
    pheno : an array of phenotypes, where each row is a sample.
    Psi : basis function values.
    time_points : a sequence of time points
    selMat : the selection matrix.

    The genotype file is assumed to have marker names as the first row,
    genotypes with sample IDs as the first column.
    """
    IDs = list(); geno = list()
    with open(flnm) as fl:
        dt = csv.reader(fl, delimiter=',', quotechar='"')
        markers = next(dt)[1:]          # not used here
        [(IDs.append(x[0]), geno.append([float(y) for y in x[1:]])) for x in dt]
    return [x.pval for x in compute_qtl([gen_model(x) for x in zip(*geno)], pheno, Psi, time_points, Regression, selMat)]
        
def compute(*args, **kwargs):
    """
    Compute normal statistics for each imputation and then calculate amended statistics.
    """
    def new_stat(x):
        lower = scoreatpercentile(x, 5.)
        higher = scoreatpercentile(x, 95.)
        tm = tmean(x, limits=(lower, higher))
        tv = tvar(x, limits=(lower, higher))
        return tm/2. + tv/8.
    res = list()
    for f in iglob("new_geno_*.csv"):
        # compute one imputation and append
        res.append(compute_one(f,*args, **kwargs))
    return np.array([2*np.log( np.exp(np.array(x)).mean() ) for x in zip(*res)])
        
