#!/usr/bin/env python
"""
Simulations comparing functional ANOVA against cross-sectional.
"""
from collections import namedtuple
import numpy as np
from matplotlib import pyplot as plt
from bspline import Bspline
from raw_regress import Regression, RegressionMissing
from simulate import DataSet, parallel_sim

from rpy2 import robjects as robj
from rpy2.robjects import r as R
import rpy2.robjects.numpy2ri

def simulate_cross_section(fcn, T, sample_size=400,is_masked=False):
    """
    Simulate additive genetic model and compare the power
    of functional ANOVA against cross-sectional means.
    """
    # overall parameters
    nruns = 1000    
    basis = Bspline
    # specify proportions of different genotypes
    prop = [0.5, 0.5]
    # generate datasets
    datasets = [DataSet(prop, fcn, sample_size, T, basis,is_masked,4,2) for i in range(nruns)]

    res_f =  parallel_sim(datasets, RegressionMissing) if is_masked else parallel_sim(datasets, Regression)

    res_c = [compute_cross_section(data.Y, data.Z.squeeze()) for data in datasets]

    return res_f, res_c

def compute_cross_section(Y,labels):
    """
    Compute ANOVA using cross-section from multivariate data.

    Parameters:
    ------------
    Y : The multivariate data where each row is one observation.
    labels : the class of samples.
    """
    
    assert len(labels) == Y.shape[0]

    means = Y.mean(axis=1)
    robj.globalEnv['y'] = robj.FloatVector(means)
    robj.globalEnv['x'] = robj.IntVector(labels)
    R("x = factor(x)")
    res = R("summary(aov(y~x-1))[[1]]$'Pr(>F)'")[0]
    return res

def compute_power(pvals, thres=0.05):
    pvals = np.asarray(pvals)
    assert pvals.ndim == 1, "pvals is not one dimensional"
    return sum(pvals < thres)/ float(pvals.shape[0])

Power = namedtuple("Power", "functional cross_section")
def simu_time(fcn, T, bottom=5,top=10, **kwargs):
    """
    Simulation of different numbers of time points.

    Parameters:
    -------------
    fcn : the function that is the genetic effect.
    T   : a function given the number of time points will generate them.
    bottom : the minimum number of time points.
    top : the maximum +1 number of time points.
    """
    powers = [list()]*(top-bottom)
    for x in range(bottom,top):
        tmp = simulate_cross_section([fcn], T(x), **kwargs)
        powers[x-bottom] = Power(compute_power(tmp[0].pvals), compute_power(tmp[1]))
    return powers

def simu_sample(fcn, T, ntm, nsamples, **kwargs):
    """
    Simulation of different sample sizes.

    Parameters:
    -------------
    ntm : a set number of time points
    nsamples : a sequence of # of samples
    """
    powers = [list()]*len(nsamples)
    for k,n in enumerate(nsamples):
        tmp = simulate_cross_section([fcn], T(ntm),n, **kwargs)
        powers[k] = Power(compute_power(tmp[0].pvals), compute_power(tmp[1]))
    return powers

def plot_sect(res, xseq, xlabel,title=""):
    """
    Plot the results of comparison between functional and cross-sectional.
    """
    first = lambda lst: [x[0] for x in lst]
    second = lambda lst: [x[1] for x in lst]
    plt.plot(xseq, first(res), linewidth=2.0, label='Functional')
    plt.plot(xseq, second(res),'--', linewidth=2.0,label='Cross-sectional')
    
    plt.legend(loc='lower right')
    plt.title(title) if title else None      # add title if set
    plt.ylabel("power")
    plt.xlabel(xlabel)
    plt.show()

if __name__ == "__main__":
    quad = lambda x: 0.001*x**2+0.1*x+2.5
    quad_T = lambda n: -50.0 + np.arange(n)*2.0
    expo = lambda x: 1.-0.1*np.exp(-0.005*x)
    expo_T = lambda n: -460. + np.arange(n)*16.0
    logit = lambda x: 1/(1+1.0*np.exp(-x))
    logit_T = lambda n: -6. + np.arange(n)*1.0
    res_quad = simu_time(quad, quad_T,top=11)
    #res_quad_missing = simu_time(quad, quad_T,top=11, is_masked=True)
    res_expo = simu_time(expo, expo_T, bottom=4)
    #res_expo_missing = simu_time(expo, expo_T, bottom=4,is_masked=True)
    res_logit = simu_time(logit, logit_T,top=9)
    #res_logit_missing = simu_time(logit, logit_T,top=9,is_masked=True)
    sam_quad = simu_sample(quad, quad_T, 7, 200+100.*np.arange(14))
    #sam_quad_missing = simu_sample(quad, quad_T, 7, 200+100.*np.arange(14),is_masked=True)
    sam_expo = simu_sample(expo, expo_T, 5, 200+100.*np.arange(5))
    #sam_expo_missing = simu_sample(expo, expo_T, 5, 200+100.*np.arange(5),is_masked=True)
    sam_logit = simu_sample(logit, logit_T, 6, 200+100.*np.arange(8))
    #sam_logit_missing = simu_sample(logit, logit_T, 6, 200+100.*np.arange(8),is_masked=True)
