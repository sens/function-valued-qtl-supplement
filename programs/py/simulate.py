#!/usr/bin/env python
"""
Simulating Stochastic processes of the form
   y(t) = b_0(t) + x_1*b_1(t)+e(t)
"""
import multiprocessing as mp
from itertools import islice
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import chi2
from pymc import gp
from pymc.gp.cov_funs import matern
from fqtl import RegressionResults, OneRegression
#from raw_regress_qr import Regression#, RegressionMissing
from bspline import Bspline
from fourier import FourierBasis
        
class StochasticProc(object):
    """
    Generate Stochastic processes in discrete form as
       Y(t) = X*B(t)+E(t)
    where t will be a vector of time-points
    """
    def __init__(self, Z, fcn_lst):
        """
        Construct Stochastic process from design matrix and a list of
        functions.

        Parameters:
        X : the design matrix.
        fcn_lst : function list.

        The dot product of the ith row of X with fcn_lst at time-points is
        the mean function for the ith Stochastic process.
        """
        self.Z = np.asarray(Z) 
        self.nsamp = self.Z.shape[0]
        self.funcs = fcn_lst
        self.cov = gp.Covariance(matern.euclidean, scale=1, amp=1, diff_degree=0.8)

        assert self.Z.shape[1] == len(fcn_lst), "X's dimensions do not match with the number of functions."

    def values(self, T):
        """
        Parameters:
        ------------
        T : a vector of time-points.

        Return values:
        ------------
        a 2-D array where rows are processes and columns samples.
        """
        result = np.empty((self.nsamp, len(T)))        
        for i in range(self.nsamp):
            mean = lambda x: np.dot(self.Z[i], np.array([f(x.squeeze()) for f in self.funcs])) # need squeeze() because pymc appends a dimension.
            result[i] = np.array( gp.Realization(gp.Mean(mean), self.cov)(T) )
        return result

def gen_genotypes(proportions, size):
    """
    Randomly generate genotypes according fixed proportion.

    Parameters:
    ---------------------
    proportions : a sequence of probabilities for genotypes.
    size : the length of generated genotypes.

    Returns:
    --------------------
    A list of length 'size' made up of 0,1,...,n-1, where n is
    the length of 'proportion.'
    """
    def encode(value):
        """
        Given a list of probabilities and a uniformly generated random values,
        generate a list of 0,1,...,n-1 according to which interval the value falls
        into.
        """
        for i in range(n)[:-1]:
            if value >= prob[i] and prob[i+1] > value:
                return i
        return n-1
    n = len(proportions)
    assert sum(proportions) == 1, "The list of probabilities do not add up to 1."
    tmp = np.random.rand(size)
    prob = np.cumsum(proportions)
    return [encode(x) for x in tmp]

def make_mask(shape):
    tmp = np.random.uniform(size=shape)
    return tmp < 0.05                   # random 5%
class DataSet(object):
    """
    A data-set for simulation.

    The DataSet class gathers into one place all the values needed for simulation
    and function calls to generate them.
    """
    def __init__(self, chances, fcn_lst, sample_size, T, basis, is_masked=False, *args, **kwargs):
        """
        Construct a random data-set for simulation.

        Parameters:
        ------------
        chances : a list of probabilities for different genotypes.
        fcn_lst : Funtion list for Stochastic process.
        T : time-points.
        basis : the basis class used for this data-set.
        args, kwargs : arguments to be passed to basis class' constructor.
        """
        self.T = T
        # generate design matrix
        geno = gen_genotypes(proportions=chances, size=sample_size)
        if len(fcn_lst) > 1:
            self.Z = np.array([np.ones(sample_size),geno]).T # add a constant term
        else:
            self.Z = np.array([geno]).T
        # generate observations
        stoc_proc = StochasticProc(self.Z, fcn_lst)
        tmp = stoc_proc.values(T)
        if is_masked:
            self.Y = np.ma.array(tmp, mask=make_mask(tmp.shape))
        else:
            self.Y = tmp
        # get Psi
        self.basis = basis(*args, **kwargs)

class SummarizeSim(object):
    """
    Summarize the result of simulations.
    """
    def __init__(self):
        pass
    def plot_hist(self):
        """
        Plot the histogram of sum of squares and super-impose
        the theorotical values of chi-square distribution.
        """
        plt.clf()
        (_,bins,_) = plt.hist(self.ssqs, bins=50, normed=True)
        x = np.linspace(0,bins[-1],100)
        y_theo = chi2.pdf(x,self.df)
        plt.plot(x,y_theo)
        #plt.savefig('Statistics-%d.png'% N, format='png')
        plt.show()
    def type1_err(self):
        crit = [chi2.ppf(1-thres, self.df) for thres in (0.1, 0.05, 0.01, 0.005, 0.001)]
        return [(self.ssqs > x).sum() / float(len(self.ssqs)) for x in crit]

class ResultsSim(RegressionResults, SummarizeSim):
    """
    """
    def __init__(self, *args, **kwargs):
        super(ResultsSim, self).__init__(*args, **kwargs)

def parallel_sim(datasets, reg_class, selMat=None):
    """
    Parallel simulation.

    Parameters:
    --------------
    datasets : a sequence of data-sets.  The member could either be objects of
               DataSet class or objects of OneRegression class.
    reg_class : the kind of regression to perform.
    """
    workers = mp.Pool()
    sim_runs = [workers.apply_async(OneRegression, (x.Y, x.basis.values(x.T), x.Z, x.T, reg_class, selMat)) for x in datasets]
    workers.close()
    workers.join()
    return ResultsSim(x.get(timeout=1) for x in sim_runs)

