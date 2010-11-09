#!/usr/bin/env python
"""
Functional QTL.
"""
import multiprocessing as mp
import numpy as np
from matplotlib import pyplot as plt
#from raw_regress_qr import Regression
        
def add_ones(seq):

    """
    Prepend a column of 1's to a matrix.
    """
    orig_mat = np.asarray(seq)        
    ones = np.ones((orig_mat.shape[0],))
    return np.vstack((ones, orig_mat)).T
           
def encode3(geno):
    """
    Given a sequence of genotypes, generate 3 encoding of it as
    the design matrix.

    Parameters:
    --------------
    geno : a sequence of genotypes.  Assumed to have 3 genotypes.

    This function will generate as many encoding as there are
    different genotypes.  Each encoding will have one genotype
    coded as 1 and all others as 0.
    """
    geno = np.asarray(geno)
    for x in set(geno):
        tmp = (geno == x)[:,np.newaxis]
        yield tmp.astype(int)
    
class OneRegression(object):
    """
    Simulate regression one times.
    """
    def __init__(self, Y, Psi, Z, T, reg_class=None, selMat=None):
        """
        Given a Stochastic process, basis, Z, and time-points
        prepare simulation.
        """
        self.Z = np.asanyarray(Z)
        self.T = np.asarray(T)
        self.sample_size = Z.shape[0]
        # prepare observations
        self.Y = Y
        # prepare Psi
        self.Psi = Psi

        # optionally run regression if it is provided
        if reg_class is not None:
            self.regress(reg_class, selMat)
    def regress(self, reg_class, selMat=None):
        """
        Make a Regression object, or an object of a descendant of
        Regression.
        Parameters:
        -----------
        reg_class: the regression class.
        """
        regression = reg_class(self.Y,self.Psi,self.Z, selMat)
        self.ssq = regression.ssq #- regression0.ssq  # to force computation early
        self.pval = regression.pval
        self.df = regression.deg_fr
        self.Yhat = regression.Yhat
        self.rss = regression.rss
        self.B = regression.B
        # for debugging purpose
        #self.cov_sq = regression.cov_sq

        return self

    def plot_mean(self):
        """
        Plot the mean curves of observations and predictions.

        This assumes that category is specified in the second column of Z.
        """
        cat_vec = np.asarray(self.Z[:,1]).squeeze()
        color_vec = ['b', 'g', 'r']
        Y = np.asanyarray(self.Y)        # Y could be a masked array, which we want to keep
        Yhat = np.asarray(self.Yhat)     # Yhat could be a matrix, which we want to change
        N = len(cat_vec)
        assert Y.shape[0] == N, "Y shape is wrong."
        assert Yhat.shape[0] == N, "Yhat shape is wrong."
    
        for cnt, cat in enumerate(set(cat_vec)):
            #plt.figure()
            Ymean = Y[cat_vec==cat,:].mean(axis=0)
            Ystd = Y[cat_vec==cat,:].std(axis=0)
            Yhatmean = Yhat[cat_vec==cat,:].mean(axis=0)
            plt.errorbar(self.T,Ymean,Ystd,fmt=color_vec[cnt]+'-')
            plt.plot(self.T, Yhatmean, color_vec[cnt]+'--')
        plt.show()

class RegressionResults(list):
    """
    Results list for Regressions whose dimension of Z and Psi are same.

    It is assumed that all regressions have same degree of freedom.

    This will be a candidate for mix-in.  For example, for simulations,
    plotting histograms and computing type1 errors can be done by taking
    results from this class and performing computation in another class.
    """
    def __init__(self,many_reg):
        """
        Collect many regression results into a list.

        Parameters:
        -------------
        many_reg : a collection of OneRegression.

        Note: internal consistency is not guaranteed if this list is
        augumented in any way after construction.
        """
        super(RegressionResults, self).__init__(many_reg)
        self.ssqs = np.array([x.ssq for x in self]) # this is already a list, and this obviates a bug when many_reg is an iterator
        self.pvals = np.array([x.pval for x in self])
        self.df = self[0].df
        assert all(self.df == x.df for x in self), "The degree of freedom is not identical."
        self.T = self[0].T
        assert all(np.all(self.T == x.T) for x in self), "The time-points are not identical."

def compute_qtl(Zs, Y, Psi, T, reg_class, selMat=None):
        """        
        Parallel computing of QTLs.

        Parameters:
        markers : the name for which QTLs are being computed
        Zs : a sequence of Z.
        The rest are passed to OneRegression.
        """        
        ## workers = mp.Pool()
        ## runs = [workers.apply_async(OneRegression,(Y,Psi,Z,T,reg_class, selMat)) for Z in Zs]
        ## regs = [x.get(timeout=10) for x in runs]
        
        regs = [OneRegression(Y,Psi,np.asarray(Z),T,reg_class, selMat)for Z in Zs]
                
        return regs 
            
    
