#!/usr/bin/env python
"""
Raw functional regression using QR for numerical stability and efficience.
"""

import numpy as np
from numpy import ma
from numpy.linalg import qr, inv, cholesky, solve
from scipy.stats import chi2

from utils import *
from sparseDiag import SparseDiag

class Regression(object):
    def __init__(self,  Y, Psi, Z, selMat=None):
        """
        Constructor for Regression.

        Parameters:
        -------------------
        Y : The raw observation array where each row is a random processes at all time points
            and each column is all random processes at certain time.  For consistence's sake, Y
            should be an array, not a matrix.
            The dimension should be Nxn, where n is the length of ts.
        Psi : A matrix whose rows are the values of basis functions at one time point and whose columns
              are one basis function at all time points.  The dimension is therefore (degree)xn, where
              n is the length of ts.
        Z : Design matrix.  Its dimension should be Nxq, where q is the length of beta(t).
        selMat : A selection matrix choosing certain rows of B so that only certain covariates are included
                 in hypothesis testing.
        """
        ## Y needs to be an array, not a matrix.  This is important in calculating covariance.
        self.N, self.q = Z.shape
        self.n, self.k = Psi.shape
        self.Y = Y
        self.Psi = np.mat(Psi)
        self.Z = np.mat(Z)
        if selMat is None:              # no selection
            self.deg_fr = self.q * self.k
            self.selMat = np.mat(np.eye(self.q))
        else:                           # only test certain covariates
            self.deg_fr = len(selMat) * self.k
            self.selMat = np.mat(selMat)

        self._pre_compute()

    def _pre_compute(self):
        """
        Pre-computing some quantities.
        """
        ZQ, ZR = qr(self.Z)
        PsiQ, PsiR = qr(self.Psi)
        # save intermediate results for later reuse
        ZinvRQ = solve(ZR, ZQ.T)
        PsiinvRQ = solve(PsiR, PsiQ.T)
        self.kprod = np.kron(self.selMat * ZinvRQ, PsiinvRQ)
        # B can be directly computed from QR decomposition
        self._B = ZinvRQ * self.Y * PsiinvRQ.T 
        # compute cov_y and cov_b
        tmp = data_cov(self.Yerr.T, self.deg_fr) #data_cov(self.Y.T)
        tmp_L = cholesky(tmp)
        cov_y = SparseDiag(tmp_L,self.N)
        tmp = cov_y.left_prod(self.kprod)
        self.tmp = tmp
        #self.cov = tmp * tmp.T
        self.cov_sq = np.mat(qr(tmp.T, mode='r'))
        self.Y = np.mat(self.Y)         # need Y to be a matrix for products below
        

    @property
    def B(self):
        return np.asarray(self._B)
    @property
    def Yhat(self):
        return np.asarray(self.Z * self._B * self.Psi.T)
    @property
    def Yerr(self):
        return np.asarray( self.Y - self.Yhat )
    @property
    def rss(self):
        return np.linalg.norm(self.Yerr, ord='fro')
    @property
    def ssq(self):
        if not hasattr(self, '_ssq'):
            tmp = np.dot(self.kprod, vec(np.asarray(self.Y.T))[:,np.newaxis]) # vec expects an array
            sq = solve(self.cov_sq.T, tmp) # double transpose in order to solve function
            self._ssq = float(sq.T * sq)
        return self._ssq
    @property
    def pval(self):
        return chi2.sf(self.ssq, self.deg_fr)

def basic_reg(Y,Psi,Z,selMat):
    """
    Basic version of regression, without QR.  For debugging purpose.
    """    
    N, q = Z.shape
    n, k = Psi.shape
    deg_fr = len(selMat)*q
    tmp = data_cov(Y.T, deg_fr)
    cov_y =  np.kron(np.eye(N), tmp)
    vecYT = vec(Y.T)
    Y,Psi, Z = [np.mat(x) for x in (Y,Psi,Z)]
    kinv = np.kron(selMat*inv(Z.T*Z), inv(Psi.T*Psi))
    kzp = np.kron(Z,Psi)
    b = kinv * kzp.T * np.mat(vecYT[:,np.newaxis])
    tmp = kinv * kzp.T
    cov = tmp * cov_y * tmp.T
    ssq = float(b.T * inv(cov) * b)
    # gather results to return
    from collections import namedtuple
    regre = namedtuple("Regre", "B,cov,ssq")
    return regre(b, cov, ssq)

def select(rows):
    """
    Generate a selection matrix so that rows specified in 'rows' are selected
    if left-multiplied by the selection-matrix.

    Note that rows are specified by Numpy's convention.
    """
    tmp = list(rows)
    tmp.sort()
    tmp_mat = np.eye(tmp[-1]+1)
    return np.mat(tmp_mat[tmp,:])
