#!/usr/bin/env python
"""
All things related to Fourier series as related to functional data analysis.
"""

import numpy as np
from numpy.linalg import lstsq
import math

class FourierBasis(object):
    """
    Implement a Fourier basis function as a sequence of basis functions.

    """
    def __init__(self, period, degree):
        """
        A constructor for finite sequence of Fourier basis function.
        """
        self.period = period
        self.degree = degree if (degree%2==1) else (degree+1) # the degree here is same as 2*degree+1
        self.omega = 2*np.pi/period
    def values(self, tpts):
        """
        Compute the values of basis functions at specified time points.

        tpts           : time points where basis functions are to evaluated at.

        return value   : a matrix where a row is the values of basis functions
                         at one time point and a column is the values of one basis
                         functions at all time points.
        """
        tpts = np.asarray(tpts)         # in case non-array is passed in
        tpts = np.squeeze(tpts)
        l = len(tpts)
        p = np.arange(2, self.degree, 2)
        val = np.empty( (l, self.degree) )
        val[:,0] = 1/math.sqrt(2)
        
        x = tpts * self.omega / 2.
        a = x[:,np.newaxis] * p[np.newaxis,:]
        val[:, 1:-1:2] = np.sin(a)
        val[:, 2::2] = np.cos(a)
        return val/math.sqrt(self.period/2.)
        
class FourierSeries(object):
    """
    A Fourier series object is a truncated Fourier series whose
    coefficients are estimated from a time series using least-square.
    """
    def __init__(self, frbs, ts):
        """
        The constructor takes a Fourier basis and time series in order
        to estimate the coefficients.
        """
        self.basis = frbs
        self.Y = ts.values            # Viewing multivariate time series as multiple univariate time series necessitates a transpose
        self.T = ts.times
        self.coef, self.resid = self._lsq()
        self.Yhat = self._predict()
        self.Yerr = self.Y - self.Yhat

    def _lsq(self):
        """
        Performing least-square using time series and fourier basis function.

        The least-square obey the following formula:
        coeffcients = (P^TP)^{-1}P^TY,
        where P is the matrix whose rows are the values of all the basis functions at one time point,
        and columns are the values of one basis function at all time points.
        Y is a column vector of observations.
        """
        tp = self.T              # get the time points
        self.bval = self.basis.values(tp)    # get the values of basis functions at all time points
        coef, resid = lstsq(self.bval, self.Y.T)[:2] # return only the coefficients
        return coef, resid
    def _predict(self):
        """
        Using estimated coefficients to predict observations using known time points.
        """
        return np.dot(self.bval, self.coef).T
        
def regress(design, Y):
    """
    Perform functional regression using basis expansion of orthogonal basis functions.

    design : a fixed, design matrix
    Y : the coefficients of basis expansion. Note that this should be the transpose of the coefficients obtained from least-square.
    """
    return lstsq(design, Y)[:2]
