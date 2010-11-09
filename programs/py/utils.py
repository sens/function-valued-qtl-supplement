#!/usr/src/env python
"""
A collection of utilities.
"""

import numpy as np
from numpy import ma

def vec(A):
    """
    vec operation.

    vec(A) is the stacking of columns of A.

    Parameters:
    ---------------------
    A : This is an array.
    """
    return A.ravel('F')

def vect(A):
    """
    vec-tranpose operation.

    vect(A) is equal to vec(A^T)^T.

    Parameters:
    ---------------------
    A : This is an array.
    """
    return A.ravel('C')

def data_cov(Y,df):
    """
    Compute the variance-covariance matrix from N independent samples.

    If there are N samples, each of which is a vector y_i, then
    we want to compute
               1/(N-n) \Sigm_i^N (y_i-mu)(y_i-mu)^T,
    where mu is a vector of means.
    
    Parameters:
    ----------------------
    Y : An array whose columns are different samples.  It is assumed that the ordering within
        each column is same, that is, the ith element of all colunms are the measurements of same
        quantity.

    Returns:
    ----------------------
    Variance-covariance matrix of the data.
    """
    Y = np.asarray(Y)
    if len(Y.shape) == 1 or Y.shape[1] == 1:
        return 0
    n, N = Y.shape               # number of columns is the number of samples
    
    return np.mat(np.cov(Y, rowvar=1))*(N-1)/(N-df)
