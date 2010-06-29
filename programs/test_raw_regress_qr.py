#!/usr/bin/env python
"""
The test module for raw_regress.py
"""

import numpy as np
from numpy import ma
from matplotlib import pylab as plt

from bspline import Bspline
from simulate import DataSet
from raw_regress_qr import Regression, basic_reg
from fqtl import OneRegression
from utils import *
    
def test_data_cov():
    def slow_data_cov(Y):
        """
        This is the literate implementation of
                  1/(N-n)\Sigm_i^N (y_i-mu)(y_i-mu)^T,
        where mu is a vector of means.
        """
        Y = np.mat(Y)
        n, N = Y.shape
        ave_Y = np.mean(Y,1)
        tmp = 0
        for i in range(N):
            diff = Y[:,i] - ave_Y
            tmp += diff * diff.T
        return tmp / (N-1)
 
    # test degenerate case
    assert data_cov(np.array([4]))==0
    assert float(data_cov(np.arange(5).reshape((5,1))))==0
    assert data_cov(np.arange(5).reshape((1,5)))== np.var(np.arange(5),ddof=1)
    # test regular case
    a = np.random.random( (30,100) )
    np.testing.assert_array_almost_equal(data_cov(a),slow_data_cov(a))
            
def test_vec_vect():
    from utils import vec, vect
    # degenerate case
    assert vec(np.array([5])) == vect(np.array([5]))
    # regular case
    A = np.random.random( (73,143) )
    assert np.all( vec(A.T) == vect(A).T )
    B = np.asarray(np.mat('3 4 9 12 90; 21 38 87 6 27; 49 52 5 18 29'))
    assert np.all( vec(B) == np.array([ [3,21,49,4,38,52,9,87,5,12,6,18,90,27,29]])[:,np.newaxis] )

def test_Regression():
    # set up a random dataset
    sample_size = 400
    T = np.linspace(0, 1, 20)   
    basis = Bspline
    propor = [0.25, 0.25, 0.5]
    fcn_lst = [lambda x: 0.2*x**2-0.2*x+0.05,
               lambda x: np.zeros((len(x),))]#, lambda x: np.ones((len(x),))] # should be equal to num of genotypes + 1
    data = DataSet(propor, fcn_lst, sample_size, T, basis, False, 3, 2)
    Psi = data.basis.values(data.T)
    # compute and compare
    desired_res = Regression(data.Y, Psi, data.Z)
    actual_res = basic_reg(data.Y,Psi,data.Z,np.eye(2))
    np.testing.assert_almost_equal(vect(desired_res.B)[:,np.newaxis], actual_res.B)
    np.testing.assert_almost_equal(desired_res.ssq, actual_res.ssq)
    cov = np.dot(desired_res.cov_sq.T, desired_res.cov_sq)
    np.testing.assert_almost_equal(cov, actual_res.cov)
