#!/usr/bin/env python
"""
Simulate to generate histograms and type I errors.
"""
import numpy as np
from bspline import Bspline
from raw_regress_qr import Regression, select
from simulate import DataSet, parallel_sim


def simu(nruns = 1000, sample_size=200,is_masked=False):    
    T = np.linspace(0, 1, 20)   
    basis = Bspline
    propor = [0.25, 0.25, 0.5]
    fcn_lst = [lambda x: 0.2*x**2-0.2*x+0.05,
               lambda x: np.zeros((len(x),))]#, lambda x: np.ones((len(x),))] # should be equal to num of genotypes + 1
    datasets = (DataSet(propor, fcn_lst, sample_size, T, basis, is_masked, 3, 2) for i in range(nruns))
    res_s = parallel_sim(datasets, Regression, select([1]))
    
    return res_s

if __name__ == "__main__":
    res = simu(5000, 300,is_masked=False)
