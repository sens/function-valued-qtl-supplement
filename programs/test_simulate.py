#!/usr/bin/env python
"""
Test various classes and functions in simulate.py.
"""

from simulate import *

def test_StochasticProc():
    """
    Testing StochasticProc.
    """
    # functions list
    f_1 = lambda x: np.cos(x)*np.exp(-x/10)
    f_0 = lambda x: 0.92 - 0.43*x - 1.42*(x**2)
    # X
    sample_size = 5
    x_1 = [0.5] * sample_size    
    X = np.array([[1]*sample_size, x_1]).T
    T = np.linspace(0,2,200)
    
    stoc_proc = StochasticProc(X, [f_0, f_1])
    Y = stoc_proc.values(T)

    # plot two Ys
    assert Y.shape == (sample_size,200), "The shape of Y is wrong."
    for i in range(sample_size):
        plt.figure()
        plt.plot(T, Y[i])
        plt.title("Y%d"%i)
    plt.show()
    #raise Exception()
