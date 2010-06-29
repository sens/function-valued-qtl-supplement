from pygsl import bspline
from pygsl import multifit
import numpy as np
from scipy import random
from matplotlib import pylab as plt

def test_bspline():
    x = np.linspace(0,15,200)
    ym = np.cos(x)*np.exp(-x/10)
    y = ym + random.normal(0,0.5,(200,))

    cubic = bspline.bspline(5, 9)
    cubic.knots(np.array([0,1,2,4,5,7,11,13,15]))
    X1 = cubic.eval_vector(x)
    from FDA.bspline import Bspline
    X = Bspline(5,9,np.array([0,1,2,4,5,7,11,13,15])).values(x)
    assert np.all(X==X1)
    c, cov, chisq = multifit.linear(X,y)
    yhat1, yerr1 = multifit.linear_est_matrix(X,c,cov)
    yhat2 = np.dot(X,c)
    #plt.plot(x,yhat2,'r.-', x, y, 'bo', x,ym, 'g.--'); plt.show()
