#!/usr/bin/env python
"""
B-spline basis functions.
"""
import numpy as np
from pygsl.bspline import bspline

class Bspline(object):
    """
    B-spline basis functions.
    """
    def __init__(self, order, nbreak, knots=None):
        """
        Construct a b-spline basis.

        Parameters:
        ------------
        order :  the order of polynomial. Note order = degree + 1.
        nbreak : the number of breakpoints.
        knots : a sequence of non-decreasing values to put knots on.

        Note that the number of coefficients, also the degree of freedom,
        is equal to number of break-points + order - 2.

        The number of knots should be at least as many as the number of
        breakpoints.
        """
        self.order = order
        self.nbreak = nbreak
        self.knots = knots

        self.basis = bspline(order, nbreak)
        
    def values(self, tpts):
        """
        Evaluate the b-spline at some time-points.

        Parameters:
        -------------
        tpts : a sequence of time-points.

        Return:
        ------------
        An array whose rows are basis functions at one time-point,
        and columns are one basis function at all time-points.
        """
        x = np.asarray(tpts)
        knots = self.knots
        if knots is None:
            self.basis.knots_uniform(x[0], x[-1])
        else:
            self.basis.knots(np.asarray(knots))
        return self.basis.eval_vector(x)
