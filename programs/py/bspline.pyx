
"""
B-spline basis functions.
"""
cimport python_exc
import numpy as np
cimport numpy as np
#from pygsl.bspline import bspline

cdef extern from "gsl/gsl_bspline.h":
    ctypedef struct gsl_vector:
        size_t size
        
    ctypedef struct gsl_bspline_workspace:
        size_t n    
    gsl_bspline_workspace* gsl_bspline_alloc(size_t k, size_t nbreak)
    void gsl_bspline_free(gsl_bspline_workspace* w)
    int gsl_bspline_knots(gsl_vector* breakpts, gsl_bspline_workspace* w)
    int gsl_bspline_knots_uniform(double a, double b, gsl_bspline_workspace* w)
    int gsl_bspline_eval(double x, gsl_vector* B, gsl_bspline_workspace* w)
#cdef extern from "gsl/gsl_vector.h":
    gsl_vector* gsl_vector_alloc(size_t n)
    void gsl_vector_free(gsl_vector* v)
    double gsl_vector_get( gsl_vector* v, size_t i)
    void gsl_vector_set(gsl_vector* v, size_t i, double x)
    

cdef class MyBspline:
    cdef gsl_bspline_workspace* w
    def __cinit__(self, int k, int nbreak): 
        self.w = gsl_bspline_alloc(k, nbreak)
        if self.w==NULL:
            python_exc.PyErr_NoMemmory()
    cpdef int knots_uniform(self,double start, double end):
        return gsl_bspline_knots_uniform(start, end, self.w)
    cpdef int knots(self,object seq):
        cdef gsl_vector* bkpts
        cdef int tmp
        cdef int n = len(seq)
        bkpts = gsl_vector_alloc(n)
        if bkpts==NULL:
            python_exc.PyErr_NoMemmory()
        for ind,x in enumerate(seq):
            gsl_vector_set(bkpts, ind, x)
        tmp = gsl_bspline_knots(bkpts, self.w)
        gsl_vector_free(bkpts)
        return tmp
    cpdef eval_vector(self,object seq):
        cdef gsl_vector* tmp
        cdef size_t m = len(seq)
        cdef size_t n = self.w.n
        cdef np.ndarray[np.float64_t, ndim=2] result = np.empty((m,n),dtype=np.float64)
        cdef double x
        
        tmp = gsl_vector_alloc(n)
        for i in range(m):
            x = seq[i]
            gsl_bspline_eval(x, tmp, self.w)
            for j in range(n):
                result[i,j] = gsl_vector_get(tmp,j)
        gsl_vector_free(tmp)
        return result
    def __dealloc__(self):
        gsl_bspline_free(self.w)

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

        self.basis = MyBspline(order, nbreak)
        
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
            self.basis.knots(knots)
        return self.basis.eval_vector(x)
