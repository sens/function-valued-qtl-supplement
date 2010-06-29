#!/usr/bin/env python
"""
A simple implementation of block diagonal matrix in quadratic product.
"""

import numpy as np

class SparseDiag(object):
    """
    A block diagonal matrix is a matrix made up of a series of
    matrices along its diagonal.
    """
    def __init__(self, data1, repeat=0):
        """
        Construct a block diagonal matrix by either a sequence of matrices
        or one matrix and a count of how many repeats.
        """
        if repeat > 0:
            self.data = np.mat(data1)
            self.repeat = repeat
            self.shape = tuple(x*repeat for x in self.data.shape)
        else:
            self.data = tuple(np.mat(x) for x in data1)
            self.shape = reduce(lambda x,y: (x[0]+y[0],x[1]+y[1]), [x.shape for x in self.data])
            # note that we are not setting repeat in this case
    def quadratic(self,a,b):
        """
        Compute a * self * b

        Return the result as a matrix.
        """
        assert a.shape[1]==self.shape[0] and self.shape[1]==b.shape[0], "Matrix dimensions do not match."
        result = np.mat( np.zeros((a.shape[0],b.shape[1])) )
        a_start = 0                     # the beginning of submatrix for a
        b_start = 0                     # similarly for b
        if hasattr(self, 'repeat'):
            seq_mat = (self.data for i in range(self.repeat))
        else:
            seq_mat = self.data
        for x in seq_mat:
            a_end = a_start + x.shape[0]
            b_end = b_start + x.shape[1]
            result += a[:,a_start:a_end] * x * b[b_start:b_end,:]
            a_start, b_start = a_end, b_end
        return result
    def left_prod(self, left):
        """
        Compute left*self
        """
        assert left.shape[1] == self.shape[0], "Incomfortable dimensions."
        result = np.mat( np.empty((left.shape[0], self.shape[1])) )
        left_start = 0
        if hasattr(self, 'repeat'):
            seq_mat = (self.data for i in range(self.repeat))
        else:
            seq_mat = self.data
        for x in seq_mat:
            left_end = left_start + x.shape[0]
            result[:, left_start:left_end] = left[:,left_start:left_end] * x
            left_start = left_end
        return result
