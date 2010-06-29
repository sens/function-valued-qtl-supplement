#!/usr/bin/env python

import numpy as np

class TimeSeries(object):
    """
    A time series class that contains both value and time.
    """

    def __init__(self, values, times):
        """
        The constructor takes values and the time corresponding to the values.

        values : values must be either a vector or a 2-dimensional
                 array.  If it is a vector then the time series is univariate.
                 If it is a 2-dimensional array, the columns are multivariate values.
        times : a sequence of times points.
        """
        self.values = np.asarray(values)
        self.times = np.asarray(times)
        # check that the dimensions and lengths match
        assert len(self.values.shape) < 3, "The values are more than 2 dimensions."
        assert len(self.times.shape)==1, "The times are not a vector."
        assert len(self.times) == self.values.shape[1], "The number of values columns and time points are not same."
        
    def __len__(self):
        return len(self.times)
