This document is related to the Python part of the package, which is an independent implementation and stands on its own.

1. Dependence
The basic dependence is Python (ver. 2.6.5), Python packages Numpy (ver. 1.5.0) and Scipy (ver. 0.8.0) for numeric computation, package matplotlib (ver. 1.0.0) for plotting in Python, and Cython (ver. 0.13) for interfacing C code; GNU GSL (ver. 1.14) is used for computing B-splines.  All the above software are usually packaged by various Linux distributions, but Windows and Mac users need to install them.  Python package PyML (ver. 0.7.4.1) is rarely packaged so it needs to be installed by hand.

For comparison against cross-sectional method we need "aov" package in R (ver. 2.12.0) and the bridge from Python to R called rpy (ver. 2.1.5).


2. Files

fourier.py  Fourier basis function, not actively used.
raw_regress_qr.py  The main functional regression implementation.
simulate.py  Various classes and utility functions for simulation.
simu_hist.py   The simulation file for type-I errors.
simu_crst.py   The simulation file for power calculation.
fqtl.py   The interface to functional regression.
sparseDiag.py   A custom implementation of sparse diagonal matrix and its product.
utils.py   A few utility functions.
impute.py   A simple example of multiple imputation using functional regression.  This file is less tested than others.


3. Running examples

To rerun simulations for type-I errors use the following lines from Python prompt for 5000 runs and 300 samples.
from simu_hist import *
res = simu(5000, 300,is_masked=False)

To plot the histograms,
res.plot_hist()
and to find out the type-I errors for 0.1, 0.05, 0.01, 0.005, 0.001
res.type1err()


For the power comparison, the script will compute all 6 comparisons and plot the figures by running
it.  Note that due to quirks of matplotlib, figures can only be seen one at a time, so the second
figure will be shown after closing the first figure.

It is recommended that those scripts be run within ipython using "%run" command.  This is comparable to running the scripts but afterward the user gets a fully functioning Python shell with access to all the variables and functions defined in the scripts.
