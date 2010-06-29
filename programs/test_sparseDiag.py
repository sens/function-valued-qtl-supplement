#!/usr/bin/env python
from sparseDiag import SparseDiag
import numpy as np

def test_repeat():
    """
    Testing SparseDiag's ability to deal with repeat block-diagonal matrices.
    """
    repeat = 2                          # two repeats
    elem = np.random.uniform(10, size=(10,12)) # random repeated-matrix
    subj = SparseDiag(elem, repeat)
    A = np.mat( np.kron(np.eye(repeat), elem) ) # set up explicitly the full block-diagonal matrix
    # set up two flanking matrices
    a = np.mat( np.random.standard_normal((2, 10*repeat)) )
    b = np.mat( np.random.standard_normal((12*repeat,3)) )
    # compute and compare
    mine = subj.quadratic(a,b)
    correct = a * A * b
    np.testing.assert_almost_equal(correct, mine)

def test_sequence():
    """
    Testing SparseDiag's ability to deal with non-repeating block-diagonal matrices.
    """
    # set up two different matrices as diagonal
    mat1 = np.random.uniform(20, size=(31,19))
    mat2 = np.random.standard_normal((61,27))
    subj = SparseDiag((mat1,mat2))
    # explicitly construct the block-diagonal matrix
    zero1 = np.zeros((31,27))
    zero2 = np.zeros((61,19))    
    A = np.vstack( (np.hstack((mat1,zero1)),np.hstack((zero2,mat2))) )
    # randomly generate two flanking matrices
    a = np.mat( np.random.standard_normal((13,A.shape[0])) )
    b = np.mat( np.random.standard_normal((A.shape[1],21)) )
    # compute and compare
    mine = subj.quadratic(a,b)
    correct = a * A * b
    np.testing.assert_almost_equal(correct, mine)

def test_left_prod():
    """
    Testing SparseDiag's left-multiplying ability.
    Only repeated block-diagonal matrices are tested.
    """
    repeat = 2                          # two repeats
    # set up block-diagonal matrix
    elem = np.random.standard_normal((52,52))
    subj = SparseDiag(elem, repeat)
    A = np.mat( np.kron(np.eye(repeat), elem) )
    # set up left flanking matrix
    a1 = np.random.standard_normal((43,52))
    a2 = np.random.standard_normal((43,52))
    a = np.hstack( (a1,a2) )
    # compute and compare
    mine = subj.left_prod(a)
    correct = a * A 
    np.testing.assert_almost_equal(correct, mine)
