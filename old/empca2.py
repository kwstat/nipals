#!/usr/bin/env python

# Modifications of the original code by Kevin Wright
# Re-naming, comments, etc.
# Emacs: M-x run-python, C-c C-l
# data[i] means ith row?
# Insert breakpoint with pdb.set_trace()
#   cont, exit

"""
Weighted Principal Component Analysis using Expectation Maximization
Classic PCA is great but it doesn't know how to handle noisy or missing
data properly.  This module provides Weighted Expectation Maximization PCA,
an iterative method for solving PCA while properly weighting data.
Missing data is simply the limit of weight=0.
Given data[nobs, nvar] and weights[nobs, nvar],
    m = empca(data, weights, options...)
Returns a Model object m, from which you can inspect the eigenvectors,
coefficients, and reconstructed model, e.g.
    pylab.plot( m.eigvec[0] )
    pylab.plot( m.data[0] )
    pylab.plot( m.model[0] )
    
For comparison, two alternate methods are also implemented which also
return a Model object:
    m = lower_rank(data, weights, options...)
    m = classic_pca(data)  #- but no weights or even options...
    
Stephen Bailey, Spring 2012
"""

from __future__ import division, print_function

import pdb # debugging
import numpy as np
import sys
from scipy.sparse import dia_matrix
import scipy.sparse.linalg
import math

# ------------------------------------------------------------------------------

class Model(object):
    """
    A wrapper class for storing data, eigenvectors, and coefficients.
    
    Returned by empca() function.  Useful member variables:
      Inputs: 
      - data   [nobs, nvar]
      - weights[nobs, nvar]
      - eigvec [ncomp, nvar]
      
      Calculated from those inputs:
      - coeff  [nobs, ncomp] - coeffs to reconstruct data using eigvec
      - model  [nobs, nvar] - reconstruction of data using eigvec,coeff
    
    Not yet implemented: eigenvalues, mean subtraction/bookkeeping
    """

    
    def __init__(self, eigvec, data, weights):
        """
        Create a Model object with eigenvectors, data, and weights.
        
        Dimensions:
        - eigvec [ncomp, nvar]  = [k, j]. Each row is unit length
        - data   [nobs, nvar]  = [i, j]
        - weights[nobs, nvar]  = [i, j]
        - coeff  [nobs, ncomp]  = [i, k]        
        """
        self.eigvec = eigvec
        self.ncomp = eigvec.shape[0]
        
        self.set_data(data, weights)

        
    def set_data(self, data, weights):
        """
        Assign a new data[nobs,nvar] and weights[nobs,nvar] to use with
        the existing eigenvectors.  Recalculates the coefficients and
        model fit.
        """
        self.data = data
        self.weights = weights

        self.nobs = data.shape[0]
        self.nvar = data.shape[1]
        self.coeff = np.zeros( (self.nobs, self.ncomp) )
        self.model = np.zeros( self.data.shape )
        
        #- Calculate degrees of freedom
        ii = np.where(self.weights>0)
        self.dof = self.data[ii].size - self.eigvec.size  - self.ncomp*self.nobs
        
        #- Cache variance of unmasked data
        self._unmasked = ii
        self._unmasked_data_var = np.var(self.data[ii])
        
        self.calc_coeff()

# ------------------------------------------------------------------------------
        
    def calc_coeff(self):
        """
        Solve for coeff[i,k] such that data[i] ~= Sum_k: coeff[i,k] * eigvec[k]
        """

        for i in range(self.nobs):
            self.coeff[i] = solve_weighted(self.eigvec.T, self.data[i], self.weights[i])

        self.calc_xhat()

# ------------------------------------------------------------------------------        
    def calc_eigvec(self):
        """
        Solve for eigvec[k,j] such that data[i] = Sum_k: coeff[i,k] eigvec[k]
        """

        #- Utility function; faster than numpy.linalg.norm()
        def norm(x):
            return np.sqrt(np.dot(x, x))
            
        #- Make copy of data so we can modify it inside this function only
        data = self.data.copy()

        #- Calculate the eigenvectors one by one, up to ncomp
        for k in range(self.ncomp):
            c = self.coeff[:, k]
            for j in range(self.nvar):
                w = self.weights[:, j]
                x = data[:, j]
                # self.eigvec[k, j] = c.dot(w*x) / c.dot(w*c)
                # self.eigvec[k, j] = w.dot(c*x) / w.dot(c*c)
                cw = c*w
                self.eigvec[k, j] = x.dot(cw) / c.dot(cw)
                                                
            #- Remove this vector from the data before continuing with next
            #? Alternate: Resolve for coefficients before subtracting?
            #- Loop replaced with equivalent np.outer(c,v) call (faster)
            # for i in range(self.nobs):
            #     data[i] -= self.coeff[i,k] * self.eigvec[k]
                                
            data -= np.outer(self.coeff[:,k], self.eigvec[k])    

        ktmp = self.eigvec.copy()
        
        #- Renormalize and re-orthogonalize eigvec. This is GramSchmidt.
        self.eigvec[0] /= norm(self.eigvec[0])
        for k in range(1, self.ncomp):
            for kx in range(0, k):
                c = np.dot(self.eigvec[k], self.eigvec[kx])
                self.eigvec[k] -=  c * self.eigvec[kx]
            self.eigvec[k] /= norm(self.eigvec[k])

        self.calc_xhat() # Recalculate model
# ------------------------------------------------------------------------------
           
    def calc_xhat(self):
        """
        Uses eigenvectors and coefficients to calculate model xhat
        """
        # xhat = P'C, can't we just use matrix multiplication?
        for i in range(self.nobs):
            self.model[i] = self.eigvec.T.dot(self.coeff[i])

# ------------------------------------------------------------------------------

    def chi2(self):
        """
        Returns sum( (model-data)^2 / weights )
        """
        # This uses [ (xhat-data) * sqrt(wt) ] ^ 2
        # Wouldn't it be easier to use (xhat-data)^2 * wt ?
        delta = (self.model - self.data) * np.sqrt(self.weights)
        return np.sum(delta**2)
        
# ------------------------------------------------------------------------------

    def rchi2(self):
        """
        Returns reduced chi2 = chi2/dof
        """
        return self.chi2() / self.dof
        
    def _model_vec(self, i):
        """Return the model using just eigvec i"""
        return np.outer(self.coeff[:, i], self.eigvec[i])
        
    def R2vec(self, i):
        """
        Return fraction of data variance which is explained by vector i.
        Notes:
          - Does *not* correct for degrees of freedom.
          - Not robust to data outliers.
        """
        
        d = self._model_vec(i) - self.data
        return 1.0 - np.var(d[self._unmasked]) / self._unmasked_data_var

    # ------------------------------------------------------------------------------

    def R2(self, ncomp=None):
        """
        Return fraction of data variance which is explained by the first
        ncomp vectors.  Default is R2 for all vectors.
        
        Notes:
          - Does *not* correct for degrees of freedom.
          - Not robust to data outliers.
        """
        if ncomp is None:
            mx = self.model
        else:            
            mx = np.zeros(self.data.shape)
            for i in range(ncomp):
                mx += self._model_vec(i)
            
        d = mx - self.data

        #- Only consider R2 for unmasked data
        return 1.0 - np.var(d[self._unmasked]) / self._unmasked_data_var
                
def _random_orthonormal(ncomp, nvar, seed=1):
    """
    Return array of random orthonormal vectors A[ncomp, nvar] 
    Doesn't protect against rare duplicate vectors leading to 0s
    """

    if seed is not None:
        np.random.seed(seed)

    A = np.random.normal(size=(ncomp, nvar))
    # fixme: temporary identity matrix
    # A =  np.array([[1.0, 0, 0, 0, 0],
    #                [0, 1, 0, 0, 0],
    #                [0, 0, 1, 0, 0],
    #                [0, 0, 0, 1, 0] ])
    for i in range(ncomp):
        A[i] /= np.linalg.norm(A[i])

    for i in range(1, ncomp):
        for j in range(0, i):
            A[i] -= np.dot(A[j], A[i]) * A[j]
            A[i] /= np.linalg.norm(A[i])

    return A # ncomp*nvar, each ROW is unit length

def solve_weighted(A, b, w):
    """
    Solve Ax = b with weights w; return x
    
    A : 2D array
    b : 1D array length A.shape[0]
    w : 1D array same length as b
    """
  
    #- Apply weights
    # nvar = len(w)
    # W = dia_matrix((w, 0), shape=(nvar, nvar))
    # bx = A.T.dot( W.dot(b) )
    # Ax = A.T.dot( W.dot(A) )

    b = A.T.dot( w*b )
    A = A.T.dot( (A.T * w).T )

    if isinstance(A, scipy.sparse.spmatrix):
        x = scipy.sparse.linalg.spsolve(A, b)
    else:
        x = np.linalg.lstsq(A, b)[0]
        
    return x

    
#-------------------------------------------------------------------------

def empca(data, weights=None, maxiter=25, ncomp=5, randseed=None, silent=False):
    """
    Iteratively solve data[i] = Sum_j: c[i,j] p[j] using weights
    
    Input:
    - data[nobs, nvar]
    - weights[nobs, nvar]
      
    Optional:
    - maxiter    : maximum number of iterations
    - ncomp     : number of model vectors
    - randseed : random number generator seed; None to not re-initialize
    
    Returns Model object
    """

    #pdb.set_trace()
    
    if weights is None:
        weights = np.ones(data.shape)

    #- Basic dimensions
    nobs, nvar = data.shape
    assert data.shape == weights.shape

    #- degrees of freedom for reduced chi2
    ii = np.where(weights > 0)
    dof = data[ii].size - ncomp*nvar - ncomp*nobs 

    #- Initial M step - starting random guess
    eigvec = _random_orthonormal(ncomp, nvar, seed=randseed)
    # store eigvec into Model object
    model = Model(eigvec, data, weights)
    
    #pdb.set_trace()
    
    if not silent:
        print("       iter        R2             rchi2")
    
    for k in range(maxiter):
        # E step
        model.calc_coeff()
        # M step
        model.calc_eigvec()
        if not silent:
            print('EMPCA %2d/%2d  %15.8f %15.8f' % \
                (k+1, maxiter, model.R2(), model.rchi2()))
            sys.stdout.flush()

    # One last time with latest coefficients
    model.calc_coeff()

    if not silent:
        print("R2:", model.R2())
    
    return model

    
def _main1():
    np.random.seed(1)
    nobs = 100
    nvar = 200
    ncomp = 3
    data = np.zeros(shape=(nobs, nvar))

    #- Generate data
    x = np.linspace(0, 2*np.pi, nvar)
    for i in range(nobs):
        for k in range(ncomp):
            c = np.random.normal()
            data[i] += 5.0*ncomp//(k+1)**2 * c * np.sin(x*(k+1))

    #- Add noise
    sigma = np.ones(shape=data.shape)
    for i in range(nobs//10):
        sigma[i] *= 5
        sigma[i, 0:nvar//4] *= 5

    weights = 1.0 / sigma**2    
    noisy_data = data + np.random.normal(scale=sigma)

    print("Testing empca")
    m0 = empca(noisy_data, weights, maxiter=20)
    print("done")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
    
def _main2():
    # We need to use 120.0 to force double-precision array instead
    # of integer
    B1 = np.array([
        [50, 67, 90, 98, 120.0],
        [55, 71, 93, 102, 129],
        [65, 76, 95, 105, 134],
        [50, 80, 102, 130, 138],
        [60, 82, 97, 135, 151],
        [65, 89, 106, 137, 153],
        [75, 95, 117, 133, 155] ])
    # Manually centered and scaled in R since I don't know how to here
    B1 = np.array([
        [-1.0954451, -1.3268069, -1.0825318, -1.2645627, -1.4934803],
        [-0.5477226, -0.9185587, -0.7577722, -1.0346422, -0.8214141],
        [ 0.5477226, -0.4082483, -0.5412659, -0.8622019, -0.4480441],
        [-1.0954451,  0.0000000,  0.2165064,  0.5748012, -0.1493480],
        [ 0.0000000,  0.2041241, -0.3247595,  0.8622019,  0.8214141],
        [ 0.5477226,  0.9185587,  0.6495191,  0.9771621,  0.9707622],
        [ 1.6431677,  1.5309311,  1.8403040,  0.7472416,  1.1201102] ])
    B1wt = np.array([
        [1, 1, 1, 1, 1.0],
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],                    
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1] ])
    
    # B2 has 0 value/weight in two cells in upper-left corner
    B2 = np.array([
        [0, 67, 90, 98, 120.0],
        [0, 71, 93, 102, 129],
        [65, 76, 95, 105, 134],
        [50, 80, 102, 130, 138],
        [60, 82, 97, 135, 151],
        [65, 89, 106, 137, 153],
        [75, 95, 117, 133, 155] ])
    B2 = np.array([
        [         0, -1.3268069, -1.0825318, -1.2645627, -1.4934803],
        [         0, -0.9185587, -0.7577722, -1.0346422, -0.8214141],
        [ 0.2201928, -0.4082483, -0.5412659, -0.8622019, -0.4480441],
        [-1.4312529,  0.0000000,  0.2165064,  0.5748012, -0.1493480],
        [-0.3302891,  0.2041241, -0.3247595,  0.8622019,  0.8214141],
        [ 0.2201928,  0.9185587,  0.6495191,  0.9771621,  0.9707622],
        [ 1.3211565,  1.5309311,  1.8403040,  0.7472416,  1.1201102]])
    B2wt = np.array([[0, 1, 1, 1, 1.0],
                     [0, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1],                    
                     [1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1],
                     [1, 1, 1, 1, 1] ])
    
    #m1 = empca(B1, B1wt, maxiter=20, ncomp=4)
    m2 = empca(B2, B2wt, maxiter=20, ncomp=4)

    np.set_printoptions(precision=3)

    float_formatter = lambda x: "%.3f" % x
    np.set_printoptions(formatter={'float_kind':float_formatter})

    # print("")
    # print("----- B1 -----")
    # print("Coeff (scores)")
    # print(m1.coeff)
    # print("Eigvec (loadings)")
    # print(m1.eigvec.transpose())
    # print("")
    # print("P'P (includes eigenvalues)")
    # print(np.dot(m1.coeff.transpose(), m1.coeff))
    print("----- B2 -----")
    print("Coeff (scores)")
    print(m2.coeff)
    print("Eigvec (loadings)")
    print(m2.eigvec.transpose())

# ------------------------------------------------------------------------------
    
if __name__ == '__main__':
    _main2()
