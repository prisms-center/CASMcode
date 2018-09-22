from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import random

def eci(individual, coef):
  """ 
  Return ECI as list of tuple containing (index, coef).
  
  Arguments
  ---------
    
    individual: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
    
    coef: array-like of shape (n_selected, 1)
      The non-zero coefficients.
  
  
  Returns
  -------
    
    eci: List[(int, float)]
      A list of tuple containing the basis function index and coefficient value
      for basis functions with non-zero coefficients: [(index, coef), ...]
  
  """
  return list(zip(indices(individual), list(coef)))
  
  
def indices(individual):
  """ 
  Convert list of bool to list of indices of True values.
  
  Arguments
  ---------
    
    individual: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
    
  
  Returns
  -------
    
    indices: List[int]
      List of indices of True values.
  
  """
  return [i for i in range(len(individual)) if individual[i] ]


def wHullDist(hull_dist, A=1.0, B=1.0, kT=1.0, **kwargs):
  """
  Returns A*np.exp(-hull_dist/kT) + B
  """
  return A*np.exp(-hull_dist/kT) + B


def wEmin(y, A=1.0, B=1.0, kT=1.0):
  """
  Returns A*np.exp(-(y - emin)/kT) + B, where emin = np.min(y)
  """
  emin = np.min(y)
  return A*np.exp(-(y - emin)/kT) + B


def wEref(y, A=1.0, B=1.0, kT=1.0, E0=0.0):
  """
  Returns A*np.exp(-(y - E0)/kT) + B, where y >= E0; 1.0 otherwise
  """
  w = A*np.exp(-(y - E0)/kT) + B
  w[list(np.where(y < E0)[0])] = 1.0
  return w

def set_sample_weight(sample_weight, y=None, X=None):
  """ 
  Calculate weighted data and weighted target values.
  
  Ordinary least squares minimizes
    (y-X*b).transpose() * (y-X*b)
  
  where 'X' is the correlation matrix of shape (Nvalue, Nbfunc), and 'y' 
  is a vector of Nvalue calculated properties, and 'b' are the fitting 
  coefficients (ECI).

  Weighted least squares minimizes
    (y-X*b).transpose() * W * (y-X*b)
  
  Using the SVD, and given that W is Hermitian:
    U * S * U.transpose() == W
  
  Define L such that:
    L.transpose() = U * sqrt(S)
  
  Then we can write the weighted least squares problem using:
    (y-X*b).transpose() * L.transpose() * L * (y-X*b)
  
  Or:
    (L*y-L*X*b).transpose() * (L*y-L*X*b)
    
  So, if weights are included, then the linear model is changed from
    X*b = y  ->  L*X*b = L*y
    
    
  Arguments
  ---------
    
    sample_weight: None, 1d array-like of shape (n_samples,1), or 2d array-like of shape (n_samples, n_samples)
      Sample weights.
      
      if sample_weight is None: (unweighted)
        W = np.matlib.eye(N) 
      if sample_weight is 1-dimensional:
        W = np.diag(sample_weight)*Nvalue/np.sum(sample_weight) 
      if sample_weight is 2-dimensional (must be Hermitian, positive-definite):
        W = sample_weight*Nvalue/np.sum(sample_weight)
    
    X: array-like of shape (n_samples, n_features), optional
      The training input samples (correlations).
    
    y: array-like of shape: (n_samples, 1), optional
      The target values (property values).
    
  Returns
  -------
    
    (weighted_y, weighted_X, W, L)
    
    weighted_y: None, or array-like of shape: (n_samples, 1)
      If y given as input, the weighted target values, weighted_y = L*y. 
    
    weighted_X: None, or array-like of shape: (n_samples, n_features)
      If X given as input, the weighted training input data, weighted_X = L*x.
    
    W: array-like of shape: (n_samples, n_samples)
      Contains sample weights. 
    
    L: array-like of shape: (n_samples, n_samples)
      Used to generate weighted_X and weighted_y, W = L * L.transpose(). 
    
      
  Notes
  -----
    
  Returns (weighted_y, weighted_X, W, L) 
  """
  # check sample_weight and convert to square matrix
  W = None
  L = None
  weighted_y = None
  weighted_X = None
  
  if sample_weight is None:
    W = np.identity(y.shape[0])
  elif len(sample_weight.shape) == 1:
    n_samples = len(sample_weight)
    W = np.diag(sample_weight)*n_samples/np.sum(sample_weight)
  elif len(sample_weight.shape) == 2:
    n_samples = len(sample_weight)
    W = sample_weight*n_samples/np.sum(sample_weight)
  else:
    raise Exception("Error in set_sample_weight: sample_weight dimension > 2")
  
  # weighted data
  U, S, V = np.linalg.svd(W)
  L = U.dot(np.diag(np.sqrt(S))).transpose()
  
  if X is not None:
    weighted_X = np.dot(L, X)
  
  if y is not None:
    weighted_y = np.dot(L, y)
  
  return (weighted_y, weighted_X, W, L)
