from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import sklearn.model_selection
import numpy as np
import pickle
from math import sqrt
from casm.learn.linear_model import LinearRegressionForLOOCV
from casm.learn.tools import indices

def LeaveOneOutForLLS(n_samples):
  """ 
  Returns train and test sets each including all n_samples.
  
  Arguments
  ---------
    
    n_samples: int
      The number of samples
  
  
  Returns
  -------
    
    cv: [(range(n_samples), range(n_samples))]
  
  """
  return [(range(n_samples), range(n_samples))] 


def cvCustom(filename=None):
  """
  Read CV generator or train / test sets from a pickle file.
  """
  if filename == None:
    raise Exception("Error using cvCustom: filename is None")
  
  with open(filename, 'rb') as f:
    return pickle.load(f)
  


def cross_val_score(estimator, X, individual, y=None, scoring=None, cv=None, penalty=0.0, fit_params=None):
  """
  Evaluate CV score for a particular individual.
  
  Arguments
  ---------
    
    estimator: estimator object implementing 'fit' and 'predict'
      The object to use to fit the data.
    
    X: array-like of shape (n_samples, n_features)
      The training input samples (correlations).
    
    individual: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
    
    y: array-like of shape: (n_samples, 1)
      The target values (property values).
    
    scoring: string, callable or None, optional, default=None
      A string or a scorer callable object / function with signature 
      scorer(estimator, X, y). The parameter for sklearn.model_valdiation.cross_val_score,
      default = None, uses estimator.score().
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
      
    penalty: float, optional, default=0.0
      The CV score is increased by 'penalty*sum(individual)'.
    
    fit_params: dict, optional,
      Parameters to pass to the fit method of the estimator.
  
  
  Returns
  -------
    
    score: float
      The cross validation score
  
  
  Notes
  -----
  
  Equivalent to:
  
    scores = cross_val_score(
      estimator,
      X[:,indices(individual)],
      y=y,
      scoring=scoring,
      cv=cv,
      fit_params=fit_params)
    return sqrt(np.mean(scores)) + penalty*sum(individual),
  """
  scores = sklearn.model_selection.cross_val_score(
    estimator,
    X[:,indices(individual)],
    y=y,
    scoring=scoring,
    cv=cv,
    fit_params=fit_params)
  return sqrt(np.mean(scores)) + penalty*sum(individual),
  
