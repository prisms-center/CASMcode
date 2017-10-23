from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import numpy as np
import sklearn

class LinearRegressionForLOOCV(sklearn.base.BaseEstimator):
  """ 
  LinearRegression estimator that fits by constructing hat matrix for faster
  LOOCV score calculations.
  
  Attributes
  ----------
    
    pinv: boolean, optional, default=True
      If true, use the psuedo-inverse in solving the least squares problem.
    
    coef_: array-like of shape (n_samples, 1)
      Estimated coefficients for the linear regression problem.
        
    H_: array-like of shape (n_samples, n_features)
      H = X*(X.transpose()*X).inverse()*X.transpose()
  """
  
  def __init__(self, pinv=True, **kwargs):
    """
    Calculates linear least squares solution for X*b = y.
    
    Arguments
    ---------
    
      pinv: boolean, optional, default=True
        If true, use the psuedo-inverse in solving the least squares problem.
      
      kwargs: keyword args
        Other keyword args are ignored
    """
    self.pinv = pinv
  
  
  def fit(self, X, y):
    """
    Calculates linear least squares solution for X*b = y.
    
    
    Arguments
    ---------
    
      X: array-like of shape (n_samples, n_features)
        The input data
      
      y: array-like of shape (n_samples, 1)
        The values
      
    
    Notes
    -----
    
      Derivation:
        X.t*X*b = X.t*y, orthogonal projection
        b = (X.t*X).inv * X.t * y, solves the least squares problem
        b = S * y,
      where
        S = ( X.t * X ).inv * X.t
      
      and 
        y_pred = X * b = H * y,
      where
        H = X * S
      
      Calculates attributes: coef_, H_
    """
    if self.pinv == False:
      S = np.linalg.inv(X.transpose().dot(X)).dot(X.transpose())
    else:
      S = np.linalg.pinv(X.transpose().dot(X)).dot(X.transpose())
    
    # coef
    self.coef_ = np.dot(S, y)
    
    # stores results for LOOCV formula
    self.H_ = np.dot(X, S)
    
  
  def predict(self, X):
    """
    Predict using the linear model.
    
    
    Arguments
    ---------
    
      X: array-like of shape (n_samples, n_features)
        The input data
    
    
    Returns
    -----
    
      y_pred: array-like of shape (n_samples, 1)
        The predicted values
    """
    return np.dot(X, self.coef_)
  
  
  def score(self, X, y):
    """
    Return the LOOCV score using the linear model.
    
    
    Arguments
    ---------
    
      X: array-like of shape (n_samples, n_features)
        The input data. (This parameter is ignored)
      
      y: array-like of shape (n_samples, 1)
        The values
      
    
    Returns
    -----
    
      loocv: float
        The LOOCV score, np.mean(((y - y_pred)/(1.0 - np.diag(self.H_)))**2)
    
    Notes
    -----
      Must already be fit. 'X' parameter is ignored.
      
    """
    y_pred = np.dot(self.H_, y)
    return np.mean(((y - y_pred)/(1.0 - np.diag(self.H_)))**2)
