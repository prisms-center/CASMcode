import sklearn.cross_validation
import numpy as np
from math import sqrt
from casm.learn.linear_model import LinearRegressionForLOOCV
from casm.learn.tools import indices

def LeaveOneOutForLLS(n_samples, **kwargs):
  """ 
  Fake cv generator yields train and test sets including all n_samples.
  
  Arguments
  ---------
    
    n_samples: int
      The number of samples
  
  
  Returns
  -------
    
    cv: [(range(n_samples), range(n_samples))]
  
  """
  return [(range(n_samples), range(n_samples))] 


#def evaluate_loocv_explicit(individual, estimator, X, y, penalty=0.0):
#  """ 
#  Evaluate LOOCV by explicitly fitting N times: 
#  
#    LOOCV = sqrt(sum_i((value_i - pred_(-i))**2)) + penalty*sum(individual)
#  
#  Arguments
#  ---------
#    
#    individual: List[bool] of length n_features
#      This is a boolean list of shape [n_features], in which an element is True 
#      iff its corresponding feature is selected for retention.
#    
#    estimator: estimator object implementing 'fit' and 'predict'
#      The object to use to fit the data.
#    
#    X: array-like of shape (n_samples, n_features)
#      The training input samples (correlations).
#    
#    y: array-like of shape: (n_samples, 1)
#      The target values (property values).
#    
#    penalty: float, optional, default=0.0
#      The CV score is increased by 'penalty*sum(individual)'.
#  
#  
#  Returns
#  -------
#    
#    loocv: float
#      The Leave-One-Out CV score, plus penalty term.
#  
#  """
#  
#  n_samples = X.shape[0]
#  
#  # explicitly calculate LOO
#  loo = sklearn.cross_validation.LeaveOneOut(n_samples)
#  res_sqr = np.zeros(n_samples)
#  for train, test in loo:
#    estimator.fit(X[train,indices(individual)], y[train])
#    t = test[0]
#    res_sqr[t] = (y[test] - estimator.predict(X[test,indices(individual)]))**2
#  cv = sqrt(np.mean(res_sqr))
#  
#  return cv + penalty*sum(individual),
#
#
#def evaluate_lls_loocv(individual, X, y, penalty=0.0):
#  """ 
#  Evaluate the linear least squares LOOCV score using LinearRegressionForLOOCV:
#  
#  Equivalent to:
#    estimator = LinearRegressionForLOOCV()
#    C = X[:,indices(individual)]
#    estimator.fit(C, y)
#    return sqrt(estimator.score(C, y)) + penalty*sum(individual)
#  """
#  estimator = LinearRegressionForLOOCV()
#  C = X[:,indices(individual)]
#  estimator.fit(C, y)
#  return sqrt(estimator.score(C, y)) + penalty*sum(individual)


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
      scorer(estimator, X, y). The parameter for sklearn.cross_validation.cross_val_score,
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
  scores = sklearn.cross_validation.cross_val_score(
    estimator,
    X[:,indices(individual)],
    y=y,
    scoring=scoring,
    cv=cv,
    fit_params=fit_params)
  
  return sqrt(np.mean(scores)) + penalty*sum(individual),

