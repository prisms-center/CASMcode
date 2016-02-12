import sklearn.cross_validation
import numpy as np
from math import sqrt
from casm.fit.linear_model import LinearRegressionForLOOCV
from casm.fit.tools import indices

def LeaveOneOutForLLS(n, **kwargs):
  """ Fake cv generator yields train and test sets including all n data points """
  return [(range(n), range(n))] 


def evaluate_loocv_explicit(individual, estimator, X, y, penalty=0.0):
  """
  Evaluate LOOCV by explicitly solving N least squares problems: 
  
    LOOCV = sqrt(sum_i((value_i - pred_(-i))**2)/sum_i(w_i)) + penalty*sum(individual)
  """
  
  Nvalue = X.shape[0]
  
  # explicitly calculate LOO
  loo = sklearn.cross_validation.LeaveOneOut(Nvalue)
  res_sqr = np.zeros(Nvalue)
  for train, test in loo:
    estimator.fit(X[train,indices(individual)], y[train])
    t = test[0]
    res_sqr[t] = (y[test] - estimator.predict(X[test,indices(individual)]))**2
  cv = sqrt(np.mean(res_sqr))
  
  return cv + penalty*sum(individual),


def evaluate_lls_loocv(individual, X, y, penalty=0.0):
  """ 
  Evaluate the linear least squares LOOCV score using LinearRegressionForLOOCV:
  
  Equivalent to:
    estimator = LinearRegressionForLOOCV()
    C = X[:,indices(individual)]
    estimator.fit(C, y)
    return sqrt(estimator.score(C, y)) + penalty*sum(individual)
  """
  estimator = LinearRegressionForLOOCV()
  C = X[:,indices(individual)]
  estimator.fit(C, y)
  return sqrt(estimator.score(C, y)) + penalty*sum(individual)


def cross_val_score(estimator, X, individual, y=None, scoring=None, cv=None, penalty=0.0, fit_params=None):
  """ 
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
  #print "begin casm.fit.cross_validation.cross_val_score"
  scores = sklearn.cross_validation.cross_val_score(
    estimator,
    X[:,indices(individual)],
    y=y,
    scoring=scoring,
    cv=cv,
    fit_params=fit_params)
  
#  print "scores:", scores
#  print "np.mean(scores):", np.mean(scores)
#  print "sqrt(np.mean(scores)):", sqrt(np.mean(scores))
#  print "sum(individual):", sum(individual)
#  print "penalty*sum(individual):", penalty*sum(individual)
#  print "total:", sqrt(np.mean(scores)) + penalty*sum(individual)
  
  return sqrt(np.mean(scores)) + penalty*sum(individual),

