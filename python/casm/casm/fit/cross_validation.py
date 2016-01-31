import sklearn.cross_validation
from math import sqrt
from casm.fit.linear_model import LinearRegressionForLOOCV


def LeaveOneOutForLLS(n):
  """ Fake cv generator yields train and test sets including all n data points """
  return [(range(n), range(n))] 


def evaluate_loocv_explicit(individual, fdata, model):
  """
  Evaluate weighted LOOCV by explicitly solving N least squares problems 
  
  weighted LOOCV = sqrt(sum_i((value_i - pred_(-i))**2)/sum_i(w_i))
  """
  
  # explicitly calculate LOO
  loo = sklearn.cross_validation.LeaveOneOut(fdata.Nvalue)
  res_sqr = np.zeros(fdata.Nvalue)
  for train, test in loo:
    model.fit(fdata.corr[train,indices(individual)], fdata.value[train])
    t = test[0]
    res_sqr[t] = (fdata.value[test] - model.predict(fdata.corr[test,indices(individual)]))**2
  cv = sqrt(np.mean(res_sqr))
  
  return cv + fdata.penalty*sum(individual),


def evaluate_lls_loocv(individual, fdata):
  """ 
  Evaluate the linear least squares LOOCV score using LinearRegressionForLOOCV
  """
  model = LinearRegressionForLOOCV()
  C = fdata.corr[:,indices(individual)]
  model.fit(C, fdata.value)
  return sqrt(model.score(C, fdata.value)) + fdata.penalty*sum(individual)


def evaluate_lls_wloocv(individual, fdata):
  """ 
  Evaluate the weighted linear least squares LOOCV score using LinearRegressionForLOOCV
  """
  model = LinearRegressionForLOOCV()
  C = fdata.wcorr[:,indices(individual)]
  model.fit(C, fdata.wvalue)
  return sqrt(model.score(C, fdata.wvalue)) + fdata.penalty*sum(individual)


def evaluate_cv(individual, fdata, model):
  """ 
  Equivalent to:
    scores = cross_val_score(
      mdata.model,
      fdata.corr[:,indices(individual)],
      y=fdata.value,
      scoring=mdata.scoring,
      cv=fdata.cv)
    return sqrt(np.mean(scores)),
  """
  scores = sklearn.cross_validation.cross_val_score(
    model,
    fdata.corr[:,indices(individual)],
    y=fdata.value,
    scoring=fdata.scoring,
    cv=fdata.cv)
  
  return sqrt(np.mean(scores)) + fdata.penalty*sum(individual),


def evaluate_wcv(individual, fdata, model):
  """ 
  Equivalent to:
    scores = cross_val_score(
      mdata.model,
      fdata.wcorr[:,indices(individual)],
      y=fdata.wvalue,
      scoring=mdata.scoring,
      cv=fdata.cv)
    return sqrt(np.mean(scores)),
  """
  scores = sklearn.cross_validation.cross_val_score(
    model,
    fdata.wcorr[:,indices(individual)],
    y=fdata.wvalue,
    scoring=fdata.scoring,
    cv=fdata.cv)
  
  return sqrt(np.mean(scores)) + fdata.penalty*sum(individual),

