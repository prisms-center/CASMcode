import numpy as np
import sklearn

class LinearRegressionForLOOCV(sklearn.base.BaseEstimator):
  """ LinearRegression estimator that fits by constructing hat matrix for faster
      LOOCV score calculations 
  """
  
  def __init__(self, penalty=0.0, **kwargs):
    self.penalty=penalty
    pass
  
  def fit(self, X, y):
    """
    Calculates linear least squares solution for:
      X*b = y,
    
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
    
    Calculates attributes:
      coef_: numpy array, length matches X.shape[1]
      H_: numpy matrix, shape: (X.shape[0], X.shape[0])
        where H = X*(X.transpose()*X).inverse()*X.transpose()
    
    """
    #print "begin casm.fit.linear_model.LinearRegressionForLOOCV.fit"
    S = np.linalg.pinv(X.transpose().dot(X)).dot(X.transpose())
    self.coef_ = np.dot(S, y)
    
    # stores results for LOOCV formula
    self.H_ = np.dot(X, S)
    
  
  def predict(self, X):
    #print "begin casm.fit.linear_model.LinearRegressionForLOOCV.predict"
    return np.dot(X, self.coef_)
  
  
  def score(self, X, y):
    """LOOCV score"""
    #print "begin casm.fit.linear_model.LinearRegressionForLOOCV.score"
    y_pred = np.dot(self.H_, y)
    return np.mean(((y - y_pred)/(1.0 - np.diag(self.H_)))**2)
