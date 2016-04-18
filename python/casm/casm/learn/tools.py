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
  return zip(indices(individual), list(coef))
  
  
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
  return [i for i in xrange(len(individual)) if individual[i] ]


def wHullDist(hull_dist, A=1.0, B=1.0, kT=1.0):
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
  
  Uses sample weights to calculate
  
    L*X * b = L*y 
  
  a weighted linear model where the weights are given by W = L * L.transpose().
  
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
    
    weighted_X: None, or array-like of shape: (n_samples, n_features)
      If X given as input, the weighted training input data, weighted_X = L*x.
    
    weighted_y: None, or array-like of shape: (n_samples, 1)
      If y given as input, the weighted target values, weighted_y = L*y. 
    
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
  L = np.linalg.cholesky(W)
  
  if X is not None:
    weighted_X = np.dot(L, X)
  
  if y is not None:
    weighted_y = np.dot(L, y)
  
  return (weighted_y, weighted_X, W, L)


def initNRandomOn(container, n_features, n_features_init):
  """ 
  Initialize a random list of bool.
  
  Arguments
  ---------
    
    container: Type
      The type of container to return
    
    n_features: int
      The total length of the container
    
    n_features_init: int
      The number of True elements in the container
  
  Returns
  -------
    
    container_: container
      A container of the specified type with 'n_features' True and the rest False.
  
  """ 
  bits = [False]*n_features
  for i in random.sample(range(n_features), n_features_init):
    bits[i] = True
  return container(bits)


class Constraints(object):
  """ 
  Holds constraints on individuals in evolutionary algorithms.
  
  Attributes
  ----------
    
    n_features_min: int
      The minimum allowed number of selected features. Must be >=1.
    
    n_features_max: int or str
      The maximum allowed number of selected features. String "all" for no limit.
    
    fix_on: 1d array-like of int
      The indices of features to fix on
    
    fix_off: 1d array-like of int
      The indices of features to fix off
  
  """ 
  def __init__(self, n_features_min=1, n_features_max="all", fix_on=np.array([], dtype=int), fix_off=np.array([], dtype=int)):
    """ 
    Arguments
    ---------
      
      n_features_min: int, optional, default=1 
        The minimum allowed number of selected features. Must be >=1.
      
      n_features_max: int or str, optional, default "all".
        The maximum allowed number of selected features. String "all" for no limit.
      
      fix_on: 1d array-like of int, optional, default=np.array([], dtype=int)
        The indices of features to fix on
      
      fix_off: 1d array-like of int, optional, default=np.array([], dtype=int)
        The indices of features to fix off
    
    """ 
    if n_features_min < 1:
      raise ValueError("n_features_min must be >= 1")
    self.n_features_min = n_features_min
    self.n_features_max = n_features_max
    self.fix_on = np.array(fix_on, dtype=int)
    self.fix_off = np.array(fix_off, dtype=int)
    

def check_constraints(constraints):
  """ 
  Returns a decorator that ensures that individuals obey constraints.
  
  The decorator will first select on/off any features that are constrained by
  'fix_on' or 'fix_off'. Then, if the number of features selected is less than
  'n_features_min', features will be randomly turned on. Finally, if the number
  of features selected is less than 'n_features_max', features will be randomly
  turned off.
  
  Arguments
  ---------
    constraints: casm.learn.Constraints
      A Constraints object containing the requested constraints
  
  
  Returns
  -------
    
    decorator: function
      A decorator that can be used to wrap evolutionary algorithm methods to 
      ensure that constraints remain satisfied.
  """
  def decorator(func):
    def wrapper(*args, **kargs):
      offspring = func(*args, **kargs)
      for child in offspring:
        
        # first ensure fix_on and fix_off
        for i in constraints.fix_on:
          child[i] = True
        for i in constraints.fix_off:
          child[i] = False
        
        # now make sure enough are 'on'
        Non = sum(child)
        if Non < constraints.n_features_min:
          # get all 'off'
          off = np.where(np.asarray(child) == False)[0]
          
          # remove those that are fix_off
          candidate = set(off).difference(set(constraints.fix_off))
          
          # select enough of the candidates to reach the minimum
          turn_on = random.sample(candidate, constraints.n_features_min-Non)
          for i in turn_on:
            child[i] = True
        
        # now make sure enough are 'off'
        if constraints.n_features_max != "all":
          if Non > constraints.n_features_max:
            # get all 'on'
            on = np.where(np.asarray(child) == True)[0]
            
            # remove those that are fix_on
            candidate = set(on).difference(set(constraints.fix_on))
            
            # select enough of the candidates to reach the maximum
            turn_off = random.sample(candidate, Non-constraints.n_features_max)
            for i in turn_off:
              child[i] = False
        
      return offspring
    return wrapper
  return decorator

