import numpy as np
import random

def eci(individual, coef):
  """ 
  Return a list of tuple: [(index, coef), ...]
  
  index: a list containing the index of non-zero eci
  coef: the eci value corresponding to the index
  
  """
  return zip(indices(individual), coef.tolist())
  
  
def indices(individual):
  """ convert bitstring to indices for slicing """
  indices = []
  for i in range(len(individual)):
    if individual[i]:
      indices.append(i)
  return indices


def wHullDist(hull_dist, A=1.0, B=1.0, kT=1.0):
  """
  Returns A*np.exp(-hull_dist/kT) + B
  """
  return A*np.exp(-hull_dist/kT) + B


def wEmin(value, A=1.0, B=1.0, kT=1.0):
  """
  Returns A*np.exp(-(value - emin)/kT) + B, where emin = np.min(value)
  """
  emin = np.min(value)
  return A*np.exp(-(value - emin)/kT) + B


def wEref(value, A=1.0, B=1.0, kT=1.0, E0=0.0):
  """
  Returns A*np.exp(-(value - E0)/kT) + B, where value >= E0; 1.0 otherwise
  """
  w = A*np.exp(-(value - E0)/kT) + B
  w[np.where(value < E0)[0].tolist()] = 1.0
  return w

def set_sample_weight(sample_weight, value=None, corr=None):
  """ Returns (wvalue, wcorr, W, L) """
  # check sample_weight and convert to square matrix
  Nvalue = len(sample_weight)
  W = None
  L = None
  wvalue = None
  wcorr = None
  
  if sample_weight is None:
    W = np.identity(self.value.shape[0])
  elif len(sample_weight.shape) == 1:
    W = np.diag(sample_weight)*Nvalue/np.sum(sample_weight)
  elif len(sample_weight.shape) == 2:
    W = sample_weight*Nvalue/np.sum(sample_weight)
  else:
    raise Exception("Error in set_sample_weight: sample_weight dimension > 2")
  
  # weighted data
  L = np.linalg.cholesky(W)
  
  if corr is not None:
    wcorr = np.dot(L, corr)
  
  if value is not None:
    wvalue = np.dot(L, value)
  
  return (wvalue, wcorr, W, L)


def initNRandomOn(container, Nbfunc, Nbfunc_init):
  """ Initialize a list of length 'Nbfunc' with 'Nbfunc_init' elements 1 and the rest 0""" 
  bits = [False]*Nbfunc
  for i in random.sample(range(Nbfunc), Nbfunc_init):
    bits[i] = True
  return container(bits)


class Constraints(object):
  def __init__(self, Nbfunc_min=1, Nbfunc_max="all", FixOn=np.array([], dtype=int), FixOff=np.array([], dtype=int)):
    if Nbfunc_min < 1:
      raise ValueError("Nbfunc_min must be >= 1")
    self.Nbfunc_min = Nbfunc_min
    self.Nbfunc_max = Nbfunc_max
    self.FixOn = np.array(FixOn, dtype=int)
    self.FixOff = np.array(FixOff, dtype=int)
    

def check_constraints(constraints):
  """ 
  Ensure that number of basis functions is in range [min, max] while also
  ensuring that basis functions in FixOn remain on, and in FixOff remain off.
  """
  def decorator(func):
    def wrapper(*args, **kargs):
      offspring = func(*args, **kargs)
      for child in offspring:
        
        # first ensure FixOn and FixOff
        for i in constraints.FixOn:
          child[i] = True
        for i in constraints.FixOff:
          child[i] = False
        
        # now make sure enough are 'on'
        Non = sum(child)
        if Non < constraints.Nbfunc_min:
          # get all 'off'
          off = np.where(np.asarray(child) == False)[0]
          
          # remove those that are FixOff
          candidate = set(off).difference(set(constraints.FixOff))
          
          # select enough of the candidates to reach the minimum
          turn_on = random.sample(candidate, constraints.Nbfunc_min-Non)
          for i in turn_on:
            child[i] = True
        
        # now make sure enough are 'off'
        if constraints.Nbfunc_max != "all":
          if Non > constraints.Nbfunc_max:
            # get all 'on'
            on = np.where(np.asarray(child) == True)[0]
            
            # remove those that are FixOn
            candidate = set(on).difference(set(constraints.FixOn))
            
            # select enough of the candidates to reach the maximum
            turn_off = random.sample(candidate, Non-constraints.Nbfunc_max)
            for i in turn_off:
              child[i] = False
        
      return offspring
    return wrapper
  return decorator

