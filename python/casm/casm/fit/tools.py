import numpy as np

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
  Returns A*np.exp(-(value - E0)/kT) + B, where value > E0; 1.0 otherwise
  """
  w = A*np.exp(-(value - E0)/kT) + B
  w[np.where(value <= E0)[0].tolist()] = 1.0
  return w


def initNRandomOn(Nbfunc, Nbfunc_init):
  """ Initialize a list of length 'Nbfunc' with 'Nbfunc_init' elements 1 and the rest 0""" 
  bits = [False]*Nbfunc
  for i in random.sample(range(Nbfunc), Nbfunc_init):
    bits[i] = True
  return bits


def check_constraints(Nbfunc_min=0, Nbfunc_max="all", FixOn=[], FixOff=[]):
  """ 
  Ensure that number of basis functions is in range [min, max] while also
  ensuring that basis functions in FixOn remain on, and in FixOff remain off.
  """
  def decorator(func):
    def wrapper(*args, **kargs):
      offspring = func(*args, **kargs)
      for child in offspring:
        
        # first ensure FixOn and FixOff
        child[FixOn] = True
        child[FixOff] = False
        
        # now make sure enough are 'on'
        Non = sum(child)
        if Non < min:
          # get all 'off'
          off = np.where(np.array(child) == False)[0].tolist()
          
          # remove those that are FixOff
          candidate = set(off).difference(set(FixOff))
          
          # select enough of the candidates to reach the minimum
          turn_on = random.sample(candidate, Nbfunc_min-Non)
          child[turn_on] = True
        
        # now make sure enough are 'off'
        if Nbfunc_max != "all":
          if Non > Nbfunc_max:
            # get all 'on'
            on = np.where(np.array(child) == True)[0].tolist()
            
            # remove those that are FixOn
            candidate = set(on).difference(set(FixOn))
            
            # select enough of the candidates to reach the maximum
            turn_off = random.sample(candidate, Non-Nbfunc_max)
            child[turn_off] = False
        
      return offspring
    return wrapper
  return decorator

