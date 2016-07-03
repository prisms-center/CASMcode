from evolve import *

class DirectSelection(BaseEstimator, SelectorMixin):
  """
  Directly specify which basis functions should be selected in an individual.
  
  The 'indiv' attribute is generated when 'fit' is called.
  
  
  Attributes
  ----------
    
    estimator:  estimator object implementing 'fit'
        The estimator specified by the input settings. Not actually used.
      
    indices: List[int]
      The indices of basis functions to be selected
    
    indiv: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
  """
  
  def __init__(self, 
               estimator, 
               indices=None, 
               bitstring=None):
    """
    Construct using either 'indices' or 'bitstring' to specify which basis
    functions should be selected.
    
    Arguments
    ---------
      
      estimator:  estimator object implementing 'fit'
        The estimator specified by the input settings.
      
      indices: List[int], optional
        The indices of basis functions
      
      bitstring: str, optional
        String consisting of '0' and '1', with '1' corresponding to selected 
        basis functions.
    
    """
    
    self.estimator = estimator
    
    if (indices is None) == (bitstring is None):
      raise Exception("Error constructing DirectSelection: Use either 'indices' "
                      "or 'bitstring' kwarg to specify selected features")
    
    if bitstring is not None:
      indices = [index for index,bit in enumerate(bitstring) if int(bit)]
    self.indices = indices
    
    self.indiv = None
  
    
  def fit(self, X, y):
    """
    Generates self.indiv, an individual of length n_features with basis functions
    selected as specified by self.indices
    
    Arguments
    ---------
    
      X: array-like of shape (n_samples, n_features)
        The input data
      
      y: array-like of shape (n_samples, 1)
        The values
    
    
    Returns
    -------
      
      self: returns an instance of self.
    """
    selected = [False] * X.shape[1]
    for i in self.indices:
      selected[i] = True
    if sum(selected) == 0:
      raise Exception("Error using DirectSelection: No basis functions selected")
    self.indiv = casm.learn.creator.Individual(selected)
    return self
  
  
  def _get_support_mask(self):
    """
    Return the individual specified at construction.
    
    Returns
    -------
      
      support: List[bool] of length n_features
        This is a boolean list of shape [n_features], in which an element is True 
        iff its corresponding feature is selected for retention.
    
    """
    return self.indiv
  

