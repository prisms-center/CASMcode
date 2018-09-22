"""CASM machine learning tools"""
from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

## This part needs to be in global scope for parallization #####################  
from deap import creator
from deap import base
from deap.tools import HallOfFame

# we'll want to minimize a cv score
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))

# each individual is a list of True or False indicating if each basis function should 
# be included in the model
creator.create("Individual", list, fitness=creator.FitnessMin, input=None)

def empty_individual(n_features):
  selected = [False] * n_features
  indiv = creator.Individual(selected)
  return indiv

class EqualIndividual(object):
  """ 
  Functor class for checking equivalence of individuals in a HallOfFame.
    
  Attributes
  ----------
    
    rel_tol: number
      Individuals are considered equal if features selected are equal and 
      abs((A.fitness.values[0] - B.fitness.values[0]) / A.fitness.values[0] ) < rel_tol
    
  """
  
  def __init__(self, rel_tol=1e-6):
    """
    Arguments
    ---------
    
      rel_tol: number, optional, default=1e-6
        Individuals are considered equal if features selected are equal and 
        abs((A.fitness.values[0] - B.fitness.values[0]) / A.fitness.values[0] ) < rel_tol
    """
    self.rel_tol = rel_tol
  
  def __call__(self, A, B):
    """
    Arguments
    ---------
    
      A, B: Individuals to be compared
        Individuals are considered equal if features selected are equal and 
        abs((A.fitness.values[0] - B.fitness.values[0]) / A.fitness.values[0] ) < rel_tol
    
    Returns
    ---------
      result: boolean
        True if A and B compare equal, False otherwise
    
    """
    if A == B:
      if abs((A.fitness.values[0] - B.fitness.values[0]) / A.fitness.values[0]) < self.rel_tol:
        return True
    return False

def create_halloffame(maxsize, rel_tol=1e-6):
  return HallOfFame(maxsize, similar=EqualIndividual(rel_tol))

################################################################################  


from casm.learn.fit import example_input_Lasso, example_input_LassoCV, example_input_RFE, \
 example_input_GeneticAlgorithm, example_input_IndividualBestFirst, \
 example_input_PopulationBestFirst, example_input_DirectSelection, \
 open_input, set_input_defaults, \
 FittingData, TrainingData, \
 print_input_help, print_individual, print_population, print_halloffame, print_eci, \
 to_json, open_halloffame, save_halloffame, \
 checkspecs, checkhull

from casm.learn.feature_selection import fit_and_select
from casm.learn.direct_selection import direct_fit

__all__ = __all__ = [
  'create_halloffame',
  'EqualIndividual',
  'empty_individual',
  'example_input_Lasso', 
  'example_input_LassoCV', 
  'example_input_RFE',
  'example_input_GeneticAlgorithm',
  'example_input_IndividualBestFirst',
  'example_input_PopulationBestFirst', 
  'example_input_DirectSelection',
  'open_input',
  'set_input_defaults',
  'FittingData', 
  'TrainingData',
  'print_input_help', 
  'print_individual', 
  'print_population', 
  'print_halloffame', 
  'print_eci',
  'to_json', 
  'open_halloffame', 
  'save_halloffame',
  'checkspecs', 
  'checkhull', 
  'fit_and_select', 
  'direct_fit'
]
