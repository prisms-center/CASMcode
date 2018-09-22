from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import copy
from casm.learn import empty_individual
from casm.learn.model_selection import cross_val_score
from casm.learn.fit import open_halloffame, make_fitting_data, make_estimator, \
  add_individual_detail, print_halloffame


def make_population(n_features, input_options):
  """
  Read direct selection population options and use to construct a population of
  individuals to be fit.
  
  Arguments
  ---------
    
    n_features: int
      The size of the individuals to be constructed
    
    input_options: iterable of dict
      Each dict contains options specifying an individual or individuals to be
      included in the population. Three forms are recognized:
      
      bitstring:
        Ex.: {"bitstring" : "01110001100..."}
        The string may have fewer than n_features digits, but not more. If fewer,
        the rest are assumed to be '0'.
      
      indices:
        Ex.: {"indices" : [1, 2, 3, 7, 8]}
        Specify which features should be 'on' by index.
      
      from_halloffame:
        Ex.: {"from_halloffame" : "my_halloffame.pkl", "individuals" : [0, 2, 5]}
        Specifies a hall of fame and particular individuals (by index) in the hall
        of fame to include in the population. "individuals" is optional, with the
        default behaviour including all individuals in the hall of fame.
  
  Returns
  -------
    
    population: List[individuals]
      A list of the individuals to be fit
    
  """
  population = []
  for opt in input_options:
    
    if "bitstring" in opt:
      indiv = empty_individual(n_features)
      for index, bit in enumerate(opt["bitstring"]):
        if int(bit):
          indiv[index] = True
      if sum(indiv) == 0:
        raise Exception("Error using making individual from bitstring: No basis functions selected")
      population.append(indiv)
    
    elif "indices" in opt:
      indiv = empty_individual(n_features)
      for index in opt["indices"]:
        indiv[index] = True
      if sum(indiv) == 0:
        raise Exception("Error using making individual from bitstring: No basis functions selected")
      population.append(indiv)
    
    elif "from_halloffame" in opt:
      hall = open_halloffame(opt["from_halloffame"])
      if "individuals" in opt:
        for index in opt["individuals"]:
          population.append(hall[index])
      else:
        for indiv in hall:
          population.append(indiv)
    
    else:
      print("Unrecognized option:\n")
      print(opt)
      raise Exception("Error making DirectSelection population: unrecognized option")
  
  return population


def direct_fit(input, save=True, verbose=True, read_existing=True, hall=None):
  """
  Fit ECI and add details for a set of individuals specified for feature_selection
  method 'DirectSelection' via the 'population' kwarg.
  
  Arguments
  ---------
    
    input: dict
      The input settings as a dict
    
    save: boolean, optional, default=True
      Save a pickle file containing the training data and scoring metric. The file
      name, which can be specified by input["fit_data_filename"], defaults to "fit_data.pkl".
    
    verbose: boolean, optional, default=True
      Print information to stdout.
    
    read_existing: boolean, optional, default=True
      If it exists, read the pickle file containing the training data and scoring 
      metric. The file name, which can be specified by input["fit_data_filename"], 
      defaults to "fit_data.pkl".
    
    hall: deap.tools.HallOfFame, optional, default=None
      A Hall Of Fame to add resulting individuals to
    
  
  Returns
  -------
    
    population: List[individuals]
      A list of the individuals that were fit
    
  """
  # construct FittingData
  fdata = make_fitting_data(input, save=True, verbose=verbose, read_existing=True)
  
  kwargs = input["feature_selection"]["kwargs"]
  
  population = make_population(fdata.n_features, kwargs["population"])
  
  for indiv_i, indiv in enumerate(population):
    
    if verbose:
      print("Begin fitting individual", indiv_i, "of", len(population))
    
    _input = copy.deepcopy(input)
    use_saved_estimator = kwargs.get("use_saved_estimator", False)
    
    # check for 'use_saved_estimator'
    if use_saved_estimator and getattr(indiv, "input", None) is not None:
      _input["estimator"] = indiv.input["estimator"]
    
    estimator = make_estimator(_input)
    
    indiv.fitness.values = cross_val_score(estimator, fdata.weighted_X, indiv, 
      y=fdata.weighted_y, scoring=fdata.scoring, cv=fdata.cv, penalty=fdata.penalty)
    add_individual_detail(indiv, estimator, fdata, _input, selector=None)
    
    if verbose:
      print("  DONE\n")
  
  if verbose:
    print("Result:")
    print_halloffame(population);
  
  if hall is not None:
    hall.update(population)
    
  return population

