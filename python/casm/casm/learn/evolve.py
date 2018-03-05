from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import copy
from operator import attrgetter
import os
import pickle
import random
import time

import deap
import deap.tools
from deap.tools import HallOfFame
import deap.algorithms
import numpy as np
import six
from sklearn.base import BaseEstimator
from sklearn.feature_selection.base import SelectorMixin

import casm.learn

def initNRandomOn(container, n_features, n_features_init):
  """ 
  Initialize a random container of bool.
  
  Arguments
  ---------
    
    container: Type
      The type of container to return. Expects container(List) to be valid.
    
    n_features: int
      The total length of the container
    
    n_features_init: int
      The number of True elements in the container
  
  Returns
  -------
    
    result: container
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
    if n_features_max != "all":
      n_features_max = int(n_features_max)
    self.n_features_min = n_features_min
    self.n_features_max = n_features_max
    self.fix_on = np.array(fix_on, dtype=int)
    self.fix_off = np.array(fix_off, dtype=int)
  
  def check(self, indiv):
    # first ensure fix_on and fix_off
    for i in self.fix_on:
      if not indiv[i]:
        return False
    
    for i in self.fix_off:
      if indiv[i]:
        return False
        
    Non = sum(indiv)
    if Non < self.n_features_min:
      return False
    if Non > self.n_features_max:
      return False
    
    return True
    

def restrict_constraints(constraints):
  """
  Returns a decorator that removes offspring that don't obey constraints.
  
  The decorator does:
    offspring = func(*args, **kargs)
    return [child for child in offspring if constraints.check(child)]
  
  Arguments
  ---------
    constraints: casm.learn.evolve.Constraints
      A Constraints object containing the requested constraints
  
  
  Returns
  -------
    
    decorator: function
      A decorator that can be used to wrap evolutionary algorithm methods that 
      return a population to ensure that constraints remain satisfied.
  """
  def decorator(func):
    def wrapper(*args, **kargs):
      offspring = func(*args, **kargs)
      return [child for child in offspring if constraints.check(child)]
    return wrapper
  return decorator
    

def enforce_constraints(constraints):
  """ 
  Returns a decorator that ensures that individuals obey constraints.
  
  The decorator will first select on/off any features that are constrained by
  'fix_on' or 'fix_off'. Then, if the number of features selected is less than
  'n_features_min', features will be randomly turned on. Finally, if the number
  of features selected is less than 'n_features_max', features will be randomly
  turned off.
  
  Arguments
  ---------
    constraints: casm.learn.evolve.Constraints
      A Constraints object containing the requested constraints
  
  
  Returns
  -------
    
    decorator: function
      A decorator that can be used to wrap evolutionary algorithm methods that 
      return a population to ensure that constraints remain satisfied.
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


def initialize_population(n_population, toolbox, filename=None, verbose=True):
  """
  if filename is not None and os.path.exists(filename):
    load initial population
    # may be List[Individual] or HallOfFame, which is converted to List[Individual]
  else:
    create random initial population of size n_population via toolbox.population
  """
  # load or construct initiali population
  
  if filename is not None and os.path.exists(filename):
    if verbose:
      print("Loading initial population:", filename)
    with open(filename, 'rb') as f:
      pop = pickle.load(f)
    if isinstance(pop, HallOfFame):
      # convert to List
      pop = [indiv for indiv in pop]
  else:
    if verbose:
      print("Constructing initial population")
    pop = toolbox.population(n_population)
  
  return pop


def save_population(pop, filename, verbose=False):
  if verbose:
    print("\nPickling population to:", filename)
  with open(filename, 'wb') as f:
      pickle.dump(pop, f, protocol=2)


def initialize_halloffame(filename=None, n_halloffame=25, verbose=False):
  hall = casm.learn.create_halloffame(n_halloffame)
  if verbose:
    print("# Hall of Fame size:", n_halloffame, "\n")
  
  if os.path.exists(filename):
    with open(filename, 'rb') as f:
      existing_hall = pickle.load(f)
    if verbose:
      print("Loading Hall of Fame:", filename)
    hall.update(existing_hall)
  return hall


def save_halloffame(hall, filename, verbose=False):
  if verbose:
    print("\nPickling Hall of Fame to:", filename)
  with open(filename, 'wb') as f:
      pickle.dump(hall, f, protocol=2)


def evaluate_all(pop, toolbox):
  """
  Evaluate and set fitness of all individual's with invalid fitness.
  
  Arguments
  ---------
    
    pop: iterable of individual
      Population to be evaluated
    
    toolbox: deap.base.Toolbox
      Contains methods used by deap during evolution. Expected to contain:
        
        toolbox.evaluate(indiv): To evaluate an individual's fitness
        
        toolbox.map(func, List[individual]): To map function executions
  
  Returns
  -------
    
    nevals: int
      Number of evaluations performed
  
  """
  # evaluate initial fitness
  invalid_ind = [ind for ind in pop if not ind.fitness.valid]
  fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
  for ind, fit in zip(invalid_ind, fitnesses):
    ind.fitness.values = fit
  return len(invalid_ind)
  

class EvolutionaryParams(object):
  """
  Holds parameters used by evolutionary algorithms.
  
  Attributes
  ----------
    
    n_population: int
      Population size.
    
    n_generation: int
      Number of generations to for each repitition. Results are saved between
      repetitions.
    
    n_repetition: int
      Number of repititions to perform.
    
    n_features_init: int
      Number of randomly selected features to initialize each individual with.
  
    pop_begin_filename: string
      Filename where the initial population is read from, if it exists.
    
    pop_end_filename: string
      Filename where the final population is saved.
    
    halloffame_filename: string
      Filename where a hall of fame is saved.
    
    n_halloffame: int
      Number of individuals to save in the hall of fame
  """
  
  def __init__(self, n_population=100, n_generation=10, n_repetition=100, n_features_init=1, 
               pop_begin_filename = "population_begin.pkl",
               pop_end_filename = "population_end.pkl",
               halloffame_filename = "evolve_halloffame.pkl",
               filename_prefix = "",
               n_halloffame = 25):
    """
    Arguments
    ---------
    
      n_population: integer, optional, default=100
        Population size.
      
      n_generation: integer, optional, default=10
        Number of generations to for each repitition. Results are saved between
        repetitions.
      
      n_repetition: integer, optional, default=100
        Number of repititions to perform.
      
      n_features_init: int, optional, default=5
        Number of randomly selected features to initialize each individual with.
      
      pop_begin_filename: string, optional, default="population_begin.pkl"
        Filename where the initial population is read from, if it exists.
      
      pop_end_filename: string, optional, default="population_end.pkl"
        Filename where the final population is saved.
      
      halloffame_filename: string, optional, default="halloffame.pkl"
        Filename where a hall of fame is saved.
      
      filename_prefix: string
        Prefix for filenames, typically taken from input file filename excluding 
        extension.
    
      n_halloffame: int, optional, default=25
        Number of individuals to save in the hall of fame
  
    """
    self.n_population = n_population
    self.n_generation = n_generation
    self.n_repetition = n_repetition
    
    self.n_features_init = n_features_init
    
    if len(filename_prefix) != 0:
      filename_prefix += "_"
    
    self.pop_begin_filename = filename_prefix + pop_begin_filename
    self.pop_end_filename = filename_prefix + pop_end_filename
    
    self.halloffame_filename = filename_prefix + halloffame_filename
    self.n_halloffame = n_halloffame


def default_stats(funcs=None):
  """
  Returns default deap.tools.Statistics object logging 'avg', 'std', 'min', 'max'.
  
  Arguments
  ---------
    
    funcs: dict, optional, default=(see below)
      A dict of functions to include in the deap.tools.Statistics object. 
      
      Default:  
        {
          "avg":np.mean, 
          "std":np.std, 
          "min":np.min, 
          "max":np.max
        }
        
  
  Returns
  -------
    stats: deap.tools.Statistics
      The statistics to be logged.
  
  """
  if funcs is None:
    funcs = {"avg":np.mean, "std":np.std, "min":np.min, "max":np.max}
  stats = deap.tools.Statistics(key=lambda ind: ind.fitness.values)
  for key, val in six.iteritems(funcs):
    stats.register(key, val)
  return stats
    

class Log(object):
  """
  Attributes
  ----------
    
    stats: deap.tools.Statistics
        The statistics to be logged.
    
    logbook: deap.tools.Logbook
    
  """
  
  def __init__(self, stats=default_stats()):
    """
    Arguments
    ---------
      
      stats: deap.tools.Statistics, optional, default=default_stats()
        The statistics to be logged.
    """
    self.stats = stats
    self.logbook = deap.tools.Logbook()
    self.logbook.header = ['gen', 'nevals'] + (self.stats.fields if self.stats else [])
  
  def record(self, pop, gen, nevals, verbose=False):
    """
    Record stats for generation in logbook.
    
    Arguments
    ---------
      
      pop: iterable of individual
        Population to be minimized
      
      gen: integer
        Current generation index
      
      nevals: integer
        Number of times 'evaluate' has been called
      
      verbose: boolean
        Print information to stdout. If true, will print logbook.stream.
    """
    record = self.stats.compile(pop) if self.stats else {}
    self.logbook.record(gen=gen, nevals=nevals, **record)
    if verbose:
      print(self.logbook.stream)


def single_flip_children(parent):
  """
  Return all children of parent that differ by one feature on/off.
  
  Arguments
  ---------
    
    parent: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
  
      
  Returns
  -------
    
    offspring: iterable of individual
      Population to be evaluated
   
  """
  offspring = []
  for index in range(len(parent)):
    child = copy.deepcopy(parent)
    del child.fitness.values
    child[index] = not child[index]
    offspring.append(child)
  return offspring


def best_child(indiv, toolbox):
  """
  Return best child of a particular individual.
  
  Arguments
  ---------
    
    indiv: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
    
    toolbox: deap.base.Toolbox
      Contains methods used by deap during evolution. Expected to contain:
      
        toolbox.children(parent): Returns offspring (List[individual]) of parent
        
        toolbox.evaluate(indiv): To evaluate an individual's fitness
        
        toolbox.map(func, List[individual]): To map function executions
  
  Returns
  -------
    
    (best, nevals)
    
    best: List[bool] of length n_features
      The child with best fitness value
    
    nevals: integer
        Number of times 'evaluate' has been called
  """
  # generate children
  offspring = toolbox.children(indiv)
  nevals = evaluate_all(offspring, toolbox)
  return max(offspring, key=lambda child: child.fitness), nevals


def eaIndividualBestFirst(pop, toolbox, n_generation=10, halloffame=None, stats=default_stats(), verbose=False):
  """
  Evolutionary algorithm that minimizes each individual in isolation by proposing
  children and replacing with the most fit child.
  
  Arguments
  ---------
    
    pop: iterable of individual
      Population to be minimized
    
    toolbox: deap.base.Toolbox
      Contains methods used during evolution. Expected to contain:
        
        toolbox.children(parent): Returns offspring (List[individual]) of parent
        
        toolbox.evaluate(indiv): To evaluate an individual's fitness
        
        toolbox.map(func, List[individual]): To map function executions
    
    n_generation: int, optional, default=10
      Number of generations to run.
    
    halloffame: deap.tools.HallOfFame, optional, default=None
      A Hall Of Fame containing the optimal solutions, as judged by CV score.
    
    stats: deap.tools.Statistics, optional, default=default_stats()
      The statistics to be logged.
    
    verbose: boolean
      Print information to stdout.
    
  
  Returns
  -------
    
    (pop, logbook)
    
    pop: iterable of individual
      The final population
    
    logbook: deap.tools.Logbook 
      Logbook with the statistics of the evolution
  """
  remaining = len(pop)
  gen = 0
  
  # Evaluate & record in logbook
  nevals = evaluate_all(pop, toolbox)
  for indiv in pop:
    indiv.minimized = False
  
  log = Log(stats)
  log.record(pop, gen=gen, nevals=nevals, verbose=verbose)
  
  # initial halloffame update
  if halloffame is not None:
    halloffame.update(pop)
  
  
  while gen < n_generation:
    
    if remaining == 0:
      pop = initialize_population(len(pop), toolbox, verbose=verbose)
      remaining = len(pop)
      nevals = evaluate_all(pop, toolbox)
      for indiv in pop:
        indiv.minimized = False
      
    else:
      nevals_sum = 0
      # for each individual in population,
      for index, indiv in enumerate(pop):
        if not indiv.minimized:
          best, nevals = best_child(indiv, toolbox)
          nevals_sum += nevals
          
          if best.fitness > indiv.fitness:
            best.minimized = False
            pop[index] = best
          else:
            indiv.minimized = True
            remaining -= 1
          
    gen += 1
    
    # Update the (local) hall of fame with the generated individuals every 
    # n_generation generations
    if halloffame is not None:
        halloffame.update(pop)
    
    # record in logbook
    log.record(pop, gen=gen, nevals=nevals_sum, verbose=verbose)
    
    
  
  return (pop, log.logbook)


def eaPopulationBestFirst(pop, toolbox, n_generation=10, halloffame=None, stats=default_stats(), verbose=False):
  """
  Evolutionary algorithm that minimizes a population by selecting the most fit
  non-parent, generating children, and updating the population with the most
  fit individuals until the entire population has been a parent.
  
  Arguments
  ---------
    
    pop: iterable of individual
      Population to be minimized
    
    toolbox: deap.base.Toolbox
      Contains methods used during evolution. Expected to contain:
        
        toolbox.children(parent): Returns offspring (List[individual]) of parent
        
        toolbox.evaluate(indiv): To evaluate an individual's fitness
        
        toolbox.map(func, List[individual]): To map function executions
    
    n_generation: int, optional, default=10
      Number of generations between saving the hall of fame. Note, this only
      controls how often the hall of fame is saved. Minimization continues until
      all individuals in the population have had children.
    
    halloffame: deap.tools.HallOfFame, optional, default=None
      A Hall Of Fame containing the optimal solutions, as judged by CV score.
    
    stats: deap.tools.Statistics, optional, default=default_stats()
      The statistics to be logged.
    
    verbose: boolean
      Print information to stdout.
    
  
  Returns
  -------
    
    (pop, logbook)
    
    pop: iterable of individual
      The final population
    
    logbook: deap.tools.Logbook 
      Logbook with the statistics of the evolution
  """
  
  gen = 0
  in_pop = pop
  
  # Evaluate & record in logbook
  nevals = evaluate_all(in_pop, toolbox)
  
  # use a HallOfFame for the population
  pop = casm.learn.create_halloffame(len(in_pop))
  pop.update(in_pop)
  
  # set as non-parents, if not specified
  for indiv in pop:
    if not hasattr(indiv, "parent"):
      indiv.parent = False
  
  # record initial population stats
  log = Log(stats)
  log.record(pop, gen=gen, nevals=nevals, verbose=verbose)
  
  while gen < n_generation:
    
    # get non-parents
    nonparents = [indiv for indiv in pop if not indiv.parent]
    
    if len(nonparents) == 0:
      new_pop = initialize_population(len(pop), toolbox, verbose=verbose)
      for indiv in new_pop:
        indiv.parent = False
      pop = casm.learn.create_halloffame(len(new_pop))
      pop.update(new_pop)
      nevals = evaluate_all(pop, toolbox)
      nonparents = [indiv for indiv in pop if not indiv.parent]
    
    # select most fit non-parent
    next_parent = max(nonparents, key=lambda indiv: indiv.fitness)
    next_parent.parent = True
    
    # generate children
    offspring = toolbox.children(next_parent)
    nevals = evaluate_all(offspring, toolbox)
    
    # set offspring as non-parents
    for indiv in offspring:
      indiv.parent = False
  
    # update population
    pop.update(offspring)
  
    gen += 1
    
    # Update the (local) hall of fame with the generated individuals every 
    # n_generation generations
    if halloffame is not None:
      halloffame.update(pop)
    
    # record in logbook
    log.record(pop, gen=gen, nevals=nevals, verbose=verbose)
  
  return (pop, log.logbook)


class EvolutionaryFeatureSelection(BaseEstimator, SelectorMixin):
  
  def __init__(self, 
               algorithm, 
               estimator, 
               scoring=None, 
               cv=None, 
               penalty=0.0,
               evolve_params_kwargs=dict(), 
               constraints_kwargs=dict(),
               alg_args=list(),
               alg_kwargs=dict(),
               stats=default_stats(),
               verbose=True):
    
    self.algorithm = algorithm
    self.alg_args = alg_args
    self.alg_kwargs = alg_kwargs
    
    self.estimator = estimator
    self.scoring = scoring
    self.cv = cv
    
    self.evolve_params = EvolutionaryParams(**evolve_params_kwargs)
    self.constraints = Constraints(**constraints_kwargs)
    
    self.penalty = penalty
    self.verbose = verbose
    
    self.toolbox = deap.base.Toolbox()
    self.stats = stats
    
  
  def _run(self):
    """
    Run the specified evolutionary algorithm.
    
    Arguments
    ---------
      
      algorithm: func,
        The evolutionary algorithm to perform
      
      alg_args: List[func]
        A list of functions used to generate positional arguments for 'algorithm'
        by calling f(self).
      
      alg_kwargs: dict
        A dict of key:function pairs used to generate keyword arguments for 'algorithm'.
        
          Example: 
            
            alg_args = [
              operator.attrgetter("pop"), 
              operatot.attrgetter("toolbox")]
              
            alg_kwargs = {
              "stats":operator.attrgetter("stats"),
              "halloffame":operatot.attrgetter("halloffame")
            }
          
          Then:
            
            in_args = [f(self) for f in alg_args]
            
            in_kwargs = dict()
            for key, f in six.iteritems(alg_kwargs):
              in_kwargs[key] = f(self)
            
            self.pop, self.logbook = algorithm(*in_args, **in_kwargs)
    
    
    Returns
    -------
    
      self: returns an instance of self.
    """
    ## read or construct hall of fame
    self.halloffame = initialize_halloffame(
      filename=self.evolve_params.halloffame_filename,
      n_halloffame=self.evolve_params.n_halloffame)
    
    
    ## read or construct initial population
    self.toolbox.decorate("population", enforce_constraints(self.constraints))
    self.pop = initialize_population(self.evolve_params.n_population, self.toolbox, 
      filename=self.evolve_params.pop_begin_filename, verbose=self.verbose)
    self.pop_begin = copy.deepcopy(self.pop)
    
    ## Run algorithm
    for rep in range(self.evolve_params.n_repetition):
      
      if self.verbose:
        print("Begin", rep+1, "of", self.evolve_params.n_repetition, "repetitions")
      
      if self.verbose:
        print("Begin", self.evolve_params.n_generation, "generations")
        t = time.clock()
      
      in_args = [f(self) for f in self.alg_args]
      in_kwargs = dict()
      for key, f in six.iteritems(self.alg_kwargs):
        in_kwargs[key] = f(self)
      
      self.pop, self.logbook = self.algorithm(*in_args, **in_kwargs)
      self.pop_end = copy.deepcopy(self.pop)
      
      if self.verbose:
        print("Runtime:", time.clock() - t, "(s)\n")
      
      ## Print end population
      if self.verbose:
        print("\nFinal population:")
        casm.learn.print_population(self.pop_end)
      
      save_population(self.pop, filename=self.evolve_params.pop_end_filename, verbose=self.verbose)
      
      ## Print hall of fame
      if self.verbose:
        print("\nHall of Fame:")
        casm.learn.print_population(self.halloffame)

      save_halloffame(self.halloffame, filename=self.evolve_params.halloffame_filename, verbose=self.verbose)
      
    return self
  
  
  def _get_support_mask(self):
    """
    Return most fit inidividual found.
    
    Returns
    -------
      
      support: List[bool] of length n_features
        This is a boolean list of shape [n_features], in which an element is True 
        iff its corresponding feature is selected for retention.
    
    """
    return self.halloffame[0]
  
  
  def get_halloffame(self):
    """
    Returns the hall of fame of most fit individuals generated by 'fit'.
    
    Returns
    -------
      
      hall: deap.tools.HallOfFame
        A Hall Of Fame containing the optimal solutions, as judged by CV score.
    
    """
    return self.halloffame     


class GeneticAlgorithm(EvolutionaryFeatureSelection):
  """
  Genetic algorithm for regression by optimizing a CV score.
  
  Implements deap.algorithms.eaSimple, using selTournament, for selection, 
  cxUniform for mating, and mutFlipBit for mutation. The probability of mating 
  and mutating is set to 1.0. 
  
  The population of solutions is saved
  
  For more details see: http://deap.gel.ulaval.ca/doc/0.8/api/algo.html#complete-algorithms
  
  Attributes
  ---------
    
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
    
    scoring: string
      A string or a scorer callable object / function with signature 
      scorer(estimator, X, y). The parameter for sklearn.model_selection.cross_val_score,
      default = None, uses estimator.score().
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
    
    penalty: float
      The CV score is increased by 'penalty*(number of selected basis function)'
    
    evolve_params: casm.learn.evolve.EvolutionaryParams
      A EvolutionaryParams object containing parameters. 
    
    constraints: casm.learn.evolve.Constraints
      A Constraints object containing constraints on allowed individuals.
    
    selTournamentSize: int
      Tournament size. A larger tournament size weeds out less fit 
      individuals more quickly, while a smaller tournament size weeds out 
      less fit individuals more gradually.

    cxUniformProb: float
      Probability of swapping bits during mating.

    mutFlipBitProb: float 
      Probability of mutating bits "constraints".
    
    CrossOverProb: 1.0
      Probability of performing mating.

    MutateProb: 1.9
      Probability of performing mutation.
    
    toolbox: deap.base.Toolbox
      Contains methods used by deap during evolution.
      
    verbose: boolean
      Print information to stdout.
    
    hall: deap.tools.HallOfFame
      A Hall Of Fame containing the optimal solutions, as judged by CV score.
    
    pop_begin: iterable of individuals
      The initial population of individual solutions.
    
    pop_end: iterable of individuals
      The final population of individual solutions.
    
    
  
  """
  
  def __init__(self, estimator, scoring=None, cv=None, penalty=0.0,
               evolve_params_kwargs=dict(), constraints_kwargs=dict(),
               selTournamentSize=3, cxUniformProb=0.5, mutFlipBitProb=0.01,
               verbose=True):
    """
    Arguments
    ---------
      
      estimator:  estimator object implementing 'fit'
        The estimator specified by the input settings.
      
      scoring: string, callable or None, optional, default=None
        A string or a scorer callable object / function with signature 
        scorer(estimator, X, y). The parameter is passed to sklearn.model_selection.cross_val_score,
        has a default=None which uses estimator.score().
      
      cv: cross-validation generator or an iterable, optional, default=None
        Provides train/test splits. The parameter is passed to sklearn.model_selection.cross_val_score,
        has a default=None which uses KFold cross-validation with k=3.
      
      penalty: float, optional, default=0.0
        The CV score is increased by 'penalty*(number of selected basis function)'
      
      evolve_params_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of EvolutionaryParams object containing
        parameters. See casm.learn.evolve.EvolutionaryParams for possibilities.
      
      constraints_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of Constraints object containing
        constraints on allowed individuals. See casm.learn.evolve.Constraints for 
        possibilities.
      
      selTournamentSize: int, optional, default=3
        Tournament size. A larger tournament size weeds out less fit 
        individuals more quickly, while a smaller tournament size weeds out 
        less fit individuals more gradually.
  
      cxUniformProb: float, optional, default=0.5
        Probability of swapping bits during mating.
  
      mutFlipBitProb: float, optional, default=0.01 
        Probability of mutating bits "constraints".
      
      verbose: boolean, optional, default=True
        Print information to stdout.
    
    """
    super(GeneticAlgorithm, self).__init__(
      deap.algorithms.eaSimple,
      estimator,
      scoring=scoring, 
      cv=cv, 
      penalty=penalty,
      evolve_params_kwargs=evolve_params_kwargs, 
      constraints_kwargs=constraints_kwargs, 
      verbose=verbose,
      alg_args=[
        attrgetter("pop"), 
        attrgetter("toolbox"), 
        attrgetter("CrossOverProb"), 
        attrgetter("MutateProb"),
        lambda obj: obj.evolve_params.n_generation], 
      alg_kwargs={
        "verbose": attrgetter("verbose"),
        "halloffame": attrgetter("halloffame"),
        "stats": attrgetter("stats")
      })
    
    self.selTournamentSize = selTournamentSize
    self.cxUniformProb = cxUniformProb
    self.mutFlipBitProb = mutFlipBitProb
    
    # currently fixed:
    self.CrossOverProb = 1.0
    self.MutateProb = 1.0
    
    ## Selection
    self.toolbox.register("select", deap.tools.selTournament, tournsize=self.selTournamentSize)
    
    ## Crossover 
    self.toolbox.register("mate", deap.tools.cxUniform, indpb=self.cxUniformProb)
    self.toolbox.decorate("mate", enforce_constraints(self.constraints))
      
    ## Mutation
    self.toolbox.register("mutate", deap.tools.mutFlipBit, indpb=self.mutFlipBitProb)
    self.toolbox.decorate("mutate", enforce_constraints(self.constraints))
    
    
  
  def fit(self, X, y):
    """
    Run genetic algorithm to generate optimal solutions for X*b = y.
    
    The best solutions generated are stored in self.halloffame.
    
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
    
    ## Register toolbox functions that depend on X,y
    self.toolbox.register("individual", initNRandomOn, 
      casm.learn.creator.Individual, X.shape[1], self.evolve_params.n_features_init)
    self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.individual)
    self.toolbox.register("evaluate", casm.learn.model_selection.cross_val_score, 
      self.estimator, X, y=y, scoring=self.scoring, cv=self.cv, penalty=self.penalty)
    
    return self._run() 
      


class IndividualBestFirst(EvolutionaryFeatureSelection):
  """
  Individual best first algorithm for regression by optimizing a CV score.
  
  Implements a best first search optimization for each individual in the initial 
  population. Each individual in the population is minimized by repeatedly begin 
  replaced by its most fit child.
  
  By default, children are generated by generating all the individual that differ
  from the parent by +/- 1 selected feature. But any generating function may be
  given that generates offspring from a single parent individual.
  
  
  Attributes
  ---------
    
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
    
    scoring: string
      A string or a scorer callable object / function with signature 
      scorer(estimator, X, y). The parameter for sklearn.model_selection.cross_val_score,
      default = None, uses estimator.score().
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
    
    penalty: float
      The CV score is increased by 'penalty*(number of selected basis function)'
    
    evolve_params: casm.learn.evolve.EvolutionaryParams
      A EvolutionaryParams object containing parameters. 
    
    constraints: casm.learn.evolve.Constraints
      A Constraints object containing constraints on allowed individuals.
    
    children: func
      A function with signature 'offspring = func(indiv)'. 
    
    toolbox: deap.base.Toolbox
      Contains methods used by deap during evolution.
      
    verbose: boolean
      Print information to stdout.
    
    hall: deap.tools.HallOfFame
      A Hall Of Fame containing the optimal solutions, as judged by CV score.
    
    pop_begin: iterable of individuals
      The initial population of individual solutions.
    
    pop_end: iterable of individuals
      The final population of individual solutions.
  """
  
  def __init__(self, estimator, scoring=None, cv=None, penalty=0.0,
               evolve_params_kwargs=dict(), constraints_kwargs=dict(),
               children=single_flip_children,
               verbose=True):
    """
    Arguments
    ---------
      
      estimator:  estimator object implementing 'fit'
        The estimator specified by the input settings.
      
      scoring: string, callable or None, optional, default=None
        A string or a scorer callable object / function with signature 
        scorer(estimator, X, y). The parameter is passed to sklearn.model_selection.cross_val_score,
        has a default=None which uses estimator.score().
      
      cv: cross-validation generator or an iterable, optional, default=None
        Provides train/test splits. The parameter is passed to sklearn.model_selection.cross_val_score,
        has a default=None which uses KFold cross-validation with k=3.
      
      penalty: float, optional, default=0.0
        The CV score is increased by 'penalty*(number of selected basis function)'
      
      evolve_params_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of EvolutionaryParams object containing
        parameters. See casm.learn.evolve.EvolutionaryParams for possibilities.
      
      constraints_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of Constraints object containing
        constraints on allowed individuals. See casm.learn.evolve.Constraints for 
        possibilities.
      
      children: func, optional, default=single_flip_children
        A function with signature 'offspring = func(indiv)'. 
    
      verbose: boolean, optional, default=True
        Print information to stdout.
    
    """
    
    # add a count of the number of individuals fully optimized
    stats = default_stats()
    stats_indiv = deap.tools.Statistics(key=lambda indiv: indiv.minimized)
    stats_indiv.register("completed", sum)
    mstats = deap.tools.MultiStatistics(fitness=stats, optimization=stats_indiv)
    
    super(IndividualBestFirst, self).__init__(
      eaIndividualBestFirst,
      estimator,
      scoring=scoring, 
      cv=cv, 
      penalty=penalty,
      evolve_params_kwargs=evolve_params_kwargs, 
      constraints_kwargs=constraints_kwargs, 
      verbose=verbose,
      alg_args=[
        attrgetter("pop"), 
        attrgetter("toolbox"), 
        lambda obj: obj.evolve_params.n_generation], 
      alg_kwargs={
        "verbose": attrgetter("verbose"),
        "halloffame": attrgetter("halloffame"),
        "stats": attrgetter("stats")
      },
      stats=mstats)
    
    self.children = children
    
    self.toolbox.register("children", self.children)
    self.toolbox.decorate("children", restrict_constraints(self.constraints))
    
  
  def fit(self, X, y):
    """
    Run individual best first algorithm to generate optimal solutions for X*b = y.
    
    The best solutions generated are stored in self.halloffame.
    
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
    ## Register toolbox functions that depend on X,y
    self.toolbox.register("individual", initNRandomOn, 
      casm.learn.creator.Individual, X.shape[1], self.evolve_params.n_features_init)
    self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.individual)
    self.toolbox.register("evaluate", casm.learn.model_selection.cross_val_score, 
      self.estimator, X, y=y, scoring=self.scoring, cv=self.cv, penalty=self.penalty)
    
    return self._run()



class PopulationBestFirst(EvolutionaryFeatureSelection):
  """
  Population best first algorithm for regression by optimizing a CV score.
  
  Implements a best first search optimization for a population of individual
  solutions. Each individual is associated with a 'status' that indicates
  whether it has had children yet or not. At each step, the most fit individual
  that hasn't had children has children and the population is updated to keep
  only the 'n_population' most fit individuals. The algorithm stops when all 
  individuals in the population have had children.
  
  By default, children are generated by generating all the individual that differ
  from the parent by +/- 1 selected feature. But any generating function may be
  given that generates offspring from a single parent individual.
  
  
  Attributes
  ----------
    
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
    
    scoring: string
      A string or a scorer callable object / function with signature 
      scorer(estimator, X, y). The parameter for sklearn.model_selection.cross_val_score,
      default = None, uses estimator.score().
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
    
    penalty: float
      The CV score is increased by 'penalty*(number of selected basis function)'
    
    evolve_params: casm.learn.evolve.EvolutionaryParams
      A EvolutionaryParams object containing parameters. 
    
    constraints: casm.learn.evolve.Constraints
      A Constraints object containing constraints on allowed individuals.
    
    children: func
      A function with signature 'offspring = func(indiv)'. 
    
    toolbox: deap.base.Toolbox
      Contains methods used by deap during evolution.
      
    verbose: boolean
      Print information to stdout.
    
    hall: deap.tools.HallOfFame
      A Hall Of Fame containing the optimal solutions, as judged by CV score.
    
    pop_begin: iterable of individuals
      The initial population of individual solutions.
    
    pop_end: iterable of individuals
      The final population of individual solutions.
  """
  
  def __init__(self, estimator, scoring=None, cv=None, penalty=0.0,
               evolve_params_kwargs=dict(), constraints_kwargs=dict(),
               children=single_flip_children,
               verbose=True):
    """
    Arguments
    ---------
      
      estimator:  estimator object implementing 'fit'
        The estimator specified by the input settings.
      
      scoring: string, callable or None, optional, default=None
        A string or a scorer callable object / function with signature 
        scorer(estimator, X, y). The parameter is passed to sklearn.model_selection.cross_val_score,
        has a default=None which uses estimator.score().
      
      cv: cross-validation generator or an iterable, optional, default=None
        Provides train/test splits. The parameter is passed to sklearn.model_selection.cross_val_score,
        has a default=None which uses KFold cross-validation with k=3.
      
      penalty: float, optional, default=0.0
        The CV score is increased by 'penalty*(number of selected basis function)'
      
      evolve_params_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of EvolutionaryParams object containing
        parameters. See casm.learn.evolve.EvolutionaryParams for possibilities.
      
      constraints_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of Constraints object containing
        constraints on allowed individuals. See casm.learn.evolve.Constraints for 
        possibilities.
      
      children: func, optional, default=single_flip_children
        A function with signature 'offspring = func(indiv)'. 
    
      verbose: boolean, optional, default=True
        Print information to stdout.
    
    """
    # add a count of the number of individuals fully optimized
    stats = default_stats()
    stats_indiv = deap.tools.Statistics(key=lambda indiv: indiv.parent)
    stats_indiv.register("parents", sum)
    mstats = deap.tools.MultiStatistics(fitness=stats, optimization=stats_indiv)
    
    super(PopulationBestFirst, self).__init__(
      eaPopulationBestFirst,
      estimator,
      scoring=scoring, 
      cv=cv, 
      penalty=penalty,
      evolve_params_kwargs=evolve_params_kwargs, 
      constraints_kwargs=constraints_kwargs, 
      verbose=verbose,
      alg_args=[
        attrgetter("pop"), 
        attrgetter("toolbox"), 
        lambda obj: obj.evolve_params.n_generation], 
      alg_kwargs={
        "verbose": attrgetter("verbose"),
        "halloffame": attrgetter("halloffame"),
        "stats": attrgetter("stats")
      },
      stats=mstats)
    
    self.children = children
    
    self.toolbox.register("children", self.children)
    self.toolbox.decorate("children", restrict_constraints(self.constraints))
    
  
  def fit(self, X, y):
    """
    Run population best first algorithm to generate optimal solutions for X*b = y.
    
    The best solutions generated are stored in self.halloffame.
    
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
    ## Register toolbox functions that depend on X,y
    self.toolbox.register("individual", initNRandomOn, 
      casm.learn.creator.Individual, X.shape[1], self.evolve_params.n_features_init)
    self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.individual)
    self.toolbox.register("evaluate", casm.learn.model_selection.cross_val_score, 
      self.estimator, X, y=y, scoring=self.scoring, cv=self.cv, penalty=self.penalty)
    
    return self._run()



