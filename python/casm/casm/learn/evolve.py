#!/usr/bin/env python
import os, sys, json, time, pickle, copy
from math import sqrt

import casm.learn
import pandas
import numpy as np
from sklearn.base import BaseEstimator
from sklearn.feature_selection.base import SelectorMixin
from sklearn.metrics import mean_squared_error
import deap
import deap.tools
import deap.algorithms

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
               halloffame_filename = "halloffame.pkl",
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
      
      n_halloffame: int, optional, default=25
        Number of individuals to save in the hall of fame
  
    """
    self.n_population = n_population
    self.n_generation = n_generation
    self.n_repetition = n_repetition
    
    self.n_features_init = n_features_init
    
    self.pop_begin_filename = pop_begin_filename
    self.pop_end_filename = pop_end_filename
    
    self.halloffame_filename = halloffame_filename
    self.n_halloffame = n_halloffame
    

class GeneticAlgorithm(BaseEstimator, SelectorMixin):
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
      scorer(estimator, X, y). The parameter for sklearn.cross_validation.cross_val_score,
      default = None, uses estimator.score().
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
    
    penalty: float
      The CV score is increased by 'penalty*(number of selected basis function)'
    
    evolve_params: casm.learn.evolve.EvolutionaryParams
      A EvolutionaryParams object containing parameters. 
    
    constraints: casm.learn.tools.Constraints
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
        scorer(estimator, X, y). The parameter is passed to sklearn.cross_validation.cross_val_score,
        has a default=None which uses estimator.score().
      
      cv: cross-validation generator or an iterable, optional, default=None
        Provides train/test splits. The parameter is passed to sklearn.cross_validation.cross_val_score,
        has a default=None which uses KFold cross-validation with k=3.
      
      penalty: float, optional, default=0.0
        The CV score is increased by 'penalty*(number of selected basis function)'
      
      evolve_params_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of EvolutionaryParams object containing
        parameters. See EvolutionaryParams for possibilities.
      
      constraints_kwargs: dict, optional, default=dict()
        Keyword arguments for construction of Constraints object containing
        constraints on allowed individuals. See Constraints for possibilities.
      
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
    
    self.toolbox = deap.base.Toolbox()
    
    self.estimator = estimator
    self.scoring = scoring
    self.cv = cv
    
    self.evolve_params = EvolutionaryParams(**evolve_params_kwargs)
    self.constraints = casm.learn.tools.Constraints(**constraints_kwargs)
    
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
    self.toolbox.decorate("mate", casm.learn.tools.check_constraints(self.constraints))
      
    ## Mutation
    self.toolbox.register("mutate", deap.tools.mutFlipBit, indpb=self.mutFlipBitProb)
    self.toolbox.decorate("mutate", casm.learn.tools.check_constraints(self.constraints))
    
    ## penalty
    self.penalty = penalty
    
    ## verbose
    print "verbose:", verbose
    self.verbose = verbose
    print "self.verbose:", self.verbose
  
  
  def fit(self, X, y):
    """
    Run genetic algorithm to generate optimal solutions for X*b = y.
    
    Each solution is stored as an individual in self.halloffame.
    
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
    
    self.toolbox.register("individual", casm.learn.tools.initNRandomOn, 
      casm.learn.creator.Individual, X.shape[1], self.evolve_params.n_features_init)
    self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.individual)
    
    ## read or construct hall of fame
    halloffame_filename = self.evolve_params.halloffame_filename
    self.halloffame = deap.tools.HallOfFame(self.evolve_params.n_halloffame)
    if self.verbose:
      print "# GeneticAlgorithm Hall of Fame size:", self.evolve_params.n_halloffame, "\n"
    
    if os.path.exists(self.evolve_params.halloffame_filename):
      existing_hall = pickle.load(open(self.evolve_params.halloffame_filename, 'rb'))
      if self.verbose:
        print "Loading Hall of Fame:", self.evolve_params.halloffame_filename
      self.halloffame.update(existing_hall)
    
    
    ## read or construct initial population
    if os.path.exists(self.evolve_params.pop_begin_filename):
      if self.verbose:
        print "Loading initial population:", self.evolve_params.pop_begin_filename
      self.pop = pickle.load(open(self.evolve_params.pop_begin_filename, 'rb'))
    else:
      if self.verbose:
        print "Constructing random initial population"
      self.pop = self.toolbox.population(self.evolve_params.n_population)
    self.pop_begin = copy.deepcopy(self.pop)
    
    ## Fitness evaluation
    self.toolbox.register("evaluate", casm.learn.cross_validation.cross_val_score, 
      self.estimator, X, y=y, scoring=self.scoring, cv=self.cv, penalty=self.penalty)
    
    ## Run algorithm
    self.stats=None
    
    for rep in range(self.evolve_params.n_repetition):
      if self.verbose:
        print "Begin", rep+1, "of", self.evolve_params.n_repetition, "repetitions"
        
        # Stats to show during run
        self.stats = deap.tools.Statistics(key=lambda ind: ind.fitness.values)
        self.stats.register("avg", np.mean)
        self.stats.register("std", np.std)
        self.stats.register("min", np.min)
        self.stats.register("max", np.max)

      if self.verbose:
        print "Begin", self.evolve_params.n_generation, "generations"
      t = time.clock()
      self.pop, self.log = deap.algorithms.eaSimple(self.pop, self.toolbox, self.CrossOverProb, self.MutateProb, 
        self.evolve_params.n_generation, verbose=self.verbose, halloffame=self.halloffame, stats=self.stats)
      if self.verbose:
        print "Runtime:", time.clock() - t, "(s)\n"
      self.pop_end = copy.deepcopy(self.pop)
      
      ## Print end population
      if self.verbose:
        print "\nFinal population:"
        casm.learn.print_population(self.pop_end)
      
      # pickle end population
      if self.verbose:
        print "\nPickling end population to:", self.evolve_params.pop_end_filename
      f = open(self.evolve_params.pop_end_filename, 'wb')
      pickle.dump(self.pop_end, f)
      f.close() 
      
      
      ## Print hall of fame
      if self.verbose:
        print "\nGeneticAlgorithm Hall of Fame:"
        casm.learn.print_population(self.halloffame)

      # pickle hall of fame
      if self.verbose:
        print "\nPickling GeneticAlgorithm Hall of Fame to:", self.evolve_params.halloffame_filename
      f = open(self.evolve_params.halloffame_filename, 'wb')
      pickle.dump(self.halloffame, f)
      f.close()
      
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

