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
    
    
  """
  
  def __init__(self, n_population=100, n_generation=10, n_repetition=100, n_features_init=1, 
               pop_begin_filename = "population_begin.pkl",
               pop_end_filename = "population_end.pkl",
               halloffame_filename = "halloffame.pkl",
               n_halloffame = 25):
    """
    Arguments:
      n_population: integer, population size
      n_generation: integer, number of generations between saving results
      n_repetition: integer, number of repetitions of n_generation generations to perform
      n_features_init: integer, number of basis functions per individual to select in 
                   the initial population
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
  
  def __init__(self, estimator, scoring=None, cv=None, penalty=0.0,
               evolve_params_kwargs=dict(), constraints_kwargs=dict(),
               selTournamentSize=3, cxUniformProb=0.5, mutFlipBitProb=0.01,
               verbose=True):
    
    self.toolbox = deap.base.Toolbox()
    
    self.estimator = estimator
    self.scoring = scoring
    self.cv = cv
    
    self.evolve_params = EvolutionaryParams(**evolve_params_kwargs)
    self.constraints = casm.learn.tools.Constraints(**constraints_kwargs)
    
    # currently fixed:
    self.CrossOverProb = 1.0
    self.MutateProb = 1.0
    
    ## Selection
    self.toolbox.register("select", deap.tools.selTournament, tournsize=selTournamentSize)
    
    ## Crossover 
    self.toolbox.register("mate", deap.tools.cxUniform, indpb=cxUniformProb)
    self.toolbox.decorate("mate", casm.learn.tools.check_constraints(self.constraints))
      
    ## Mutation
    self.toolbox.register("mutate", deap.tools.mutFlipBit, indpb=mutFlipBitProb)
    self.toolbox.decorate("mutate", casm.learn.tools.check_constraints(self.constraints))
    
    ## penalty
    self.penalty = penalty
    
    ## verbose
    self.verbose = verbose
  
  
  def fit(self, X, y):
    
    self.toolbox.register("individual", casm.learn.tools.initNRandomOn, 
      casm.learn.creator.Individual, X.shape[1], self.evolve_params.n_features_init)
    self.toolbox.register("population", deap.tools.initRepeat, list, self.toolbox.individual)
    
    ## read or construct hall of fame
    halloffame_filename = self.evolve_params.halloffame_filename
    self.halloffame = deap.tools.HallOfFame(self.evolve_params.n_halloffame)
    print "# GeneticAlgorithm Hall of Fame size:", self.evolve_params.n_halloffame, "\n"
    
    if os.path.exists(self.evolve_params.halloffame_filename):
      existing_hall = pickle.load(open(self.evolve_params.halloffame_filename, 'rb'))
      print "Loading Hall of Fame:", self.evolve_params.halloffame_filename
      self.halloffame.update(existing_hall)
    
    
    ## read or construct initial population
    if os.path.exists(self.evolve_params.pop_begin_filename):
      print "Loading initial population:", self.evolve_params.pop_begin_filename
      self.pop = pickle.load(open(self.evolve_params.pop_begin_filename, 'rb'))
    else:
      print "Constructing random initial population"
      self.pop = self.toolbox.population(self.evolve_params.n_population)
    self.pop_begin = copy.deepcopy(self.pop)
    
    ## Fitness evaluation
    self.toolbox.register("evaluate", casm.learn.cross_validation.cross_val_score, 
      self.estimator, X, y=y, scoring=self.scoring, cv=self.cv, penalty=self.penalty)
    
    ## Run algorithm
    
    for rep in range(self.evolve_params.n_repetition):
      print "Begin", rep+1, "of", self.evolve_params.n_repetition, "repetitions"
      # Stats to show during run
      self.stats = deap.tools.Statistics(key=lambda ind: ind.fitness.values)
      self.stats.register("avg", np.mean)
      self.stats.register("std", np.std)
      self.stats.register("min", np.min)
      self.stats.register("max", np.max)

      print "Begin", self.evolve_params.n_generation, "generations"
      t = time.clock()
      self.pop, self.log = deap.algorithms.eaSimple(self.pop, self.toolbox, self.CrossOverProb, self.MutateProb, 
        self.evolve_params.n_generation, verbose=True, halloffame=self.halloffame, stats=self.stats)
      print "Runtime:", time.clock() - t, "(s)\n"
      self.pop_final = copy.deepcopy(self.pop)
      
      ## Print end population
      if self.verbose:
        print "\nFinal population:"
        casm.learn.print_population(self.pop_final)
      
      # pickle end population
      print "\nPickling end population to:", self.evolve_params.pop_end_filename
      f = open(self.evolve_params.pop_end_filename, 'wb')
      pickle.dump(self.pop_final, f)
      f.close() 
      
      
      ## Print hall of fame
      if self.verbose:
        print "\nGeneticAlgorithm Hall of Fame:"
        casm.learn.print_population(self.halloffame)

      # pickle hall of fame
      print "\nPickling GeneticAlgorithm Hall of Fame to:", self.evolve_params.halloffame_filename
      f = open(self.evolve_params.halloffame_filename, 'wb')
      pickle.dump(self.halloffame, f)
      f.close()

  
  def _get_support_mask(self):
    """Return most fit inidividual found"""
    return self.halloffame[0]
    
      
#      ## Save the overall best
#      print "\nSaving best ECI:"
#      casm.learn.print_eci(hall[0].eci)
#      write_eci(proj, hall[0].eci) 
#      
    
#    
#    print "\nUse casm.learn.plot to analyze cluster expansion predictions"
#
#
#def main(proj, input):
#  
#  ## Filenames, for now fixed
#  
#  # selection of calculated configurations to train on
#  training_selection_filename = "train"
#  
#  # store cross validation sets so that the same ones are used to compare all 'individuals'
#  cv_filename = "cv.pkl"
#  
#  # if it exists, begin with this population, else generate randomly
#  population_begin_filename = "population_begin.pkl"
#  
#  # final population
#  population_end_filename = "population_end.pkl"
#  
#  # most optimal individuals of all time, if it exists will be read in and 
#  # continuely updated with addition runs
#  hall_of_fame_filename = "hall_of_fame.pkl"
#  
#  
#  # expect 'casm.learn input.json'
#  print input
#  
#  # construct FittingData
#  fdata = make_fitting_data(input)
#  
#  # construct ModelData
#  mdata = make_estimator(input)
#  
#  # construct hall of fame
#  if os.path.exists(hall_of_fame_filename):
#    hall = pickle.load(open(hall_of_fame_filename, 'rb'))
#  else:
#    hall = deap.tools.HallOfFame(input["hall_of_fame_size"])
#  
#  # feature selection & fit method
#  
#  # store results
#
#if __name__ == "__main__":
#  
#  parser = argparse.ArgumentParser(description = 'Fit cluster expansion coefficients (ECI)')
#  parser.add_argument('filename', nargs='?', help='Input file', type=str)
#  parser.add_argument('--path', help='Path to CASM project. Default assumes the current directory is in the CASM project.', type=str, default=os.getcwd())
#  parser.add_argument('--format', help='Print input file description', action="store_true")
#  parser.add_argument('--example_input', help='Print example input file', action="store_true")
#  args = parser.parse_args()
#  
#  if args.format:
#    casm.learn.print_input_help()
#    exit()
#  
#  if args.example_input:
#    print json.dumps(casm.learn.example_input(), indent=2)
#    exit()
#  
#  # for now, assume being run
#  proj = Project(args.path)
#  input = json.load(open(args.filename, 'r'))
#  
#  main(proj, input)
#


