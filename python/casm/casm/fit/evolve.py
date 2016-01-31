#!/usr/bin/env python
import os, sys, json, time, pickle
from math import sqrt

import casm.fit
import pandas
import numpy as np
from sklearn.metrics import mean_squared_error

def default_constraints():
  return dict({"Nbfunc_min":0, "Nbfunc_max":"all", "FixOn":[], "FixOff":[]})

class GeneticAlgorithm(object):
  
  def __init__(self, model, constraints=default_constraints(), Npop=100, Ngen=10, Nrep=100, Nbfunc_init=0)
  
    self.toolbox = deap.base.Toolbox()
    
    ## specify how to initialize the Population
    self.toolbox.register("indices", init_method, **init_kwargs)
    self.toolbox.register("individual", init_method, creator.Individual, **init_kwargs)
    self.toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    
    mods = [casm.fit.tools, deap.tools]
    
    init_method = find_method(mods, options["init"]["method"])
    select_method = find_method(mods, options["select"]["method"])
    mate_method = find_method(mods, options["mate"]["method"])
    mutate_method = find_method(mods, options["mutate"]["method"])
    
    ## Selection
    self.toolbox.register("select", select_method, *options["select"]["args"])
    
    ## Crossover 
    self.toolbox.register("mate", mate_method, *options["mate"]["args"])
    self.toolbox.decorate("mate", check_constraints, **options["constraints"]["kwargs"])
      
    ## Mutation
    self.toolbox.register("mutate", mutate_method, *options["mutate"]["args"])
    self.toolbox.decorate("mutate", check_constraints, **options["constraints"]["kwargs"])
    
    
    ## Fitness evaluation
    
    # specify how to evaluate the fitness of an individual
    self.toolbox.register("evaluate", casm.fit.evaluate_wcv, fdata=fdata, model=self.model)
    
    
    ## Data collection
    
    # Store the 100 best individuals ever encountered
    
    # load existing hall of fame if exists, else construct new empty hall of fame
    if os.path.exists(hall_of_fame_filename):
      hall = pickle.load(open(hall_of_fame_filename,'rb'))
    else:
      hall = tools.HallOfFame(25)
    
    
    ## Run algorithm
    
    # Number of generations to run per repitition
    Ngen = 10
    
    # Number of repetitions
    Nrep = 100
    
    for rep in range(Nrep):
      # Stats to show during run
      stats = tools.Statistics(key=lambda ind: ind.fitness.values)
      stats.register("avg", np.mean)
      stats.register("std", np.std)
      stats.register("min", np.min)
      stats.register("max", np.max)

      print "Begin", Ngen, "generations"
      t = time.clock()
      pop_final, log = algorithms.eaSimple(pop, toolbox, CrossOverProb, MutateProb, 
        Ngen, verbose=True, halloffame=hall, stats=stats)
      print "Runtime:", time.clock() - t, "(s)"
      
      ## For individuals in the hall of fame calculate the requsted CV score, ECI, wrms, rms 
      for indiv in hall:
        
        indiv.cv = casm.fit.evaluate_wcv(indiv, fdata)
        
        fdata.model.fit(fdata.wcorr[:,casm.fit.indices(indiv)], fdata.wvalue)
        indiv.eci = casm.fit.eci(indiv, fdata.model.coef_)
        indiv.rms = sqrt(mean_squared_error(fdata.value, fdata.model.predict(fdata.corr[:,casm.fit.indices(indiv)])))
        indiv.wrms = sqrt(mean_squared_error(fdata.wvalue, fdata.model.predict(fdata.wcorr[:,casm.fit.indices(indiv)])))
      
      
      ## Print hall of fame
      print "\nHall of Fame:"
      casm.fit.print_hall_of_fame(hall)
      
      ## Save the overall best
      print "\nSaving best ECI:"
      casm.fit.print_eci(hall[0].eci)
      write_eci(proj, hall[0].eci) 
      
      # pickle end population
      print "\nPickling end population to:", population_end_filename
      f = open(population_end_filename, 'wb')
      pickle.dump(pop_final, f)
      f.close() 
      
      # pickle hall of fame
      print "\nPickling hall of fame to:", hall_of_fame_filename
      f = open(hall_of_fame_filename, 'wb')
      pickle.dump(hall, f)
      f.close()
    
    
    print "\nUse casm.fit.plot to analyze cluster expansion predictions"


def main(proj, input):
  
  ## Filenames, for now fixed
  
  # selection of calculated configurations to train on
  training_selection_filename = "train"
  
  # store cross validation sets so that the same ones are used to compare all 'individuals'
  cv_filename = "cv.pkl"
  
  # if it exists, begin with this population, else generate randomly
  population_begin_filename = "population_begin.pkl"
  
  # final population
  population_end_filename = "population_end.pkl"
  
  # most optimal individuals of all time, if it exists will be read in and 
  # continuely updated with addition runs
  hall_of_fame_filename = "hall_of_fame.pkl"
  
  
  # expect 'casm.fit input.json'
  print input
  
  # construct FittingData
  fdata = make_fitting_data(input)
  
  # construct ModelData
  mdata = make_model_data(input)
  
  # construct hall of fame
  if os.path.exists(hall_of_fame_filename):
    hall = pickle.load(open(hall_of_fame_filename, 'rb'))
  else
    hall = deap.tools.HallOfFame(input["hall_of_fame_size"])
  
  # feature selection & fit method
  
  # store results

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description = 'Fit cluster expansion coefficients (ECI)')
  parser.add_argument('filename', nargs='?', help='Input file', type=str)
  parser.add_argument('--path', help='Path to CASM project. Default assumes the current directory is in the CASM project.', type=str, default=os.getcwd())
  parser.add_argument('--format', help='Print input file description', action="store_true")
  parser.add_argument('--example_input', help='Print example input file', action="store_true")
  args = parser.parse_args()
  
  if args.format:
    casm.fit.print_input_help()
    exit()
  
  if args.example_input:
    print json.dumps(casm.fit.example_input(), indent=2)
    exit()
  
  # for now, assume being run
  proj = Project(args.path)
  input = json.load(open(args.filename, 'r'))
  
  main(proj, input)



