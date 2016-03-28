#!/usr/bin/env python
import os, json, time, pickle
from casm.project import Selection, Project, query, write_eci
from math import sqrt
from sklearn.metrics import mean_squared_error
import sklearn.feature_selection
import casm.fit
import argparse
import deap.tools


def main(sel, input, format=None, verbose=True):
  
  # construct FittingData
  fdata = casm.fit.make_fitting_data(sel, input, verbose=verbose)
  
  # construct model used for fitting
  estimator = casm.fit.make_estimator(input, verbose=verbose)
  
  # feature selection & fit method
  selector = casm.fit.make_selector(input, estimator, 
    scoring=fdata.scoring, cv=fdata.cv, penalty=fdata.penalty, verbose=verbose)
  
  # fit
  if verbose:
    print "# Fit and select..."
  t = time.clock()
  selector.fit(fdata.wcorr, fdata.wvalue)
  if verbose:
    # custom for SelectFromModel
    if hasattr(selector, "estimator_") and hasattr(selector.estimator_, "n_iter_"):
      print "#   Iterations:", selector.estimator_.n_iter_
    if hasattr(selector, "threshold"):
      print "#   Feature selection threshold:", selector.threshold
    print "#   DONE  Runtime:", time.clock() - t, "(s)\n"
  
  
  # construct hall of fame
  halloffame_filename = input.get("halloffame_filename", "halloffame.pkl")
  halloffame_size = input.get("halloffame_size", 25)
  hall = deap.tools.HallOfFame(halloffame_size)
  print "# Hall of Fame size:", halloffame_size, "\n"
  
  if os.path.exists(halloffame_filename):
    existing_hall = pickle.load(open(halloffame_filename, 'rb'))
    if verbose:
      print "Loading Hall of Fame:", halloffame_filename
    hall.update(existing_hall)
    
      
  # store results: Assume selector either has a 'halloffame' attribute, or 'get_support()' member
  if hasattr(selector, "halloffame"):
    print "Adding statistics..."
    for i in xrange(len(selector.halloffame)):
      casm.fit.add_individual_detail(selector.halloffame[i], estimator, fdata, selector, input)
    print "  DONE"
    
    print "Result:"
    casm.fit.print_halloffame(selector.halloffame)
    
    hall.update(selector.halloffame)
  else:
    print [1 if x else 0 for x in selector.get_support()]
    indiv = casm.fit.creator.Individual(selector.get_support())
    print "Adding statistics..."
    indiv.fitness.values = casm.fit.cross_validation.cross_val_score(
      estimator, fdata.wcorr, indiv, 
      y=fdata.wvalue, scoring=fdata.scoring, cv=fdata.cv, penalty=fdata.penalty)
    casm.fit.add_individual_detail(indiv, estimator, fdata, selector, input)
    print "  DONE"
    
    print "Result:"
    casm.fit.print_halloffame([indiv])
    
    hall.update([indiv])
  
  # print hall of fame
  #print "\nHall of Fame:"
  #casm.fit.print_halloffame(hall, format=format)
  
  # pickle hall of fame
  print "\nPickling Hall of Fame to:", halloffame_filename
  f = open(halloffame_filename, 'wb')
  pickle.dump(hall, f)
  f.close()

  

if __name__ == "__main__":
  
  parser = argparse.ArgumentParser(description = 'Fit cluster expansion coefficients (ECI)')
  parser.add_argument('-s', '--settings', nargs=1, help='Settings input filename', type=str)
  parser.add_argument('--format', help='Hall of fame print format', type=str, default=None)
  parser.add_argument('--path', help='Path to CASM project. Default assumes the current directory is in the CASM project.', type=str, default=os.getcwd())
  parser.add_argument('--settings_format', help='Print input file description', action="store_true")
  parser.add_argument('--exLasso', help='Print example input file using Lasso', action="store_true")
  parser.add_argument('--exLassoCV', help='Print example input file using LassoCV', action="store_true")
  parser.add_argument('--exRFE', help='Print example input file using Recursive Feature Elimination (RFE)', action="store_true")
  parser.add_argument('--exGeneticAlgorithm', help='Print example input file using GeneticAlgorithm', action="store_true")
  parser.add_argument('--hall', help='Print hall of fame summary', action="store_true")
  parser.add_argument('--indiv', nargs='+', help='Print individual summary. Expects index of individual in hall of fame', type=int)
  parser.add_argument('--select', nargs=1, help='Select individual to use', type=int)
  args = parser.parse_args()
  
  if args.settings_format:
    casm.fit.print_input_help()
    exit()
  
  if args.exLasso:
    print json.dumps(casm.fit.example_input_Lasso(), indent=2)
    exit()
  elif args.exLassoCV:
    print json.dumps(casm.fit.example_input_LassoCV(), indent=2)
    exit()
  elif args.exRFE:
    print json.dumps(casm.fit.example_input_RFE(), indent=2)
    exit()
  elif args.exGeneticAlgorithm:
    print json.dumps(casm.fit.example_input_GeneticAlgorithm(), indent=2)
    exit()
  
  if args.settings:
    
    print "Loading", args.settings[0]
    
    proj = Project(args.path)
    with open(args.settings[0], 'r') as f:
      try:
        input = json.load(f)
      except Exception as e:
        print "Error parsing JSON in", args.settings[0]
        raise e
    
    if args.hall or args.indiv:
      
      halloffame_filename = input.get("halloffame_filename", "halloffame.pkl")
      # print Hall of Fame summary
      print "Loading Hall of Fame:", halloffame_filename
      with open(halloffame_filename, 'rb') as f:
        existing_hall = pickle.load(f)
      
      if args.hall:
        casm.fit.print_halloffame(existing_hall, format=args.format)
      elif args.indiv:
        casm.fit.print_individual(existing_hall, args.indiv, format=args.format)
  
    elif args.select:
      
      halloffame_filename = input.get("halloffame_filename", "halloffame.pkl")
      # print Hall of Fame summary
      print "Loading Hall of Fame:", halloffame_filename
      with open(halloffame_filename, 'rb') as f:
        existing_hall = pickle.load(f)
      
      write_eci(proj, existing_hall[args.select[0]].eci, verbose=True)
      
    else:
      # run fitting
      
      # read training set
      sel = Selection(proj, input.get("train_filename", "train"))
      
      main(sel, input, verbose=True, format=args.format)
  
  else:
    
    print \
    """
    Fitting is performed in four steps:
    
    1) Select training data.
      Create a selection of configurations to include in the regression problem.
    
    2) Select scoring metric.
      Add sample weights to configurations if desired and select a cross validation
      method.
    
    3) Select estimator.
      Choose how to solve for ECI from calculated property and correlations. For
      instance: LinearRegression, Lasso, or Ridge regression.
    
    4) Select features.
      Select which basis functions to include in the cluster expansion. For instance,
      SelectFromModel along with a l-1 norm minimizing estimator. Or GeneticAlgorithm.
    
    
    To proceed, create a fitting directory: 
      mkdir fit1; cd fit1
    
    
    
    """
    
    parser.print_help()



