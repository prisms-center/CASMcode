#!/usr/bin/env python
import os, json, time, pickle
from casm.project import Selection, Project, query, write_eci
from math import sqrt
from sklearn.metrics import mean_squared_error
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
  parser.add_argument('--settings_example', help='Print example input file', action="store_true")
  args = parser.parse_args()
  
  if args.settings_format:
    casm.fit.print_input_help()
    exit()
  
  if args.settings_example:
    print json.dumps(casm.fit.example_input(), indent=2)
    exit()
  
  if args.settings:
    
    print "Loading", args.settings[0]
    
    # for now, assume being run
    proj = Project(args.path)
    input = json.load(open(args.settings[0], 'r'))
    
    # read training set
    sel = Selection(proj, input.get("train_filename", "train"))
    
    main(sel, input, verbose=True, format=args.format)
  
  else:
    
    parser.print_help()



