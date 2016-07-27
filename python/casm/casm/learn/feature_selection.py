from evolve import *
from casm.learn.fit import make_fitting_data, make_estimator, make_selector, \
  add_individual_detail, print_halloffame
import casm.learn.cross_validation

def fit_and_select(input, save=True, verbose=True, read_existing=True, hall=None, input_filename=None):
  """
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
    
    input_filename: str, optional
      Used to determine the name of the 'specs' file saved.
  
  
  Returns
  -------
    
    (fdata, estimator, selector)
    
    fdata: casm.learn.FittingData
      A FittingData instance containing the problem data.
    
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
    
    selector:  selector object implementing 'fit' and having either a 
               'get_support()' or 'get_halloffame()' member
      The feature selector specified by the input settings.
      
  
  """
  # construct FittingData
  if verbose:
    print "# Get fitting data..."
  fdata = make_fitting_data(input, save=True, verbose=verbose, read_existing=True, input_filename=input_filename)
    
  # construct model used for fitting
  if verbose:
    print "# Construct estimator..."
  estimator = make_estimator(input, verbose=verbose)
  
  # feature selection
  if verbose:
    print "# Construct feature selector..."
  selector = make_selector(input, estimator, 
    scoring=fdata.scoring, cv=fdata.cv, penalty=fdata.penalty, verbose=verbose)
  
  if not hasattr(selector, "get_halloffame") and not hasattr(selector, "get_support"):
    raise Exception("Selector has neither 'get_halloffame' nor 'get_support'")
  
  # fit
  if verbose:
    print "# Fit and select..."
  t = time.clock()
  selector.fit(fdata.weighted_X, fdata.weighted_y)
  
  # print calculation properties that may be of interest
  if verbose:
    # custom for SelectFromModel
    if hasattr(selector, "estimator_") and hasattr(selector.estimator_, "n_iter_"):
      print "#   Iterations:", selector.estimator_.n_iter_
    if hasattr(selector, "estimator_") and hasattr(selector.estimator_, "alpha_"):
      print "#   Alpha:", selector.estimator_.alpha_
    if hasattr(selector, "threshold"):
      print "#   Feature selection threshold:", selector.threshold
    print "#   DONE  Runtime:", time.clock() - t, "(s)\n"
  
  if hall is not None:
    # store results: Assume selector either has a 'get_halloffame()' attribute, or 'get_support()' member
    if hasattr(selector, "get_halloffame"):
      if verbose:
        print "Adding statistics..."
      selector_hall = selector.get_halloffame()
      for i in xrange(len(selector_hall)):
        add_individual_detail(selector_hall[i], estimator, fdata, input, selector=selector)
      if verbose:
        print "  DONE\n"
      
      if verbose:
        print "Result:"
      print_halloffame(selector_hall)
      
      hall.update(selector_hall)
    
    elif hasattr(selector, "get_support"):
      indiv = casm.learn.creator.Individual(selector.get_support())
      if verbose:
        print "Adding statistics..."
      indiv.fitness.values = casm.learn.cross_validation.cross_val_score(
        estimator, fdata.weighted_X, indiv, 
        y=fdata.weighted_y, scoring=fdata.scoring, cv=fdata.cv, penalty=fdata.penalty)
      add_individual_detail(indiv, estimator, fdata, input, selector=selector)
      if verbose:
        print "  DONE\n"
      
      if verbose:
        print "Result:"
      print_halloffame([indiv])
      
      hall.update([indiv])
    
  return (fdata, estimator, selector)


