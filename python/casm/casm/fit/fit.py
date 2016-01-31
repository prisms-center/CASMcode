import sklearn.linear_model
import sklearn.cross_validation
import random
import numpy as np
from math import sqrt

## This part needs to be in global scope for parallization #####################  
from deap import creator
from deap import base

# we'll want to minimize a cv score
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))

# each individual is a list of True or False indicating if each basis function should 
# be included in the model
creator.create("Individual", list, fitness=creator.FitnessMin, input=None)
################################################################################  


def find_method(mods, attrname):
  for m in mods:
    if hasattr(m, attrname):
      return getattr(m, attrname)
  print "ERROR: Could not find a method named:", attrname
  print "Tried:", mods
  raise AttributeError("Could not find: " + attrname)

def example_input():
  input = dict()

  # regression model
  input["model"] = dict()
  input["model"]["method"] = "LinearRegression"
  input["model"]["kwargs"] = None
  
  # property begin fit
  input["property"] = "formation_energy"
  
  # sample weighting
  input["weight"] = dict()
  input["weight"]["method"] = None
  input["weight"]["kwargs"] = None
  
  # cross validation
  input["cv"] = dict()
  input["cv"]["method"] = "LeaveOneOut"
  input["cv"]["kwargs"] = None
  input["cv"]["penalty"] = 0.0
  
  # hall of fame
  input["hall_of_fame_size"] = 25
  
  
  return input


def print_input_help():
  
  print \
  """
  Input files:
    Settings file: A JSON file containing settings describing how to perform the fit.
    'train': A configuration selection file used to fit the cluster expansion. 
    'population_begin.pkl': An optional input file providing an initial set of
      candidate ECI sets.
  
  Generated files:
    'fit_data.pkl': Stores cross validation sets, training data, model weights and
      other information that can be used when running repeatedly.
    'hall_of_fame.pkl': Stores the best ECI sets found, as determined by the CV
      score.
    'population_end.pkl': Stores the results of the most recent optimization. Can
      be renamed 'population_begin.pkl' to use as the initial state of a new run.
  
  Settings file description:
  ------------------------------------------------------------------------------
  {
  
  # A scikit-learn linear model name and keyword args used to construct the model 
  # object. Options include: 'LinearRegression', 'Ridge', 'Lasso', etc. 
  # See: http://scikit-learn.org/stable/modules/linear_model.html
  # By default, the kwarg "fit_intercept" is set to False.
    "model": {
      "method": "LinearRegression", 
      "kwargs": null
    },
  
  # Method to use for weighting training data. 
  #
  # If weights are included, then the linear model is changed from
  #   corr*ECI = property  ->  L*corr*ECI = L*property, 
  #
  # where 'corr' is the correlation matrix of shape (Nvalue, Nbfunc),
  # and 'property' is a vector of Nvalue calculated properties, and 
  # W = L*L.transpose() is the weight matrix.
  #
  # By default, W = np.matlib.eye(Nvalue) (unweighted).
  #
  # If the weighting method provides 1-dimensional input (this is typical), in
  # a numpy array called 'w':
  #   W = np.diag(w)*Nvalue/np.sum(w)
  #
  # If the 'custom2d' method is used, the input W_in must by Hermitian, 
  # positive-definite and is normalized by:
  #   W = W_in*Nvalue/np.sum(W_in)
  #
  # The weighting methods are:
  #   'wHullDist': Weight according to w_i = A*exp(-hull_dist/kT) + B, where A, B, 
  #     and kT are user-defined kwargs parameters, and hull_dist is the distance 
  #     from the convex hull of the training data
  #   'wEmin': Weight according to w_i = A*exp(-dist_from_minE/kT) + B,
  #     where A, B, and kT are user-defined kwargs parameters, and dist_from_minE 
  #     is calculated from the training data
  #   'wEref': weight according to w_i = A*exp(-(formation_energy - E0)/kT) + B, for
  #     (formation_energy - E0) > 0.0; and w_i = 1.0 if (formation_energy - E0) <= 0.0.
  #     where A, B, E0, and kT are user-defined kwargs parameters.
  #   'wCustom': Weights are read from a column titled 'weight' in the training data 
  #     selection file.
  #   'wCustom2d': Weights are read from columns in the training data selection file,
  #     which are expected to be titled 'weight(0)' ... 'weight(Nvalue-1)'  
    "weight": {
      "method": null, 
      "kwargs": null
    }
    
  # Name of property to be fit, as used for input to 'casm query -k'
    "property": "formation_energy", 
  
  # Hall of fame size, the number of best sets of ECI to store in 'hall_of_fame.pkl',
  # as determined by CV score.
    "hall_of_fame_size": 25, 
  
  # A scikit-learn cross validation method to use to generate cross validation sets.
  #
  # Options include 'KFold', 'ShuffleSplit', 'LeaveOneOut', etc.
  # See: http://scikit-learn.org/stable/modules/cross_validation.html
  #
  # The cv score reported is:
  #
  #   cv = sqrt(np.mean(scores)) + (Number of non-zero ECI)*penalty, 
  #
  # where 'scores' is an array containing the mean squared error calculated for  
  # each training/testing set, '(Number of non-zero ECI)' is the number of basis 
  # functions with non-zero ECI, and 'penalty' is the user-input penalty per basis 
  # function (default=0.0).
  # By default, the kwarg "shuffle" is set to True.
    "cv": {
      "method": "LeaveOneOut", 
      "kwargs": null
      "penalty": 0.0
    }, 
  
  # Feature selection method to use:
  #
  # Options include classes in casm.fit.feature_selection and sklearn.feature_selection:
  # Evolutionary algorithms, from casm.fit.feature_selection, are implemented
  # using deap: http://deap.readthedocs.org/en/master/index.html
  #   "GeneticAlgorithm": implements deap.algorithms.eaSimple, using selTournament,
  #     for selection, cxUniform for mating, and mutFlipBit for mutation. The
  #     probabilty of mating and mutating is set to 1.0.
  #     Options for "kwargs":
  #       "Npop": int, (default 100) Population size. This many random initial 
  #         starting individuals are created.
  #       "Ngen": int, (default 10) Number of generations between saving the hall 
  #         of fame.
  #       "Nrep": int, (default 100) Number of repetitions of Ngen generations. 
  #         Each repetition begins with the existing final population.
  #       "Nbfunc_init: int or "all", (default 0) Number of randomly selected 
  #          basis functions to initialize each individual with.
  #       "selTournamentSize": int, (default 3). Tournament size. A larger 
  #          tournament size weeds out less fit individuals more quickly, while
  #          a smaller tournament size weeds out less fit individuals more
  #          gradually.
  #       "cxUniformProb": number, (default 0.5) Probability of swapping bits 
  #         during mating.
  #       "mutFlipBitProb": number, (default 0.01) Probability of mutating bits
  #       "constraints": See below.
  #   "IndividualBestFirst": Best first search optimization for each individual 
  #     in the initial population. At each step, all the 'children' that differ
  #     by +/- 1 selected basis function are evaluated and the most fit child
  #     of each child is chosen to replace it's parent, until the CV score is 
  #     minimized.
  #     Options for "kwargs":
  #       "Npop": int, (default 100) Population size. This many random initial 
  #         starting individuals are minimized and the results saved in the hall 
  #         of fame.
  #       "Ngen": int, (default 10) Number of generations between saving the hall 
  #         of fame.
  #       "Nrep": int, (default 10) Number of repetitions for minimizing Npop 
  #         individuals. Each repetition begins with a new population of random 
  #         individuals.
  #       "Nbfunc_init: (default 5) Number of randomly selected basis functions 
  #          to initialize each individual with.
  #       "constraints": See below.
  #   "PopulationBestFirst": Each individual is associated with a 'status' that 
  #     is '1' if that individuals children have been evaluated, and '0' if they  
  #     have not been evaluated. At each step, the children of the most fit 
  #     individual with status '0' are evaluated and the population is updated to 
  #     keep only the 'Npop' most fit individuals. The algorithm stops when all 
  #     individuals in the population have status '1'.
  #       "Npop": int, (default 100) Population size. This many random initial 
  #         starting individuals are minimized and the results saved in the hall 
  #         of fame.
  #       "Ngen": int, (default 10) Number of generations between saving the hall 
  #         of fame.
  #       "Nrep": int, (default 10) Number of repetitions for minimizing Npop 
  #         individuals. Each repetition begins with a new population of random 
  #         individuals.
  #       "Nbfunc_init: Number of randomly selected basis functions to initialize
  #          each individual with.
  #       "constraints": See below.
  #
  # The evolutionary algorithms have an optional set of "constraints" parameters
  # that may restrict the number of basis functions selected to some range, or
  # enforce some basis functions to have or not have coefficients:
  #   "Nbfunc_min": (integer) At least Nbfunc_min basis functions must be selected at all times
  #   "Nbfunc_max": (integer or "all") No more than Nbfunc_max basis functions may be selected at 
  #     any time. Default is "all".
  #   "FixOn": An array of indices of basis functions that must be included.
  #   "FixOff": An array of indices of basis functions that may not be included.
  #
  # From sklearn.feature_selection, see: http://scikit-learn.org/stable/modules/feature_selection.html
  #   "SelectFromModel": Directly fit the chosen model using all basis functions  
  #     and select only basis functions with coefficents smaller than a "threshold"
  #     (default None). 
  #   "RFE", "RFECV": Recursive feature selection
    "feature_selection" : {
      "method": "GeneticAlgorithm",
      "kwargs": {
        "Npop": 100,
        "Ngen": 10,
        "Nrep": 100,
        "Nbunc_init": 0,
        "selTournamentSize": 3,
        "cxUniformProb": 0.5,
        "mutFlipBitProb": 0.01,
        "constraints": {
          "Nbfunc_min": 0,
          "Nbfunc_max": "all",
          "FixOn": [],
          "FixOff": []
        }
      }
    }
  }
  ------------------------------------------------------------------------------
  """
  
  


class FittingData(object):
  """ 
  FittingData holds correlations, property values, sample weights, etc. used
  to solve:
  
    wcorr * ECI = wvalue 
  
  a weighted linear model where the weights are given by W = L * L.transpose(),
  and 
    
    wcorr = L * corr
    wvalue = L * value
    
  Attributes:
    corr: numpy array of correlations, shape: (Nvalue, Nbfunc)
    value: numpy array of calculated values to fit to, shape: (Nvalue, 1)
    cv: a scikit-learn type cross validation iterator
    Nvalue: number of property values to be fit to
    Nbfunc: number of basis functions calculated
    W: numpy array containing sample weights, shape: (Nvalue, Nvalue)
    L: wvalue = L*value, wcorr = L*corr, where W = L * L.transpose(), 
    wcorr: numpy array of weighted correlations, shape: (Nvalue, Nbfunc)
      wcorr = L*corr, where W = L * L.tranpose()
    wvalue: numpy array of calculated values to fit to, shape: (Nvalue, 1)
      wvalue = W'*value, where W = W' * W'.transpose()
    scoring: parameter for sklearn.cross_validation.cross_val_score
        default = None, uses model.score()
  """
  
  def __init__(self, corr, value,cv, sample_weight=[], scoring=None):
    """
    Arguments:
      corr: numpy array of correlations, shape: (Nvalue, Nbfunc)
      value: numpy array of calculated values to fit to, shape: (Nvalue, 1)
      model: scikit-learn type model used for fitting
      cv: a scikit-learn type cross validation iterator
      sample_weight: numpy array of sample weights, shape: (Nvalue,)
        if sample_weight == None: (default, unweighted)
          W = np.matlib.eye(N) 
        if sample_weight is 1-dimensional:
          W = np.diag(sample_weight)*Nvalue/np.sum(sample_weight) 
        if sample_weight is 2-dimensional (must be Hermitian, positive-definite):
          W = sample_weight*Nvalue/np.sum(sample_weight) 
      scoring: parameter for sklearn.cross_validation.cross_val_score
        default = None, uses model.score()
    """
    self.corr = corr
    self.value = value
    
    # Number of configurations and basis functions
    self.Nvalue, self.Nbfunc = self.corr.shape
    
    # check sample_weight and convert to square matrix
    if sample_weight == None:
      self.W = np.identity(self.value.shape[0])
    elif len(sample_weight.shape) == 1:
      self.W = np.diag(sample_weight)*self.Nvalue/np.sum(sample_weight)
    elif len(sample_weight.shape) == 2:
      self.W = sample_weight*self.Nvalue/np.sum(sample_weight)
    else:
      raise Exception("Error in LinearRegressionForLOOCV.fit: sample_weight dimension > 2")
    
    # weighted data
    self.L = np.linalg.cholesky(self.W)
    self.wcorr = np.dot(self.L, self.corr)
    self.wvalue = np.dot(self.L, self.value)
    
    # cv sets
    self.cv = cv
    
    # scoring
    self.scoring = scoring


def make_fitting_data(proj, input, verbose = True):
  """ 
  Construct a FittingData instance, either by reading existing 'fit_data.pkl',
  or from an input file settings.
  
  Arguments:
    proj: a casm.fit.Project
    input: input file as dict
      
  """
  
  # property, weight, and cv inputs should remain constant
  # model and feature_selection might change
  
  filename = "fit_data.pkl"
  
  if os.path.exists(filename):
    print "Reading existing fitting data from fit_data.pkl..."
    fdata = pickle.load(open(filename, 'rb'))
    print "  DONE"
    
    if fdata.input["property"] != input["property"]:
      print \
      """
      ERROR: Input file and stored data differ. Input 'property' has changed.
      """
      print "Stored data:", fdata.input["property"]
      print "Input:", input["property"]
    if fdata.input["cv"] != input["cv"]:
      print \
      """
      ERROR: Input file and stored data differ. Input 'cv' has changed.
      """
      print "Stored data:", fdata.input["cv"]
      print "Input:", input["cv"]
    if fdata.input["weight"] != input["weight"]:
      print \
      """
      ERROR: Input file and stored data differ. Input 'weight' has changed.
      """
      print "Stored data:", fdata.input["weight"]
      print "Input:", input["weight"]
    
  else:
  
    # get current training selection
    sel = Selection(proj, path="train")
    
    ## property
    
    # get property name (required)
    property = input["property"]
    if verbose:
      print "Property:", property
    
    ## query data
    
    # values to query
    columns = ["corr", property, "hull_dist(train,atom_frac)"]
    
    # perform query
    df = query(proj, columns, sel)

    i = df.columns.tolist().index(property)
    
    # columns of interest, as numpy arrays
    corr = df.iloc[:,:i].values
    value = df.iloc[:,i].values
    hull_dist = df.iloc[:,i+1].values
    
    Nvalue = value.shape[0]
    
    if verbose:
      print "# Training configurations:", Nvalue
    if verbose:
      print "# Basis functions:", i
    
    ## weight (optional)
    
    # set defaults if not provided
    sample_weight = None
    if "weight" not in input:
      input["weight"] = dict()
    if "kwargs" not in input["weight"]:
      input["weight"]["kwargs"] = dict()
    if "method" not in input["weight"]:
      input["weight"] = None
    
    # get kwargs
    kwargs = input["weight"]["kwargs"]
          
    # use method to get weights
    if input["weight"]["method"] == "wCustom":
      if verbose:
        print "Reading custom weights"
      df = sel.df()
      sample_weight = df["weight"].values
    elif input["weight"]["method"] == "wCustom2d":
      if verbose:
        print "Reading custom2d weights"
      cols = ["weight(" + str(i) + ")" for i in range(len(Nvalue))]
      sample_weight = df.iloc[:,cols].values
    elif input["weight"]["method"] == "wHullDist":
      sample_weight = wHullDist(hull_dist, **kwargs)
    elif input["weight"]["method"] == "wEmin":
      sample_weight = wEmin(value, **kwargs)
    elif input["weight"]["method"] == "wEref":
      sample_weight = wEref(value, **kwargs)
          
    if verbose:
      print "Weighting:", input["weight"]["method"]
      print "  kwargs:", json.dumps(kwargs)
    
    
    ## cv
    
    # set kwargs (shuffle=True) and penalty (0.0) defaults
    if "kwargs" not in input["cv"]:
      input["cv"]["kwargs"] = dict()
    if "shuffle" not in input["cv"]["kwargs"]:
      input["cv"]["kwargs"]["shuffle"] = True
    kwargs = input["cv"]["kwargs"]
    if "penalty" not in input["cv"]:
      input["cv"]["penalty"] = 0.0
    
    
    # get cv method (required)
    cv_method = None
    if input["model"]["method"] == "LinearRegression" and input["cv"]["method"] == "LeaveOneOut":
      cv_method = casm.fit.LeaveOneOutForLLS
    else:
      cv_method = getattr(sklearn.cross_validation, input["cv"]["method"])
    cv = cv_method(Nvalue, **kwargs)
    if verbose:
      print "CV:", input["cv"]["method"]
      print "  kwargs:", json.dumps(kwargs) 
    
    # set model (required) and scoring
    scoring = sklearn.metrics.make_scorer(mean_squared_error, greater_is_better=True)
    if input["model"]["method"] == "LinearRegression" and input["cv"]["method"] == "LeaveOneOut":
      model_method = casm.fit.LinearRegressionForLOOCV
      scoring = None
    else:
      model_method = getattr(sklearn.cross_validation, input["model"]["method"])
    model = mode_method(**kwargs)
    
    # get penalty
    penalty = input["cv"]["penalty"]
    if verbose:
      print "  penalty:", penalty 
    
    fdata = casm.fit.FittingData(corr, value, cv, sample_weight=sample_weight)
    
    fdata.input = dict()
    fdata.input["cv"] = input["cv"]
    fdata.input["weight"] = input["weight"]
    fdata.input["property"] = input["property"]
    
    pickle.dump(fdata, open(filename, 'wb'))
    
    

#class ModelData(object):
#  """ 
#  FittingData holds a scikit-learn type model and scoring method
#  
#  Attributes:
#    model: scikit-learn type model used for fitting
#    
#  """
#  
#  def __init__(self, model, scoring=None):
#    """
#    Arguments:
#    
#      model: a scikit-learn type linear model
#      
#    """
#     
#    # model used for fitting
#    if model == None:
#      model = casm.fit.standard_model()
#    self.model = model
    


def make_model(input, verbose = True):
  """
  Construct model object from input file settings.
  
  Arguments:
    input: input file as dict
  
  """""
  
  ## model
  
  # get kwargs (default: fit_intercept=False)
  kwargs = dict()
  if "kwargs" in input["model"]:
    if input["model"]["kwargs"]:
      kwargs = input["model"]["kwargs"]
  if "fit_intercept" not in kwargs:
    kwargs["fit_intercept"] = False
  
  if input["model"]["method"] == "LinearRegression" and input["cv"]["method"] == "LeaveOneOut":
    model_method = casm.fit.LinearRegressionForLOOCV
  else:
    model_method = getattr(sklearn.cross_validation, input["model"]["method"])
  model = mode_method(**kwargs)
  
  if verbose:
    print "Model:", input["model"]["method"]
    print "  kwargs:", json.dumps(kwargs)  
  
  return model

  
def eci(individual, coef):
  """ 
  Return a list of tuple: [(index, coef), ...]
  
  index: a list containing the index of non-zero eci
  coef: the eci value corresponding to the index
  
  """
  return zip(indices(individual), coef.tolist())
  
  
def indices(individual):
  """ convert bitstring to indices for slicing """
  indices = []
  for i in range(len(individual)):
    if individual[i]:
      indices.append(i)
  return indices


def print_hall_of_fame(hall):
  """ 
  Print all individual in hall of fame.
  
  Format:
      0: 0111000100000000000000100001000000000100...   Nbfunc: 14      CV:  0.021686393   LOOCV:  0.020368764   Note: This was the best
      1: 1111100100001100001010000000100001001000...   Nbfunc: 23      CV:  0.021858133   LOOCV:  0.018797794
      2: 1111000100000000001000100000101001000010...   Nbfunc: 18      CV:  0.021908411   LOOCV:  0.019696178
  ...
  """
  form_str = "{0:5}: {1}   Nbfunc: {2:<5}   CV: {3:< 12.8g}   WRMS: {4:< 12.8g}   RMS: {5:< 12.8g}"
  for i in range(len(hall)):
    bitstr = ""
    for j in range(min(len(hall[i]),40)):
      if hall[i][j]:
        bitstr += '1'
      else:
        bitstr += '0'
    if len(hall[i]) > 40:
      bitstr += "..."  
    if hasattr(hall[i], 'note'):
      print (form_str + "   Note: {5}").format(i, bitstr, sum(hall[i]), hall[i].fitness.values[0], hall[i].wrms, hall[i].rms, hall[i].note)
    else:
      print form_str.format(i, bitstr, sum(hall[i]), hall[i].fitness.values[0], hall[i].wrms, hall[i].rms)


def print_eci(eci):
  """
  Print ECI.
  
  Format:
      1: -1.53460558686
      2:  0.574571156376
      3:  1.04379783648
  ...
  """
  for bfunc in eci:
    print "{index:>5}: {value:< .12g}".format(index=bfunc[0], value=bfunc[1])


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
  model = make_model(input)
  
  # construct hall of fame
  if os.path.exists(hall_of_fame_filename):
    hall = pickle.load(open(hall_of_fame_filename, 'rb'))
  else
    hall = deap.tools.HallOfFame(input["hall_of_fame_size"])
  
  # feature selection & fit method
  selector = make_selector(model, input, fdata, hall)
  selector.fit(model.wcorr, model.wvalue)
  
  # store and analyze results 
  if not hasattr(selector, "halloffame"):
    support = selector.get_support()
    indiv = creator.Individual(input=input)
    for i in range(len(suport)):
      indiv[i] = support[i]
    
  # print hall of fame
  
  

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

