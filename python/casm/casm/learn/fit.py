import sklearn.linear_model
import sklearn.cross_validation
import sklearn.metrics
import random, re, time, os, types, json, pickle, copy, uuid
from os.path import splitext, basename
import numpy as np
from math import sqrt
from casm.project import Project, Selection, query
import casm.learn.linear_model
import casm.learn.tools
import pandas

def _find_method(mods, attrname):
  for m in mods:
    if hasattr(m, attrname):
      return getattr(m, attrname)
  print "ERROR: Could not find a method named:", attrname
  print "Tried:", mods
  raise AttributeError("Could not find: " + attrname)


def example_input_Lasso():
  input = dict()
  
  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_folds"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "Lasso"
  input["estimator"]["kwargs"] = dict()
  input["estimator"]["kwargs"]["alpha"] = 1e-4
  input["estimator"]["kwargs"]["max_iter"] = 1e6
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "SelectFromModel"
  input["feature_selection"]["kwargs"] = None
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input


def example_input_LassoCV():
  input = dict()

  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_folds"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "LassoCV"
  input["estimator"]["kwargs"] = dict()
  input["estimator"]["kwargs"]["eps"] = 1e-6
  input["estimator"]["kwargs"]["n_alphas"] = 100
  input["estimator"]["kwargs"]["max_iter"] = 1e6
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "SelectFromModel"
  input["feature_selection"]["kwargs"] = None
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input


def example_input_RFE():
  input = dict()

  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "LeaveOneOut"
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "LinearRegression"
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "RFE"
  input["feature_selection"]["kwargs"] = dict()
  input["feature_selection"]["kwargs"]["n_features_to_select"] = 25
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input


def example_input_GeneticAlgorithm():
  input = dict()
  
  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "LeaveOneOut"
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "LinearRegression"
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "GeneticAlgorithm"
  d = {
    "constraints_kwargs": { 
      "n_features_max": "all", 
      "n_features_min": 5, 
      "fix_off": [], 
      "fix_on": []
    }, 
    "selTournamentSize": 3, 
    "mutFlipBitProb": 0.01, 
    "evolve_params_kwargs": {
      "n_generation": 10, 
      "n_repetition": 25, 
      "n_features_init": 5, 
      "n_population": 100, 
      "halloffame_filename": "ga_halloffame.pkl", 
      "n_halloffame": 50
    }, 
    "cxUniformProb": 0.5
  }
  input["feature_selection"]["kwargs"] = d
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input


def example_input_IndividualBestFirst():
  input = dict()
  
  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "LeaveOneOut"
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "LinearRegression"
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "IndividualBestFirst"
  d = {
    "constraints_kwargs": { 
      "n_features_max": "all", 
      "n_features_min": 5, 
      "fix_off": [], 
      "fix_on": []
    }, 
    "evolve_params_kwargs": {
      "n_generation": 10, 
      "n_repetition": 25, 
      "n_features_init": 5, 
      "n_population": 10, 
      "halloffame_filename": "indiv_bestfirst_halloffame.pkl", 
      "n_halloffame": 50
    }
  }
  input["feature_selection"]["kwargs"] = d
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input


def example_input_PopulationBestFirst():
  input = dict()
  
  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "LeaveOneOut"
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "LinearRegression"
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "PopulationBestFirst"
  d = {
    "constraints_kwargs": { 
      "n_features_max": "all", 
      "n_features_min": 5, 
      "fix_off": [], 
      "fix_on": []
    }, 
    "evolve_params_kwargs": {
      "n_generation": 10, 
      "n_repetition": 25, 
      "n_features_init": 5, 
      "n_population": 50, 
      "halloffame_filename": "pop_bestfirst_halloffame.pkl", 
      "n_halloffame": 50
    }
  }
  input["feature_selection"]["kwargs"] = d
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input


def example_input_DirectSelection():
  input = dict()
  
  specs = dict()
  
  # data
  specs["data"] = dict()
  specs["data"]["filename"] = "train"
  specs["data"]["type"] = "selection"
  specs["data"]["X"] = "corr"
  specs["data"]["y"] = "formation_energy"
  specs["data"]["kwargs"] = None
  
  # sample weighting
  specs["weight"] = dict()
  specs["weight"]["method"] = "wHullDist"
  specs["weight"]["kwargs"] = dict()
  specs["weight"]["kwargs"]["A"] = 0.0
  specs["weight"]["kwargs"]["B"] = 1.0
  specs["weight"]["kwargs"]["kT"] = 0.01
  
  # cross validation
  specs["cv"] = dict()
  specs["cv"]["method"] = "LeaveOneOut"
  specs["cv"]["penalty"] = 0.0
  
  input["problem_specs"] = specs
  
  # regression estimator
  input["estimator"] = dict()
  input["estimator"]["method"] = "LinearRegression"
  
  # feature selection
  input["feature_selection"] = dict()
  input["feature_selection"]["method"] = "DirectSelection"
  d = {
    "use_saved_estimator" : False,
    "population" : [
      {"from_halloffame" : "my_halloffame.pkl", "individuals" : casm.NoIndent([0, 2, 5]) },
      {"bitstring": "111000"},
      {"indices" : casm.NoIndent([1, 2, 3, 12, 13, 21])}
    ]
  }
  input["feature_selection"]["kwargs"] = d
  
  # hall of fame
  input["n_halloffame"] = 25
  
  return input
  

def print_input_help():
  
  print \
  """
  Settings file description:
  ------------------------------------------------------------------------------
  {
  
  # The "problem_specs" specify which data to use and how to score candidate
  # solutions. It consists primarily of "data", "weight", and "cv" settings.
  
    "problem_specs" : {
  
  # Specifies the data to use for learning
  #
  #   A filename and filetype describing where to find data to use for learning. 
  #   Also includes the labels of the sample ('X') and target ('y') data to use 
  #   and kwargs containing additional options.
  #
  #
  # Object attributes
  # -----------------
  #
  # filename: string, optional, default="train"
  #   The path to a file containing the training data
  #
  # filetype: string, optional, default="selection"
  #   Options:
  #     "selection": path to a CASM selection file. Only the selected 
  #        configurations included in the file will be used for training. The
  #        data, typically correlations and a property, can be queried separately
  #        and do not need to be included. If they do exist, the data in the file
  #        will be used.
  #     "csv": path to a CSV file
  #     "json": path to a JSON file 
  #
  # X: string, optional, default="corr"
  #   The name of sample data. Expected to take the form "X(0)", "X(1)", etc...
  #
  # y: string, optional, default="formation_energy"
  #   The name of the target value to train with.
  #
  # kwargs: dict or null, optional, default=dict()
  #   Additional parameters to be used to get training data. 
  #  
  #   Options for 'filetype' "selection":
  #     "project_path": indicate the path to a CASM project. Default null uses
  #       the CASM project containing the current working directory.
  #   
  #   Options for 'filetype' "csv":
  #     Any options to pass to pandas.read_csv
  #
  #   Options for 'filetype' "json":
  #     Any options to pass to pandas.read_json
  
      "data" : {
        "filename": "train",
        "filetype": "selection",
        "X": "corr",
        "y": "formation_energy",
        "kwargs": null
      },
  
  # Method to use for weighting training data. 
  #
  #   If weights are included, then the linear model is changed from
  #     X*b = y  ->  L*X*b = L*y, 
  #
  #   where 'X' is the correlation matrix of shape (Nvalue, Nbfunc),
  #   and 'property' is a vector of Nvalue calculated properties, and 
  #   W = L*L.transpose() is the weight matrix.
  #
  #   By default, W = np.matlib.eye(Nvalue) (unweighted).
  #
  #   If the weighting method provides 1-dimensional input (this is typical), in
  #   a numpy array called 'w':
  #     W = np.diag(w)*Nvalue/np.sum(w)
  #
  #   If the 'custom2d' method is used, the input W_in must by Hermitian, 
  #   positive-definite and is normalized by:
  #     W = W_in*Nvalue/np.sum(W_in)
  #
  #
  # Object attributes
  # -----------------
  #
  # method: string, optional, default=null
  #   The weighting method to use
  #  
  #   Options:
  #   'wHullDist': Weight according to w_i = A*exp(-hull_dist/kT) + B, where A, B, 
  #     and kT are user-defined kwargs parameters, and hull_dist is the distance 
  #     from the convex hull of the training data, using the kwarg "hull_selection"
  #     to determine which selection of configurations to use to find the hull.
  #     The default hull_selection is "CALCULATED".
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
  #
  #
  # kwargs: dict or null, optional, default=dict()
  #   Additional parameters to be used to construct the estimator 
  #
  #   Options: (as described above)
  #     "A": float
  #     "B": float
  #     "kT": float
  #     "E0": float
  #     "hull_selection": string
    
      "weight": {
        "method": null, 
        "kwargs": null
      },
  
  # A scikit-learn cross validation method to use to generate cross validation sets.
  #
  #   The cv score reported is:
  #
  #     cv = sqrt(np.mean(scores)) + (Number of non-zero ECI)*penalty, 
  #
  #   where 'scores' is an array containing the mean squared error calculated for  
  #   each training/testing set, '(Number of non-zero ECI)' is the number of basis 
  #   functions with non-zero ECI, and 'penalty' is the user-input penalty per basis 
  #   function (default=0.0).
  #
  #
  # Object attributes
  # -----------------
  #
  # method: string
  #   A scikit-learn cross validation method. 
  #  
  #   Options include 'KFold', 'ShuffleSplit', 'LeaveOneOut', etc.
  #     See: http://scikit-learn.org/stable/modules/cross_validation.html
  #
  #     Note: The 'LinearRegression' estimator is implemented using 
  #     casm.learn.linear_model.LinearRegressionForLOOCV', which solves X*b=y using:
  #       S = np.linalg.pinv(X.transpose().dot(X)).dot(X.transpose())
  #       b = np.dot(S, y)
  #       H = np.dot(X, S)
  #       y_pred = np.dot(H, y)
  #
  #     Note: When the estimator is 'LinearRegression', the 'LeaveOneOut' 
  #     cross-validation score is calculated via:
  #
  #       LOOCV = np.mean(((y - y_pred)/(1.0 - np.diag(H)))**2)
  #
  # kwargs: dict or null, optional, default=dict()
  #   Additional parameters to be used to construct the cross-validation method constructor.
  #
  # penalty: float, optional, default=0.0
  #   The CV score is increased by 'penalty*(number of selected basis function)'
  #
  #
  
      "cv": {
        "method": "LeaveOneOut", 
        "kwargs": null,
        "penalty": 0.0
      }, 
  
  # Problem specs filename.
  #
  # Optional. Name to use for file storing the training data and CV train/test 
  # sets. The default is determined from the input filename, for example, 
  # 'my_input_specs.pkl' is used if the input file is named 'my_input.json'.
  
      "specs_filename": "problem_specs.pkl"
    
    },
  
  # A scikit-learn linear model estimator.
  #
  #
  # Object attributes
  # -----------------
  #
  # method: string
  #   A scikit-learn linear model estimator. 
  #  
  #   Options: 'LinearRegression', 'Ridge', 'Lasso', 'LassoCV', etc.
  #     See: http://scikit-learn.org/stable/modules/linear_model.html
  #
  #   Note: The 'LinearRegression' estimator is implemented using 
  #   casm.learn.linear_model.LinearRegressionForLOOCV', which solves X*b=y using:
  #     b = np.dot(S, y)
  #     S = np.linalg.pinv(X.transpose().dot(X)).dot(X.transpose())
  #     y_pred = np.dot(H, y)
  #     H = np.dot(X, S)
  #
  # kwargs: dict or null, optional, default=dict()
  #   Additional parameters to be used to construct the estimator 
  #
  #   Options for "LinearRegression":
  #     "pinv": bool, optional, default=True
  #       If True, use the pseudo-inverse via np.linalg.pinv; else use np.linalg.inv.
  #   
  #   Options for other methods:
  #     Any options to pass to the estimator construtor.
  #
  #   By default, the kwarg "fit_intercept" is set to False.
  #
  
    "estimator": {
      "method": "LinearRegression", 
      "kwargs": null
    },
  
    
  # A scikit-learn or casm.feature_selection feature selection method.
  #
  #
  # Object attributes
  # -----------------
  #
  # method: string
  #   A scikit-learn or casm.feature_selection feature selection method. 
  #  
  #   Options from sklearn.feature_selection: "SelectFromModel", "RFE", etc.
  #     See: http://scikit-learn.org/stable/modules/feature_selection.html
  #
  #   Options from casm.feature_selection: 
  #
  #   "DirectSelection": Allows directly specifying which basis functions should
  #     be included.
  #
  #     Options for "kwargs":
  #
  #       "population": List[dict]
  #          Contains a list of options specifying which individuals to fit.
  #          Options are:
  #
  #            "bitstring": Ex.: {"bitstring" : "01110001100"}
  #              String consisting of '0' and '1', with '1' corresponding to 
  #              selected basis functions. May be shorter than the total number 
  #              of possible basis functions, in which case '0' are effectively 
  #              padded to the end.
  #
  #            "indices": Ex.: {"indices" : [1, 2, 3, 7, 8]}
  #              List of indices of basis functions to be selected.
  #
  #            "from_halloffame": Ex.: { "from_halloffame" : "my_halloffame.pkl", 
  #                                      "individuals" : [0, 2, 5]}
  #               Specifies a hall of fame .pkl file and particular individuals 
  #               in the hall (by index) to include in the population. The 
  #               "individuals" list is optional, with the default behaviour 
  #               including all individuals in the hall of fame.
  #
  #       "use_saved_estimator": boolean, optional, default=False
  #          If True, and individuals in the input population come from a HallOfFame
  #          the estimator method stored in the individual's saved input file will
  #          be used instead of the estimator specified in the current input file.
  #
  #
  #   Evolutionary algorithms, from casm.learn.feature_selection, are implemented
  #   using deap: http://deap.readthedocs.org/en/master/index.html
  #   
  #   "GeneticAlgorithm": Implements deap.algorithms.eaSimple, using selTournament,
  #     for selection, cxUniform for mating, and mutFlipBit for mutation. The
  #     probability of mating and mutating is set to 1.0.
  #
  #     Options for "kwargs":
  #
  #       "n_population": int, optional, default=100
  #          Population size. This many random initial starting individuals are 
  #          created.
  #       
  #       "n_halloffame": int, optional, default=25
  #          Maxsize of the hall of fame which holds the best individuals 
  #          encountered in any generation. Upon completion, the individuals in  
  #          this hall of fame are into your overall casm-learn hall of fame to
  #          be compared to results obtained from other fitting or feature 
  #          selection methods.
  #          
  #       "n_generation": int, optional, default=10
  #          Number of generations between saving the hall of fame.
  #
  #       "n_repetition": int, optional, default=100
  #          Number of repetitions of n_generation generations. Each repetition 
  #          begins with the existing final population.
  #
  #       "n_features_init: int or "all", optional, default=0
  #          Number of randomly selected features to initialize each individual 
  #          with.
  #
  #       "selTournamentSize": int, optional, default=3
  #          Tournament size. A larger tournament size weeds out less fit 
  #          individuals more quickly, while a smaller tournament size weeds out 
  #          less fit individuals more gradually.
  #
  #       "cxUniformProb": float, optional, default=0.5
  #          Probability of swapping bits during mating.
  #
  #       "mutFlipBitProb": float, optional, default=0.01 
  #          Probability of mutating bits
  #
  #       "constraints": dict, optional, default=dict()
  #          Keyword arguments for setting constraints on allowed individuals. 
  #          See below for options.
  #
  #
  #   "IndividualBestFirst": 
  #     Implements a best first search optimization for each individual in the initial 
  #     population. Each individual in the population is minimized by repeatedly begin 
  #     replaced by its most fit child.
  #
  #     Children are generated by generating all the individual that differ from
  #     the parent by +/- 1 selected feature. 
  #
  #
  #     Options for "kwargs":
  #
  #       "n_population": int, optional, default=100 
  #          Population size. This many random initial starting individuals are 
  #          minimized and the results saved in the hall of fame.
  #
  #       "n_halloffame": int, optional, default=25
  #          Maxsize of the hall of fame which holds the best individuals 
  #          encountered in any generation. Upon completion, the individuals in  
  #          this hall of fame are into your overall casm-learn hall of fame to
  #          be compared to results obtained from other fitting or feature 
  #          selection methods.
  #
  #       "n_generation": int, optional, default=10
  #          Number of generations between saving the hall of fame.
  #
  #       "n_repetition": int, optional, default=100
  #          Number of repetitions of n_generation generations. Each repetition 
  #          begins with the existing final population.
  #
  #       "n_features_init: int, optional, default=5
  #          Number of randomly selected features to initialize each individual 
  #          with.
  #
  #       "constraints": dict, optional, default=dict()
  #          Keyword arguments for setting constraints on allowed individuals. 
  #          See below for options.
  #
  #
  #   "PopulationBestFirst": 
  #     Implements a best first search optimization for a population of individual
  #     solutions. Each individual is associated with a 'status' that indicates
  #     whether it has had children yet or not. At each step, the most fit individual
  #     that hasn't had children has children and the population is updated to keep
  #     only the 'n_population' most fit individuals. The algorithm stops when all 
  #     individuals in the population have had children.
  #
  #     Children are generated by generating all the individual that differ from
  #     the parent by +/- 1 selected feature. 
  #     
  #
  #     Options for "kwargs":
  #
  #       "n_population": int, optional, default=100 
  #          Population size. This many random initial starting individuals are 
  #          included in the starting population, which is minimized, and the 
  #          results are saved in the hall of fame.
  #
  #       "n_halloffame": int, optional, default=25
  #          Maxsize of the hall of fame which holds the best individuals 
  #          encountered in any generation. Upon completion, the individuals in  
  #          this hall of fame are into your overall casm-learn hall of fame to
  #          be compared to results obtained from other fitting or feature 
  #          selection methods.
  #
  #       "n_generation": int, optional, default=10
  #          Number of generations between saving the hall of fame.
  #
  #       "n_repetition": int, optional, default=100
  #          Number of repetitions of n_generation generations. Each repetition 
  #          begins with the existing final population.
  #
  #       "n_features_init: int, optional, default=5
  #          Number of randomly selected features to initialize each individual 
  #          with.
  #
  #       "constraints": dict, optional, default=dict()
  #          Keyword arguments for setting constraints on allowed individuals. 
  #          See below for options.
  #
  #   The evolutionary algorithms have an optional set of "constraints" parameters
  #   that may restrict the number of basis functions selected to some range, or
  #   enforce some basis functions to have or not have coefficients:
  #
  #   Options for "constraints":
  #     "n_features_min": int, optional, default=1
  #        The minimum allowed number of selected features. Must be >=1.
  #
  #     "n_features_max": int or str, optional, default="all"
  #        The maximum allowed number of selected features. String "all" for no limit.
  #
  #     "fix_on": 1d array-like of int, optional, default=[]
  #        The indices of features to fix on
  #  
  #     "fix_off": 1d array-like of int, optional, default=[]
  #        The indices of features to fix off
  
    "feature_selection" : {
      "method": "GeneticAlgorithm",
      "kwargs": {
        "n_population": 100,
        "n_generation": 10,
        "n_repetition": 100,
        "Nbunc_init": 0,
        "selTournamentSize": 3,
        "cxUniformProb": 0.5,
        "mutFlipBitProb": 0.01,
        "constraints": {
          "n_features_min": 1,
          "n_features_max": "all",
          "fix_on": [],
          "fix_off": []
        }
      }
    },
  
  # Hall of fame size. 
  #
  # The number of individuals to store in 'halloffame.pkl',
  # as determined by CV score. Default=25.
  
    "n_halloffame": 25
  
  }
  ------------------------------------------------------------------------------
  """
  
  
class FittingData(object):
  """ 
  FittingData holds feature values, target values, sample weights, etc. used
  to solve:
  
    L*X * b = L*y 
  
  a weighted linear model where the weights are given by W = L * L.transpose().
    
  Attributes
  ----------
    
    X: array-like of shape (n_samples, n_features)
      The training input samples (correlations).
    
    y: array-like of shape: (n_samples, 1)
      The target values (property values).
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
    
    n_samples: int
      The number of samples / target values (number of rows in X)
    
    n_features: int
      The number of features (number of columns in X)
    
    W: array-like of shape: (n_samples, n_samples)
      Contains sample weights. 
    
    L: array-like of shape: (n_samples, n_samples)
      Used to generate weighted_X and weighted_y, W = L * L.transpose(). 
    
    weighted_X: array-like of shape: (n_samples, n_features)
      Weighted training input data, weighted_X = L*x.
    
    weighted_y: array-like of shape: (n_samples, 1)
      Weighted target values, weighted_y = L*y. 
      
    scoring: string, callable or None, optional, default: None
      A string or a scorer callable object / function with signature 
      scorer(estimator, X, y). The parameter for sklearn.cross_validation.cross_val_score,
      default = None, uses estimator.score().
    
    penalty: float, optional, default=0.0
      The CV score is increased by 'penalty*(number of selected basis function)'
  """
  
  def __init__(self, X, y, cv, sample_weight=[], scoring=None, penalty=0.0):
    """
    Arguments
    ---------
    
      X: array-like of shape (n_samples, n_features)
        The training input samples (correlations).
      
      y: array-like of shape: (n_samples, 1)
        The target values (property values).
      
      cv: cross-validation generator or an iterable
        Provides train/test splits
      
      sample_weight: None, 1-d array-like of shape: (n_samples, 1), or 2-d array-like of shape: (n_samples, n_samples)
        Sample weights.
        
        if sample_weight is None: (default, unweighted)
          W = np.matlib.eye(N) 
        if sample_weight is 1-dimensional:
          W = np.diag(sample_weight)*Nvalue/np.sum(sample_weight) 
        if sample_weight is 2-dimensional (must be Hermitian, positive-definite):
          W = sample_weight*Nvalue/np.sum(sample_weight) 
      
      scoring: string, callable or None, optional, default=None
        A string or a scorer callable object / function with signature 
        scorer(estimator, X, y). The parameter for sklearn.cross_validation.cross_val_score,
        default = None, uses estimator.score().
        
      penalty: float, optional, default=0.0
        The CV score is increased by 'penalty*(number of selected basis function)'
    """
    self.X = X
    self.y = y
    
    # Number of configurations and basis functions
    self.n_samples, self.n_features = self.X.shape
    
    # weight
    self.weighted_y, self.weighted_X, self.W, self.L = casm.learn.tools.set_sample_weight(
      sample_weight, X=self.X, y=self.y)
    
    # cv sets
    self.cv = cv
    
    # scoring
    self.scoring = scoring
    
    # penalty
    self.penalty = penalty


def make_fitting_data(input, save=True, verbose=True, read_existing=True, input_filename=None):
  """ 
  Construct a FittingData instance, either by reading existing 'fit_data.pkl',
  or from an input settings.
  
  Arguments
  ---------
    
    input: dict
      The input settings as a dict
    
    save: boolean, optional, default=True
      Save a pickle file containing the training data and scoring metric. The file
      name, which can be specified by input["fit_data_filename"], defaults to 
      splitext(basename(input_filename))[0] + "_specs.pkl", i.e. 'my_input_specs.pkl'
      if the input file is named 'my_input.json'. If no input_filename is given,
      then the default is "problem_specs.pkl".
    
    verbose: boolean, optional, default=True
      Print information to stdout.
    
    read_existing: boolean, optional, default=True
      If it exists, read the pickle file containing the training data and scoring 
      metric. The file name, which can be specified by input["fit_data_filename"], 
      defaults to "fit_data.pkl".
    
    input_filename: str, optional
      Used to determine the name of the 'specs' file saved.
  
  
  Returns
  -------
    
    fdata: casm.learn.FittingData
      A FittingData instance constructed based on the input parameters.
      
  """
  if "problem_specs" not in input:
    input["problem_specs"] = dict()
  
  specs = input["problem_specs"]
  
  # set data defaults if not provided
  if "data" not in specs:
    specs["data"] = dict()
  if "kwargs" not in specs["data"] or specs["data"]["kwargs"] is None:
    specs["data"]["kwargs"] = dict()
  
  # set weight defaults if not provided
  sample_weight = None
  if "weight" not in specs:
    specs["weight"] = dict()
  if "kwargs" not in specs["weight"] or specs["weight"]["kwargs"] is None:
    specs["weight"]["kwargs"] = dict()
  if "method" not in specs["weight"]:
    specs["weight"]["method"] = None
  
  # set cv kwargs defaults
  if "kwargs" not in specs["cv"] or specs["cv"]["kwargs"] is None:
    specs["cv"]["kwargs"] = dict()
  
  # property, weight, and cv inputs should remain constant
  # estimator and feature_selection might change
  default_filename = "problem_specs.pkl"
  if input_filename is not None:
    default_filename = splitext(basename(input_filename))[0] + "_specs.pkl"
  fit_data_filename = specs.get("specs_filename", default_filename)
  
  if read_existing and os.path.exists(fit_data_filename):
    if verbose:
      print "# Reading existing fitting data from:", fit_data_filename
    fdata = pickle.load(open(fit_data_filename, 'rb'))
    if verbose:
      print "#   DONE\n"
    
    s = "Problem specifications have changed.\n\n" + \
        "To proceed with the existing specs adjust your input settings \"problem_specs\" to match.\n" + \
        "To proceed with the new specs run in a new directory or delete '" + fit_data_filename + "'."
    
    def check_input(name):
      if fdata.input["problem_specs"][name] != specs[name]:
        print "ERROR: Input file and stored data differ. Input '" + name + "' has changed."
        print "Stored data:\n", json.dumps(fdata.input["problem_specs"][name], indent=2)
        print "Input:\n", json.dumps(specs[name], indent=2)
        print s
        exit()
    
    for name in ["data", "cv", "weight"]:
      check_input(name)
    
    
  else:
    
    ## get data ####
    
    filename = specs["data"].get("filename", "train")
    data_type = specs["data"].get("type", "selection").lower()
    X_name = specs["data"].get("X", "corr")
    y_name = specs["data"].get("y", "formation_energy")
    hull_dist_name = "hull_dist"
    
    if data_type == "selection":
      
      # read training set
      proj = Project(specs["data"]["kwargs"].get("project_path", None), verbose=verbose)
      
      sel = Selection(proj, filename, all=False)
      
      # get property name (required)
      property = specs["data"].get("filename", "formation_energy")
      
      ## if necessary, query data
      columns = [X_name, y_name, 'is_calculated']
      if specs["weight"]["method"] == "wHullDist":
        hull_selection = specs["weight"]["kwargs"].get("hull_selection", "CALCULATED")
        hull_dist_name = "hull_dist(" + hull_selection + ",atom_frac)"
        if verbose:
          print "# wHullDist: Will calculate hull distance:", hull_dist_name
        columns.append(hull_dist_name)
      
      # perform query
      if verbose:
        print "# Querying CASM:", columns
      sel.query(columns)
      if verbose:
        print "#   DONE\n"
      
      data = sel.data
      
    elif data_type.lower() == "csv":
      # populate from csv file
      data = pandas.read_csv(f, **specs["data"]["kwargs"])
    
    elif data_type.lower() == "json":
      # populate from json file
      data = pandas.read_json(self.path, **specs["data"]["kwargs"])
    
    # columns of interest, as numpy arrays
    X = data.loc[:,[x for x in sel.data.columns if re.match(X_name + "\([0-9]*\)", x)]].values
    y = data.loc[:,y_name].values
    if specs["weight"]["method"] == "wHullDist":
      hull_dist = data.loc[:,hull_dist_name]
    
    n_samples = X.shape[0]
    n_features = X.shape[1]
    
    if verbose:
      print "# Target:", y_name
      print "# Training samples:", n_samples
      print "# Features:", n_features
    
    ## weight (optional)
    
    # get kwargs
    weight_kwargs = copy.deepcopy(specs["weight"]["kwargs"])
          
    # use method to get weights
    if specs["weight"]["method"] == "wCustom":
      if verbose:
        print "Reading custom weights"
      sample_weight = data["weight"].values
    elif specs["weight"]["method"] == "wCustom2d":
      if verbose:
        print "Reading custom2d weights"
      cols = ["weight(" + str(i) + ")" for i in xrange(n_samples)]
      sample_weight = data.loc[:,cols].values
    elif specs["weight"]["method"] == "wHullDist":
      sample_weight = casm.learn.tools.wHullDist(hull_dist, **weight_kwargs)
    elif specs["weight"]["method"] == "wEmin":
      sample_weight = casm.learn.tools.wEmin(y, **weight_kwargs)
    elif specs["weight"]["method"] == "wEref":
      sample_weight = casm.learn.tools.wEref(y, **weight_kwargs)
          
    if verbose:
      print "# Weighting:"
      print json.dumps(specs["weight"], indent=2), "\n"
    
    
    ## cv
    cv_kwargs = copy.deepcopy(specs["cv"]["kwargs"])
    
    # get cv method (required user input) 
    cv_method = _find_method([sklearn.cross_validation], specs["cv"]["method"])
    cv = cv_method(n_samples, **cv_kwargs)
    
    if verbose:
      print "# CV:"
      print json.dumps(specs["cv"], indent=2), "\n"
    
    ## scoring
    scoring = sklearn.metrics.make_scorer(sklearn.metrics.mean_squared_error, greater_is_better=True)
    
    ## penalty
    penalty = specs["cv"].get("penalty", 0.0)
    
    fdata = casm.learn.FittingData(X, y, cv, 
      sample_weight=sample_weight, scoring=scoring, penalty=penalty)
    
    fdata.input = dict()
    fdata.input["problem_specs"] = specs
    
    if save == True:
      pickle.dump(fdata, open(fit_data_filename, 'wb'))
  
  # during runtime only, if LinearRegression and LeaveOneOut, update fdata.cv and fdata.scoring
  # to use optimized LOOCV score method
  if input["estimator"].get("method", "LinearRegression") == "LinearRegression" and specs["cv"]["method"] == "LeaveOneOut":
    fdata.scoring = None
    fdata.cv = casm.learn.cross_validation.LeaveOneOutForLLS(fdata.weighted_y.shape[0])
  
  return fdata 
  

def make_estimator(input, verbose = True):
  """
  Construct estimator object from input settings.
  
  Arguments
  ---------
    
    input: dict
      The input settings as a dict
    
    verbose: boolean, optional, default=True
      Print information to stdout.
    
  
  Returns
  -------
    
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
  
  """""
  
  ## estimator
  
  # get kwargs (default: fit_intercept=False)
  kwargs = copy.deepcopy(input["estimator"].get("kwargs", dict()))
  if "fit_intercept" not in kwargs:
    kwargs["fit_intercept"] = False
  
  # Use casm version of LinearRegression
  if input["estimator"]["method"] == "LinearRegression":
    estimator_method = casm.learn.linear_model.LinearRegressionForLOOCV
  else:
    estimator_method = _find_method([sklearn.linear_model], input["estimator"]["method"])
  estimator = estimator_method(**kwargs)
  
  if verbose:
    print "# Estimator:"
    print json.dumps(input["estimator"], indent=2), "\n"
    
  
  return estimator


def make_selector(input, estimator, scoring=None, cv=None, penalty=0.0, verbose=True):
  """ 
  Construct selector object from input settings
  
  Arguments
  ---------
    
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
  
    scoring: string, callable or None, optional, default: None
      A string or a scorer callable object / function with signature 
      scorer(estimator, X, y). The parameter for sklearn.cross_validation.cross_val_score,
      default = None, uses estimator.score().
    
    cv: cross-validation generator or an iterable
      Provides train/test splits
    
    penalty: float, optional, default=0.0
      The CV score is increased by 'penalty*(number of selected basis function)'
    
    verbose: boolean, optional, default=True
      Print information to stdout.
    
  
  Returns
  -------
    
    selector:  selector object implementing 'fit' and having either a 
               'get_support()' or 'get_halloffame()' member
      The feature selector specified by the input settings.
    
  """
  
  # read input, construct and return a feature selector
  #
  # The feature selector should act like a sklearn.feature_selection class and
  # inherit from sklearn.base.BaseEstimator and sklearn.feature_selection.SelectorMixin,
  
  kwargs = copy.deepcopy(input["feature_selection"].get("kwargs", dict()))
  if kwargs is None:
    kwargs = dict()
  if "evolve_params_kwargs" in kwargs:
    if "halloffame_filename" not in kwargs["evolve_params_kwargs"]:
      kwargs["evolve_params_kwargs"]["halloffame_filename"] = input.get("halloffame_filename", "halloffame.pkl")
    if "n_halloffame" not in kwargs["evolve_params_kwargs"]:
      kwargs["evolve_params_kwargs"]["n_halloffame"] = input.get("n_halloffame", 25)
  
  if verbose:
    print "# Feature Selection:"
    print json.dumps(input["feature_selection"], indent=2), "\n"
  
  mods = [casm.learn.feature_selection, sklearn.feature_selection]
  
  selector_method = _find_method(mods, input["feature_selection"]["method"])
  
  # check if 'cv', 'scoring', 'penalty' are allowed kwargs 
  arg_count = selector_method.__init__.func_code.co_argcount
  allowed_kwargs = selector_method.__init__.func_code.co_varnames[:arg_count]
  
  if "cv" in allowed_kwargs:
    kwargs["cv"] = cv
  if "scoring" in allowed_kwargs:
    kwargs["scoring"] = scoring
  if "penalty" in allowed_kwargs:
    kwargs["penalty"] = penalty
  if "verbose" in allowed_kwargs:
    kwargs["verbose"] = verbose
  
  selector = selector_method(estimator, **kwargs)
  
  return selector


def add_individual_detail(indiv, estimator, fdata, input, selector=None):
  """
  Adds attributes to an individual describing the details of the method used 
  calculate it and the it's prediction ability. 
  
  Adds the attributes:
  
    eci: List[(int, float)]
      A list of tuple containing the basis function index and coefficient value
      for basis functions with non-zero coefficients: [(index, coef), ...]
    
    rms: float
      The root mean square prediction error of the unweighted problem
    
    wrms: float
      The root mean square prediction error of the weighted problem
    
    estimator_method: string
      The estimator method class name
    
    feature_selection_method: string,
      The feature_selection method class name
    
    note: string
      Additional notes describing the individual
    
    input: dict
      The input settings
    
    id: uuid.UUID
      Unique id, generated by uuid.uuid4()
  
  
  Arguments
  ---------
    
    indiv: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
     
    estimator:  estimator object implementing 'fit'
      The estimator specified by the input settings.
    
    fdata: casm.learn.FittingData
      A FittingData instance containing the problem data.
    
    input: dict
      The input settings
    
    selector:  selector object implementing 'fit' and having either a 
               'get_support()' or 'get_halloffame()' member, optional, default=None
      The feature selector specified by the input settings.
    
    
  Note
  ----
    Individuals should already have a 'fitness.values' attribute with the cv
    score. As is the convention in the 'deap' package, the 'fitness.values' 
    attribute is a tuple with the first element being the cv score.
    
  """
  # eci
  estimator.fit(fdata.weighted_X[:,casm.learn.tools.indices(indiv)], fdata.weighted_y)
  indiv.eci = casm.learn.tools.eci(indiv, estimator.coef_)
  
  # rms and wrms
  indiv.rms = sqrt(sklearn.metrics.mean_squared_error(
    fdata.y, estimator.predict(fdata.X[:,casm.learn.tools.indices(indiv)])))
  indiv.wrms = sqrt(sklearn.metrics.mean_squared_error(
    fdata.weighted_y, estimator.predict(fdata.weighted_X[:,casm.learn.tools.indices(indiv)])))
  
  # estimator (hide implementation detail)
  indiv.estimator_method = type(estimator).__name__
  if indiv.estimator_method == "LinearRegressionForLOOCV":
    indiv.estimator_method = "LinearRegression"
  
  # feature_selection
  if selector is not None:
    indiv.feature_selection_method = type(selector).__name__
  elif not hasattr(indiv, "feature_selection_method"):
    indiv.feature_selection_method = "DirectSelection"
  
  # note
  indiv.note = input.get("note", "")
  
  # input settings
  indiv.input = input
  
  # unique id
  indiv.id = uuid.uuid4()
  
  return indiv


def bitstr(indiv, n_bits_max=None):
  if n_bits_max is None:
    n_bits_max = len(indiv)
  bitstr = ""
  for j in range(min(len(indiv),n_bits_max)):
    if indiv[j]:
      bitstr += '1'
    else:
      bitstr += '0'
  if len(indiv) > n_bits_max:
    bitstr += "..." 
  return bitstr
  

def print_population(pop):
  """ 
  Print all individual in a population.
  
  Example:
    
    Index: Selected                                    #Selected    CV           
    ----------------------------------------------------------------------------------------------------
        0: 0111011110000111000001001100000100010000... 25           0.015609282  
        1: 0111011110000111000001001101000100010000... 25           0.015611913  
        2: 0111011110000111000001001100000100010000... 24           0.015619745  
    ...
  
  
  Arguments
  ---------
    
    pop: List-like of List[bool] of length n_features
      A population, a list-like container of individuals. Each individual is a 
      boolean list of shape [n_features], in which an element is True iff its 
      corresponding feature is selected for retention.
  """
  print "{0:5}: {1:43} {2:<12} {3:<12}".format("Index", "Selected", "#Selected", "CV")
  print "-"*100
  form_str = "{0:5}: {1} {2:<12} {3:<12.8g}"
  for i in range(len(pop)):
    print form_str.format(i, bitstr(pop[i], 40), sum(pop[i]), pop[i].fitness.values[0])


def to_json(index, indiv):
  """
  Serialize an individual to JSON records. 
  
  Keys in each record: 
    "index": index of individual in hall of fame
    "id": uuid str for individual
    "selected": str of 0's and 1's indicating selected features
    "n_selected": number of selected features
    "cv": CV score
    "rms": root-mean-square error
    "wrms": weighted root-mean-square error
    "estimator_method": name of estimator method
    "features_selection_method": name of feature selection method
    "note": a descriptive note
    "eci": List of (feature index, coefficient value) pairs
    "input": input settings dict
    
  Arguments
  ---------
    
    index: int
      The index in hall of fame of the individual
      
    indiv: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True 
      iff its corresponding feature is selected for retention.
  
  
  Note
  ----
  
    ECI are serialized using cls=casm.NoIndent, so when writing with json.dump or 
    json.dumps, include 'cls=casm.NoIndentEncoder'.
  
  """
  d = dict()
  d["selected"] = bitstr(indiv)
  d["index"] = index
  d["id"] = str(indiv.id)
  d["n_selected"] = sum(indiv)
  d["cv"] = indiv.fitness.values[0]
  d["rms"] = indiv.rms
  d["wrms"] = indiv.wrms
  d["estimator_method"] = indiv.estimator_method
  d["feature_selection_method"] = indiv.feature_selection_method
  d["note"] = indiv.note
  d["eci"] = []
  for bfunc in indiv.eci:
    d["eci"].append(casm.NoIndent(bfunc))
  d["input"] = indiv.input
  return d


def to_dataframe(indices, hall):
  """
  Convert hall of fame data to pandas.DataFrame. 
  
  Columns: 
    "index": index of individual in hall of fame
    "id": uuid str for individual
    "selected": str of 0's and 1's indicating selected features
    "n_selected": number of selected features
    "cv": CV score
    "rms": root-mean-square error
    "wrms": weighted root-mean-square error
    "estimator_method": name of estimator method
    "features_selection_method": name of feature selection method
    "note": a descriptive note
    "eci": JSON string of a List of (feature index, coefficient value) pairs
    "input": input settings dict
  
  
  Arguments
  ---------
    
    indices: List[int]
      The indices in hall of fame of the individuals to include
      
    hall: deap.tools.HallOfFame
      A Hall Of Fame of ECI sets
    
  """
  data = [to_json(i, hall[i]) for i in indices]
  for d in data:
    d["eci"] = json.dumps(d["eci"], cls=casm.NoIndentEncoder)
  return pandas.DataFrame.from_records(data)
    


def _print_individual(index, indiv, format=None):
  """ 
  Print all individual in hall of fame.
  
  Index: Selected                                    #Selected    CV          wRMS        
  ----------------------------------------------------------------------------------------------------
      0: 0111011110000111000001001100000100010000... 25           0.015609282  0.014073401 
      1: 0111011110000111000001001101000100010000... 25           0.015611913  0.014045382 
      2: 0111011110000111000001001100000100010000... 24           0.015619745  0.01411583  
  ...
  """
  if format is None:
    if len(indiv) > 40:
      bitstr_len = 43
    else:
      bitstr_len = len(indiv)
    form_str = "{0:5}: {1:<" + str(max([12, bitstr_len])) + "} {2:<12} {3:<12.8g} {4:<12.8g} {5:<12.8g} {6:<24} {7:<24} {8}"
    print form_str.format(index, bitstr(indiv,40), sum(indiv), indiv.fitness.values[0], 
      indiv.rms, indiv.wrms, indiv.estimator_method, indiv.feature_selection_method, indiv.note)
    return
    
  elif format.lower() == "json":
    print json.dumps(to_json(index,indiv), indent=2, cls=casm.NoIndentEncoder)
    return
    
  elif format.lower() == "details":
    print "##"
    print "Index:", index
    print "ID:", indiv.id
    print "Selected:", bitstr(indiv)
    print "#Selected:", sum(indiv)
    print "CV:", indiv.fitness.values[0]
    print "RMS:", indiv.rms
    print "wRMS:", indiv.wrms
    print "Estimator:", indiv.estimator_method
    print "FeatureSelection:", indiv.feature_selection_method
    print "Note:", indiv.note
    print "ECI:\n"
    print_eci(indiv.eci)
    print "Input:\n", json.dumps(indiv.input, indent=2)
    return


def _print_halloffame_header(hall):
  """ 
  Print header for hall of fame.
  """
  if len(hall[0]) > 40:
    bitstr_len = 43
  else:
    bitstr_len = len(hall[0])
  print ("{0:5}: {1:<" + str(max([12, bitstr_len])) + "} {2:<12} {3:<12} {4:<12} {5:<12} {6:<24} {7:<24} {8}").format(
    "Index", "Selected", "#Selected", "CV", "RMS", "wRMS", "Estimator", "FeatureSelection", "Note")
  print "-"*(6+max([12, bitstr_len])+13*5+25*2)


def print_individual(hall, indices, format=None):
  """ 
  Print selected individuals from hall of fame.
  
  Arguments:
  ----------
  
    hall: deap.tools.HallOfFame
      A Hall Of Fame of ECI sets
    
    index: List[int]
      Indices of individual in hall to be printed
    
    format: str, optional, default=None
      Options: 
        None:      to print summary only
        "details": to print more
        "csv":     to print as CSV 
        "json":    to print as JSON  
  """
  if format is None:
    _print_halloffame_header(hall)
    for index in indices:
      _print_individual(index, hall[index], format=format)
    return
    
  elif format.lower() == "json":
    h = []
    for index in indices:
      d = to_json(index, hall[index])
      h.append(d)
    print json.dumps(h, indent=2, cls=casm.NoIndentEncoder)
  
  elif format.lower() == "csv":
    df = to_dataframe(range(len(hall)), hall) 
    print df.to_csv()
    
  elif format.lower() == "details":
    for index in indices:
      _print_individual(index, hall[index], format=format)
    print ""
  return


def print_halloffame(hall, format=None):
  """ 
  Print all individual in hall of fame.
  
  Arguments:
  ----------
  
    hall: deap.tools.HallOfFame
      A Hall Of Fame of ECI sets
    
    format: str, optional, default=None
      Options: 
        None:      to print summary only
        "details": to print more
        "csv":     to print as CSV 
        "json":    to print as JSON 
  """
  if format is None:
    _print_halloffame_header(hall)
    for index, indiv in enumerate(hall):
      _print_individual(index, indiv, format=format)
      
    return
    
  elif format.lower() == "json":
    h = []
    for index, indiv in enumerate(hall):
      d = to_json(index, indiv)
      h.append(d)
    print json.dumps(h, indent=2, cls=casm.NoIndentEncoder)
      
    
  elif format.lower() == "csv":
    df = to_dataframe(range(len(hall)), hall) 
    print df.to_csv()
    
  elif format.lower() == "details":
    for index, indiv in enumerate(hall):
      _print_individual(index, indiv, format=format)
      print ""
    
    return


def print_eci(eci):
  """
  Print ECI.
  
  Format:
      1: -1.53460558686
      2:  0.574571156376
      3:  1.04379783648
  ...
  
  Arguments
  ---------
    
    eci: List[(int, float)]
      A list of tuple containing the basis function index and coefficient value
      for basis functions with non-zero coefficients: [(index, coef), ...]
  
  """
  for bfunc in eci:
    print "{index:>5}: {value:< .12g}".format(index=bfunc[0], value=bfunc[1])


def open_halloffame(halloffame_filename, verbose=False):
  """
  Open hall of fame from .pkl file.
  
  Arguments
  ---------
    
    halloffame_filename: string
      Name of .pkl file containing the hall of fame
    
    verbose: boolean, optional, default=True
      Print information to stdout.
  
  Returns
  -------
    
    hall: deap.tools.HallOfFame
      A Hall Of Fame of ECI sets
    
  """
  if verbose:
    print "Loading Hall of Fame:", halloffame_filename
        
  with open(halloffame_filename, 'rb') as f:
    hall = pickle.load(f)
  
  # for backwards compatibility, add id if not existing
  for indiv in hall:
    if not hasattr(indiv, "id"):
      indiv.id = uuid.uuid4()
  
  return hall


def save_halloffame(hall, halloffame_filename, verbose=False):
  """
  Save hall of fame as .pkl file.
  
  Arguments
  ---------
    
    hall: deap.tools.HallOfFame
      A Hall Of Fame of ECI sets
    
    halloffame_filename: string
      Name of .pkl file to write the hall of fame
    
    verbose: boolean, optional, default=True
      Print information to stdout.
    
  """
  if verbose:
    print "\nPickling Hall of Fame to:", halloffame_filename
  f = open(halloffame_filename, 'wb')
  pickle.dump(hall, f)
  f.close()


