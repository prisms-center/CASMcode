from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

# conda's current version of pandas raises these warnings, but they are safe
# see: https://stackoverflow.com/questions/40845304
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import sklearn.linear_model
import sklearn.model_selection
import sklearn.metrics
import random, re, time, os, types, json, pickle, copy, uuid, shutil, tempfile
import numpy as np
import pandas
import six

if six.PY2:
    from funcsigs import signature
else:
    from inspect import signature

from math import sqrt
from os.path import splitext, basename, join

import casm.learn.linear_model
import casm.learn.tools
import casm.learn.selection_wrapper
from casm.misc import noindent
from casm.project import Project, Selection, query, write_eci

def _find_method(mods, attrname):
  for m in mods:
    if hasattr(m, attrname):
      return getattr(m, attrname)
  print("ERROR: Could not find a method named:", attrname)
  print("Tried:", mods)
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
  specs["cv"]["kwargs"]["n_splits"] = 10
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
  specs["cv"]["kwargs"]["n_splits"] = 10
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
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_splits"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
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
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_splits"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
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
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_splits"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
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
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_splits"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
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
  specs["cv"]["method"] = "KFold"
  specs["cv"]["kwargs"] = dict()
  specs["cv"]["kwargs"]["n_splits"] = 10
  specs["cv"]["kwargs"]["shuffle"] = True
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
      {"from_halloffame" : "my_halloffame.pkl", "individuals" : noindent.NoIndent([0, 2, 5]) },
      {"bitstring": "111000"},
      {"indices" : noindent.NoIndent([1, 2, 3, 12, 13, 21])}
    ]
  }
  input["feature_selection"]["kwargs"] = d

  # hall of fame
  input["n_halloffame"] = 25

  return input


def print_input_help():

  print("""
  Settings file description:
  ------------------------------------------------------------------------------
  {

  # The "problem_specs" options specify which data to use and how to score
  # candidate solutions. It consists primarily of "data", "weight", and "cv"
  # settings.

    "problem_specs" : {

  # The "problem_specs"/"data" options specify the data to use for learning
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
  #       the CASM project containing the current working directory. This option
  #       is not currently implemented, but included as a placeholder. Currently
  #       'casm-learn' must be run from inside a CASM project.
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

  # The "problem_specs"/"weight" options specify the method to use for weighting
  # training data.
  #
  #   Ordinary least squares minimizes
  #     (y-X*b).transpose() * (y-X*b)
  #
  #   where 'X' is the correlation matrix of shape (Nvalue, Nbfunc), and 'y'
  #   is a vector of Nvalue calculated properties, and 'b' are the fitting
  #   coefficients (ECI).
  #
  #   Weighted least squares minimizes
  #     (y-X*b).transpose() * W * (y-X*b)
  #
  #   Using the SVD, and given that W is Hermitian:
  #     U * S * U.transpose() == W
  #
  #   Define L such that:
  #     L.transpose() = U * sqrt(S)
  #
  #   Then we can write the weighted least squares problem using:
  #     (y-X*b).transpose() * L.transpose() * L * (y-X*b)
  #
  #   Or:
  #     (L*y-L*X*b).transpose() * (L*y-L*X*b)
  #
  #   So, if weights are included, then the linear model is changed from
  #     X*b = y  ->  L*X*b = L*y
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

  # The "problem_specs"/"cv" options specify how to generate cross validation sets.
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
  #   A scikit-learn or casm cross validation method.
  #
  #   Options include 'KFold', 'ShuffleSplit', 'LeaveOneOut', etc.
  #     See: http://scikit-learn.org/stable/modules/model_selection.html
  #
  #   CASM also provides the following method:
  #    'cvCustom': Read a scikit-learn type 'cv' generator or training/test sets
  #       from a pickle file. This can be used to load the 'cv' data written by
  #       'casm-learn --checkspecs' by using the required kwarg 'filename'.
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
  #   Additional parameters to be used to construct the cross-validation method
  #   constructor.
  #
  # penalty: float, optional, default=0.0
  #   The CV score is increased by 'penalty*(number of selected basis function)'

      "cv": {
        "method": "KFold",
        "kwargs": {
          "n_splits": 10,
          "shuffle": true
        },
        "penalty": 0.0
      },

  # The "problem_specs"/"specs_filename" option:
  #
  # Optional. Name to use for file storing the training data and CV train/test
  # sets. The default is determined from the input filename, for example,
  # 'my_input_specs.pkl' is used if the input file is named 'my_input.json'.

      "specs_filename": "problem_specs.pkl"

    },

  # The "estimator" option specifies a linear model estimator.
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

    "estimator": {
      "method": "LinearRegression",
      "kwargs": null
    },


  # The "feature_selection" option specifies a feature selection method.
  #
  # Object attributes
  # -----------------
  #
  # method: string
  #   A scikit-learn or casm.learn.feature_selection feature selection method.
  #
  #   Options from sklearn.feature_selection: "SelectFromModel", "RFE", etc.
  #     See: http://scikit-learn.org/stable/modules/feature_selection.html
  #
  #   Options from casm.learn.feature_selection:
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
  #       "constraints_kwargs": dict, optional, default=dict()
  #          Keyword arguments for setting constraints on allowed individuals.
  #          See below for options.
  #
  #       "evolve_params_kwargs": dict, optional, default=dict()
  #          Keyword arguments for controlling how long the algorithm runs, how
  #          new random individuals are initialized, when restart files are
  #          written, and the names of the files. See below for options.
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
  #       "constraints_kwargs": dict, optional, default=dict()
  #          Keyword arguments for setting constraints on allowed individuals.
  #          See below for options.
  #
  #       "evolve_params_kwargs": dict, optional, default=dict()
  #          Keyword arguments for controlling how long the algorithm runs, how
  #          new random individuals are initialized, when restart files are
  #          written, and the names of the files. See below for options.
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
  #     Options for "kwargs":
  #
  #       "constraints_kwargs": dict, optional, default=dict()
  #          Keyword arguments for setting constraints on allowed individuals.
  #          See below for options.
  #
  #       "evolve_params_kwargs": dict, optional, default=dict()
  #          Keyword arguments for controlling how long the algorithm runs, how
  #          new random individuals are initialized, when restart files are
  #          written, and the names of the files. See below for options.
  #
  #
  #   The evolutionary algorithms have an optional set of "constraints_kwargs"
  #   parameters that may restrict the number of basis functions selected to some
  #    range, or enforce some basis functions to have or not have coefficients:
  #
  #   Options for "constraints_kwargs":
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
  #
  #
  #   The evolutionary algorithms share an optional set of "evolve_params_kwargs"
  #   parameters that control how long the algorithm runs, how new random
  #   individuals are initialized, when restart files are written, and the names
  #   of the files:
  #
  #   Options for "evolve_params_kwargs":
  #
  #     "n_population": int, optional, default=100
  #        Initial population size. This many random initial starting individuals
  #         are created if no "pop_begin_filename" file exists.
  #
  #     "n_halloffame": int, optional, default=25
  #        Maxsize of the hall of fame which holds the best individuals
  #        encountered in any generation. Upon completion, the individuals in
  #        this hall of fame are into your overall casm-learn hall of fame to
  #        be compared to results obtained from other fitting or feature
  #        selection methods.
  #
  #     "n_generation": int, optional, default=10
  #        Number of generations between saving the hall of fame.
  #
  #     "n_repetition": int, optional, default=100
  #        Number of repetitions of n_generation generations. Each repetition
  #        begins with the existing final population.
  #
  #     "n_features_init: int or "all", optional, default=0
  #        Number of randomly selected features to initialize each individual
  #        with.
  #
  #     "pop_begin_filename": string, optional, default="population_begin.pkl"
  #        Filename suffix where the initial population is read from, if it
  #        exists. For example, if "filename_prefix" is "Ef_kfold10" and
  #        "pop_begin_filename" is "population_begin.pkl", then the initial
  #        population is read from the file "Ef_kfold10_population_begin.pkl".
  #
  #        The population file may contain either a  list of individual,
  #        as written to the "population_end.pkl" file, or a HallOfFame
  #        instance, as written to either an "evolve_halloffame.pkl" file or
  #        overall casm-learn "halloffame.pkl" file.
  #
  #     "pop_end_filename": string, optional, default="population_end.pkl"
  #        Filename where the final population is saved. For example, if
  #        "filename_prefix" is "Ef_kfold10" and "pop_end_filename" is
  #        "population_end.pkl", then the final population is saved to the
  #        file "Ef_kfold10_population_end.pkl".
  #
  #     "halloffame_filename": string, optional, default="evolve_halloffame.pkl"
  #        Filename where a hall of fame is saved holding the best individuals
  #        encountered in any generation. For example, if "filename_prefix" is
  #        "Ef_kfold10" and "halloffame_filename" is "evolve_halloffame.pkl",
  #        then it is saved to the file "Ef_kfold10_evolve_halloffame.pkl".
  #
  #     "filename_prefix": string, optional
  #        Prefix for filenames, default uses input file filename excluding
  #        extension. For example, if input file is named "Ef_kfold10.json", then
  #        "Ef_kfold10_population_begin.pkl", "Ef_kfold10_population_end.pkl", and
  #        "Ef_kfold10_evolve_halloffame.pkl" are used.

    "feature_selection" : {
      "method": "GeneticAlgorithm",
      "kwargs": {
        "selTournamentSize": 3,
        "cxUniformProb": 0.5,
        "mutFlipBitProb": 0.01,
        "constraints": {
          "n_features_min": 1,
          "n_features_max": "all",
          "fix_on": [],
          "fix_off": []
        },
        "evolve_params_kwargs": {
          "n_population": 100,
          "n_generation": 10,
          "n_repetition": 100,
          "n_features_init": 0
        }
      }
    },

  # The "halloffame_filename" option:
  #
  # Optional. Default = "halloffame.pkl"
  # Name to use for file storing the best results obtained to date, as determined
  # by the CV score. This enables comparison of the results of various estimator
  # or feature selection methods.

      "halloffame_filename": "halloffame.pkl"

    },

  # The "n_halloffame" option:
  #
  # Optional. Default = 25
  # The number of individuals to store in the hall of fame.

    "n_halloffame": 25

  # The "checkspecs" option:
  #
  #   Currently, these settings are used with the '--checkspecs' option to control
  #   the output files containing training data (including calculated weights) and
  #   cv generators or training / testing sets.
  #
  #
  # Object attributes
  # -----------------
  #
  # data_filename: string
  #   The path to the file where the training data (including weights) should be
  #   written
  #
  # data_kwargs: dict or null, optional, default=dict()
  #   Additional parameters to be used to write training data.
  #
  #   Options for input 'filetype' "selection":
  #     None.
  #
  #   Options for input 'filetype' "csv":
  #     Any options to pass to pandas.to_csv
  #
  #   Options for input 'filetype' "json":
  #     Any options to pass to pandas.to_json
  #
  # cv_filename: string
  #   The path to the file where the cv generator or train/tests sets should be
  #   written as pickle file
  #
  # cv_kwargs: dict or null, optional, default=dict()
  #   Additional parameters to be used to write cv data using pickle.dump.

    "checkspecs" : {
      "data_filename": "check_train",
      "data_kwargs": null,
      "cv_filename": "check_cv.pkl",
      "cv_kwargs": null
    },

  # The "checkhull" option:
  #
  #   Currently, these settings are used with the '--checkhull' option to
  #   calculate convex hull properties.
  #
  #
  # Object attributes
  # -----------------
  #
  # selection: str, optional, default="ALL"
  #   A CASM selection (either 'casm select' output filename or one of the
  #   standard selections: "MASTER", "CALCULATED", or "ALL") containing all the
  #   configurations to be considered. The DFT convex hull is generated from
  #   the subset of this selection for which 'is_calculated' is true.
  #
  # write_results: bool, optional, default=False
  #   If True, write CASM selection files containing the output data. Output
  #   selection files are named "checkhull_(problem_specs_prefix)_(i)_(selname)",
  #   where 'problem_specs_prefix' is input["problem_specs_prefix"], 'i' is the
  #   index of the individual in the hall of fame, and 'selname' is one of:
  #     "dft_gs" : DFT calculated ground states
  #     "clex_gs" : predicted ground states
  #     "gs_missing" : DFT ground states that are not predicted ground states
  #     "gs_spurious" : Predicted ground states that are not DFT ground states
  #     "uncalculated" : Predicted ground states and near ground states that have not been calculated
  #     "below_hull" : All configurations predicted below the prediction of the DFT hull
  #
  # primitive_only: bool, optional, default=True
  #   If True, only use primitive configurations to construct the convex hull,
  #   else use all selected configurations.
  #
  # uncalculated_range: number, optional, default=0.0
  #   Include all configurations with clex_hull_dist less than this value (+hull_tol)
  #   in the "uncalculated" configurations results. Default only includes predicted
  #   ground states.
  #
  # ranged_rms: List[number], optional, default=[0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
  #   Calculates the root-mean-square error for DFT calculated configurations
  #   within a particular range (in eV/unitcell) of the DFT hull. The list
  #   provides all the ranges for which the RMSE is requested.
  #
  # composition: str, optional, default="atom_frac"
  #   Composition argument use for 'casm query' properties 'hull_dist' and
  #   'clex_hull_dist'. For thermodynamic ground states, use "atom_frac".
  #
  # hull_tol: number, optional, default=proj.settings.data["lin_alg_tol"]
  #   Tolerance used for identify hull states
  #
  # dim_tol: number, optional, default=1e-8
  #   Tolerance for detecting composition dimensionality
  #
  # bottom_tol: number, optional, default=1e-8
  #   Tolerance for detecting which facets form the convex hull bottom

    "checkhull" : {
      "selection": "ALL",
      "write_results": true,
      "primitive_only": true,
      "uncalculated_range": 1e-8,
      "ranged_rms": [0.001, 0.005, 0.01, 0.05, 0.1, 0.5],
      "composition": "atom_frac",
      "hull_tol": 1e-8,
      "dim_tol": 1e-8,
      "bottom_tol": 1e-8
    }

  }
  ------------------------------------------------------------------------------
  """)


def default_filename(prefix, default, suffix):
  """
  Make a default filename from the input file filename.

  Implements:
    if prefix is not None:
      return prefix + suffix
    else:
      return default

  Arguments
  ---------

    prefix: str or None
      The prefix used in determining the default filename.

    default: str
      Filename if input_filename is None

    suffix: str
      If input_filename is not None, append this suffix to the input_filename
      (excluding extension) to make the filename.


  Returns
  --------

    filename: str
      The generated default filename

  """
  if prefix is not None:
    return prefix + suffix
  else:
    return default


def set_input_defaults(input, input_filename=None):
  """
  Set common input defaults. Currently, includes everything except "checkspecs"
  and "checkhull" defaults.

  Arguments
  ---------

    input: dict
      The input settings as a dict

    input_filename: str, optional, default=None
      The input settings filename, which is used in determining the default
      problem specs filename.


  Returns
  ---------

    input: dict
      The input settings as a dict, with defaults added

  """
  if "problem_specs" not in input:
    input["problem_specs"] = dict()

  specs = input["problem_specs"]

  if "problem_specs_prefix" not in specs:
    if input_filename is None:
      specs["problem_specs_prefix"] = ""
    else:
      specs["problem_specs_prefix"] = splitext(basename(input_filename))[0]

  if "specs_filename" not in specs:
    specs["specs_filename"] = default_filename(specs["problem_specs_prefix"], "problem_specs.pkl", "_specs.pkl")

  # set data defaults if not provided
  if "data" not in specs:
    specs["data"] = dict()

  defaults = {
    "filename":"train",
    "filetype":"selection",
    "X":"corr",
    "y":"formation_energy",
  }

  for key, val in six.iteritems(defaults):
    if key not in specs["data"]:
      specs["data"][key] = val

  if "kwargs" not in specs["data"] or specs["data"]["kwargs"] is None:
    specs["data"]["kwargs"] = dict()
  if specs["data"]["filetype"] == "selection":
    if "project_path" not in specs["data"]["kwargs"]:
      specs["data"]["kwargs"]["project_path"] = None

  # set weight defaults if not provided
  sample_weight = None
  if "weight" not in specs:
    specs["weight"] = dict()
  if "kwargs" not in specs["weight"] or specs["weight"]["kwargs"] is None:
    specs["weight"]["kwargs"] = dict()
  if "method" not in specs["weight"]:
    specs["weight"]["method"] = None
  if specs["weight"]["method"] == "wHullDist":
    if "hull_selection" not in specs["weight"]["kwargs"]:
      specs["weight"]["kwargs"]["hull_selection"] = "CALCULATED"

  # set cv defaults
  if "kwargs" not in specs["cv"] or specs["cv"]["kwargs"] is None:
    specs["cv"]["kwargs"] = dict()
  if "penalty" not in specs["cv"]:
    specs["cv"]["penalty"] = 0.0

  # set estimator defaults
  if "method" not in input["estimator"]:
    input["estimator"]["method"] = "LinearRegression"
  if "kwargs" not in input["estimator"] or input["estimator"]["kwargs"] is None:
    input["estimator"]["kwargs"] = dict()
  if "fit_intercept" not in input["estimator"]["kwargs"]:
    input["estimator"]["kwargs"]["fit_intercept"] = False

  # set feature_selection defaults
  if "kwargs" not in input["feature_selection"] or input["feature_selection"]["kwargs"] is None:
    input["feature_selection"]["kwargs"] = dict()
  kwargs = input["feature_selection"]["kwargs"]
  if "evolve_params_kwargs" in kwargs:
    evolve_kwargs = kwargs["evolve_params_kwargs"]
    if "halloffame_filename" not in evolve_kwargs:
      evolve_kwargs["halloffame_filename"] = "evolve_halloffame.pkl"
    if "n_halloffame" not in evolve_kwargs:
      evolve_kwargs["n_halloffame"] = 25
    if "filename_prefix" not in evolve_kwargs and input_filename is not None:
      evolve_kwargs["filename_prefix"] = specs["problem_specs_prefix"]


  # hall of fame
  if "halloffame_filename" not in input:
    input["halloffame_filename"] = default_filename(specs["problem_specs_prefix"], "halloffame.pkl", "_halloffame.pkl")
  if "n_halloffame" not in input:
    input["n_halloffame"] = 25

  return input


def open_input(input_filename):
  """
  Read casm-learn input file into a dict

  Arguments
  ---------

    input_filename: str
      The path to the input file

  Returns
  -------
    input: dict
      The result of reading the input file and running it through
      casm.learn.set_input_defaults
  """
  # open input and always set input defaults before doing anything else
  with open(input_filename, 'rb') as f:
    try:
      input = set_input_defaults(json.loads(f.read().decode('utf-8')), input_filename)
    except Exception as e:
      print("Error parsing JSON in", args.settings[0])
      raise e
  return input

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
      scorer(estimator, X, y). The parameter for sklearn.model_selection.cross_val_score,
      default = None, uses estimator.score().

    penalty: float, optional, default=0.0
      The CV score is increased by 'penalty*(number of selected basis function)'

    data: pandas.DataFrame, optional, default=None
        Optionally, store TrainingData.data with weighted_X and weighted_y data
        added. No checks are made for consistency of tdata.X, tdata.y and X and
        y or other parameters.

  """

  def __init__(self, X, y, cv, sample_weight=[], scoring=None, penalty=0.0, tdata=None):
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
        scorer(estimator, X, y). The parameter for sklearn.model_selection.cross_val_score,
        default = None, uses estimator.score().

      penalty: float, optional, default=0.0
        The CV score is increased by 'penalty*(number of selected basis function)'

      tdata: TrainingData instance, optional, default=None
        Optionally, store TrainingData.data with weighted_X and weighted_y data
        added. No checks are made for consistency of tdata.X, tdata.y and X and
        y or other parameters.
    """
    self.X = X
    self.y = y

    # Number of configurations and basis functions
    self.n_samples, self.n_features = self.X.shape

    # weight
    self.sample_weight = sample_weight
    self.weighted_y, self.weighted_X, self.W, self.L = casm.learn.tools.set_sample_weight(
      sample_weight, X=self.X, y=self.y)

    # cv sets
    self.cv = cv

    # scoring
    self.scoring = scoring

    # penalty
    self.penalty = penalty

    # data
    if tdata is not None:
      self.data = tdata.data.copy()
      for i in range(self.n_features):
        self.data.loc[:,"weighted_" + tdata.X_name + "(" + str(i) + ")"] = self.weighted_X[:,i]
      for i in range(self.n_features):
        self.data.loc[:,"weighted_" + tdata.y_name] = self.weighted_y


class TrainingData(object):
  """
  TrainingData is a data structure used to collect data from the training
  data file.


  Attributes
  ----------

    filename: str
      The name of the training data file

    filetype: str, one of ["selection", "csv", "json"]
      The data format of the training data file

    X_name: str
      The name of the X columns in the training data file

    X: array-like of shape (n_samples, n_features)
      The training input samples (correlations).

    y_name: str
      The name of the y column in the training data file

    y: array-like of shape: (n_samples, 1)
      The target values (property values).

    n_samples: int
      The number of samples / target values (number of rows in X)

    n_features: int
      The number of features (number of columns in X)

    sel: casm.Selection, (exists if filetype=="selection")
      The selection specifying the training data

    data: pandas.DataFrame
      Contains the X and y data

    hull_dist_name: str, (exists if weight method=="wHullDist")
      The name of the hull_dist column in the training data file

    hull_dist: array-like of shape: (n_samples, 1), (exists if weight method=="wHullDist")
      The hull distance values.


  """
  def __init__(self, input, verbose=True):
    """
    Arguments
    ---------

      input: dict
        The input settings as a dict

      verbose: boolean, optional, default=True
        Print information to stdout.

    """
    specs = input["problem_specs"]

    self.filename = specs["data"]["filename"]
    self.filetype = specs["data"]["filetype"].lower()
    self.X_name = specs["data"]["X"]
    self.y_name = specs["data"]["y"]
    hull_dist_name = "hull_dist"

    if self.filetype == "selection":

      # read training set
      proj = Project(specs["data"]["kwargs"]["project_path"], verbose=verbose)

      sel = Selection(proj, self.filename, all=False)

      # get property name (required)
      property = specs["data"]["y"]

      ## if necessary, query data
      columns = [x for x in [self.y_name, "is_calculated"] if x not in sel.data.columns]
      if len([x for x in sel.data.columns if re.match(self.X_name + "\([0-9]*\)", x)]) == 0:
        columns.append(self.X_name)
      if specs["weight"]["method"] == "wHullDist":
        hull_selection = specs["weight"]["kwargs"]["hull_selection"]
        hull_dist_name = "hull_dist(" + hull_selection + ",atom_frac)"
        if verbose:
          print("# wHullDist: Will calculate hull distance:", hull_dist_name)
        columns.append(hull_dist_name)

      # perform query
      if len(columns):
        sel.query(columns, verbose=verbose)

      data = sel.data
      self.sel = sel

    elif self.filetype.lower() == "csv":
      # populate from csv file
      data = pandas.read_csv(self.filename, **specs["data"]["kwargs"])

    elif self.filetype.lower() == "json":
      # populate from json file
      data = pandas.read_json(self.filename, **specs["data"]["kwargs"])


    # columns of interest, as numpy arrays
    X = data.loc[:,[x for x in data.columns if re.match(self.X_name + "\([0-9]*\)", x)]].values
    y = data.loc[:,self.y_name].values
    if specs["weight"]["method"] == "wHullDist":
      self.hull_dist_name = hull_dist_name
      self.hull_dist = data.loc[:,hull_dist_name]

    self.X = X
    self.y = y
    self.data = data

    self.n_samples = X.shape[0]
    self.n_features = X.shape[1]


def read_sample_weight(input, tdata, verbose=True):
  """
  Read input file and read or calculate sample weights

  Arguments
  ---------

    input: dict
      The input settings as a dict

    tdata: TrainingData instance
      tdata.data is populated with a "weight" column, or the "weight" column
      is used to generate "sample_weight" output, depending on the weight method

    verbose: boolean, optional, default=True
      Print information to stdout.


  Returns
  ---------

    sample_weight: None, 1d array-like of shape (n_samples,1), or 2d array-like of shape (n_samples, n_samples)
      Sample weights.

      if sample_weight is None: (unweighted)
        W = np.matlib.eye(N)
      if sample_weight is 1-dimensional:
        W = np.diag(sample_weight)*Nvalue/np.sum(sample_weight)
      if sample_weight is 2-dimensional (must be Hermitian, positive-definite):
        W = sample_weight*Nvalue/np.sum(sample_weight)
  """

  specs = input["problem_specs"]

  if verbose:
    print("# Weighting:")
    print(json.dumps(specs["weight"], indent=2))

  # get kwargs
  weight_kwargs = copy.deepcopy(specs["weight"]["kwargs"])

  # use method to get weights
  if specs["weight"]["method"] == "wCustom":
    if verbose:
      print("# Reading custom weights")
    sample_weight = tdata.data["weight"].values
  elif specs["weight"]["method"] == "wCustom2d":
    if verbose:
      print("# Reading custom2d weights")
    cols = ["weight(" + str(i) + ")" for i in range(tdata.n_samples)]
    sample_weight = tdata.data.loc[:,cols].values
  elif specs["weight"]["method"] == "wHullDist":
    sample_weight = casm.learn.tools.wHullDist(tdata.hull_dist, **weight_kwargs)
    tdata.data.loc[:,"weight"] = sample_weight
  elif specs["weight"]["method"] == "wEmin":
    sample_weight = casm.learn.tools.wEmin(tdata.y, **weight_kwargs)
    tdata.data.loc[:,"weight"] = sample_weight
  elif specs["weight"]["method"] == "wEref":
    sample_weight = casm.learn.tools.wEref(tdata.y, **weight_kwargs)
    tdata.data.loc[:,"weight"] = sample_weight
  else:
    sample_weight = None

  if verbose:
    print("")

  return sample_weight


def make_fitting_data(input, save=True, verbose=True, read_existing=True):
  """
  Construct a FittingData instance, either by reading existing 'fit_data.pkl',
  or from an input settings.

  Arguments
  ---------

    input: dict
      The input settings as a dict

    save: boolean, optional, default=True
      Save a pickle file containing the training data and scoring metric. The file
      name is specified by input["problem_specs"]["specs_filename"].  See/use
      set_input_defaults for default values.

    verbose: boolean, optional, default=True
      Print information to stdout.

    read_existing: boolean, optional, default=True
      If it exists, read the pickle file containing the training data and scoring
      metric. The file name, which can be specified by input["fit_data_filename"],
      defaults to "fit_data.pkl".


  Returns
  -------

    fdata: casm.learn.FittingData
      A FittingData instance constructed based on the input parameters.

  """
  if verbose:
    print("# Get problem data...")

  specs = input["problem_specs"]

  # property, weight, and cv inputs should remain constant
  # estimator and feature_selection might change
  fit_data_filename = specs["specs_filename"]

  if read_existing and os.path.exists(fit_data_filename):
    if verbose:
      print("# Reading existing problem specs from:", fit_data_filename)
    fdata = pickle.load(open(fit_data_filename, 'rb'))
    if verbose:
      print("#   DONE\n")

    s = "Problem specifications have changed.\n\n" + \
        "To proceed with the existing specs adjust your input settings \"problem_specs\" to match.\n" + \
        "To proceed with the new specs run in a new directory or delete '" + fit_data_filename + "'."

    def check_input(name):
      if fdata.input["problem_specs"][name] != specs[name]:
        print("ERROR: Input file and stored data differ. Input '" + name + "' has changed.")
        print("Stored data:\n", json.dumps(fdata.input["problem_specs"][name], indent=2))
        print("Input:\n", json.dumps(specs[name], indent=2))
        print(s)
        exit()

    for name in ["data", "cv", "weight"]:
      check_input(name)

  else:

    ## get data ####

    tdata = TrainingData(input, verbose=verbose)

    if verbose:
      print("# Target:", tdata.y_name)
      print("# Training samples:", tdata.n_samples)
      print("# Features:", tdata.n_features, "\n")

    ## weight (optional)
    sample_weight = read_sample_weight(input, tdata, verbose=verbose)

    ## cv
    cv_kwargs = copy.deepcopy(specs["cv"]["kwargs"])

    # get cv method (required user input)
    cv_method = _find_method([sklearn.model_selection, casm.learn.model_selection], specs["cv"]["method"])
    cv = cv_method(**cv_kwargs)

    if verbose:
      print("# CV:")
      print(json.dumps(specs["cv"], indent=2), "\n")

    ## scoring
    scoring = sklearn.metrics.make_scorer(sklearn.metrics.mean_squared_error, greater_is_better=True)

    ## penalty
    penalty = specs["cv"]["penalty"]

    fdata = casm.learn.FittingData(tdata.X, tdata.y, cv,
      sample_weight=sample_weight, scoring=scoring, penalty=penalty, tdata=tdata)

    fdata.input = dict()
    fdata.input["problem_specs"] = specs

    if save == True:
      with open(fit_data_filename, 'wb') as f:
          pickle.dump(fdata, f, protocol=2)

    if verbose:
      print("# Writing problem specs to:", fit_data_filename)
      print("# To inspect or customize the problem specs further, use the '--checkspecs' method\n")


  # during runtime only, if LinearRegression and LeaveOneOut, update fdata.cv and fdata.scoring
  # to use optimized LOOCV score method
  if input["estimator"]["method"] == "LinearRegression" and specs["cv"]["method"] == "LeaveOneOut":
    fdata.scoring = None
    fdata.cv = casm.learn.model_selection.LeaveOneOutForLLS(fdata.weighted_y.shape[0])

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

  if verbose:
    print("# Estimator:")
    print(json.dumps(input["estimator"], indent=2), "\n")

  ## estimator

  # get kwargs (default: fit_intercept=False)
  kwargs = copy.deepcopy(input["estimator"]["kwargs"])

  # Use casm version of LinearRegression
  if input["estimator"]["method"] == "LinearRegression":
    estimator_method = casm.learn.linear_model.LinearRegressionForLOOCV
  else:
    estimator_method = _find_method([sklearn.linear_model], input["estimator"]["method"])
  estimator = estimator_method(**kwargs)



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
      scorer(estimator, X, y). The parameter for sklearn.model_selection.cross_val_score,
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

  kwargs = copy.deepcopy(input["feature_selection"]["kwargs"])

  if verbose:
    print("# Feature Selection:")
    print(json.dumps(input["feature_selection"], indent=2), "\n")

  mods = [casm.learn.feature_selection, sklearn.feature_selection, casm.learn.selection_wrapper]

  selector_method = _find_method(mods, input["feature_selection"]["method"])

  # check if 'cv', 'scoring', 'penalty' are allowed kwargs

#  arg_count = selector_method.__init__.func_code.co_argcount
#  allowed_kwargs = selector_method.__init__.func_code.co_varnames[:arg_count]

#  print("allowed_kwargs:", str(allowed_kwargs))

#  if "cv" in allowed_kwargs:
#    kwargs["cv"] = cv
#  if "scoring" in allowed_kwargs:
#    kwargs["scoring"] = scoring
#  if "penalty" in allowed_kwargs:
#    kwargs["penalty"] = penalty
#  if "verbose" in allowed_kwargs:
#    kwargs["verbose"] = verbose

  sig = signature(selector_method.__init__)
#  print("sig.parameters.keys():", str(sig.parameters.keys()))

  if "cv" in sig.parameters:
    kwargs["cv"] = cv
  if "scoring" in sig.parameters:
    kwargs["scoring"] = scoring
  if "penalty" in sig.parameters:
    kwargs["penalty"] = penalty
  if "verbose" in sig.parameters:
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

    mean_absolute_error: float
      The mean absolute prediction error of the unweighted problem

    wmean_absolute_error: float
      The mean absolute prediction error of the weighted problem

    max_absolute_error: float
      The maximum absolute prediction error of the unweighted problem

    wmax_absolute_error: float
      The maximum absolute prediction error of the weighted problem

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

  y_pred = estimator.predict(fdata.X[:,casm.learn.tools.indices(indiv)])
  weighted_y_pred = estimator.predict(fdata.weighted_X[:,casm.learn.tools.indices(indiv)])

  # rms and wrms
  indiv.rms = sqrt(sklearn.metrics.mean_squared_error(fdata.y, y_pred))
  indiv.wrms = sqrt(sklearn.metrics.mean_squared_error(fdata.weighted_y, weighted_y_pred))

  # mean_absolute_error
  indiv.mean_absolute_error = sklearn.metrics.mean_absolute_error(fdata.y, y_pred)
  indiv.wmean_absolute_error = sklearn.metrics.mean_absolute_error(fdata.weighted_y, weighted_y_pred)

  # max_absolute_error
  indiv.max_absolute_error = np.absolute(fdata.y - y_pred).max()
  indiv.wmax_absolute_error = np.absolute(fdata.weighted_y - weighted_y_pred).max()

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


def checkhull(input, hall, indices=None, verbose=True):
  """
  Check for hull properties and add pandas.DataFrame attributes to each individual:

    "dft_gs" : DFT calculated ground states
    "clex_gs" : predicted ground states
    "gs_missing" : DFT ground states that are not predicted ground states
    "gs_spurious" : Predicted ground states that are not DFT ground states
    "uncalculated" : Predicted ground states and near ground states that have not been calculated
    "below_hull" : All configurations predicted below the prediction of the DFT hull
    "ranged_rms": root-mean-square error calculated for the subset of configurations
      whose DFT formation energy lies within some range of the DFT convex hull.
      Currently calculated for ranges 0.001, 0.005, 0.01, 0.05, 0.1, 0.5 eV/unit cell

  Arguments
  ---------

    input: dict
      The input settings as a dict. This is expected to have the object
      input["checkspecs"], with the following attributes:

        selection: str, optional, default="ALL"
          A CASM selection (either 'casm select' output filename or one of the
          standard selection "MASTER", "CALCULATED", or "ALL") containing all the
          configurations to be considered. The DFT convex hull is generated from
          the subset of this selection for which 'is_calculated' is true.

        write_results: bool, optional, default=False
          If True, write casm selection files containing the output data. Output
          selection files are named "checkhull_(problem_specs_prefix)_(i)_(selname)",
          where 'problem_specs_prefix' is input["problem_specs_prefix"], 'i' is the
          index of the individual in the hall of fame, and 'selname' is one of:
            "dft_gs" : DFT calculated ground states
            "clex_gs" : predicted ground states
            "gs_missing" : DFT ground states that are not predicted ground states
            "gs_spurious" : Predicted ground states that are not DFT ground states
            "uncalculated" : Predicted ground states and near ground states that have not been calculated
            "below_hull" : All configurations predicted below the prediction of the DFT hull

        primitive_only: bool, optional, default=True
          If True, only use primitive configurations to construct the convex hull,
          else use all selected configurations.

        uncalculated_range: number, optional, default=0.0
          Include all configurations with clex_hull_dist less than this value (+hull_tol)
          in the "uncalculated" configurations. Default only includes predicted
          ground states.

        ranged_rms: List[number], optional, default=[0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
          Calculates the root-mean-square error for DFT calculated configurations
          within a particular range (in eV/unitcell) of the DFT hull. The list
          provides all the ranges for which the RMSE is requested.

        composition: str, optional, default="atom_frac"
          Composition argument use for 'casm query' properties 'hull_dist' and
          'clex_hull_dist'. For thermodynamic ground states, use "atom_frac".

        hull_tol: number, optional, default=proj.settings.data["lin_alg_tol"]
          Tolerance used for identify hull states

        dim_tol: number, optional, default=1e-8
          Tolerance for detecting composition dimensionality

        bottom_tol: number, optional, default=1e-8
          Tolerance for detecting which facets form the convex hull bottom


        Example:

          "checkhull" : {
            "selection": "ALL",
            "write_results": True
            "primitive_only": true,
            "uncalculated_range": 0.0,
            "ranged_rms": [0.001, 0.005, 0.01, 0.05, 0.1, 0.5],
            "composition": "atom_frac",
            "hull_tol": 1e-8,
            "dim_tol": 1e-8,
            "bottom_tol": 1e-8,
          }


    hall: deap.tools.HallOfFame, optional, default=None
      A Hall Of Fame containing individuals to check convex hull properties

    selection: str
      The name of the selection to use for calculating the convex hull

    write_results: bool, optional, default=True


    indices: List[int], optional, default=None
      If given, only check hull properties for specified individuals in the hall
      of fame.

    verbose: boolean, optional, default=True
      Print information to stdout.


  """

  # set checkhull default settings
  if "checkhull" not in input:
    input["checkhull"] = dict()

  opt = {
    "selection":"ALL",
    "write_results":False,
    "primitive_only":True,
    "uncalculated_range":0.0,
    "hull_tol":1e-8,
    "composition":"atom_frac",
    "dim_tol":1e-8,
    "bottom_tol":1e-8,
    "ranged_rms":[0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
  }

  for key, val in six.iteritems(opt):
    if key not in input["checkhull"] or input["checkhull"][key] is None:
      input["checkhull"][key] = val

  # construct FittingData (to check consistency with input file)
  fdata = make_fitting_data(input, save=True, verbose=verbose, read_existing=True)

  # open Project to work on
  proj = Project(input["problem_specs"]["data"]["kwargs"]["project_path"], verbose=verbose)

  # tolerance for finding configurations 'on_hull'
  d = input["checkhull"]
  selection = d["selection"]
  composition = d["composition"]
  write_results = d["write_results"]
  hull_tol = d["hull_tol"]
  uncalculated_range = d["uncalculated_range"]
  dim_tol = d["dim_tol"]
  bottom_tol = d["bottom_tol"]
  primitive_only = d["primitive_only"]

  # save current default clex
  orig_clex = proj.settings.default_clex
  clex = copy.deepcopy(orig_clex)
  all_eci = proj.dir.all_eci(clex.property, clex.calctype, clex.ref, clex.bset)

  # not sure of all the edge cases, enforcing these seems simpler for now...
  if clex.name != "formation_energy":
    print("default clex:", clex.name)
    print("use 'casm settings --set-default-clex formation_energy' to change the default clex")
    raise Exception("Error using checkhull: default clex must be 'formation_energy'")

  if input["problem_specs"]["data"]["y"] != "formation_energy":
    print("property:", input["problem_specs"]["data"]["y"])
    raise Exception("Error using checkhull: property must be 'formation_energy'")

  # load hull selection
  sel = Selection(proj, selection, all=False)

  # create a temporary file containing the subset of the selection that is calculated
  tmp_dir = join(proj.dir.casm_dir(), "tmp")
  if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)
  dft_selection = join(tmp_dir, selection + "_dft")
  dft_hull_selection = join(tmp_dir, selection + "_dft_hull")
  proj.capture("select -c " + selection + " --set-off 'not(is_calculated)' -f -o " + dft_selection)

  # properties to query
  selected = "selected"
  is_primitive = "is_primitive"
  is_calculated = "is_calculated"
  configname = "configname"
  dft_hull_dist_long = "hull_dist(" + dft_selection + "," + composition + ")"
  dft_hull_dist = "dft_hull_dist"
  clex_hull_dist_long = "clex_hull_dist(" + selection + "," + composition + ")"
  clex_hull_dist = "clex_hull_dist"
  clex_dft_hull_dist_long = "clex_hull_dist(" + dft_hull_selection + "," + composition + ")"
  clex_dft_hull_dist = "clex_dft_hull_dist"
  comp = "comp"
  dft_Eform = "formation_energy"
  clex_Eform = "clex(formation_energy)"

  query_cols = [comp, is_calculated, configname, dft_hull_dist_long, dft_Eform]
  if primitive_only:
    query_cols.append(is_primitive)
  sel.query(query_cols)


  compcol = []
  for col in sel.data.columns:
    if len(col) >= 5 and col[:5] == "comp(":
      compcol.append(col)
  compcol.sort()

  if indices is None:
    indices = range(len(hall))

  # write eci.json
  eci = "__tmp"
  clex.eci = eci
  if eci not in all_eci:
    proj.capture("settings --new-eci " + eci)
  else:
    proj.capture("settings --set-eci " + eci)

  # for each individual specified...
  for indiv_i in indices:

    if verbose:
      print("-- Check: individual", indiv_i, " --")
      print_individual(hall, [indiv_i])
      print("")

    # write ECI to use
    indiv = hall[indiv_i]
    write_eci(proj, indiv.eci, fit_details=casm.learn.to_json(indiv_i, indiv), clex=clex, verbose=verbose)

    if primitive_only:
      df = sel.data[sel.data.loc[:,is_primitive] == 1].sort_values(compcol)
    else:
      df = sel.data.sort_values(compcol)
    df_calc = df[df.loc[:,is_calculated] == 1].apply(pandas.to_numeric, errors='ignore')
    dft_gs = df_calc[df_calc.loc[:,dft_hull_dist_long] < hull_tol]

    # create dft hull selection
    dft_hull_sel = Selection(proj, dft_hull_selection)
    dft_hull_sel.save(data=dft_gs, force=True)

    # query:
    sel.query([clex_hull_dist_long, clex_Eform, clex_dft_hull_dist_long], force=True)

    if primitive_only:
      df = sel.data[sel.data.loc[:,is_primitive] == 1].sort_values(compcol)
    else:
      df = sel.data.sort_values(compcol)
    df.rename(
      inplace=True,
      columns={
        dft_hull_dist_long : dft_hull_dist,
        clex_hull_dist_long : clex_hull_dist,
        clex_dft_hull_dist_long : clex_dft_hull_dist
      })
    df_calc = df[df.loc[:,is_calculated] == 1].apply(pandas.to_numeric, errors='ignore')

    clex_gs = df[df.loc[:,clex_hull_dist] < hull_tol]
    dft_gs = df_calc[df_calc.loc[:,dft_hull_dist] < hull_tol]

    uncalculated = df[(df.loc[:,clex_hull_dist] < uncalculated_range + hull_tol) & (df.loc[:,is_calculated] == 0)]
    gs_spurious = df_calc[(df_calc.loc[:,clex_hull_dist] < hull_tol) & ~(df_calc.loc[:,dft_hull_dist] < hull_tol)]
    gs_missing = df_calc[~(df_calc.loc[:,clex_hull_dist] < hull_tol) & (df_calc.loc[:,dft_hull_dist] < hull_tol)]
    below_hull = df[df.loc[:,clex_dft_hull_dist] < -hull_tol]

    indiv.clex_gs = clex_gs
    indiv.dft_gs = dft_gs
    indiv.uncalculated = uncalculated
    indiv.gs_spurious = gs_spurious
    indiv.gs_missing = gs_missing
    indiv.below_hull = below_hull

    to_drop = [selected,is_primitive,is_calculated]

    def printer(attr, title):
      prefix = input["problem_specs"]["problem_specs_prefix"]
      output_name = "checkhull"
      if len(prefix) != 0:
         output_name += "_" + prefix
      output_name += "_" + str(indiv_i) + "_" + attr

      df = getattr(indiv, attr)
      kwargs = {"index":False}
      if df.shape[0]:
        if verbose:
          print(title + ":")
          print(df.drop(to_drop, axis=1, errors='ignore').to_string(**kwargs))
          if write_results:
            print("write:", output_name, "\n")
          else:
            print("")
        if write_results:
          output_sel = Selection(proj, output_name)
          output_sel.save(data=df, force=True)

    printer("dft_gs", "DFT ground states")
    printer("clex_gs", "Predicted ground states")
    printer("gs_missing", "DFT ground states that are not predicted ground states")
    printer("gs_spurious", "Predicted ground states that are not DFT ground states")
    printer("uncalculated", "Predicted ground states and near ground states that have not been calculated")
    printer("below_hull", "All configurations predicted below the prediction of the DFT hull")

    # calculated ranged rms
    def calc_ranged_rms(r):
      ranged_df = df_calc[df_calc.loc[:,dft_hull_dist] < r + hull_tol].copy()
      n = ranged_df.shape[0]
      ranged_df.loc[:,"diff"] = ranged_df.loc[:,dft_Eform] - ranged_df.loc[:,clex_Eform]
      ranged_df.loc[:,"diff_sqr"] = ranged_df.loc[:,"diff"]*ranged_df.loc[:,"diff"]
      rms = sqrt(sklearn.metrics.mean_squared_error(ranged_df.loc[:,dft_Eform], ranged_df.loc[:,clex_Eform]))
      return (n, rms)

    ranged_rms = []
    for r in input["checkhull"]["ranged_rms"]:
      res = calc_ranged_rms(r)
      ranged_rms.append({"range": r, "n_config": res[0], "rms": res[1] })
    ranged_rms.append({"range": "all", "n_config":df_calc.shape[0], "rms":indiv.rms})
    if verbose:
      print("ranged_rms:")
      for d in ranged_rms:
        print("  RMS error for", d["n_config"], "calculated configurations within", end=' ')
        print(d["range"], "eV/unitcell of the DFT hull:", d["rms"])
      print("")

    sel.data.drop([clex_dft_hull_dist_long], axis=1, inplace=True)

    indiv.ranged_rms = ranged_rms
    indiv.checkhull_settings = input["checkhull"]

    if verbose:
      print("\n")

  # reset eci setting
  proj.capture("settings --set-eci " + orig_clex.eci)

  # remove __tmp eci and tmp selections
  shutil.rmtree(proj.dir.eci_dir(clex))
  os.remove(dft_selection)
  os.remove(dft_hull_selection)


def checkspecs(input, verbose=True):
  """
  Output data and cv files containing the current problem specs.

  Arguments
  ---------

    input: dict
      The input settings as a dict. This is expected to have the object
      input["checkspecs"], with the following attributes:

        data_filename: string
          The path to the file where the training data (including weights) should be
          written

        data_kwargs: dict or null, optional, default=dict()
          Additional parameters to be used to write training data.

          Options for input 'filetype' "selection":
            None.

          Options for input 'filetype' "csv":
            Any options to pass to pandas.to_csv

          Options for input 'filetype' "json":
            Any options to pass to pandas.to_json

        cv_filename: string
          The path to the file where the cv generator or train/tests sets should be
          written as pickle file

        cv_kwargs: dict or null, optional, default=dict()
          Additional parameters to be used to write cv data using pickle.dump.

        Example:

          "checkspecs" : {
            "data_filename": "check_train",
            "data_kwargs": null,
            "cv_filename": "check_cv.pkl",
            "cv_kwargs": null
          }

    verbose: boolean, optional, default=True
      Print information to stdout.

  """
  specs = input["problem_specs"]

  # set checkspecs default settings
  if "checkspecs" not in input:
    input["checkspecs"] = dict()
  if "data_filename" not in input["checkspecs"]:
    suffix = "_data" + splitext(basename(specs["data"]["filename"]))[1]
    input["checkspecs"]["data_filename"] = default_filename(specs["problem_specs_prefix"], "check_data", suffix)
  if "data_kwargs" not in input["checkspecs"] or input["checkspecs"]["data_kwargs"] is None:
    input["checkspecs"]["data_kwargs"] = dict()
  if "cv_filename" not in input["checkspecs"]:
    input["checkspecs"]["cv_filename"] = default_filename(specs["problem_specs_prefix"], "check_cv.pkl", "_cv.pkl")
  if "cv_kwargs" not in input["checkspecs"] or input["checkspecs"]["cv_kwargs"] is None:
    input["checkspecs"]["cv_kwargs"] = dict()

  specs = input["problem_specs"]

  # construct or load problem specs
  fdata = casm.learn.fit.make_fitting_data(input, save=True, verbose=verbose, read_existing=True)

  if verbose:
    print("# Check problem specs:")
    print(json.dumps(input["checkspecs"], indent=2), "\n")

  # print training data
  filetype = specs["data"]["filetype"]
  filename = input["checkspecs"]["data_filename"]
  kwargs = input["checkspecs"]["data_kwargs"]

  if filetype == "selection":
    proj = Project(specs["data"]["kwargs"]["project_path"], verbose=verbose)
    sel = Selection(proj, filename, all=False)
    sel.save(data=fdata.data, force=True)
  elif filetype == "csv":
    fdata.data.to_csv(filename, **kwargs)
  elif filetype == "json":
    fdata.data.to_json(filename, **kwargs)

  if verbose:
    print("# Wrote training data (including weights) to:", filename)

  # print cv
  cv_filename = input["checkspecs"]["cv_filename"]
  cv_kwargs = input["checkspecs"]["cv_kwargs"]
  with open(cv_filename, 'wb') as f:
    pickle.dump(fdata.cv, f, protocol=2)

  if verbose:
    print("# Wrote CV data to:", cv_filename, "\n")


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
  print("{0:5}: {1:43} {2:<12} {3:<12}".format("Index", "Selected", "#Selected", "CV"))
  print("-"*100)
  form_str = "{0:5}: {1} {2:<12} {3:<12.8g}"
  for i in range(len(pop)):
    print(form_str.format(i, bitstr(pop[i], 40), sum(pop[i]), pop[i].fitness.values[0]))


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
    "wrms": weighted problem root-mean-square error
    "mean_absolute_error": mean absolute error
    "wmean_absolute_error": weighted problem mean absolute error
    "max_absolute_error": maximum absolute error
    "wmax_absolute_error": weighted problem maximum absolute error
    "estimator_method": name of estimator method
    "features_selection_method": name of feature selection method
    "note": a descriptive note
    "eci": List of (feature index, coefficient value) pairs
    "input": input settings dict

  Keys added by --checkhull, if y="formation_energy": pandas.DataFrame
    "dft_gs" : DFT calculated ground states
    "clex_gs" : predicted ground states
    "gs_missing" : DFT ground states that are not predicted ground states
    "gs_spurious" : Predicted ground states that are not DFT ground states
    "uncalculated" : Predicted ground states and near ground states that have not been calculated
    "below_hull" : All configurations predicted below the prediction of the DFT hull
    "ranged_rms": root-mean-square error calculated for the subset of configurations
      whose DFT formation energy lies within some range of the DFT convex hull


  Arguments
  ---------

    index: int
      The index in hall of fame of the individual

    indiv: List[bool] of length n_features
      This is a boolean list of shape [n_features], in which an element is True
      iff its corresponding feature is selected for retention.


  Note
  ----

    ECI are serialized using cls=casm.misc.noident.noindent.NoIndent, so when writing with json.dump or
    json.dumps, include 'cls=casm.misc.noindent.noindent.NoIndentEncoder'.

  """
  d = dict()
  d["index"] = index
  d["id"] = str(indiv.id)
  d["selected"] = bitstr(indiv)
  d["n_selected"] = int(sum(indiv))
  d["cv"] = indiv.fitness.values[0]

  attr = [
    "rms", "wrms", "ranged_rms",
    "mean_absolute_error", "wmean_absolute_error",
    "max_absolute_error", "wmax_absolute_error",
    "estimator_method", "feature_selection_method",
    "note",
    "checkhull_settings"
  ]
  for a in attr:
    if hasattr(indiv, a):
      d[a] = getattr(indiv,a)

  for attr in ["clex_gs", "dft_gs", "uncalculated", "gs_spurious", "gs_missing", "below_hull"]:
    if hasattr(indiv, attr):
      d[attr] = json.loads(getattr(indiv, attr).to_json(orient='records'))
  d["eci"] = []
  for bfunc in indiv.eci:
    d["eci"].append(noindent.NoIndent(bfunc))

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
    "mean_absolute_error": mean absolute error
    "wmean_absolute_error": weighted problem mean absolute error
    "max_absolute_error": maximum absolute error
    "wmax_absolute_error": weighted problem maximum absolute error
    "estimator_method": name of estimator method
    "features_selection_method": name of feature selection method
    "note": a descriptive note
    "eci": JSON string of a List of (feature index, coefficient value) pairs
    "input": input settings dict

  Keys added by --checkhull, if y="formation_energy": pandas.DataFrame
    "dft_gs" : DFT calculated ground states
    "clex_gs" : predicted ground states
    "gs_missing" : DFT ground states that are not predicted ground states
    "gs_spurious" : Predicted ground states that are not DFT ground states
    "uncalculated" : Predicted ground states and near ground states that have not been calculated
    "below_hull" : All configurations predicted below the prediction of the DFT hull
    "ranged_rms": root-mean-square error calculated for the subset of configurations
      whose DFT formation energy lies within some range of the DFT convex hull


  Arguments
  ---------

    indices: List[int]
      The indices in hall of fame of the individuals to include

    hall: deap.tools.HallOfFame
      A Hall Of Fame of ECI sets

  """
  data = [to_json(i, hall[i]) for i in indices]
  for d in data:
    d["eci"] = json.dumps(d["eci"], cls=noindent.NoIndentEncoder)
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
    print(form_str.format(index, bitstr(indiv,40), sum(indiv), indiv.fitness.values[0],
      indiv.rms, indiv.wrms, indiv.estimator_method, indiv.feature_selection_method, indiv.note))
    return

  elif format.lower() == "json":
    print(json.dumps(to_json(index,indiv), indent=2, cls=noindent.NoIndentEncoder))
    return

  elif format.lower() == "details":
    print("##")
    print("Index:", index)
    print("ID:", indiv.id)
    print("Selected:", bitstr(indiv))
    print("#Selected:", sum(indiv))
    print("CV:", indiv.fitness.values[0])

    attr = [
      "rms", "wrms",
      "mean_absolute_error", "wmean_absolute_error",
      "max_absolute_error", "wmax_absolute_error",
      "estimator_method", "feature_selection_method",
      "note",
    ]
    for a in attr:
      if hasattr(indiv, a):
        print(a + ":", getattr(indiv, a))
    if hasattr(indiv, "ranged_rms"):
      print("ranged_rms:\n", json.dumps(indiv.ranged_rms, indent=2))


    print("eci:\n")
    print_eci(indiv.eci)

    for attr in ["clex_gs", "dft_gs", "uncalculated", "gs_spurious", "gs_missing", "below_hull"]:
      if hasattr(indiv, attr):
        print(attr + ":")
        print(getattr(indiv, attr).to_string(index=False))

    print("Input:\n", json.dumps(indiv.input, indent=2))
    if hasattr(indiv, "checkhull_settings"):
      print("Checkhull settings:\n", json.dumps(indiv.checkhull_settings, indent=2))
    return


def _print_halloffame_header(hall):
  """
  Print header for hall of fame.
  """
  if len(hall[0]) > 40:
    bitstr_len = 43
  else:
    bitstr_len = len(hall[0])
  print(("{0:5}: {1:<" + str(max([12, bitstr_len])) + "} {2:<12} {3:<12} {4:<12} {5:<12} {6:<24} {7:<24} {8}").format(
    "Index", "Selected", "#Selected", "CV", "RMS", "wRMS", "Estimator", "FeatureSelection", "Note"))
  print("-"*(6+max([12, bitstr_len])+13*5+25*2))


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
    print(json.dumps(h, indent=2, cls=noindent.NoIndentEncoder))

  elif format.lower() == "csv":
    df = to_dataframe(range(len(hall)), hall)
    print(df.to_csv())

  elif format.lower() == "details":
    for index in indices:
      _print_individual(index, hall[index], format=format)
    print("")
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
    print(json.dumps(h, indent=2, cls=noindent.NoIndentEncoder))


  elif format.lower() == "csv":
    df = to_dataframe(range(len(hall)), hall)
    print(df.to_csv())

  elif format.lower() == "details":
    for index, indiv in enumerate(hall):
      _print_individual(index, indiv, format=format)
      print("")

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
    print("{index:>5}: {value:< .12g}".format(index=bfunc[0], value=bfunc[1]))


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
    print("Loading Hall of Fame:", halloffame_filename)

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
    print("\nPickling Hall of Fame to:", halloffame_filename)
  with open(halloffame_filename, 'wb') as f:
      pickle.dump(hall, f, protocol=2)
