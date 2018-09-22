from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import argparse
import deap.tools
import json
import os
import six
import sys

import casm.learn
from casm.project import Project, Selection, write_eci
  

def main(argv = None):
  if argv is None:
    argv = sys.argv[1:]
  parser = argparse.ArgumentParser(description = 'Fit cluster expansion coefficients (ECI)')
  parser.add_argument('--desc', help='Print extended usage description', action="store_true")
  parser.add_argument('-s', '--settings', nargs=1, help='Settings input filename', type=str)
  parser.add_argument('--format', help='Hall of fame print format. Options are "details", "json", or "csv".', type=str, default=None)
  #parser.add_argument('--path', help='Path to CASM project. Default assumes the current directory is in the CASM project.', type=str, default=os.getcwd())
  parser.add_argument('--settings-format', help='Print input file description', action="store_true")
  parser.add_argument('--exLasso', help='Print example input file using Lasso', action="store_true")
  parser.add_argument('--exLassoCV', help='Print example input file using LassoCV', action="store_true")
  parser.add_argument('--exRFE', help='Print example input file using Recursive Feature Elimination (RFE)', action="store_true")
  parser.add_argument('--exGeneticAlgorithm', help='Print example input file using GeneticAlgorithm', action="store_true")
  parser.add_argument('--exIndividualBestFirst', help='Print example input file using IndividualBestFirst', action="store_true")
  parser.add_argument('--exPopulationBestFirst', help='Print example input file using PopulationBestFirst', action="store_true")
  parser.add_argument('--exDirectSelection', help='Print example input file using DirectSelection', action="store_true")
  parser.add_argument('--hall', help='Print hall of fame summary', action="store_true")
  parser.add_argument('--indiv', nargs='+', help='Specify particular individuals by index in hall of fame', type=int)
  parser.add_argument('--select', nargs=1, help='Select individual to use', type=int)
  parser.add_argument('--checkhull', help='Check convex hull properties using the provided selection', action="store_true", default=False)
  parser.add_argument('--checkspecs', help='Output data and cv files containing the current problem specs', action="store_true", default=False)
  parser.add_argument('-q','--quiet', help='Quiet output', action="store_true", default=False)
  args = parser.parse_args(argv)
  
  args.verbose = not args.quiet
  
  if args.settings_format:
    casm.learn.print_input_help()
    return
  
  if args.exLasso:
    print(six.u(json.dumps(casm.learn.example_input_Lasso(), indent=2)))
    return
  elif args.exLassoCV:
    print(six.u(json.dumps(casm.learn.example_input_LassoCV(), indent=2)))
    return
  elif args.exRFE:
    print(six.u(json.dumps(casm.learn.example_input_RFE(), indent=2)))
    return
  elif args.exGeneticAlgorithm:
    print(six.u(json.dumps(casm.learn.example_input_GeneticAlgorithm(), indent=2)))
    return
  elif args.exIndividualBestFirst:
    print(six.u(json.dumps(casm.learn.example_input_IndividualBestFirst(), indent=2)))
    return
  elif args.exPopulationBestFirst:
    print(six.u(json.dumps(casm.learn.example_input_PopulationBestFirst(), indent=2)))
    return
  elif args.exDirectSelection:
    print(six.u(json.dumps(casm.learn.example_input_DirectSelection(), indent=2, cls=casm.NoIndentEncoder)))
    return
  
  if args.settings:
    
    if args.verbose:
      print("Loading", args.settings[0])
    
    input = casm.learn.open_input(args.settings[0])
    
    if args.hall:
      
      halloffame_filename = input["halloffame_filename"]
      # print Hall of Fame summary
      existing_hall = casm.learn.open_halloffame(halloffame_filename, args.verbose)
      
      if args.indiv:
        casm.learn.print_individual(existing_hall, args.indiv, format=args.format)
      elif args.hall:
        casm.learn.print_halloffame(existing_hall, format=args.format)
      
    elif args.checkhull:
      
      halloffame_filename = input["halloffame_filename"]
      # print Hall of Fame summary
      existing_hall = casm.learn.open_halloffame(halloffame_filename, args.verbose)
      
      casm.learn.fit.checkhull(input, existing_hall, indices=args.indiv, verbose=args.verbose)
      
      # pickle hall of fame
      casm.learn.save_halloffame(existing_hall, halloffame_filename, args.verbose)
      
    elif args.select:
      
      specs = input["problem_specs"]
      halloffame_filename = input["halloffame_filename"]
      # print Hall of Fame summary
      existing_hall = casm.learn.open_halloffame(halloffame_filename, args.verbose)
      
      proj = Project(specs["data"]["kwargs"]["project_path"])
      
      index = args.select[0]
      indiv = existing_hall[index]
      
      write_eci(proj, indiv.eci, casm.learn.to_json(index, indiv), verbose=args.verbose)
      
    elif args.checkspecs:
      
      casm.learn.checkspecs(input, verbose=args.verbose)
           
    else:
      
      # construct hall of fame
      specs = input["problem_specs"]
      halloffame_filename = input["halloffame_filename"]
      halloffame_size = input["n_halloffame"]
      hall = casm.learn.create_halloffame(halloffame_size)
      if args.verbose:
        print("# Hall of Fame size:", halloffame_size, "\n")
      
      if os.path.exists(halloffame_filename):
        existing_hall = casm.learn.open_halloffame(halloffame_filename, args.verbose)
        hall.update(existing_hall)
      
      # run fitting
      if input["feature_selection"]["method"] == "DirectSelection":
        casm.learn.direct_fit(input, verbose=args.verbose, hall=hall)
      else:
        casm.learn.fit_and_select(input, verbose=args.verbose, hall=hall)
      
      # pickle hall of fame
      casm.learn.save_halloffame(hall, halloffame_filename, args.verbose)
    
  elif args.desc:
    
    print("""
    
    1) Specify the problem:
      
      'casm-learn' helps solve the problem:
      
        X*b = y,
      
      where:
        
        X: 2d matrix of shape (n_samples, n_features)
          The correlation matrix, holding the evaluated basis functions. The 
          entry X[config, bfunc] holds the average value of the 'bfunc' cluster
          basis function for configuration 'config'. The number of configurations
          is 'n_samples' and the number of cluster basis functions is 'n_features'.
        
        y: 1d matrix of shape (n_samples, 1)
          The calculated properties being fit to. The most common case is that
          y[config] holds the formation energy calculated for configuration 
          'config'.
        
        b: 1d matrix of shape (n_features, 1)
          The effective cluster interactions (ECI) being solved for.
      
      To specify this problem, the 'casm-learn' input file specifies which
      configurations to fit to (the training data), how to weight the 
      configurations, and how to compare solutions via cross-validation. 
            
      
      Training data may be input via a 'casm select' output file. The default
      name expected is 'train'. So to use all calculated configurations, you 
      could create a directory in your CASM project where you will perform
      fitting and generate a 'train' file:
        
        cd /my/casm/project
        mkdir fit_1 && cd fit_1
        casm select --set is_calculated -o train
      
      
      Example 'casm-learn' JSON input files can be output by the
      'casm-learn --exMethodName' options:
        
        casm-learn --exGeneticAlgorithm > fit_1_ga.json
        casm-learn --exRFE > fit_1_rfe.json
        ...etc..
      
      By default, these settings files are prepared for fitting formation_energy,
      using the 'train' configuration selection. Edit the file as needed, and
      see 'casm-learn --settings-format' for help.
      
      
      When weighting configurations, the problem is transformed:
      
        X*b = y  ->  L*X*b = L*y, 
  
      where, W = L*L.tranpose(): 
        
        W: 2d matrix of shape (n_samples, n_samples)
          The weight matrix is specified in the casm-learn input file. If the 
          weighting method provides 1-dimensional input (this is typical, i.e.
          a weight for each configuration), in an array called 'w', then:
          
            W = diag(w)*n_samples/sum(w),
            
          diag(w) being the diagonal matrix with 'w' along the diagonal.
      
      
      A cross-validation score is used for comparing generated ECI. The cv score
      reported is:
  
        cv = sqrt(mean(scores)) + N_nonzero_eci*penalty, 
  
      where:
      
        scores: 1d array of shape (number of train/test sets)
          The mean squared error calculated for each training/testing set
        
        N_nonzero_eci: int
          The number of basis functions with non-zero ECI
          
        penalty: number, optional, default=0.0
          Is the user-input penalty per basis function that can be used to
          favor solutions with a small number of non-zero ECI
      
      See 'casm-learn --settings-format' for help specifying the cross-validation
      training and test sets using options from scikit-learn. It is usually
      important to use the 'shuffle'=true option so that configurations are
      randomly added to train/test sets and not ordered by supercell size.
      
      
      When you run 'casm-learn' with a new problem specification the first time, 
      it generates a "problem specs" file that stores the training data, weights, 
      and cross-validation train/test sets. Then, when running subsequent times,
      the data can be loaded more quickly, and the cross-validation can be 
      performed using the same train/test sets. 'casm-learn' will attempt to
      prevent you from re-running with a different problem specification so that
      solutions can be compared via their cv score in an "apples-to-apples" 
      manner. The default name for the "specs" file is determined from the input
      filename. For example, 'my_input_specs.pkl' is used if the input file is 
      named 'my_input.json'. See 'casm-learn --settings-format' for more help.
      
      
      The '--checkspecs' option can be used to write output files with the 
      generated problem specs data. Amont other things, this can be used to 
      adjust weights manually or save and re-use train/test sets. See 
      'casm-learn --settings-format' for more help.
      
    
    2) Select estimator and feature selection methods
      
      The "estimator" option specifies a linear model estimator that determines
      how to solve the linear problem L*X*b = L*b, for b.
      
      The "feature_selection" option specifies a feature selection method that
      determines which features (ECI) should be considered for the solution. The
      remaining are effectively set to 0.0 when calculating the cluster 
      expansion. Generally there is a tradeoff: By limiting the number of 
      features included in the cluster expansion Monte Carlo calculations can be
      more efficient, but at a possible loss of accuracy. Be careful to avoid
      overfitting however. If your cross validation scheme does not provide 
      enough testing data, you may fit your training data very well, but not
      have an accurate extrapolation to other configurations.
      
      See 'casm-learn --settings-format' for help specifying the estimator and
      feature selection methods. Assuming you are using the GeneticAlgorithm and
      have named your input file 'fit_1_ga.json', run:
        
        casm-learn -s fit_1_ga.json
      
      'casm-learn' will run and eventually store its results.  For a single 
      problem specification (step 1, the settings in "problem_specs"
      in the 'casm-learn' input file), you may try many different estimation
      and feature selection methods and use the cv score to compare results. All
      the results for a single problem specification can be stored in a 'Hall Of
      Fame' that collects the N individual solutions with the best cv scores. To
      view these results use:
        
        casm-learn -s fit_1_ga.json --hall
      
      For more details, or to output the results for further analysis in JSON or 
      CSV format, there is a '--format' option. To view only particular 
      individuals in the hall of fame, there is a '--indiv' option.
    
    
    3) Analyze results
      
      The above steps (1) and (2) may be repeated many times as you attempt to
      optimize your ECI. Solutions for different problems (i.e. different 
      weighting schemes, re-calculating with more training data) may be compared 
      based on scientific knowledge, for instance, which predicts the 0K ground 
      state configurations correctly, or from analysis of Monte Carlo results.
      
      The '--checkhull' option provides a simple way to check the 0K ground 
      states and can create 'casm select' style output files with enumerated but
      uncalculated configurations that are predicted to be low energy. These can
      then be used to generate more training data and re-fit the ECI.
      
      When you have generated ECI that you wish to use in Monte Carlo 
      calculations, use the '--select' option to write an 'eci.json' file into
      your CASM project for the currently selected cluster expansion (as listed
      by 'casm settings -l).
    
    
    4) Use results
    
      Once an 'eci.json' file has been written, you can run Monte Carlo 
      calculations. See 'casm monte -h' and 'casm format --monte' for help. 
      
    """)
    
  else:
    
    parser.print_help()

if __name__ == "__main__":
  main()

