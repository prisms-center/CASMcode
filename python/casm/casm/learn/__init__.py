"""CASM machine learning tools"""

## This part needs to be in global scope for parallization #####################  
from deap import creator
from deap import base

# we'll want to minimize a cv score
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))

# each individual is a list of True or False indicating if each basis function should 
# be included in the model
creator.create("Individual", list, fitness=creator.FitnessMin, input=None)
################################################################################  


from fit import example_input_Lasso, example_input_LassoCV, example_input_RFE, \
 example_input_GeneticAlgorithm, example_input_IndividualBestFirst, \
 example_input_PopulationBestFirst, print_input_help, FittingData, \
 fit_and_select, print_individual, print_population, print_halloffame, print_eci

__all__ = dir()