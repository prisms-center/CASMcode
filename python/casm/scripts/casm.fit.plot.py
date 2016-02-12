#!/usr/bin/env python
from casm.project import Selection, Project, query
import os
import pandas
import numpy as np
from string import ascii_lowercase

import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt

def plot(proj):
  """ 
  Plot calculated and predicted formation energy and hulls
  """
  

class ClexAnalysis(object):

  def __init__(self, proj=None, train=None, predict=None):
    """
    Arguments:
      proj: a casm.project.Project, or None. If None, defaults to Project()
      train: a casm.project.Selection, or None. If None, defaults to Selection(proj, path="train").
      predict a casm.project.Selection, or None. If None, defaults to Selection(proj, path="predict").
    """
    
    if proj==None:
      proj = Project()
    self.proj = proj
    
    if train==None:
      train = Selection(proj, path="train")
    self.train = train
    
    if predict == None:
      predict = Selection(proj, path="predict")
    self.predict = predict
    
  
  def check_hull(self, tol=1e-8):
    
    configname = 'configname'
    comp = 'comp'
    
    calc_prop = 'formation_energy'
    pred_prop = 'clex(formation_energy)'
    
    long_calc_hull = 'hull_dist({0},atom_frac)'.format(self.train.path)
    calc_hull = 'hull_dist(train,atom_frac)'
    
    long_pred_hull = 'clex_hull_dist({0},atom_frac)'.format(self.predict.path)
    pred_hull = 'clex_hull_dist(predict,atom_frac)'
    
    long_pred_calc_hull = 'clex_hull_dist({0},atom_frac)'.format(self.train.path)
    pred_calc_hull = 'clex_hull_dist(train,atom_frac)'
    
    replace = [(long_calc_hull, calc_hull), (long_pred_hull, pred_hull), (long_pred_calc_hull, pred_calc_hull)]
    
    is_calculated = 'is_calculated'
    
    hull_tol = tol
    
    
    df_train = query(self.proj, [
      configname,
      comp, 
      calc_prop, 
      long_calc_hull,
      pred_prop], 
      self.train)
    
    # shorten column names
    for from_name, to_name in replace:
      df_train.columns = [n.replace(from_name, to_name) for n in df_train.columns]
    
    comp_axes = []
    for c in ascii_lowercase:
      name = 'comp(' + c + ')'
      if name in df_train.columns:
        comp_axes.append(name)
    print comp_axes
    
    # sort by hull_dist, followed by comp to make plotting the calculated hull easy
    sort_order = [calc_hull] + comp_axes
    df_train.sort_values(sort_order,inplace=True)
    
    # select the calculated hull 
    df_calc_hull = df_train[ df_train[calc_hull] < hull_tol ]
    
    print "\nCalculated ground states (among the training set):"
    print df_calc_hull.to_string(index=False)
    
    
    df_predict = query(self.proj, [
      configname,
      comp, 
      calc_prop,
      pred_prop, 
      long_calc_hull,
      long_pred_hull,
      long_pred_calc_hull,
      is_calculated], 
      self.predict)
    
    # shorten column names
    for from_name, to_name in replace:
      df_predict.columns = [n.replace(from_name, to_name) for n in df_predict.columns]
    
    # sort by clex_hull_dist, followed by comp to make plotting the predicted hull easy
    sort_order = [pred_hull] + comp_axes
    df_predict.sort_values(sort_order,inplace=True)
    
    
    # configurations on the predicted hull
    df_pred_hull = df_predict[ df_predict[pred_hull] < hull_tol]
    
    print "\nPredicted ground states:"
    print df_pred_hull.to_string(index=False)
    
    
    # configurations on the predicted hull that are not calculated
    df_pred_gs_uncalculated = df_pred_hull[ df_pred_hull[is_calculated] == 0 ]
    
    if len(df_pred_gs_uncalculated.index):
      print "\nUncalculated predicted ground states:"
      print df_pred_gs_uncalculated.to_string(index=False)
    
    # configurations on the predicted hull that are calculated
    df_is_calculated = df_pred_hull[ df_pred_hull[is_calculated] == 1 ]
    
    # configurations on the predicted hull that are calculated and not on the calculated hull
    df_pred_gs_spurious = df_is_calculated[ df_is_calculated[calc_hull].astype(float) > hull_tol ]
    
    if len(df_pred_gs_spurious.index):
      print "\nSpuriously predicted ground states:"
      print df_pred_gs_spurious.to_string(index=False)
    
    self.df_train = df_train
    self.df_predict = df_predict
    self.df_calc_hull = df_calc_hull
    self.df_pred_hull = df_pred_hull
    self.df_pred_gs_uncalculated = df_pred_gs_uncalculated
    self.df_pred_gs_spurious = df_pred_gs_spurious
    

  def plot_hull(self, comp_axis='comp(a)', fig_filename='hull.eps', all_predicted=False):
    """ 
    Use results from 'check_hull' to plot training values, predicted values, hull
    """          
    df_train = self.df_train
    df_predict = self.df_predict
    df_calc_hull = self.df_calc_hull
    df_pred_hull = self.df_pred_hull
    calc_prop = 'formation_energy'
    pred_prop = 'clex(formation_energy)'
    
    xlabel = comp_axis
    ylabel = calc_prop
    
    
    fig1 = plt.figure()
    
    # plot calculated formation_energy
    plt.scatter(df_train[comp_axis], df_train[calc_prop], 
      s=40, marker='o', facecolors='none', edgecolors='r', label='Calculated')
    
    # plot calculated hull
    plt.plot(df_calc_hull[comp_axis], df_calc_hull[calc_prop],
      'r--', label='Calculated hull')
    
    # plot prediction of calculated hull states
    plt.plot(df_calc_hull[comp_axis], df_calc_hull[pred_prop],
      'gx--', label='Predicted calculated hull')
    
    if all_predicted:
      # plot predicted formation_energy (all predicted)
      plt.scatter(df_predict[comp_axis], df_predict[pred_prop],
        c='b', s=20, marker='o', edgecolors='none', label='Predicted')
    else:
      # plot predicted formation_energy (only for training set)
      plt.scatter(df_train[comp_axis], df_train[pred_prop],
        c='b', s=20, marker='o', edgecolors='none', label='Predicted')
    
    # plot predicted hull 
    plt.plot(df_pred_hull[comp_axis], df_pred_hull[pred_prop], 
      'b--', label='Predicted hull')
    
    
    # label axes
    xl = plt.xlabel(xlabel)
    yl = plt.ylabel(ylabel)
    
    # legend
    lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
      
    # save image
    print "\nSaving plot:", fig_filename
    plt.savefig(fig_filename, bbox_inches='tight', bbox_extra_artist=[xl, yl, lgd])
    
    if 'DISPLAY' in os.environ:
      plt.show()
  

if __name__ == "__main__":
  
  # select a project (with path=None, default uses project containing current directory)
  proj = Project(casm_exe='/Users/bpuchala/Work/codes/CASMcode/bin/casm')
  
  ## Check and plot hull using current eci 
  checker = ClexAnalysis(proj)
  
  checker.check_hull()
  
  checker.plot_hull()



