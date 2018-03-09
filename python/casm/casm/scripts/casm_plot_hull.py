from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import bokeh.plotting
import bokeh.models

import casm.plotting

usage_desc = """
Plot the convex hull of calculated formation energies

If you have 'casm view' setup, then clicking on a configuration in the plot
will attempt to use 'casm view' to view that configuration.

Input file attributes:

  figure_kwargs: JSON object (optional)
    Input arguments for bokeh.models.Figure
  
  series: JSON array of one JSON object
    An array with one JSON object:

    project: str (optional, default=os.getcwd())
      Path to CASM project to get data from
    
    selection: str (optional, default="MASTER")
      Path to selection to use
    
    x: str
      Query to use for x-values. Should be one of 'comp(x)', 'comp_n(X)', or
      'atom_frac(X)'.
    
    y: str
      Query to use for y-values. Should be one of 'formation_energy' or 
      'formation_energy_per_atom'.
    
    tooltips: JSON array of str (optional, default=[])
      Additional properties to query and include in 'tooltips' info that appears
      when hovering over a data point. Set to 'null' to disable.
    
    dft_style: JSON object (optional, default=casm.plotting.dft_hull_style({})
      Visual styling attribute values for calculated energies. Defaults are 
      determined casm.plotting.dft_hull_style. The 'marker' value should be one 
      of the methods of bokeh.models.Figure
    
    clex_style: JSON object (optional, default=casm.plotting.dft_hull_style({})
      Visual styling attribute values for predicted energies. Defaults are 
      determined casm.plotting.clex_hull_style. The 'marker' value should be one 
      of the methods of bokeh.models.Figure
"""    

input_example = {
  "figure_kwargs": {
    "plot_height": 400,
    "plot_width": 800,
    "tools": "crosshair,pan,reset,box_zoom,wheel_zoom,save"
  },
  "series": [
    {
      "project": None,
      "selection": "MASTER",
      "hull_selection":"MASTER",
      "x": "comp(a)",
      "y": "formation_energy",
      "tooltips": [
        "scel_size"
      ]
    }
  ]
}

style_example = {
  "hull_line_dash": "", 
  "hover_alpha": 0.7, 
  "fill_alpha": 0.5, 
  "selected": {
    "color": "red", 
    "radii": 10
  }, 
  "hull_line_width": 1.0, 
  "on_hull": {
    "line_color": "red", 
    "line_width": 2.0, 
    "line_alpha": 0.8
  }, 
  "marker": "circle", 
  "hover_color": "orange", 
  "unselected": {
    "color": "gray", 
    "radii": 7
  }, 
  "off_hull": {
    "line_color": "blue", 
    "line_width": 0.0, 
    "line_alpha": 0.0
  },
  'clex_of_dft_hull_line_dash': "4 4"
}

class PlotHullCommand(casm.plotting.PlotTypeCommand):
    
    @classmethod
    def name(self):
        return "hull"
    
    @classmethod
    def short_desc(self):
        return "Plot convex hull"
    
    @classmethod
    def long_desc(self):
        return usage_desc
    
    @classmethod
    def style_example(self):
        return style_example
    
    @classmethod
    def input_example(self):
        return input_example
    
    @classmethod
    def plot(self, doc, args):
        with open(args.input, 'rb') as f:
            input = json.loads(f.read().decode('utf-8'))
        
        data = casm.plotting.PlottingData()
        figure_kwargs = input.get('figure_kwargs', casm.plotting.default_figure_kwargs)
        fig = bokeh.plotting.Figure(**figure_kwargs)
        tap_action = casm.plotting.TapAction(data)
        
        hullplot = casm.plotting.ConvexHullPlot(data, **input['series'][0])
        hullplot.query()
        hullplot.plot(fig, tap_action)
        
        # add tools
        fig.add_tools(tap_action.tool())
        fig.add_tools(bokeh.models.BoxSelectTool(renderers=hullplot.renderers))
        fig.add_tools(bokeh.models.LassoSelectTool(renderers=hullplot.renderers))
        
        doc.add_root(fig)

