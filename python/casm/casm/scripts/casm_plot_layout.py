from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import copy
import json
import bokeh.layouts
import bokeh.plotting
import bokeh.models

import casm.plotting

usage_desc = """
Layout several plots of CASM query output

If you have 'casm view' setup, then clicking on a configuration in the plot
will attempt to use 'casm view' to view that configuration.

Input file attributes:

  n_plots_each_row: JSON array
    Number of plots in each row.  Ex: [1, 2] results in 1 plot in the first row
    and two plots in the second row
  
  subplots: JSON array of JSON objects
    Each element of 'subplots' specifies the position of the subplot and 
    provides the required data for plotting.
  
    pos: array of two int
      Row and index of the subplot, as: [row, index]. Index from 0.

    figure_kwargs: JSON object (optional)
      Input arguments for bokeh.models.Figure
  
    series: JSON array of JSON objects
      An array of JSON objects, each containing information for one series to be
      plotted. See 'casm.plot.hull' and 'casm.plot.scatter' for format. Additionally,
      the 'type' is required:

      type: str
        One of 'histogram', 'hull', 'rankplot', 'scatter'
"""

input_example = {
  "n_plots_each_row": [1, 2],
  "subplots": [
    {
      "pos": [0, 0],
      "figure_kwargs": {
        "plot_height": 400,
        "plot_width": 800,
        "tools": "crosshair,pan,reset,box_zoom,wheel_zoom,save"
      },
      "series": [
        {
          "type": "hull",
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
    },
    {
      "pos": [1, 0],
      "figure_kwargs": {
        "plot_height": 400,
        "plot_width": 400,
        "tools": "crosshair,pan,reset,box_zoom,wheel_zoom,save"
      },
      "series": [
        {
          "type": "scatter",
          "project": None,
          "selection": "MASTER",
          "x": "basis_deformation",
          "y": "formation_energy",
          "tooltips": [
            "scel_size",
            "volume_relaxation"
          ],
          "legend":"formation_energy"
        }
      ]
    },
    {
      "pos": [1, 1],
      "figure_kwargs": {
        "plot_height": 400,
        "plot_width": 400,
        "tools": "crosshair,pan,reset,box_zoom,wheel_zoom,save"
      },
      "series": [
        {
          "type": "scatter",
          "project": None,
          "selection": "MASTER",
          "x": "lattice_deformation",
          "y": "formation_energy",
          "tooltips": [
            "scel_size",
            "volume_relaxation"
          ],
          "legend":"formation_energy"
        }
      ]
    }
  ]
}

def make_layout(val, n_plots_each_row):
    layout = []
    for n in n_plots_each_row:
        layout.append([copy.deepcopy(val) for i in range(n)])
    return layout

class PlotLayoutCommand(casm.plotting.PlotTypeCommand):
    
    @classmethod
    def name(cls):
        return "layout"
    
    @classmethod
    def short_desc(cls):
        return "Layout several plots"
    
    @classmethod
    def long_desc(cls):
        return usage_desc
    
    @classmethod
    def style_example(cls):
        return style_example
    
    @classmethod
    def input_example(cls):
        return input_example
    
    @classmethod
    def plot(cls, doc, args):
        with open(args.input, 'rb') as f:
            layout_input = json.loads(f.read().decode('utf-8'))
        
        data = casm.plotting.PlottingData()
        layout = make_layout(None, layout_input['n_plots_each_row'])
        
        options = {
          'hull':casm.plotting.ConvexHullPlot,
          'scatter':casm.plotting.Scatter,
          'rankplot':casm.plotting.RankPlot,
          'histogram':casm.plotting.Histogram
        }
        for subplot in layout_input['subplots']:
            for index, series in enumerate(subplot['series']):
                series_name = str(subplot['pos'][0]) + "." + str(subplot['pos'][1]) + "." + str(index)
                series['self'] = options[series['type']](data=data, index=index, series_name=series_name, **series)
        
        # first query data necessary for all series in every figure
        for subplot in layout_input['subplots']:
            for series in subplot['series']:
                series['self'].query()
            
        # next create plots
        for subplot in layout_input['subplots']:
        
            figure_kwargs = subplot.get('figure_kwargs', casm.plotting.default_figure_kwargs)
        
            r = subplot['pos'][0]
            c = subplot['pos'][1]
            
            ## Construct a figure
            fig = bokeh.plotting.Figure(**figure_kwargs)
            tap_action = casm.plotting.TapAction(data)
            renderers = []
            for series in subplot['series']:
                series['self'].plot(fig, tap_action=tap_action)
                renderers += series['self'].renderers
            
            # add tools
            fig.add_tools(tap_action.tool())
            fig.add_tools(bokeh.models.BoxSelectTool(renderers=renderers))
            fig.add_tools(bokeh.models.LassoSelectTool(renderers=renderers))
            layout[r][c] = fig
        
        gplot = bokeh.layouts.layout(layout, plot_width=400, plot_height=400)
        
        doc.add_root(gplot)
