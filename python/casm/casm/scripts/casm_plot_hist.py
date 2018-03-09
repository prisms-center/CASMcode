from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import casm.plotting
import json
import bokeh.plotting
import bokeh.models

usage_desc = """
Histogram of CASM query output

If you have 'casm view' setup, then clicking on a configuration in the plot
will attempt to use 'casm view' to view that configuration.

Input file attributes:

  figure_kwargs: JSON object (optional)
    Input arguments for bokeh.models.Figure
  
  series: JSON array of JSON objects
    An array with one JSON object:

    project: str (optional, default=os.getcwd())
      Path to CASM project to get data from
    
    selection: str (optional, default="MASTER")
      Path to selection to use
    
    x: str
      Query to use for x-values
    
    hist_kwargs: JSON object (optional, default={})
      Arguments to pass to numpy.histogram to generate the histogram.
    
    style: JSON object (optional, default=casm.plotting.scatter_series_style(<series index>, {})
      Visual styling attribute values. Defaults are determined casm.plotting.scatter_series_style.
      The 'marker' value should be one of the methods of bokeh.models.Figure
    
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
      "x": "comp(a)",
      "hist_kwargs": {
        "bins": 10
      }
    }
  ]
}

style_example = {
  "selected": {
    "color": "blue", 
  }
}

class PlotHistCommand(casm.plotting.PlotTypeCommand):
    
    @classmethod
    def name(cls):
        return "hist"
    
    @classmethod
    def short_desc(cls):
        return "Histogram plot"
    
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
            input = json.loads(f.read().decode('utf-8'))
        
        data = casm.plotting.PlottingData()
        figure_kwargs = input.get('figure_kwargs', casm.plotting.default_figure_kwargs)
        fig = bokeh.plotting.Figure(**figure_kwargs)
        tap_action = casm.plotting.TapAction(data)
        renderers = []
        
        for index, series in enumerate(input['series']):
            series['self'] = casm.plotting.Histogram(data=data, index=index, **series)
        
        # first query data necessary for all series
        for series in input['series']:
            series['self'].query()
        
        for series in input['series']:
            series['self'].plot(fig, tap_action)
            renderers += series['self'].renderers
        
        # add tools
        fig.add_tools(tap_action.tool())
        fig.add_tools(bokeh.models.BoxSelectTool(renderers=renderers))
        fig.add_tools(bokeh.models.LassoSelectTool(renderers=renderers))
        
        doc.add_root(fig)

    