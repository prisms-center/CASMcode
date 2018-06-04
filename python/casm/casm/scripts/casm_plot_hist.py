from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

from casm.project import Project, Selection
import casm.plotting
import argparse
import os
import sys
import json
from bokeh.io import curdoc
from bokeh.client import push_session
import bokeh.plotting
import bokeh.models

input_help = "Input file"
desc_help = "Print extended usage description"
usage_desc = """
Histogram of CASM query output

Before running you must start the bokeh server:
- install bokeh with: 'pip install bokeh'
- start server: 'bokeh serve'

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
    

Example input file:

{
  "figure_kwargs": {
    "plot_height": 400,
    "plot_width": 800,
    "tools": "crosshair,pan,reset,resize,box_zoom"
  },
  "series": [
    {
      "project": null,
      "selection": "MASTER",
      "x": "comp(a)",
      "hist_kwargs": {
        "bins": 10
      }
    }
  ]
}

'style' example:

{
  "selected": {
    "color": "blue", 
  }
}

"""

def main():
    parser = argparse.ArgumentParser(description = 'Histogram plot')
    parser.add_argument('--desc', help=desc_help, default=False, action="store_true")
    parser.add_argument('input', nargs='?', help=input_help, type=str)
    
    
    # ignore 'mc casm'
    args = parser.parse_args(sys.argv[1:])
    
    if args.desc:
        print(usage_desc)
        sys.exit(1)
    
    with open(args.input, 'r') as f:
        input = json.load(f)
    
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
    
    # set up session
    session = push_session(curdoc())
    curdoc().add_root(fig)
    session.show() # open the document in a browser
    session.loop_until_closed() # run forever

if __name__ == "__main__":
    main()
