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
from bokeh.server.server import Server
import bokeh.plotting
import bokeh.models

input_help = "Input file"
desc_help = "Print extended usage description"
usage_desc = """
Rank plot of CASM query output

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
    
    to_query: List[str] (optional, default=[])
      List of queries to be made, for instance if necessary for scoring.
    
    scoring_query: str
      Query to use for scoring. If not in 'query', will be added automatically
    
    scoring_module: str (optional, default=None)
      If 'scoring_query' is None, will attempt load the 'scoring_function' from
      the 'scoring_module' and use that for scoring.
    
    scoring_function: str (optional, default=None)
      If 'scoring_query' is None, will attempt load the 'scoring_function' from
      the 'scoring_module' and use that for scoring. Expected signature is:
        pandas.Series = f(casm.project.Selection)
    
    tooltips: JSON array of str (optional, default=[])
      Additional properties to query and include in 'tooltips' info that appears
      when hovering over a data point. Set to 'null' to disable.
    
    style: JSON object (optional, default=casm.plotting.scatter_series_style(<series index>, {})
      Visual styling attribute values. Defaults are determined casm.plotting.scatter_series_style.
      The 'marker' value should be one of the methods of bokeh.models.Figure
    

Example input file:

{
  "figure_kwargs": {
    "plot_height": 400,
    "plot_width": 800,
    "tools": "crosshair,pan,reset,box_zoom,wheel_zoom,save"
  },
  "series": [
    {
      "project": null,
      "selection": "MASTER",
      "scoring_query": "hull_dist(CALCULATED,atom_frac)",
      "to_query": [
        "configname",
        "basis_deformation",
        "lattice_deformation"
      ]
    }
  ]
}

'style' example:

{
  "hover_alpha": 0.7, 
  "hover_color": "orange", 
  "selected": {
    "line_color": "blue", 
    "line_alpha": 0.0, 
    "line_width": 0.0, 
    "color": "blue", 
    "fill_alpha": 0.5, 
    "radii": 10
  }, 
  "unselected": {
    "line_color": "gray", 
    "line_alpha": 0.0, 
    "line_width": 0.0, 
    "color": "gray", 
    "fill_alpha": 0.3, 
    "radii": 7
  }, 
  "marker": "circle"
}

"""

def plot(doc, args):
    with open(args.input, 'r') as f:
        input = json.load(f)
    
    data = casm.plotting.PlottingData()
    figure_kwargs = input.get('figure_kwargs', casm.plotting.default_figure_kwargs)
    fig = bokeh.plotting.Figure(**figure_kwargs)
    tap_action = casm.plotting.TapAction(data)
    renderers = []
    
    for index, series in enumerate(input['series']):
        series['self'] = casm.plotting.RankPlot(data=data, index=index, **series)
    
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

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(description = 'Rank plot')
    parser.add_argument('--desc', help=desc_help, default=False, action="store_true")
    parser.add_argument('input', nargs='?', help=input_help, type=str)
    args = parser.parse_args(argv)
    
    if args.desc:
        print(usage_desc)
        return
    elif args.input is not None:
        def f(doc):
            plot(doc, args)
        server = Server({'/': f}, num_procs=1)
        server.start()
        
        print('Opening on http://localhost:5006/')
        print('Enter Ctrl+C to stop')
        try:
            server.io_loop.add_callback(server.show, "/")
            server.io_loop.start()
        except KeyboardInterrupt as e:
            print('\nStopping...')
            pass
    else:
        parser.print_help()
        return

