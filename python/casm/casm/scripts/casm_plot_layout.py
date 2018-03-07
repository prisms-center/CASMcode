from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

from casm.project import Project, Selection
import casm.plotting
import argparse
import os
import sys
import json
import copy
from bokeh.io import curdoc
from bokeh.client import push_session
from bokeh.server.server import Server
import bokeh.plotting
import bokeh.models
import bokeh.layouts

input_help = "Input file"
desc_help = "Print extended usage description"
usage_desc = """
Layout plots of CASM query output

Before running you must start the bokeh server:
- install bokeh with: 'pip install bokeh'
- start server: 'bokeh serve'

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
        One of 'hull', 'scatter'.

Example input file:

{
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
          ... hull plot input ...
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
          ... scatter plot input ...
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
          ... scatter plot input ...
        }
      ]
    }
  ]
}

"""

def make_layout(val, n_plots_each_row):
    layout = []
    for n in n_plots_each_row:
        layout.append([copy.deepcopy(val) for i in range(n)])
    return layout

def plot(doc, args):
    with open(args.input, 'r') as f:
        layout_input = json.load(f)
    
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

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(description = 'Layout plot')
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
