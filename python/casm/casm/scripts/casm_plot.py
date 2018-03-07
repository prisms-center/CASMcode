from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import argparse
import subprocess
import sys
from bokeh.server.server import Server

import casm.scripts.casm_plot_hist
import casm.scripts.casm_plot_hull
import casm.scripts.casm_plot_layout
import casm.scripts.casm_plot_rankplot
import casm.scripts.casm_plot_scatter

plot_types = {
    "hist": casm.scripts.casm_plot_hist.main,
    "hull": casm.scripts.casm_plot_hull.main,
    "layout": casm.scripts.casm_plot_layout.main,
    "rankplot": casm.scripts.casm_plot_rankplot.main,
    "scatter": casm.scripts.casm_plot_scatter.main,
    }

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(description = 'Make useful plots')
    parser.add_argument('type', help="Plot type", type=str, default="", nargs="?", metavar="<type>")
    parser.add_argument('args', help="Plotting arguments", nargs=argparse.REMAINDER, metavar="...")
    parser.add_argument('--desc', help="Print plot type list", action="store_true", default=False)
    #parser.add_argument('--serve', help="Start bokeh server for plots", action="store_true", default=False)
    args = parser.parse_args(argv)

#     if args.serve:
#         subprocess.Popen(['bokeh', 'serve'])
#         return
    
    if args.type in plot_types:
        plot_types[args.type](argv[1:])
    elif args.desc:
        parser.print_help()
        print("")
        print("available plot types:")
        for type in plot_types:
            print("  " + type)
        print("")
        print("For help using a plot type: 'casm plot <type> --help'\n")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()