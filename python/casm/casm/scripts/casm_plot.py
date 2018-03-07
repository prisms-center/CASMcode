from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import argparse
import subprocess
import sys

import casm.scripts.casm_plot_hist as cp_hist
import casm.scripts.casm_plot_hull as cp_hull
import casm.scripts.casm_plot_layout as cp_layout
import casm.scripts.casm_plot_rankplot as cp_rankplot
import casm.scripts.casm_plot_scatter as cp_scatter

plot_type_cls = [
    cp_hist.PlotHistCommand,
    cp_hull.PlotHullCommand,
    cp_layout.PlotLayoutCommand,
    cp_rankplot.PlotRankplotCommand,
    cp_scatter.PlotScatterCommand]

plot_types = {cls.name():cls.run for cls in plot_type_cls}

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    parser = argparse.ArgumentParser(description = 'Make useful plots')
    parser.add_argument('type', help="Plot type", type=str, default="", nargs="?", metavar="<type>")
    parser.add_argument('args', help="Plotting arguments", nargs=argparse.REMAINDER, metavar="...")
    parser.add_argument('--desc', help="Print plot type list", action="store_true", default=False)
    args = parser.parse_args(argv)
    
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
