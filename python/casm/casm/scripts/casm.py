from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

import argparse
import json
import os
import pty
import subprocess
import sys

import casm.scripts.casm_calc
import casm.scripts.casm_learn
import casm.scripts.casm_plot
from casm.project import API, Project

def _exec(argv=None):
    subprocess.Popen(argv)

def _shell(argv=None):
    orig = os.environ['PS1']
    os.environ['PS1'] = "<shell> " + os.environ['PS1']
    pty.spawn('/bin/bash')
    os.environ['PS1'] = orig
        

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    try:
        libcasm_commands = json.loads(API().command_list())
        libcasm_commands.remove('--version')
        found_libcasm = True
    except:
        libcasm_commands = []
        found_libcasm = False
    
    python_commands = {
        'calc': casm.scripts.casm_calc.main,
        'learn': casm.scripts.casm_learn.main,
        'plot': casm.scripts.casm_plot.main,
        'shell': _shell,
        'exec': _exec}
    
    commands = sorted(libcasm_commands + list(python_commands.keys()))
    
    if not found_libcasm:
        print("Could not find libcasm. Please check your installation.")
        return
    
    parser = argparse.ArgumentParser(description = 'CASM: First-principles based statistical mechanics')
    parser.add_argument('command', help="CASM command to execute", type=str, default="", nargs="?", metavar="<command>")
    parser.add_argument('args', help="CASM command arguments", nargs=argparse.REMAINDER, metavar="...")
    parser.add_argument('--desc', help="Print command list", action="store_true", default=False)
    parser.add_argument('--version', help="Print casm version", action="store_true", default=False)
    parser.add_argument('--path', help="Path to project. Default uses project containing current working directory.", type=str, default="")
    
    args = parser.parse_args(argv)
    
    if args.version:
        args.command = "version"
        args.args = []
    if args.desc:
        parser.print_help()
        print("")
        print("available commands:")
        for cmd in commands:
            print("  " + cmd)
        print("")
        print("For help using a command: 'casm <command> --help'\n")
        print("For step by step help use: 'casm status -n'\n")
    elif args.command in python_commands:
        python_commands[args.command](argv[1:])
    elif args.command in libcasm_commands:
        _api = API()
        primclex = _api.primclex_null()
        log = _api.stdout()
        debug_log = log
        err_log = log
        args.args = ["\'" + x + "\'" for x in args.args]
        res = _api(' '.join([args.command] + args.args), primclex, args.path, log, debug_log, err_log)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()