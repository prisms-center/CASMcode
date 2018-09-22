import abc

import argparse
import json
from bokeh.server.server import Server

class PlotTypeCommand(abc.ABC):
    
    @classmethod
    @abc.abstractmethod
    def name(cls):
        return "Derived types need to implement name"
    
    @classmethod
    @abc.abstractmethod
    def short_desc(cls):
        return "Derived types need to implement short_desc, for argparse description"
    
    @classmethod
    @abc.abstractmethod
    def long_desc(cls):
        return "Derived types need to implement short_desc, for --desc output"
    
    @classmethod
    @abc.abstractmethod
    def style_example(cls):
        return "Derived types need to implement style_example for --desc output"
    
    @classmethod
    @abc.abstractmethod
    def input_example(cls):
        return "Derived types need to implement example_input, for --example output"
    
    @classmethod
    @abc.abstractmethod
    def plot(cls, doc, args):
        """Create the bokeh Document to be shown"""
        return
    
    @classmethod
    def run(cls, argv=None, port=5006):
        if argv is None:
            argv = sys.argv[1:]
        parser = argparse.ArgumentParser(description = cls.short_desc())
        parser.add_argument('input', nargs="?", help="Input file", type=str)
        parser.add_argument('--example', help="Print example input file", default=False, action="store_true")
        parser.add_argument('--desc', help="Print extended usage description", default=False, action="store_true")
        args = parser.parse_args(argv)
        
        if args.desc:
            print(cls.long_desc())
            print('')
            print("casm plot " + cls.name() + " 'style' example:")
            print(json.dumps(cls.style_example(), indent=2))
            return
        elif args.example:
            print(json.dumps(cls.input_example(), indent=2))
            return
        elif args.input is not None:
            def f(doc):
                cls.plot(doc, args)
            server = Server({'/': f}, num_procs=1, port=port)
            server.start()
            
            print('Opening on http://localhost:' + str(port) + '/')
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
