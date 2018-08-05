from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

import json
import os
import sys
from os.path import join, exists, dirname

if 'MC_API_KEY' not in os.environ:
    print("write_mc_config.py failed: No MC_API_KEY")
    sys.exit(0)

config = {
  "apikey": os.environ['MC_API_KEY'],
  "mcurl": "https://materialscommons.org/api"
}
filename = join(os.environ['HOME'],'.materialscommons','config.json')
if not exists(dirname(filename)):
    os.makedirs(dirname(filename))
with open(filename,'w') as f:
    json.dump(config, f, indent=2)
