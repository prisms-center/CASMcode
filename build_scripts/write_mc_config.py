import json
import os
from os.path import join, exists, dirname

config = {
  "apikey": os.environ['MC_API_KEY'],
  "mcurl": "https://materialscommons.org/api"
}
filename = join(os.environ['HOME'],'.materialscommons','config.json')
if not exists(dirname(filename)):
    os.makedirs(dirname(filename))
with open(filename,'w') as f:
    json.dump(config, f, indent=2)
