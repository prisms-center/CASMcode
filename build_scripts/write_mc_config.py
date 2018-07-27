import json
import os
from os.path import join

config = {
  "apikey": os.environ['MC_API_KEY'],
  "mcurl": "https://materialscommons.org/api"
}

with open(join(os.environ['HOME'],'.materialscommons','config.json'),'w') as f:
  json.dump(config, f)
