#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas
from os.path import join, exists

index = 0
while True:
  path = 'path.' + str(index)
  resultsfile = join(path,'results.json')
  
  if not exists(resultsfile):
    break
  
  ### open Monte Carlo results.json file
  with open(resultsfile, 'r') as f:
    res = pandas.read_json(f)
  
  index += 1
  
  plt.plot(res['<comp(a)>'], res['T'], 'bo-')
plt.show()
