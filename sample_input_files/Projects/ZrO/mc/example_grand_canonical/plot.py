#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas

### open Monte Carlo results.json file
with open('results.json', 'r') as f:
  res = pandas.read_json(f)

### show data
for col in res.columns:
  print col

### basic plot
plt.plot(res['param_chem_pot(a)'], res['<comp(a)>'])
plt.show()
