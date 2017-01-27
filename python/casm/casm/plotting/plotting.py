import bokeh.plotting
import bokeh.models
import numpy as np
import pandas
import uuid
import casm.project

import os
import casm.learn
import json
import pickle
import re

import bokeh.client
import bokeh.io

from bokeh.plotting import hplot, vplot

class Session(object):
  def __init__(self, doc=None):
    if doc is None:
      doc = bokeh.io.curdoc()
    self.doc = doc
    self.session = bokeh.client.push_session(self.doc)
  
  def add_root(self, layout):
    self.doc.add_root(layout)
  
  def begin(self):
    self.session.show()
  
  def begin_interactive(self):
    self.session.show()
    self.session.loop_until_closed()
  
  def close(self):
    self.session.close()

class Messages(object):
  """TextInput for writing messages"""
  def __init__(self):
    self.widget = bokeh.models.TextInput(value='...', title="Messages:")
  
  @property
  def value(self):
    return float(self.widget.value)
  
  @value.setter
  def value(self, _value):
    self.widget.value = str(_value)

def _selected(sel):
  """Selected configurations"""
  if 'to_be_selected' in sel.data.columns:
    return sel.data['to_be_selected']==True
  else:
    return sel.data['selected']==True

def _unselected(sel):
  """Unselected configurations"""
  if 'to_be_selected' in sel.data.columns:
    return sel.data['to_be_selected']==False
  else:
    return sel.data['selected']==False
  

def hull_dist(path="CALCULATED"):
  """
  Equivalent to 'hull_dist(' + path + ',comp)'
  
  Arguments:
    path: a CASM Selection path, default "CALCULATED"
  
  """
  return 'hull_dist(' + path + ',comp)'

def hull_dist_per_atom(path="CALCULATED"):
  """
  Equivalent to 'hull_dist(' + path + ',atom_frac)'
  
  Arguments:
    path: a CASM Selection path, default "CALCULATED"
  
  """
  return 'hull_dist(' + path + ',atom_frac)'

def clex_hull_dist(path="CALCULATED"):
  """
  Equivalent to 'clex_hull_dist(' + path + ',comp)'
  
  Arguments:
    path: a CASM Selection path, default "CALCULATED"
  
  """
  return 'clex_' + hull_dist(path)

def clex_hull_dist_per_atom(path="CALCULATED"):
  """
  Equivalent to 'hull_dist(' + path + ',atom_frac)'
  
  Arguments:
    path: a CASM Selection path, default "CALCULATED"
  
  """
  return 'clex_' + hull_dist_per_atom(path)

def add_src_data(sel, name, data, force=False):
  if sel.src is None:
    sel.src = bokeh.models.ColumnDataSource(data={name:data})
    sel.selection_callbacks = []
    def callbacks():
      for f in sel.selection_callbacks:
        f()
    sel.update_selection = callbacks
  else:
    if name not in sel.src.data.keys() or force == True:
      sel.src.data[name] = data


def view_on_tap(sel, attrname, old, new):
  """
  Execute 'casm view' on tapped configuration
  """
  if len(sel.src.selected['1d']['indices']) == 1:
    index = sel.src.selected['1d']['indices'][0]
    configname = sel.src.data['configname'][index]
    args = "view " + configname
    sel.proj.command(args=args)


class ConfigurationTapAction(object):
  """
  Enable execution of a python function when a glyph is 'tapped' in a bokeh plot.
  
  Example:
    p = ... a Bokeh plot
    circles = p.circ(... A set of circles) 
    my_tap_action = ConfigurationTapAction(sel)
    def f(sel, attrname, old, new):
      if len(sel.src.selected['1d']['indices']) == 1:
        # the index of the tapped configuration:
        index = sel.src.selected['1d']['indices'][0]
        ... code to execute ...
        ... basically can ignore attrname, old, new ...
      
    # call f when a circle is tapped
    p.add_tools(my_tap_action.tool(f, [circles]))

  """
  
  def __init__(self, sel):
    self.tap_indicator_src = bokeh.models.ColumnDataSource(data={"x":[0]})
    self.js_callback = bokeh.models.CustomJS(args={'src':self.tap_indicator_src}, code="""
      var data = src.get('data');
      data['x'][0] = data['x'][0] + 1;
      src.set('data', data);
      """)
    self.sel = sel
    self.py_callback = None
  
  def __call__(self, attrname, old, new):
    self.py_callback(self.sel, attrname, old, new)
  
  def tool(self, callback=None, renderers=None):
    self.py_callback = callback
    self.tap_indicator_src.on_change('data', self)
    return bokeh.models.TapTool(callback=self.js_callback, renderers=renderers)


dft_style = {
  
  'marker': 'circle',
  
  'selected_color': "blue",
  'selected_radii': 10,
  'unselected_color': "green",
  'unselected_radii': 7,
  'fill_alpha': 0.5,
  
  'cutoff_color': "red",
  'cutoff_line_width': 1,
  'cutoff_alpha': 1.0,

  'off_hull_line_color': "blue",
  'on_hull_line_color': "#FF4105",
  'off_hull_line_width': 0.0,
  'on_hull_line_width': 2.0,
  'on_hull_line_alpha': 0.8,
  'off_hull_line_alpha': 0.0,
  
  'hull_line_width': 1.0,
  'hull_line_dash': "",

  'hover_color': "orange",
  'hover_alpha': 0.7
}

clex_style = {
  
  'marker': 'x',
  
  'selected_color': "#05ff58",
  'selected_radii': 10,
  'unselected_color': "blue",
  'unselected_radii': 7,
  'fill_alpha': 0.5,
  
  'cutoff_color': "red",
  'cutoff_line_width': 1,
  'cutoff_alpha': 1.0,

  'off_hull_line_color': "#05ff58",
  'on_hull_line_color': "#05ff58",
  'off_hull_line_width': 2.0,
  'on_hull_line_width': 2.0,
  'on_hull_line_alpha': 1.0,
  'off_hull_line_alpha': 1.0,
  
  'hull_line_width': 1.0,
  'hull_line_dash': "",
  'clex_of_dft_hull_line_dash': "4 4",

  'hover_color': "orange",
  'hover_alpha': 0.7
}

def update_dft_glyphs(sel, selected=None):
  """
  Use the selected and on_hull iterables to update glyph styles based on sel.dft_style.
  """
  
  if selected is None:
    selected = sel.data['selected']
  
  style = sel.dft_style
  
  ## dft points:
  # colors to plot selected (True), and unselected (False) configurations
  colormap = dict({True:style['selected_color'], False:style['unselected_color']})
  
  # radii to plot selected (True), and unselected (False) configurations
  radiimap = dict({True:style['selected_radii'], False:style['unselected_radii']})
  
  sel.src.data['color'] = map(lambda x: colormap[x], selected)
  sel.src.data['radii'] = map(lambda x: radiimap[x], selected)
  
  
  ## on dft hull:
  # colors to plot on_hull (True), and not on_hull (False) configurations
  linecolormap = dict({True:style['on_hull_line_color'], False:style['off_hull_line_color']})
  
  # width to plot on_hull (True), and not on_hull (False) configurations
  linewidthmap = dict({True:style['on_hull_line_width'], False:style['off_hull_line_width']})
  
  # width to plot on_hull (True), and not on_hull (False) configurations
  linealphamap = dict({True:style['on_hull_line_alpha'], False:style['off_hull_line_alpha']})
  
  on_hull = None
  if 'on_dft_hull' in sel.data.columns:
    on_hull = sel.data['on_dft_hull']
  
  if on_hull is not None:
    sel.src.data['line_color'] = map(lambda x: linecolormap[x], on_hull)
    sel.src.data['line_width'] = map(lambda x: linewidthmap[x], on_hull)
    sel.src.data['line_alpha'] = map(lambda x: linealphamap[x], on_hull)
  
def update_clex_glyphs(sel, selected=None):
  """
  Use the selected and on_hull iterables to update glyph styles based on sel.clex_style.
  """
  
  if selected is None:
    selected = sel.data['selected']
  
  ## clex points:
  
  style = sel.clex_style
  
  # colors to plot selected (True), and unselected (False) configurations
  colormap = dict({True:style['selected_color'], False:style['unselected_color']})
  
  # radii to plot selected (True), and unselected (False) configurations
  radiimap = dict({True:style['selected_radii'], False:style['unselected_radii']})
  
  sel.src.data['clex_color'] = map(lambda x: colormap[x], selected)
  sel.src.data['clex_radii'] = map(lambda x: radiimap[x], selected)
  
  
  ## on clex hull:
  # colors to plot on_hull (True), and not on_hull (False) configurations
  linecolormap = dict({True:style['on_hull_line_color'], False:style['off_hull_line_color']})
  
  # width to plot on_hull (True), and not on_hull (False) configurations
  linewidthmap = dict({True:style['on_hull_line_width'], False:style['off_hull_line_width']})
  
  # width to plot on_hull (True), and not on_hull (False) configurations
  linealphamap = dict({True:style['on_hull_line_alpha'], False:style['off_hull_line_alpha']})
  
  on_clex_hull = None
  if 'on_clex_hull' in sel.data.columns:
    on_clex_hull = sel.data['on_clex_hull']
  
  if on_clex_hull is not None:
    sel.src.data['clex_line_color'] = map(lambda x: linecolormap[x], on_clex_hull)
    sel.src.data['clex_line_width'] = map(lambda x: linewidthmap[x], on_clex_hull)
    sel.src.data['clex_line_alpha'] = map(lambda x: linealphamap[x], on_clex_hull)

def update_glyphs(sel, selected=None):
  update_dft_glyphs(sel, selected)


class ConvexHullPlot(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
    y: the value of the y-axis, 'formation_energy' by default
  """
  
  def __init__(self, sel, hull_sel=None, hull_tol=1e-8, x='comp(a)',
      y='formation_energy', dft_style=dft_style, clex_style=clex_style):
    """
    Arguments:
      sel: A CASM Selection
      x: (str) column name to use as x axis. Default='comp(a)'. Must be 'comp(x)', 'comp_n(X)', or 'atom_frac(X)'.
      y: (str) column name to use as y axis. Default='formation_energy'. Must be 'formation_energy' or 'formation_energy_per_atom'.
      hull_sel: (Selection or None) a CASM Selection, to use for the hull. Default uses sel.
      hull_tol: (number, default 1e-8) tolerance to decide if configuration is on the hull
    """
    self.sel = sel
    if hull_sel is None:
      hull_sel = sel
    self.hull_sel = hull_sel
    self.hull_sel_calculated = casm.project.Selection(sel.proj, hull_sel.path + ".calculated")
    self.dft_hull_sel = casm.project.Selection(sel.proj, hull_sel.path + ".dft_hull")
    
    self.is_comp = False
    self.is_comp_n = False
    self.is_atom_frac = False
    if re.match('\s*comp\(.*\)\s*', x):
      self.is_comp = True
      self.x_all = "comp"
    elif re.match('\s*comp_n\(.*\)\s*', x):
      self.is_comp_n = True
      self.x_all = "comp_n"
    elif re.match('\s*atom_frac\(.*\)\s*', x):
      self.is_atom_frac = True
      self.x_all = "atom_frac"
    else:
        raise ValueError("ConvexHullPlot x: '" + x + "' is not allowed.")
    self.x = x
    
    if y == 'formation_energy':
        self.Ef_per_atom = False
        self.y_clex = 'clex(formation_energy)'
    elif y == 'formation_energy_per_atom':
        self.Ef_per_atom = True
        self.y_clex = 'clex(formation_energy,per_species)'
    else:
        raise ValueError("ConvexHullPlot y: '" + y + "' is not allowed.")
    self.y = y
    
    self.hull_tol = hull_tol
    
    self.sel.dft_style = dft_style
    self.sel.clex_style = clex_style
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    
    #add data to sel.src if necessary
    self._query()
    
    # will store the plot object in self.p_
    self._plot()
    
    # store plot in layout
    self.layout = self.p
    
    # callbacks
    #self.sel.selection_callbacks.append(self.update)
    
    
  
  @property
  def p(self):
    return self.p_
  
  def _plot(self):
    
    update_dft_glyphs(self.sel)
    update_clex_glyphs(self.sel)
    
    # plot formation energy vs comp, with convex hull states
    
    # dft
    style = self.sel.dft_style
    self.p_ = bokeh.plotting.Figure(plot_width=800, plot_height=400, tools=self.tools)
    cr = getattr(self.p_, style['marker'])(self.x, self.y, source=self.sel.src, 
      size='radii', fill_color='color', fill_alpha=style['fill_alpha'], 
      line_color='line_color', line_width='line_width', line_alpha='line_alpha',
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'],
      legend="dft")
    self.hull_line = self.p_.line(self.x, self.y, source=self.calc_hull_src, 
      line_color=style['on_hull_line_color'], line_width=style['hull_line_width'],
      line_dash=style['hull_line_dash'],
      legend="dft")
    
    # clex
    style = self.sel.clex_style
    cr_clex = getattr(self.p_, style['marker'])(self.x, self.y_clex, source=self.sel.src, 
      size='clex_radii', fill_color='clex_color', fill_alpha=style['fill_alpha'], 
      line_color='clex_line_color', line_width='clex_line_width', line_alpha='clex_line_alpha',
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'],
      legend="clex")
    self.clex_hull_line = self.p_.line(self.x, self.y_clex, source=self.clex_hull_src, 
      line_color=style['on_hull_line_color'], line_width=style['hull_line_width'],
      line_dash=style['hull_line_dash'],
      legend="clex")
    
    # clex of dft hull
    style = self.sel.clex_style
    self.clex_of_dft_hull_line = self.p_.line(self.x, self.y_clex, source=self.clex_of_dft_hull_src, 
      line_color=style['on_hull_line_color'], line_width=style['hull_line_width'],
      line_dash=style['clex_of_dft_hull_line_dash'],
      legend="clex_of_dft_hull")

    self.p_.xaxis.axis_label = self.x
    self.p_.yaxis.axis_label = self.y

    # hover over a point to see 'configname', 'Ef', and 'comp(a)'
    tooltips = [
        ("configname","@configname")
    ]
    
    # tooltips for energy
    if self.Ef_per_atom:
      tooltips.append(("Ef_per_atom","@formation_energy_per_atom{1.1111}"))
      tooltips.append(("clex_Ef_per_atom","@{clex(formation_energy,per_species)}{1.1111}"))
      tooltips.append(("hull_dist_per_atom", "@{" + hull_dist_per_atom(self.hull_sel_calculated.path) + "}{1.1111}"))
      tooltips.append(("clex_hull_dist_per_atom", "@{" + clex_hull_dist_per_atom(self.hull_sel.path) + "}{1.1111}"))
      tooltips.append(("clex_of_dft_hull_dist_per_atom", "@{" + clex_hull_dist_per_atom(self.dft_hull_sel.path) + "}{1.1111}"))
    else:
      tooltips.append(("Ef","@formation_energy{1.1111}"))
      tooltips.append(("clex_Ef","@{clex(formation_energy)}{1.1111}"))
      tooltips.append(("hull_dist", "@{" + hull_dist(self.hull_sel_calculated.path) + "}{1.1111}"))
      tooltips.append(("clex_hull_dist", "@{" + clex_hull_dist(self.hull_sel.path) + "}{1.1111}"))
      tooltips.append(("clex_of_dft_hull_dist", "@{" + clex_hull_dist(self.dft_hull_sel.path) + "}{1.1111}"))
    
    # tooltips for composition
    if self.is_comp:
      for key in self.sel.src.data.keys():
        if re.match('\s*comp\(.*\)\s*', key):
          tooltips.append((key, "@{" + key + "}{1.1111}"))
    elif self.is_comp_n:
      for key in self.sel.src.data.keys():
        if re.match('\s*comp_n\(.*\)\s*', key):
          tooltips.append((key, "@{" + key + "}{1.1111}"))
    elif self.is_atom_frac:
      for key in self.sel.src.data.keys():
        if re.match('\s*atom_frac\(.*\)\s*', key):
          tooltips.append((key, "@{" + key + "}{1.1111}"))
    
    self.tap_action = ConfigurationTapAction(self.sel)

    self.p_.add_tools(self.tap_action.tool(view_on_tap, [cr, cr_clex]))
    self.p_.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[cr, cr_clex]))
    self.p_.add_tools(bokeh.models.BoxSelectTool(renderers=[cr, cr_clex]))
    self.p_.add_tools(bokeh.models.LassoSelectTool(renderers=[cr, cr_clex]))
    
  
  def _set_on_hull(self, sel, on_hull_label, hull_dist_label):
    """
    Use the 'hull_dist_label' column, and self.hull_tol to create the boolean
    column 'on_hull_label' indicating which configurations are on the hull.
    """
    sel.data.loc[_unselected(sel),on_hull_label] = False
    on_hull = map(lambda x: abs(x) < self.hull_tol, sel.data.loc[_selected(sel), hull_dist_label])
    sel.data.loc[_selected(sel),on_hull_label] = on_hull
  
  
  def _query(self):
    
    force = True
    
    # check which hull_selection configurations are calculated
    self.hull_sel.add_data('is_calculated')
    
    # store hull_selection.calculated
    self.hull_sel_calculated.save(
      data=self.hull_sel.data[self.hull_sel.data.loc[:,'is_calculated']==True],
      force=force)
    
    
    # for selection query:
    # - to plot: x_all, y, y_clex
    # - for tooltips / dft hull: configname, x_all, hull_dist, clex_hull_dist
    self.sel.query(['is_calculated', self.x_all, self.y, self.y_clex,
      hull_dist(self.hull_sel_calculated.path),
      hull_dist_per_atom(self.hull_sel_calculated.path),
      clex_hull_dist(self.hull_sel.path),
      clex_hull_dist_per_atom(self.hull_sel.path)],
      force=force)
    
    # convert to numeric (NaN) for configurations without DFT calculations 
    self.sel.data.loc[:,hull_dist(self.hull_sel_calculated.path)] = \
      pandas.to_numeric(self.sel.data[hull_dist(self.hull_sel_calculated.path)],
        errors='coerce')
    self.sel.data.loc[:,hull_dist_per_atom(self.hull_sel_calculated.path)] = \
      pandas.to_numeric(self.sel.data[hull_dist_per_atom(self.hull_sel_calculated.path)],
        errors='coerce')
    
    # find hull points
    self._set_on_hull(self.sel, 'on_dft_hull', hull_dist_per_atom(self.hull_sel_calculated.path))
    self._set_on_hull(self.sel, 'on_clex_hull', clex_hull_dist_per_atom(self.hull_sel.path))
    
    # create dft hull selection
    self.dft_hull_sel.save(data=self.sel.data[self.sel.data.loc[:,'on_dft_hull']==True], force=force)
    
    # for selection query:
    #  - for tooltips: clex_of_dft_hull_dist
    self.sel.query([
      clex_hull_dist(self.dft_hull_sel.path),
      clex_hull_dist_per_atom(self.dft_hull_sel.path)],
      force=force)
    
    # find clex_of_dft_hull line
    self._set_on_hull(self.sel, 'on_clex_of_dft_hull', clex_hull_dist_per_atom(self.dft_hull_sel.path))
    
    
    # source for points, tooltips
    for col in self.sel.data.columns:
       add_src_data(self.sel, col, self.sel.data.loc[:,col])
    
    # source for hull data
    self.calc_hull_src = bokeh.models.ColumnDataSource(data=self.sort_dft_hull_line())
    
    # source for clex hull data
    self.clex_hull_src = bokeh.models.ColumnDataSource(data=self.sort_clex_hull_line())
    
    # source for clex of dft hull data
    self.clex_of_dft_hull_src = bokeh.models.ColumnDataSource(data=self.sort_clex_of_dft_hull_line())
    
  
#  def update(self):
#    """
#    Enable update of y values
#    """
#    
#    # update y-data
#    add_hull_data(self.sel, self.hull_sel, self.dft_hull_sel, self.hull_tol, force=True)
#    add_src_data(self.sel, hull_dist(self.hull_sel.path), self.sel.data.loc[:,hull_dist(self.hull_sel.path)])
#    
#    self.update_hull_line()
#    self.update_clex_hull_line()
#    slef.update_clex_of_dft_hull_line()

  def _sort_hull_line(self, sel, on_hull_label):
    """
    Return x,y data where 'on_hull_label'==True, sorted by x 
    """
    # sort hull line data
    on_hull = sel.data.loc[:,on_hull_label] == True
    return sel.data.loc[on_hull, [self.x, self.y, self.y_clex]].sort_values([self.x])
  
  def sort_dft_hull_line(self):
    """ sort dft hull line data """
    self.sorted_dft_hull = self._sort_hull_line(self.sel, 'on_dft_hull')
    return self.sorted_dft_hull
  
  def sort_clex_hull_line(self):
    # sort clex hull line data
    self.sorted_clex_hull = self._sort_hull_line(self.sel, 'on_clex_hull')
    return self.sorted_clex_hull
  
  def sort_clex_of_dft_hull_line(self):
    # set clex of dft hull line data
    self.sorted_clex_of_dft_hull = self._sort_hull_line(self.sel, 'on_clex_of_dft_hull')
    return self.sorted_clex_of_dft_hull
  
#  def update_hull_line(self):
#    self.sort_hull_line()
#    
#    # source for hull data
#    self.calc_hull_src.data[self.x] = self.sel.data.loc[self.sorted_dft_hull.index, self.x]
#    self.calc_hull_src.data[self.y] = self.sel.data.loc[self.sorted_dft_hull.index, self.y]
#
#  def update_clex_hull_line(self):
#    self.sort_clex_hull_line()
#    
#    # source for clex hull data
#    self.clex_hull_src.data[self.x] = self.sel.data.loc[self.sorted_clex_hull.index, self.x]
#    self.clex_hull_src.data[self.y] = self.sel.data.loc[self.sorted_clex_hull.index, self.y]
#
#  def update_hull_line(self):
#    self.sort_clex_of_dft_hull_line()
#    
#    # source for clex of dft hull data
#    self.clex_of_dft_hull_src.data[self.x] = self.sel.data.loc[self.sorted_clex_of_dft_hull.index, self.x]
#    self.clex_of_dft_hull_src.data[self.y] = self.sel.data.loc[self.sorted_clex_of_dft_hull.index, self.y]


class Scatter(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
    y: the value of the y-axis, 'comp(a)' by default
  """
  
  def __init__(self, sel, x='comp(a)', y='comp(a)'):
    """
    Arguments:
      sel: A CASM Selection
      x: (str) property name to use as x axis, 'comp(a)' by default.
      y: (str) property name to use as y axis, 'comp(a)' by default.
    """
    
    self.sel = sel
    
    self.x = x
    self.sel.add_data('configname')
    add_src_data(self.sel, 'configname', self.sel.data['configname'])
    
    self.sel.add_data(x)
    add_src_data(self.sel, self.x, self.sel.data.loc[:,self.x])
    
    self.y = y
    self.sel.add_data(y)
    add_src_data(self.sel, self.y, self.sel.data.loc[:,self.y])
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    
    self._plot()
    
    # store plot in layout
    self.layout = self.p
  
  @property
  def p(self):
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'dft_style'):
      self.sel.dft_style = dft_style
      update_glyphs(self.sel)
    
    style = self.sel.dft_style
    
    # plot formation energy vs comp, with convex hull states
    if True:
      self.p_ = bokeh.plotting.Figure(plot_width=400, plot_height=400, tools=self.tools)
      
      cr = getattr(self.p_, style['marker'])(self.x, self.y, source=self.sel.src, 
        size='radii', fill_color='color', fill_alpha=style['fill_alpha'],
        line_color='line_color', line_width='line_width', line_alpha='line_alpha',
        hover_alpha=style['hover_alpha'], hover_color=style['hover_color'])
    
    else:
      # rug plot
      x_max = max(self.sel.src.data[self.x])
      x_min = min(self.sel.src.data[self.x])
      y_max = max(self.sel.src.data[self.y])
      y_min = min(self.sel.src.data[self.y])
      x_rug_length = (y_max - y_min)/20.0
      y_rug_length = (x_max - x_min)/20.0
      
      self.p_ = bokeh.plotting.Figure(
        plot_width=400, plot_height=400, tools=self.tools,
        x_range=(x_min-y_rug_length, x_max+2*y_rug_length), y_range=(y_min-x_rug_length, y_max+2*x_rug_length))
      
      cr = self.p_.circle(self.x, self.y, source=self.sel.src, 
        size='radii', fill_color='color', fill_alpha=style['fill_alpha'],
        line_color='line_color', line_width='line_width', line_alpha='line_alpha',
        hover_alpha=style['hover_alpha'], hover_color=style['hover_color'])
      
      x_rug = self.p_.segment(
        x0=self.x, x1=self.x,
        y0=y_max + x_rug_length, y1=y_max + 2.0*x_rug_length,
        source=self.sel.src,
        line_color='line_color', line_width=1.0)
      
      y_rug = self.p_.segment(
        x0=x_max + y_rug_length, x1=x_max + 2.0*y_rug_length,
        y0=self.y, y1=self.y,
        source=self.sel.src,
        line_color='line_color', line_width=1.0)
    
    
    self.p_.xaxis.axis_label = self.x
    self.p_.yaxis.axis_label = self.y

    # hover over a point to see 'configname', 'x', and 'y'
    tooltips = [
        ("configname","@configname"), 
        (self.x,"@{" + self.x + "}{1.1111}"),
        (self.y,"@{" + self.y + "}{1.1111}")
    ]
    self.tap_action = ConfigurationTapAction(self.sel)
    self.p_.add_tools(self.tap_action.tool(view_on_tap, [cr]))
    self.p_.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[cr]))
    self.p_.add_tools(bokeh.models.BoxSelectTool(renderers=[cr]))
    self.p_.add_tools(bokeh.models.LassoSelectTool(renderers=[cr]))


class Histogram(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
  """
  
  def __init__(self, sel, x='comp(a)', **kwargs):
    """
    Arguments:
      sel: A CASM Selection
      x: (str) property name to use as x axis, 'comp(a)' by default.
      kwargs: kwargs to pass to numpy.histogram
    """
    
    self.sel = sel
    
    self.x = x
    self.sel.add_data('configname')
    add_src_data(self.sel, 'configname', self.sel.data['configname'])
    self.sel.add_data(x)
    
    if 'to_be_selected' not in self.sel.data.columns:
      self.sel.data.loc[:,'to_be_selected'] = self.sel.data.loc[:,'selected']
    self.kwargs = kwargs
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    
    # for storing histogram freq/bins 
    self.src = None
    
    self._plot()
    
    # ensure that changes in 'sel.src' will update the histogram
    self.sel.selection_callbacks.append(self.update)
    
    # store plot in layout
    self.layout = self.p
  
  @property
  def p(self):
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'dft_style'):
      self.sel.dft_style = dft_style
      update_glyphs(self.sel)
    
    style = self.sel.dft_style
    
    self.p_ = bokeh.plotting.Figure(plot_width=400, plot_height=400, tools=self.tools)
    
    # set data in 'self.src'
    self.update()
    
    qd = self.p_.quad(top='top', bottom='bottom', left='left', right='right', source=self.src, 
      fill_alpha=0.0, line_color='black', line_width=1)
    
    qd_sel = self.p_.quad(top='top_sel', bottom='bottom_sel', left='left', right='right', source=self.src, 
      fill_alpha=style['fill_alpha'], fill_color=style['selected_color'],
      line_alpha=0.0, line_width=0.0)
        
    self.p_.xaxis.axis_label = self.x
    self.p_.yaxis.axis_label = self.y_label

    # hover over a point to see 'x', and 'Freq.'
    tooltips = [
        (self.x,"@{left}{1.1111}"),
        (self.y_label,"@{top}{1.1111}")
    ]
    
    self.p_.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[qd]))
  
  def update(self):
    """
    Update plot data if sel.data['selected'] or sel.data['x'] change.
    """
    
    # generate histogram and plot
    tmp_kwargs = dict(self.kwargs)
    tmp_kwargs['density'] = False
    hist, edges = np.histogram(self.sel.data.loc[:,self.x], **tmp_kwargs)
    
    sel_kwargs = dict(self.kwargs)
    sel_kwargs['bins'] = edges
    sel_kwargs['density'] = False
    
    _df = self.sel.data
    hist_sel, edges_sel = np.histogram(_df.loc[_df['to_be_selected']==True, self.x], **sel_kwargs)
    
    # need to do 'density' of selected by hand
    if 'density' in self.kwargs and self.kwargs['density']:
       hist_freq = np.array(hist)
       width = edges[1:] - edges[:-1]
       area = hist_freq.dot(width)
       hist = hist_freq*width/area
       hist_sel = hist_sel*width/area
    
    if self.src is None:
      self.src = bokeh.models.ColumnDataSource()
      self.src.data["left"] = edges[:-1]
      self.src.data["right"] = edges[1:]
    
    self.src.data["top"] = hist
    self.src.data["bottom"] = [0]*len(hist)
      
    self.src.data["top_sel"] = hist_sel
    self.src.data["bottom_sel"] = [0]*len(hist_sel)
          
    if 'density' not in self.kwargs.keys() or self.kwargs['density'] == False:
      self.y_label = "Frequency"
    else:
      self.y_label = "Density"


class GridPlot(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
  """
  
  def __init__(self, sel, 
    x=['comp(a)', 'formation_energy'], 
    y=['comp(a)', 'formation_energy'], 
    type=None,
    axis=0,
    kwargs=None):
    """
    Arguments:
      sel: A CASM Selection
      x: (list of str) property names to use as x axis
      y: (list of str) property names to use as x axis
      type = (list of class names, either Scatter or Histogram)
      axis = (0 or 1) axis that type refer to
      kwargs: (list of dict) kwargs to pass to constructors
    """
    
    self.sel = sel
    
    self.x = x
    self.y = y
    self.axis = axis
    if type is None:
      if axis == 1:
        type = [Scatter]*len(self.x)
      else:
        type = [Scatter]*len(self.y)
    self.type = type
    if kwargs is None:
      kwargs = [dict()]*len(self.type)
    self.kwargs = kwargs
    
    self._plot()
    
    # ensure that changes in 'sel.src' will update the histogram
    self.sel.selection_callbacks.append(self.update)
    
    # layout
    self.layout = bokeh.models.GridPlot(children=self.p_)
    
  
  @property
  def p(self):
    """A list of list of bokeh plots"""
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # construct plots
    self.casm_p_ = []
    self.p_ = []
    for i, _y in enumerate(self.y):
      casm_row = []
      row = []
      for j, _x in enumerate(self.x):
        if self.axis == 1:
          _index = j
        else:
          _index = i
          
        _type = self.type[_index]
        _kwargs = self.kwargs[_index]
        
        if _type is Histogram:
          casm_row.append(_type(self.sel, x=_x, **_kwargs))
          row.append(casm_row[-1].p)
        elif _type is Scatter:
          casm_row.append(_type(self.sel, x=_x, y=_y, **_kwargs))
          row.append(casm_row[-1].p)
        else:
          print "Unknown or unsupported plot type:", _type.name
          casm_row.append(None)
          row.append(None)
      
      self.casm_p_.append(casm_row)
      self.p_.append(row)
    
    # linked panning
    for r in range(len(self.y)):
      for c in range(len(self.x)):
        x_range = self.p_[0][c].x_range if (r != 0) else None
        y_range = self.p_[r][0].y_range if (c != 0) else None
        
        #self.casm_p_[r][c].layout = None
        if x_range is not None:
          self.p_[r][c].x_range = x_range
        if y_range is not None:
          self.p_[r][c].y_range = y_range
  
  def update(self):
    for row in self.casm_p_:
      for casm_p_ in row:
        if hasattr(casm_p_, 'update'):
          casm_p_.update()


class RankPlot(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    scoring: a scoring function to operate on a DataFrame and return a Series
    score_id: a UUID to use for the 'score' column name in the sel.src ColumnDataSource
    rank_id: a UUID to use for the 'rank' column name in the sel.src ColumnDataSource
    max_score: the maximum score
    min_score: the minimum score
  """
  
  def __init__(self, sel, scoring, name="Score", mode='set', tooltips=None, on_tap=view_on_tap):
    """
    Arguments:
      sel: A CASM Selection
      scoring: (function) a scoring function to operate on a row of sel.data
      name: y-axis label to describe the scoring function
      tooltips: Customize tooltips. Default includes configname, score, and rank.
      on_tap: Customize callback that occurs upon tapping a point. Default is view_on_tap.
    """
    self.sel = sel
    
    self.scoring = scoring
    self.name = name
    
    self.score_id = str(uuid.uuid4())
    self.rank_id = str(uuid.uuid4())
    
    self._max_score = None
    self._min_score = None
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    self.tooltips = tooltips
    self.on_tap = on_tap
    
    self.score()
    
    self._plot()
    
    # store plot in layout
    self.layout = self.p
    
  
  @property
  def p(self):
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'dft_style'):
      self.sel.dft_style = dft_style
      update_glyphs(self.sel)
    
    style = self.sel.dft_style
    
    _tools = [self.tools]
    
    self.p_ = bokeh.plotting.Figure(plot_width=800, plot_height=400, tools=_tools)
    p_stem = self.p_.segment(self.rank_id, 0, self.rank_id, self.score_id, source=self.sel.src, 
      line_color='color', line_width=1.0,
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'])
    
    p_circ = self.p_.circle(self.rank_id, self.score_id, source=self.sel.src, 
      size='radii', fill_color='color', fill_alpha=style['fill_alpha'],
      line_color='line_color', line_width='line_width', line_alpha='line_alpha',
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'])

    self.p_.xaxis.axis_label = "Rank"
    self.p_.yaxis.axis_label = self.name
    
    if self.tooltips is None:
      # hover over a point to see 'configname', 'score', Ef', and 'comp(a)'
      self.tooltips = [
          ("configname","@configname"), 
          ("score","@{" + self.score_id + "}{1.1111}"), 
          ("rank","@{" + self.rank_id + "}{1.1111}")
      ]
    
    if self.on_tap is not None:
      self.tap_action = ConfigurationTapAction(self.sel)
      self.p_.add_tools(self.tap_action.tool(self.on_tap, [p_circ]))
    
    self.p_.add_tools(bokeh.models.HoverTool(tooltips=self.tooltips, renderers=[p_circ]))
    self.p_.add_tools(bokeh.models.BoxSelectTool(renderers=[p_circ]))
    self.p_.add_tools(bokeh.models.LassoSelectTool(renderers=[p_circ]))
  
  @property
  def max_score(self):
    if self._max_score is None:
      self._max_score = max(self.sel.src.data[self.score_id])
    return self._max_score
  
  @property
  def min_score(self):
    if self._min_score is None:
      self._min_score = min(self.sel.src.data[self.score_id])
    return self._min_score
  
  def score(self, scoring=None):
    """
    Calculate 'score' and 'rank', but do not change selection. May also change the
    'scoring' function.
    """
    if scoring is not None:
      self.scoring = scoring
    
    try:
      # score
      self.df = pandas.DataFrame(self.scoring(self.sel), columns=['score'])
      self.df.loc[:,'selected'] = self.sel.data.loc[:,'selected']
    except:
      print "Error applying scoring function"
      raise
      
    # add 'score' to src, as self.score_id
    add_src_data(self.sel, self.score_id, self.df['score'], force=True)
    
    # rank
    self.df.loc[self.df.sort_values(['score']).index, 'rank'] = range(self.df.shape[0])
    
    # add 'rank' to src as self.rank_id
    add_src_data(self.sel, self.rank_id, self.df['rank'], force=True)
    

class RankSelect(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    scoring: a scoring function to operate on a row of sel.data
    score_id: a UUID to use for the 'score' column name in the sel.src ColumnDataSource
    rank_id: a UUID to use for the 'rank' column name in the sel.src ColumnDataSource
    max_score: the maximum score
    min_score: the minimum score
  """
  
  def __init__(self, sel, scoring, name="Score", mode='set'):
    """
    Arguments:
      sel: A CASM Selection
      scoring: (function) a scoring function to operate on a row of sel.data
      name: y-axis label to describe the scoring function
      mode: (str) 'set', 'set_on', 'set_off'
    """
    self.sel = sel
    
    self.scoring = scoring
    self.name = name
    
    if mode not in ['set', 'set_on', 'set_off']:
      raise ValueError("mode: " + mode + " not allowed. Expected 'set', 'set_on', or 'set_off'")
    self.mode = mode
    
    self.score_id = str(uuid.uuid4())
    self.rank_id = str(uuid.uuid4())
    
    self._max_score = None
    self._min_score = None
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    self.tooltips = None
    
    # bokeh DataSource for cutoff line
    self.cutoff_line_src = None
    
    # construct pre-score widgets
    self.input = dict()
    self._pre_score_input()
    
    self.score()
    
    # construct post-score widgets
    self._post_score_input()
    
    self._plot()
    
    # store plot in layout with widgets
    _hplot = hplot(self.input['cutoff'].widget, self.select_mode, self.select_action)
    self.layout = vplot(self.p, _hplot, self.msg.widget, width=self.p.plot_width)
    
  
  @property
  def p(self):
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'dft_style'):
      self.sel.dft_style = dft_style
      update_glyphs(self.sel)
    
    style = self.sel.dft_style
    
    _tools = [self.mouse, self.tools]
    
    self.p_ = bokeh.plotting.Figure(plot_width=800, plot_height=400, tools=_tools)
    p_stem = self.p_.segment(self.rank_id, 0, self.rank_id, self.score_id, source=self.sel.src, 
      line_color='color', line_width=1.0,
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'])
    
    p_circ = self.p_.circle(self.rank_id, self.score_id, source=self.sel.src, 
      size='radii', fill_color='color', fill_alpha=style['fill_alpha'],
      line_color='line_color', line_width='line_width', line_alpha='line_alpha',
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'])

    self.p_.xaxis.axis_label = "Rank"
    self.p_.yaxis.axis_label = "Score"

    if self.tooltips is None:
      # hover over a point to see 'configname', 'score', Ef', and 'comp(a)'
      tooltips = [
          ("configname","@configname"), 
          ("score","@{" + self.score_id + "}{1.1111}"), 
          ("rank","@{" + self.rank_id + "}{1.1111}")
      ]
    
    if self.cutoff_line_src is None:
      self.cutoff_line_src = bokeh.models.ColumnDataSource(
        data={
          'x':[0, self.sel.data.shape[0]], 
          'y':[self.input['cutoff'].value, self.input['cutoff'].value]
        })
      
      def f(attrname, old, new):
        self.input['cutoff'].value = new['y'][0]
      
      self.cutoff_line_src.on_change('data', f)
    
    
    self.p_.select(type=bokeh.models.CrosshairTool).dimensions = ['width']
    
    p_click=bokeh.models.CustomJS(args={'src':self.cutoff_line_src, 'loc': self.mouse_location_src}, code="""
      var data = src.get('data');
      data['y'][0] = loc.get('data')['y'][0];
      data['y'][1] = loc.get('data')['y'][0];
      src.set('data', data);
      """)
    self.p_.select(type=bokeh.models.CrosshairTool).dimensions = ['width']
    self.p_.add_tools(bokeh.models.TapTool(callback=p_click))
    
    p_line = self.p_.line(x='x',y='y', source=self.cutoff_line_src,
      line_width=style['cutoff_line_width'], line_alpha=style['cutoff_alpha'], color=style['cutoff_color'])

    self.p_.add_tools(bokeh.models.HoverTool(tooltips=self.tooltips, renderers=[p_circ]))
    self.p_.add_tools(bokeh.models.BoxSelectTool(renderers=[p_circ]))
  
  @property
  def max_score(self):
    if self._max_score is None:
      self._max_score = max(self.sel.src.data[self.score_id])
    return self._max_score
  
  @property
  def min_score(self):
    if self._min_score is None:
      self._min_score = min(self.sel.src.data[self.score_id])
    return self._min_score
  
  def score(self, scoring=None):
    """
    Calculate 'score' and 'rank', but do not change selection. May also change the
    'scoring' function.
    """
    if scoring is not None:
      self.scoring = scoring
    
    try:
      # score
      self.df = pandas.DataFrame(self.scoring(self.sel), columns=['score'])
      self.df.loc[:,'selected'] = self.sel.data.loc[:,'selected']
    except:
      print "Error applying scoring function"
      raise
      
    # add 'score' to src, as self.score_id
    add_src_data(self.sel, self.score_id, self.df['score'], force=True)
    
    # rank
    self.df.loc[self.df.sort_values(['score']).index, 'rank'] = range(self.df.shape[0])
    
    # add 'rank' to src as self.rank_id
    add_src_data(self.sel, self.rank_id, self.df['rank'], force=True)
  
  def _pre_score_input(self):
    """
    Construct gui widgets
    """
    ### mouse location
    self.mouse_location_src = bokeh.models.ColumnDataSource(data={'x':[0.0], 'y':[0.0]})
    self.mouse = bokeh.models.HoverTool(
      tooltips=None,
      callback=bokeh.models.CustomJS(args={'src': self.mouse_location_src}, code="""
        src.get('data')['y'][0] = cb_data.geometry.y;
        src.get('data')['y'][1] = cb_data.geometry.y;
        """))
    
    self.msg = Messages()
    
    self.select_mode = bokeh.models.Select(value="set", options=["set", "set_on", "set_off"])
    def select_mode_f(attrname, old, new):
      self.mode = new
      self.update_cutoff(None, None, None)
    self.select_mode.on_change('value', select_mode_f)
    
    self.select_action = bokeh.models.Select(
      value="Select action", 
      options=["Select action", "Revert", "Apply", "Apply and Save"])
    def select_action_f(attrname, old, new):
      if new == "Revert":
        self.sel.data.loc[:,'to_be_selected'] = self.sel.data.loc[:,'selected']
        self.sel.update_selection()
        update_glyphs(self.sel, selected=self.sel.data['to_be_selected'])
        self.msg.value = "Reverted to last applied selection"
      elif new == "Apply":
        self.sel.data.loc[:,'selected'] = self.sel.data.loc[:,'to_be_selected']
        self.sel.update_selection()
        update_glyphs(self.sel, selected=self.sel.data['to_be_selected'])
        self.msg.value = "Applied selection"
      elif new == "Apply and Save":
        try:
          self.sel.data.loc[:,'selected'] = self.sel.data.loc[:,'to_be_selected']
          self.sel.save(force=True)
          self.sel.update_selection()
          update_glyphs(self.sel, selected=self.sel.data['to_be_selected'])
          self.msg.value = "Applied selection and saved: " + self.sel.path
        except Exception, e:
          self.msg.value = str(e)
      self.select_action.value = "Select action"
    self.select_action.on_change('value', select_action_f)
  
  def _post_score_input(self):
    """
    Construct gui widgets
    """
    self.input['cutoff'] = FloatInput(self.max_score, title="Cutoff:")
    self.input['cutoff'].update.append(self.update_cutoff)
  
  def update_cutoff(self, attrname, old, new):
    """
    Update self.sel.data['selected'], update_glyphs(self.sel), update cutoff line.
    Ignores 'attrname', 'old', and 'new'.
    """
    
    f = None
    if self.mode == 'set':
      f = self._set_f
    elif self.mode == 'set_on':
      f = self._set_on_f
    elif self.mode == 'set_off':
      f = self._set_off_f
    
    self.sel.data.loc[:,'to_be_selected'] = self.df.apply(f, axis='columns')
    self.cutoff_line_src.data['y'] = [self.input['cutoff'].value, self.input['cutoff'].value]
    self.sel.update_selection()
    update_glyphs(self.sel, selected=self.sel.data['to_be_selected'])
  
  
  def _set_f(self, r):
    return r['score'] < self.input['cutoff'].value
    
  def _set_on_f(self, r):
    if r['score'] < self.input['cutoff'].value:
      return True
    else:
      return r['selected']
  
  def _set_off_f(self, r):
    if r['score'] < self.input['cutoff'].value:
      return False
    else:
      return r['selected']


class WeightSelect(object):
  
  def __init__(self, sel, input_filename="fit_input.json",
    hullplot_kwargs=None, A_params=None, B_params=None, kT_params=None, Eref_params=None,
    show_weightplot=True):
    """
    Interactively set weights, via wHullDist.
    
    Weights for selected configurations are applied using: A*np.exp(-hull_dist/kT) + B
    Weights for unselected configurations are set to 1.0
    
    Arguments:
      sel: A CASM Selection to use for the training set
      input_filename: The name of a fitting input file. Defaults to "fit_input.json". 
        If input_filename does not exist a default with that name is created.
      hullplot_kwargs: (dict) kwargs to pass to convex hull plots of formation energy
        and weighted formation energy
      A_params, B_params, kT_params: (dict) Settings for input sliders. Options are:
        'title', 'value', 'start', 'end', and 'step'. Units on kT are meV.
      show_weightplot: (boolean) Show a plot of sample weight values
    """
    self.sel = sel
    self.input_filename = input_filename
    self.hullplot_kwargs = hullplot_kwargs
    if not os.path.exists(input_filename):
      self.fit_input = casm.learn.example_input()
      with open(self.input_filename, 'w') as f:
        json.dump(self.fit_input, f, indent=2)
    else:
      self.fit_input = json.load(open(self.input_filename, 'r'))
    
    # create select box with weighting method options
    self.select_method = SelectInput(
      value="wHullDist", 
      options=["wHullDist", "wEmin", "wEref"], 
      title="Select Method:")
    for method in self.select_method.options:
      self.select_method.update[method] = self._update_method
    
    # create sliders for input of parameters, first set ranges
    weight_kwargs = self.fit_input["weight"].get("kwargs", dict())
    if A_params is None:
      A_params = {"title":"A", "value":0.0, "start":0.0, "end":10.0}
      if 'A' in weight_kwargs:
        A_params['value'] = weight_kwargs['A']
    self.A_params = A_params
    
    if B_params is None:
      B_params = {"title":"B", "value":1.0, "start":0.0, "end":10.0}
      if 'B' in weight_kwargs:
        B_params['value'] = weight_kwargs['B']
    self.B_params = B_params
    
    if kT_params is None:
      kT_params = {"title":"kT (meV)", "value":10.0, "start":1.0, "end":100.0, "step":1.0}
      if 'kT' in weight_kwargs:
        kT_params['value'] = weight_kwargs['kT']
    self.kT_params = kT_params
    
    self.Ef_min = self.sel.data['formation_energy'].min()
    self.Ef_max = self.sel.data['formation_energy'].max()
    delta = self.Ef_max - self.Ef_min
    if Eref_params is None:
      Eref_params = {
        "title":"Eref (meV)", 
        "value":self.Ef_min, 
        "start":(self.Ef_min - delta*1.0), 
        "end":self.Ef_max, 
        "step":0.001}
      if 'Eref' in weight_kwargs:
        Eref_params['value'] = weight_kwargs['Eref']
    self.Eref_params = Eref_params
    
    # create input widgets
    self.input = dict()
    self.input['A'] = bokeh.models.Slider(**self.A_params)
    self.input['B'] = bokeh.models.Slider(**self.B_params)
    self.input['kT'] = bokeh.models.Slider(**self.kT_params)
    self.input['Eref'] = bokeh.models.Slider(**self.Eref_params)
    
    for key, widget in self.input.iteritems():
      widget.on_change('value', self._update_wvalue)
    
    # create select action input
    self.msg = Messages()
    self.select_action = bokeh.models.Select(
      value="Select action", 
      options=["Select action", "Revert", "Apply", "Apply and Save"])
    def select_action_f(attrname, old, new):
      if new == "Revert":
        self._revert()
      elif new == "Apply":
        self._apply()
      elif new == "Apply and Save":
        self._apply_and_save()
      self.select_action.value = "Select action"
    self.select_action.on_change('value', select_action_f)
    
    # callbacks
    self.sel.selection_callbacks.append(self._update_selection)
    
    # add weighted property value to self.sel.src
    self.wvalue_id = str(uuid.uuid4())
    
    # add convex hull plot
    self.hullplot = None
    self.sel.data.loc[:,'weight'] = np.zeros((self.sel.data.shape[0],1))
    self.show_weightplot = show_weightplot
    self.weightplot = None
    self.whullplot = None
    self.layout = None
    self.ref_line = None
    self.plot_width = 400
    self._update_selection()
    self.msg.value = "..."
    
  
  def _apply(self):
    self.A_params['value'] = self.input['A'].value
    self.B_params['value'] = self.input['B'].value
    self.kT_params['value'] = self.input['kT'].value
    self.Eref_params['value'] = self.input['Eref'].value
    self.msg.value = "Applied sample weights"
  
  def _revert(self):
    self.input['A'].value = self.A_params['value']
    self.input['B'].value = self.B_params['value']
    self.input['kT'].value = self.kT_params['value']
    self.input['Eref'].value = self.Eref_params['value']
    self.msg.value = "Reverted to last applied sample weights"
  
  def _save(self):
    try:
      self.fit_input["weight"]["method"] = self.select_method.value
      self.fit_input["weight"]["kwargs"]["A"] = self.input['A'].value
      self.fit_input["weight"]["kwargs"]["B"] = self.input['B'].value
      self.fit_input["weight"]["kwargs"]["kT"] = self.input['kT'].value
      if self.select_method.value == ["wEref"]:
        self.fit_input["weight"]["kwargs"]["Eref"] = self.input['Eref'].value
      with open(self.input_filename, 'w') as f:
        json.dump(self.fit_input, f, indent=2)
      self.msg.value = "Saved input file: " + self.input_filename
    except Exception, e:
      self.msg.value = str(e)
  
  def _apply_and_save(self):
    try:
      self._apply()
      self._save()
      self.msg.value = "Applied sample weights and saved input file: " + self.input_filename
    except Exception, e:
      self.msg.value = str(e)
  
  def _set_value(self):
    """Set wvalue based on weighting parameters"""
    self.sel.data.loc[:,self.wvalue_id] = self.sel.data.loc[:,'formation_energy']
    add_src_data(self.sel, self.wvalue_id, self.sel.data.loc[:,self.wvalue_id], force=True)
  
  def _set_wvalue(self):
    """Set wvalue based on weighting parameters"""
    # calculate weights
    if self.select_method.value == "wHullDist":
      w = casm.learn.tools.wHullDist(
            self.hull_dist_values[np.where(self._selected)], 
            self.input['A'].value, 
            self.input['B'].value, 
            self.input['kT'].value*0.001)
    elif self.select_method.value == "wEmin":
      w = casm.learn.tools.wEmin(
            self.sel.data['formation_energy'].values[np.where(self._selected)], 
            self.input['A'].value, 
            self.input['B'].value, 
            self.input['kT'].value*0.001)
    elif self.select_method.value == "wEref":
      w = casm.learn.tools.wEref(
            self.sel.data['formation_energy'].values[np.where(self._selected)], 
            self.input['A'].value, 
            self.input['B'].value, 
            self.input['kT'].value*0.001,
            self.input['Eref'].value)
    
    self.sel.data.loc[:,'weight'] = 0.0
    self.sel.data.loc[self._selected,'weight'] = w
    add_src_data(self.sel, 'weight', self.sel.data.loc[:,'weight'], force=True)
    
    # update data (will update hullplot scatter points, but not convex hull line)
    self.sel.data.loc[self._selected,self.wvalue_id] = casm.learn.tools.set_sample_weight(
      w, y=self.sel.data.loc[self._selected,'formation_energy'].values)[0]
    self.sel.data.loc[self._unselected,self.wvalue_id] = self.sel.data.loc[self._unselected,'formation_energy']
    add_src_data(self.sel, self.wvalue_id, self.sel.data.loc[:,self.wvalue_id], force=True)
  
  def _update_wvalue(self, attrname, old, new):
    """Update due to changes in weighting paramaters, but no change in selection"""
    self._set_wvalue()
    
    if self.ref_line is None:
      if self.select_method.value in ["wEmin", "wEref"] :
        self.ref_line_src = bokeh.models.ColumnDataSource(
          data={'x':[self.x_min, self.x_max], 'y':[self.input['Eref'].value, self.input['Eref'].value]})
        self.ref_line = self.hullplot.p.line('x', 'y', source=self.ref_line_src, line_color='red', line_width=1.0)
    if self.ref_line is not None:
      if self.select_method.value == "wEref":
        self.ref_line_src.data['x'] = [self.x_min, self.x_max]
        self.ref_line_src.data['y'] = [self.input['Eref'].value, self.input['Eref'].value]
      elif self.select_method.value == "wEmin":
        self.ref_line_src.data['x'] = [self.x_min, self.x_max]
        self.ref_line_src.data['y'] = [self.Ef_min, self.Ef_min]
      elif self.select_method.value == "wHullDist":
        self.ref_line_src.data['x'] = []
        self.ref_line_src.data['y'] = []
    if self.select_method.value == "wEref":
      self.ref_line_src.data['y'] = [self.input['Eref'].value, self.input['Eref'].value]
    
    self.whullplot.update_hull_line()
  
  def _update_selection(self):
    """Update to reflect change in selection"""
    
    # use the unweighted values to initially plot the hull
    self._set_value()
    
    if self.hullplot_kwargs is None:
      self.hullplot_kwargs=dict()
    hullplot_kwargs = self.hullplot_kwargs
    
    # add convex hull plot, using unweighted formation energies
    if self.hullplot is None:
      hullplot_kwargs['y'] = 'formation_energy'
      self.hullplot = ConvexHullPlot(self.sel, **self.hullplot_kwargs)
      self.hullplot.p.plot_width = self.plot_width
      self.hullplot.p.yaxis.axis_label = "Formation Energy"
      self.x_min = self.sel.data[self.hullplot.x].min()
      self.x_max = self.sel.data[self.hullplot.x].max()
    else:
      self.whullplot.update()
    
    if self.show_weightplot and self.weightplot is None:
      self.weightplot = Scatter(self.sel, self.hullplot.x, 'weight')
      self.weightplot.p.plot_width = self.plot_width
      self.weightplot.p.yaxis.axis_label = "Weight"
    
    # add convex hull plot, using unweighted formation energies
    if self.whullplot is None:
      hullplot_kwargs['y'] = self.wvalue_id
      self.whullplot = ConvexHullPlot(self.sel, **self.hullplot_kwargs)
      self.whullplot.p.plot_width = self.plot_width
      self.whullplot.p.yaxis.axis_label = "Weighted Formation Energy"
    else:
      self.whullplot.update()
    
    # set weights and update the hullplot
    self.hull_dist_values = self.sel.data.loc[:, hull_dist(self.sel.path)].values
    self._selected = _selected(self.sel)
    self._unselected = _unselected(self.sel)
    self._update_wvalue(None, None, None)
    self.msg.value = "Save selection to re-calculate the convex hull"
    
    # Set up layouts and add to document
    if self.layout is None:
      self.layout = self._layout()
  
  def _update_method(self, attrname, old, new):
    self.layout.children[0].children[0].children = self._widget_layout()
    self._update_wvalue(None, None, None)
  
  def _layout(self):
    if self.show_weightplot:
      children = [self.hullplot.layout, self.whullplot.layout, self.weightplot.layout]
    else:
      children = [self.hullplot.layout, self.whullplot.layout]
    grid = bokeh.models.GridPlot(children=[children])
    row = hplot(vplot(*self._widget_layout()), grid)
    return vplot(row, self.msg.widget, width=self.whullplot.p.plot_width*3 + 300)
  
  def _widget_layout(self):
    if self.select_method.value in ["wHullDist", "wEmin"]:
      return [self.select_method.widget,
              self.input['kT'], 
              self.input['A'], 
              self.input['B'], 
              self.select_action]
    else:
      return [self.select_method.widget,
              self.input['kT'], 
              self.input['A'], 
              self.input['B'],
              self.input['Eref'], 
              self.select_action]


class ECISelection(object):
  def __init__(self, proj, input_filename="fit_input.json", cwd=None):
    """
    Arguments
    ---------
      proj: a casm.Project instance
      
      input_filename: str, optional, default "fit_input.json"
        The name of the fitting input file
      
      cwd: str, optional, default current working directory
        The directory containing input_filename
    """
    if cwd is None:
      cwd = os.getcwd()
    self.input_filename = os.path.join(cwd, input_filename,)
    self.fit_input = json.load(open(self.input_filename, 'r'))
    halloffame_filename = os.path.join(cwd, self.fit_input.get("halloffame_filename", "halloffame.pkl"))
    self.hall = pickle.load(open(halloffame_filename, 'rb'))
    
    data = []
    for index, indiv in enumerate(self.hall):
      d = {
        "selected":True, 
        "index":index, 
        "Nbfunc":sum(indiv),
        "cv":indiv.fitness.values[0], 
        "rms":indiv.rms, 
        "wrms:":indiv.wrms,
        "estimator_method":indiv.estimator_method,
        "feature_selection_method":indiv.feature_selection_method,
        "note":indiv.note}
      data.append(d)
    self.data = pandas.DataFrame(data)
    
    # shape: (#bfunc, #indiv)
    eci = np.zeros((len(indiv), len(self.hall)))
    for index, indiv in enumerate(self.hall):
      for bfunc in indiv.eci:
        eci[bfunc[0], index] = bfunc[1]
    self.eci = pandas.DataFrame(eci, columns=[str(i) for i in xrange(len(self.hall))])
    
    self.src = None



class ECISelect(object):
  
  def __init__(self, proj, input_filename="fit_input.json", cwd=None):
    """
    Arguments
    ---------
      proj: a casm.Project instance
      
      input_filename: str, optional, default "fit_input.json"
        The name of the fitting input file
      
      cwd: str, optional, default current working directory
        The directory containing input_filename
    """
    
    self.proj = proj
    
    if cwd is None:
      cwd = os.getcwd()
    self.sel = ECISelection(proj, input_filename, cwd)
    
    # rankplot: showing cv
    # scatter plot or median plot of all indiv
    # indiv plot: value vs bfunc index, colored by branch
    self._cv_rankplot()
    self._eci_stemplot(0)
    
    self.layout = hplot(self.cv_rankplot.layout, self.eci_stemplot)
    
  
  def _cv_rankplot(self):
    tooltips = [
        ("Index","@{index}"),
        ("CV score","@{cv}{1.1111}"),
        ("#Basis Functions", "@{Nbfunc}"),
        ("Estimator","@{estimator_method}"),
        ("Selection","@{feature_selection_method}"),
        ("Note","@{note}")
    ]
    add_src_data(self.sel, "index", self.sel.data["index"])
    add_src_data(self.sel, "cv", self.sel.data["cv"])
    add_src_data(self.sel, "estimator_method", self.sel.data["estimator_method"])
    add_src_data(self.sel, "feature_selection_method", self.sel.data["feature_selection_method"])
    add_src_data(self.sel, "note", self.sel.data["note"])
    add_src_data(self.sel, "Nbfunc", self.sel.data["Nbfunc"])
    
    def score(sel):
      return sel.data.loc[:,"cv"].values
    
    def _update_eci_stemplot(sel, attrname, old, new):
      """
      Change viewed individual on tap
      """
      print "tap!"
      if len(sel.src.selected['1d']['indices']) == 1:
        print "update!"
        index = sel.src.selected['1d']['indices'][0]
        sel.eci_src.data["current"] = sel.eci_src.data[str(index)]
    
    # self, sel, scoring, name="Score", mode='set', tooltips=None, on_tap=view_on_tap)
    self.cv_rankplot = RankPlot( self.sel, score, name="CV Score", 
      tooltips=tooltips, on_tap=_update_eci_stemplot)
    
    self.cv_rankplot.p.xaxis.axis_label = "Individual"
    self.cv_rankplot.p.plot_width = 400
  
  def _eci_stemplot(self, index):
    
    self.sel.eci_src = bokeh.models.ColumnDataSource(data=self.sel.eci)
    
    with open(self.proj.dir.basis(self.proj.settings.default_clex), 'r') as f:
      self.basis = json.load(f)
    
    self.sel.eci_src.data["index"] = range(self.sel.eci.shape[0])
    self.sel.eci_src.data["branch"] = [j["orbit"][0] for j in self.basis["cluster_functions"] ]
    self.sel.eci_src.data["mult"] = [j["mult"] for j in self.basis["cluster_functions"] ]
    self.sel.eci_src.data["max_length"] = [j["prototype"]["max_length"] for j in self.basis["cluster_functions"] ]
    self.sel.eci_src.data["min_length"] = [j["prototype"]["min_length"] for j in self.basis["cluster_functions"] ]
    colormap = ["black", "blue", "red", "green", "cyan", "magenta", "yellow", "orange"]
    self.sel.eci_src.data["color"] = map(lambda x: colormap[x % len(colormap)], self.sel.eci_src.data["branch"])
    self.sel.eci_src.data["current"] = self.sel.eci_src.data[str(index)]
    
    # TOOLS in bokeh plot
    _tools = ["crosshair,pan,reset,resize,box_zoom"]
    
    self.eci_stemplot = bokeh.plotting.Figure(plot_width=800, plot_height=400, 
      tools=_tools, y_range=(self.sel.eci.min().min()*1.1, self.sel.eci.max().max()*1.1))
    p_stem = self.eci_stemplot.segment("index", 0, "index", "current", source=self.sel.eci_src, 
      line_color='color', line_width=1.0)
    
    p_circ = self.eci_stemplot.circle("index", "current", source=self.sel.eci_src, 
      size=10, fill_color='color')

    self.eci_stemplot.xaxis.axis_label = "Basis Function Index"
    self.eci_stemplot.yaxis.axis_label = "ECI Value"
    
    # hover over a point to see 'ECI', other info on basis function/cluster
    tooltips = [
        ("Index","@{index}"), 
        ("ECI","@{current}{1.1111}"), 
        ("Cluster size","@{branch}"),
        ("Multiplicity", "@{mult}"),
        ("Max length", "@{max_length}"),
        ("Min length", "@{min_length}")
    ]
    
    self.eci_stemplot.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[p_circ]))
    self.eci_stemplot.add_tools(bokeh.models.BoxSelectTool(renderers=[p_circ]))
    self.eci_stemplot.add_tools(bokeh.models.LassoSelectTool(renderers=[p_circ]))
    

class FloatInput(object):
  """
  Input a float, and run actions on change of value.
  
  Attributes:
    self.widget: a bokeh.models.TextInput containing the value
    self.update: A list of callables that should be called on changes of the cutoff value.
      Signature should be f(attrname, old, new), where 'old' and 'new' are float.
  """
  
  def __init__(self, value="0.0", title="Float input:"):
    """
    Arguments:
      value: Value to use for cutoff.
    """
    
    self.widget = bokeh.models.TextInput(value=str(value), title=title)
    self.update = []
    self.widget.on_change('value', self.on_change)
  
  def on_change(self, attrname, old, new):
    for f in self.update:
      f(attrname, float(old), float(new))
  
  @property
  def value(self):
    return float(self.widget.value)
  
  @value.setter
  def value(self, _value):
    self.widget.value = str(_value)


class SelectInput(object):
  """
  A dropdown selection menu.
  
  Attributes:
    self.widget: a bokeh.models.Select containing the cutoff value
    self.update: a dict of option name str -> list of callables that should be 
      called on change of value
  """
  
  def __init__(self, value=None, options=[], title="Select:"):
    """
    Arguments:
      value: (str) initial value
      options: (list of str) value options
      title: display name
    """
    
    self.widget = bokeh.models.Select(value=value, title=title, options=options)
    self.update = dict()
    self.widget.on_change('value', self.on_change)
  
  def on_change(self, attrname, old, new):
    if new in self.update.keys():
      self.update[new](attrname, old, new)
  
  @property
  def value(self):
    return self.widget.value
  
  @value.setter
  def value(self, _value):
    self.widget.value = _value
  
  @property
  def options(self):
    return self.widget.options
  
  @options.setter
  def options(self, _value):
    self.widget.options = _value


class StrInput(object):
  """
  Input a str, and run actions on change of value.
  
  Attributes:
    self.widget: a bokeh.models.TextInput containing the value
    self.update: A list of callables that should be called on changes of the cutoff value.
      Signature should be f(attrname, old, new), where 'old' and 'new' are str.
  """
  
  def __init__(self, value="text here", title="Str input:"):
    """
    Arguments:
      value: Value to use for cutoff.
    """
    
    self.widget = bokeh.models.TextInput(value=value, title=title)
    self.update = []
    self.widget.on_change('value', self.on_change)
  
  def on_change(self, attrname, old, new):
    for f in self.update:
      f(attrname, old, new)
  
  @property
  def value(self):
    return self.widget.value
  
  @value.setter
  def value(self, _value):
    self.widget.value = _value


class GUIOptions(object):
  """
  Attributes:
    self.widget: a bokeh.models.TextInput containing the value
    self.layout: layout to include in page
    self.update: A list of callables that should be called on changes of the cutoff value.
      Signature should be f(attrname, old, new), where 'old' and 'new' are str.
    self.session: a casm.plotting.Session
  """
  
  def __init__(self, session):
    """
    Arguments:
      session: the bokeh push_session
    """
    self.session = session
    self.widget = bokeh.models.Select(
      value="Options", 
      options=["Options", "Close (without saving)"])
    self.update = []
    
    def select_action_f(attrname, old, new):
      if new == "Close (without saving)":
        self.widget.options = ["Session is closed"]
        self.widget.value = "Session is closed"
        session.close()
    
    self.update.append(select_action_f)
    
    self.widget.on_change('value', self.on_change)
    self.layout = hplot(self.widget, width=200)
  
  def on_change(self, attrname, old, new):
    for f in self.update:
      f(attrname, old, new)
  
  @property
  def value(self):
    return self.widget.value
  
  @value.setter
  def value(self, _value):
    self.widget.value = _value
 
 

