from __future__ import (absolute_import, division, print_function, unicode_literals)
from builtins import *

# conda's current version of pandas raises these warnings, but they are safe
# see: https://stackoverflow.com/questions/40845304
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import imp
import copy
import json
import os
import pickle
import re
import uuid

import bokeh.client
import bokeh.io
import bokeh.models
import bokeh.plotting
import bokeh.layouts as bk_layouts
import numpy as np
import pandas
import six

import casm
import casm.project
import casm.learn
from casm.misc import compat

int_dtypes = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
float_dtypes = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
numerics = int_dtypes + float_dtypes

class Session(object):
  def __init__(self, doc=None):
    if doc is None:
      doc = bokeh.io.curdoc()
    self.doc = doc
    self.session = bokeh.client.push_session(self.doc)

  def add_root(self, layout):
    self.doc.add_root(layout)

  def begin(self):
    self.session.show(self.doc)

  def begin_interactive(self):
    self.session.show(self.doc)
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

def set_src_data(src, name, data, force=False):
  """Add/modify data column in a ColumnDataSource, creating the column if necessary"""
  if name not in src.column_names:
    src.add(data, name)
  elif force == True:
    src.patch({name:(slice(len(data)),data)})

def add_src_data(sel, name, data, force=False):
  """Add data to Selection's ColumnDataSource, creating if necessary"""
  if sel.src is None:
    sel.src = bokeh.models.ColumnDataSource(data={name:data})
    sel.selection_callbacks = []
    def callbacks():
      for f in sel.selection_callbacks:
        f()
    sel.update_selection = callbacks
  else:
    set_src_data(sel.src, name, data, force=force)

def view_on_tap(sel, attrname, old, new):
  """
  Execute 'casm view' on tapped configuration
  """
  if len(sel.src.selected['1d']['indices']) == 1:
    index = sel.src.selected['1d']['indices'][0]
    configname = sel.src.data['configname'][index]
    args = "view " + configname
    sel.proj.capture(args=args)


class PlottingData(object):
    def __init__(self):
        self.data_ = dict()
        self.tap_action = TapAction(self)

    def add_project(self, proj_list):
        for proj in proj_list:
            if proj.path not in self.data_:
                self.data_[proj.path] = {'project':proj, 'selections':{}}
            else:
                self.data_[proj.path]['project'] = proj

    def add_selection(self, sel_list):
        for sel in sel_list:
            if sel.proj.path not in self.data_:
                self.data_[sel.proj.path] = {'project':sel.proj, 'selections':{sel.path:sel}}
            else:
                self.data_[sel.proj.path]['selections'][sel.path] = sel

    def project(self, proj_path=None, kwargs=dict()):
        if proj_path is None:
            proj_path = casm.project.project_path()
        else:
            proj_path = os.path.abspath(proj_path)
        if proj_path not in self.data_:
            proj = casm.project.Project(proj_path, **kwargs)
            self.add_project([proj])
        return self.data_[proj_path]['project']

    def selection(self, proj_path=None, sel_path="MASTER", kwargs=dict(), proj_kwargs=dict()):
        proj = self.project(proj_path, kwargs=proj_kwargs)
        if sel_path is None:
            sel_path = "MASTER"
        if sel_path not in ["MASTER", "ALL", "CALCULATED"]:
            sel_path = os.path.abspath(sel_path)
        if sel_path not in self.data_[proj.path]['selections']:
            sel = casm.project.Selection(self.project(proj.path), sel_path, **kwargs)
            self.add_selection([sel])
        return self.data_[proj.path]['selections'][sel_path]





class TapAction(object):
  """
  Enable execution of a python function when a glyph is 'tapped' in a bokeh plot.

  Example:
    p = ... a Bokeh plot
    circles = p.circ(... A set of circles)
    my_tap_action = TapAction(sel)
    def f(sel, attrname, old, new):
      if len(sel.src.selected['1d']['indices']) == 1:
        # the index of the tapped configuration:
        index = sel.src.selected['1d']['indices'][0]
        ... code to execute ...
        ... basically can ignore attrname, old, new ...

    # call f when a circle is tapped
    p.add_tools(my_tap_action.tool(f, [circles]))

  Attributes:

      data: PlottingData instance
          Use to store associated Project and Selection

      py_callback: {<sel_name>:<Python function>}
          Dict of callback functions with signature: f(Selection, attrname, old, new).
          Set the callback function using the 'tool' method. See 'view_on_tap'
          for an example.

  """

  def __init__(self, data):
    """
    Arguments:

        data: PlottingData instance
    """
    self.tap_indicator_src = bokeh.models.ColumnDataSource(data={"x":[0], "project":['None'], "selection":['None']})
    self.js_callback = bokeh.models.CustomJS(args={'src':self.tap_indicator_src}, code="""
      var data = src.data;
      data['x'][0] = data['x'][0] + 1;
      hovered_data = cb_data.source.data;
      data['project'][0] = hovered_data['project'][0];
      data['selection'][0] = hovered_data['selection'][0];
      src.data = data;
      """)
    self.data = data
    self.py_callback = dict()
    self.renderers = []

  def __call__(self, attrname, old, new):
    proj_path = new['project'][0]
    sel_path = new['selection'][0]
    sel = self.data.selection(proj_path, sel_path, proj_kwargs={"verbose":False})
    self.py_callback[sel_path](sel, attrname, old, new)

  def add_callback(self, sel, callback=None, renderers=None):
    self.py_callback[sel.path] = callback
    self.renderers += renderers

  def tool(self):
    self.tap_indicator_src.on_change('data', self)
    return bokeh.models.TapTool(callback=self.js_callback, renderers=self.renderers)


# convex hull plot

default_dft_hull_style = {
  "marker": "circle",
  "hover_color": "orange",
  "hover_alpha": 0.7,
  "fill_alpha": 0.5,
  "selected": {
    "color": "red",
    "radii": 10
  },
  "unselected": {
    "color": "gray",
    "radii": 7
  },
  "on_hull": {
    "line_color": "red",
    "line_width": 2.0,
    "line_alpha": 0.8
  },
  "off_hull": {
    "line_color": "blue",
    "line_width": 0.0,
    "line_alpha": 0.0
  },
  "hull_line_width": 1.0,
  "hull_line_dash": ""
}

default_clex_hull_style = {
  "marker": "x",
  "hover_color": "orange",
  "hover_alpha": 0.7,
  "fill_alpha": 0.5,
  "selected": {
    "color": "blue",
    "radii": 10
  },
  "unselected": {
    "color": "gray",
    "radii": 7
  },
  "on_hull": {
    "line_color": "blue",
    "line_width": 2.0,
    "line_alpha": 0.8
  },
  "off_hull": {
    "line_color": "blue",
    "line_width": 1.5,
    "line_alpha": 0.5
  },
  "hull_line_width": 1.0,
  "hull_line_dash": "",
  "clex_of_dft_hull_line_dash": "4 4"
}

def _hull_style(input, default):
    for val in ['marker', 'hover_color', 'hover_alpha', 'selected',
        'hull_line_width', 'hull_line_dash', 'fill_alpha',
        'clex_of_dft_hull_line_dash']:
        if val in default and val not in input:
            input[val] = default[val]
    for tmp in ['selected', 'unselected']:
        if tmp not in input:
            input[tmp] = dict()
        for val in ['color', 'radii']:
            if val not in input[tmp]:
                input[tmp][val] = default[tmp][val]
    for tmp in ['on_hull', 'off_hull']:
        if tmp not in input:
            input[tmp] = dict()
        for val in ['line_color', 'line_width', 'line_alpha']:
            if val not in input[tmp]:
                input[tmp][val] = default[tmp][val]
    return input

def dft_hull_style(input):
    return _hull_style(input, default_dft_hull_style)

def clex_hull_style(input):
    return _hull_style(input, default_clex_hull_style)


rank_select_cutoff_style = {
  "color": "red",
  "line_width": 1,
  "alpha": 1.0
}

default_scatter_style = {
  "marker": "circle",
  "hover_color": "orange",
  "hover_alpha": 0.7,
  "selected": {
    "color": "blue",
    "radii": 10,
    "line_color": "blue",
    "line_width": 0.0,
    "line_alpha": 0.0,
    "fill_alpha": 0.5
  },
  "unselected": {
    "color": "green",
    "radii": 7,
    "line_color": "gray",
    "line_width": 0.0,
    "line_alpha": 0.0,
    "fill_alpha": 0.3
  }
}

default_figure_kwargs = {
    "plot_height": 400,
    "plot_width": 800,
    "tools": "crosshair,pan,reset,box_zoom,wheel_zoom,save"
}


def scatter_series_style(index, input):
  selected_color = ["blue", "red", "green", "black", "cyan", "magenta", "yellow", "orange"]
  unselected_color = ["gray", "gray", "gray", "gray", "gray", "gray", "gray", "gray"]
  marker = ["circle", "triangle", "square", "diamond", "x", "cross"]

  for val in ['hover_color', 'hover_alpha']:
      if val not in input:
          input[val] = default_scatter_style[val]

  if 'marker' not in input:
      input['marker'] = marker[index % 6]

  if 'selected' not in input:
      input['selected'] = copy.deepcopy(default_scatter_style['selected'])
      input['selected']['color'] = selected_color[index % 8]
      input['selected']['line_color'] = selected_color[index % 8]
  else:
      for key in ['radii', 'line_width', 'line_alpha', 'fill_alpha']:
          if key not in input['selected']:
              input['selected'][key] = default_scatter_style['selected'][key]
      if 'color' not in input['selected']:
          input['selected']['color'] = selected_color[index % 8]
      if 'line_color' not in input['selected']:
          input['selected']['line_color'] = selected_color[index % 8]

  if 'unselected' not in input:
      input['unselected'] = copy.deepcopy(default_scatter_style['unselected'])
      input['unselected']['color'] = unselected_color[index % 8]
      input['unselected']['line_color'] = unselected_color[index % 8]
  else:
      for key in ['radii', 'line_width', 'line_alpha', 'fill_alpha']:
          if key not in input['unselected']:
              input['unselected'][key] = default_scatter_style['unselected'][key]
      if 'color' not in input['unselected']:
          input['unselected']['color'] = unselected_color[index % 8]
      if 'line_color' not in input['unselected']:
          input['unselected']['line_color'] = unselected_color[index % 8]

  return input

def rankplot_style(input):
    index = 0

    if 'selected' not in input:
      input['selected'] = dict()
    if 'line_width' not in input['selected']:
      input['selected']['line_width'] = 1.0
    if 'line_alpha' not in input['selected']:
      input['selected']['line_alpha'] = 1.0

    if 'unselected' not in input:
      input['unselected'] = dict()
    if 'line_width' not in input['unselected']:
      input['unselected']['line_width'] = 1.0
    if 'line_alpha' not in input['selected']:
      input['unselected']['line_alpha'] = 1.0

    return scatter_series_style(index, input)


def update_dft_hull_glyphs(sel, style, selected=None):
  """
  Use the selected and on_hull iterables to update glyph styles based on sel.dft_style.
  """
  if selected is None:
    selected = sel.data['selected']

  ## clex points:

  for value in ['color', 'radii']:
    cmap = dict({True:style['selected'][value], False:style['unselected'][value]})
    set_src_data(sel.src, value, list(map(lambda x: cmap[x], selected)), force=True)

  on_dft_hull = None
  if 'on_dft_hull' in sel.data.columns:
    on_dft_hull = sel.data['on_dft_hull']

  if on_dft_hull is not None:
    for value in ['line_color', 'line_width', 'line_alpha']:
      cmap = dict({True:style['on_hull'][value], False:style['off_hull'][value]})
      set_src_data(sel.src, value, list(map(lambda x: cmap[x], on_dft_hull)), force=True)


def update_clex_hull_glyphs(sel, style, selected=None):
  """
  Use the selected and on_hull iterables to update glyph styles based on sel.clex_style.
  """

  if selected is None:
    selected = sel.data['selected']

  ## clex points:

  for value in ['color', 'radii']:
    cmap = dict({True:style['selected'][value], False:style['unselected'][value]})
    set_src_data(sel.src, 'clex_' + value, list(map(lambda x: cmap[x], selected)), force=True)

  on_clex_hull = None
  if 'on_clex_hull' in sel.data.columns:
    on_clex_hull = sel.data['on_clex_hull']

  if on_clex_hull is not None:
    for value in ['line_color', 'line_width', 'line_alpha']:
      cmap = dict({True:style['on_hull'][value], False:style['off_hull'][value]})
      set_src_data(sel.src, 'clex_' + value, list(map(lambda x: cmap[x], on_clex_hull)), force=True)


def update_scatter_glyphs(sel, style, id, selected=None):

  if selected is None:
    selected = sel.data['selected']

  ## clex points:

  for value in ['color', 'radii', 'line_color', 'line_width', 'line_alpha', 'fill_alpha']:
    cmap = dict({True:style['selected'][value], False:style['unselected'][value]})
    set_src_data(sel.src, value + '.' + id, list(map(lambda x: cmap[x], selected)), force=True)


class ConvexHullPlot(object):
  """
  Attributes:
      data
      sel
      hull_sel
      hull_sel_calculated
      dft_hull_sel
      is_comp, is_comp_n, is_atom_frac
      x, x_all, y, y_clex
      hull_tol
      dft_style, clex_style
      tooltips, tooltips_exclude
      r_dft, r_dft_hull_line
      r_clex, r_clex_hull_line, r_clex_of_dft_hull_line
      renderers
      sorted_dft_hull
      sorted_clex_hull
      sorted_clex_of_dft_hull
  """

  def __init__(self, data=None, project=None, selection="MASTER", hull_selection=None,
    x='comp(a)', y='formation_energy', tooltips=[],
    dft_style={}, clex_style={}, hull_tol=1e-8, index=0, series_name=None, type=None):
    """
    Arguments:
      data: PlottingData object
      project: path to CASM project (optional, default=casm.project.project_path())
      selection: path for CASM Selection
      hull_selection: path for CASM Selection to use for the hull. Default uses selection.
      x: (str) column name to use as x axis. Default='comp(a)'. Must be 'comp(x)', 'comp_n(X)', or 'atom_frac(X)'.
      y: (str) column name to use as y axis. Default='formation_energy'. Must be 'formation_energy' or 'formation_energy_per_atom'.
      tooltips:
      dft_style:
      clex_style:
      hull_tol: (number, default 1e-8) tolerance to decide if configuration is on the hull
    """
    if data is None:
        self.data = PlottingData()
    self.data = data

    if project is None:
        project = casm.project.project_path()
    proj = self.data.project(project, kwargs={"verbose":False})

    prim = proj.prim
    if prim.n_independent_compositions != 1:
        print("in project:", proj.path)
        print("n_independent_compositions:", prim.n_independent_compositions)
        raise Exception("Currently ConvexHullPlot only works for binary alloys")

    if selection is None:
      selection = "MASTER"
    self.sel = self.data.selection(project, selection, proj_kwargs={"verbose":False})

    if hull_selection is None:
      hull_selection = selection
    self.hull_sel = self.data.selection(project, hull_selection, proj_kwargs={"verbose":False})

    self.hull_sel_calculated = casm.project.Selection(self.sel.proj, self.hull_sel.path + ".calculated")
    self.dft_hull_sel = casm.project.Selection(self.sel.proj, self.hull_sel.path + ".dft_hull")

    self.data.add_selection([self.hull_sel_calculated, self.dft_hull_sel])

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

    self.dft_style = casm.plotting.dft_hull_style(copy.deepcopy(dft_style))
    self.clex_style = casm.plotting.clex_hull_style(copy.deepcopy(clex_style))

    self.tooltips = tooltips


  def plot(self, fig=None, tap_action=None):

    if fig is None:
        fig = bokeh.models.Figure(**default_figure_kwargs)

    update_dft_hull_glyphs(self.sel, self.dft_style)
    update_clex_hull_glyphs(self.sel, self.clex_style)

    # plot formation energy vs comp, with convex hull states
    self.renderers = []

    # dft
    style = self.dft_style

    self.r_dft = getattr(fig, style['marker'])(self.x, self.y, source=self.sel.src,
      size='radii', fill_color='color', fill_alpha=style['fill_alpha'],
      line_color='line_color', line_width='line_width', line_alpha='line_alpha',
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'],
      legend="dft")
    self.r_dft_hull_line = fig.line(self.x, self.y, source=self.calc_hull_src,
      line_color=style['on_hull']['line_color'], line_width=style['hull_line_width'],
      line_dash=style['hull_line_dash'],
      legend="dft")

    # clex
    style = self.clex_style
    self.r_clex = getattr(fig, style['marker'])(self.x, self.y_clex, source=self.sel.src,
      size='clex_radii', fill_color='clex_color', fill_alpha=style['fill_alpha'],
      line_color='clex_line_color', line_width='clex_line_width', line_alpha='clex_line_alpha',
      hover_alpha=style['hover_alpha'], hover_color=style['hover_color'],
      legend="clex")
    self.r_clex_hull_line = fig.line(self.x, self.y_clex, source=self.clex_hull_src,
      line_color=style['on_hull']['line_color'], line_width=style['hull_line_width'],
      line_dash=style['hull_line_dash'],
      legend="clex")

    # clex of dft hull
    style = self.clex_style
    self.r_clex_of_dft_hull_line = fig.line(self.x, self.y_clex, source=self.clex_of_dft_hull_src,
      line_color=style['on_hull']['line_color'], line_width=style['hull_line_width'],
      line_dash=style['clex_of_dft_hull_line_dash'],
      legend="clex_of_dft_hull")

    self.renderers = [self.r_dft,
      self.r_dft_hull_line, self.r_clex,
      self.r_clex_hull_line, self.r_clex_of_dft_hull_line]

    fig.xaxis.axis_label = self.x
    fig.yaxis.axis_label = self.y

    if tap_action is not None:
        tap_action.add_callback(self.sel, view_on_tap, [self.r_dft, self.r_clex])

    # hover over a point to see 'configname', 'Ef', and 'comp(a)'
    tooltips = [
        ("configname","@configname")
    ]

    if self.tooltips is not None:
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

        for col in self.tooltips:

            if col in float_dtypes:
                tooltips.append((col,"@{" + col + "}{1.1111}"))
            else:
                tooltips.append((col,"@{" + col + "}"))

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

        fig.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[self.r_dft, self.r_clex]))
    else:
        fig.add_tools(bokeh.models.HoverTool(tooltips=None, renderers=[self.r_dft, self.r_clex]))


  def _set_on_hull(self, sel, on_hull_label, hull_dist_label):
    """
    Use the 'hull_dist_label' column, and self.hull_tol to create the boolean
    column 'on_hull_label' indicating which configurations are on the hull.
    """
    sel.data.loc[_unselected(sel), on_hull_label] = False
    on_hull = list(map(lambda x: abs(x) < self.hull_tol, sel.data.loc[_selected(sel), hull_dist_label]))
    sel.data.loc[_selected(sel), on_hull_label] = on_hull


  def query(self):

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
      clex_hull_dist_per_atom(self.hull_sel.path)] + self.tooltips,
      force=force)

    # add project and selection path
    N = self.sel.data.shape[0]
    self.sel.add_data('project', [self.sel.proj.path]*N)
    self.sel.add_data('selection', [self.sel.path]*N)


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


class Scatter(object):
  """
  Attributes:
    data:
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
    y: the value of the y-axis, 'comp(a)' by default
    legend:
    index:
    series_name:
    style:
    tooltips:
  """

  def __init__(self, data=None, project=None, selection="MASTER",
    x='comp(a)', y='formation_energy', tooltips=[],
    legend=None, style={}, index=0, series_name=None, type=None):
    """
    Arguments:
        data: PlottingData object
        project: path to CASM project (optional, default=casm.project.project_path())
        selection: path for CASM Selection
        x: (str) column name to use as x axis. Default='comp(a)'
        y: (str) column name to use as y axis. Default='formation_energy'
        tooltips:
        legend:
        style:
        index:
        series_name:
        type:
    """
    if data is None:
        self.data = PlottingData()
    self.data = data

    if project is None:
        project = casm.project.project_path()

    if selection is None:
      selection = "MASTER"
    self.sel = self.data.selection(project, selection, proj_kwargs={"verbose":False})

    self.x = x
    self.y = y

    if legend is None:
        legend = y
    elif legend == "off" or legend == "none":
        legend = None
    self.legend = legend
    self.index = index
    if series_name is None:
        series_name = str(self.index)
    self.series_name = series_name
    self.style = scatter_series_style(self.index, copy.deepcopy(style))

    self.tooltips = tooltips


  def query(self):
      columns = [self.x, self.y]
      if self.tooltips is not None:
          columns += self.tooltips
      self.sel.query(columns)
      for col in self.sel.data.columns:
          add_src_data(self.sel, col, self.sel.data.loc[:,col])

      # add project and selection path
      N = self.sel.data.shape[0]
      self.sel.add_data('project', [self.sel.proj.path]*N)
      self.sel.add_data('selection', [self.sel.path]*N)


  def plot(self, fig=None, tap_action=None):

      update_scatter_glyphs(self.sel, self.style, self.series_name)

      # avoid having a legend name same as a column name
      while self.legend in self.sel.src.column_names:
          self.legend += ' '

      self.r = getattr(fig, self.style['marker'])(
          self.x,
          self.y,
          source=self.sel.src,
          size='radii.' + self.series_name,
          fill_color='color.' + self.series_name,
          fill_alpha='fill_alpha.' + self.series_name,
          line_color='line_color.' + self.series_name,
          line_width='line_width.' + self.series_name,
          line_alpha='line_alpha.' + self.series_name,
          hover_alpha=self.style['hover_alpha'],
          hover_color=self.style['hover_color'],
          legend=self.legend)

      self.renderers = [self.r]

      if self.index == 0:
          fig.xaxis.axis_label = self.x
          fig.yaxis.axis_label = self.y

      if tap_action is not None:
          tap_action.add_callback(self.sel, casm.plotting.view_on_tap, [self.r])


      if self.tooltips is not None:
          tooltips = []

          for col in self.tooltips + [self.x, self.y]:
              if col in float_dtypes:
                  tooltips.append((col,"@{" + col + "}{1.1111}"))
              else:
                  tooltips.append((col,"@{" + col + "}"))

          fig.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[self.r]))
      else:
          fig.add_tools(bokeh.models.HoverTool(tooltips=None, renderers=[self.r]))


class Histogram(object):
  """
  Attributes:
    data:
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
    y: the value of the y-axis, 'comp(a)' by default
    legend:
    index:
    series_name:
    style:
    tooltips:

    src:
    y_label:
  """

  def __init__(self, data=None, project=None, selection="MASTER",
    x='comp(a)', hist_kwargs={}, tooltips=[],
    style={}, index=0, series_name=None, type=None):
    """
    Arguments:
        data: PlottingData object
        project: path to CASM project (optional, default=casm.project.project_path())
        selection: path for CASM Selection
        x: (str) column name to use as x axis. Default='comp(a)'
        hist_kwargs:
        tooltips:
        legend:
        style:
        index:
        series_name:
        type:
    """
    if data is None:
        self.data = PlottingData()
    self.data = data

    if project is None:
        project = casm.project.project_path()

    if selection is None:
      selection = "MASTER"
    self.sel = self.data.selection(project, selection, proj_kwargs={"verbose":False})

    self.x = x
    self.hist_kwargs = hist_kwargs

    self.index = index
    if series_name is None:
        series_name = str(self.index)
    self.series_name = series_name
    self.style = scatter_series_style(self.index, copy.deepcopy(style))

    self.tooltips = tooltips

    self.src = None

  def query(self):
      columns = [self.x]
      self.sel.query(columns)
      for col in self.sel.data.columns:
          add_src_data(self.sel, col, self.sel.data.loc[:,col])

      if 'to_be_selected' not in self.sel.data.columns:
          self.sel.data.loc[:,'to_be_selected'] = self.sel.data.loc[:,'selected']

      # add project and selection path
      N = self.sel.data.shape[0]
      self.sel.add_data('project', [self.sel.proj.path]*N)
      self.sel.add_data('selection', [self.sel.path]*N)

      # create histogram data
      tmp_kwargs = dict(self.hist_kwargs)
      tmp_kwargs['density'] = False
      hist, edges = np.histogram(self.sel.data.loc[:,self.x], **tmp_kwargs)

      sel_kwargs = dict(self.hist_kwargs)
      sel_kwargs['bins'] = edges
      sel_kwargs['density'] = False

      hist_sel, edges_sel = np.histogram(
          self.sel.data.loc[self.sel.data['to_be_selected']==True, self.x],
          **sel_kwargs)

      # need to do 'density' of selected by hand
      if 'density' in self.hist_kwargs and self.hist_kwargs['density']:
           hist_freq = np.array(hist)
           width = edges[1:] - edges[:-1]
           area = hist_freq.dot(width)
           hist = hist_freq*width/area
           hist_sel = hist_sel*width/area

      if self.src is None:
          self.src = bokeh.models.ColumnDataSource()
          set_src_data(self.src, "left", edges[:-1], force=True)
          set_src_data(self.src, "right", edges[1:], force=True)

      set_src_data(self.src, "top", hist, force=True)
      set_src_data(self.src, "bottom", [0]*len(hist), force=True)

      set_src_data(self.src, "top_sel", hist_sel, force=True)
      set_src_data(self.src, "top_unsel", hist - hist_sel, force=True)
      set_src_data(self.src, "bottom_sel", [0]*len(hist_sel), force=True)

      if 'density' not in self.hist_kwargs.keys() or self.hist_kwargs['density'] == False:
          self.y_label = "Frequency"
      else:
          self.y_label = "Density"


  def plot(self, fig=None, tap_action=None):

      update_scatter_glyphs(self.sel, self.style, self.series_name)

      self.r_quad = fig.quad(
          top='top',
          bottom='bottom',
          left='left',
          right='right',
          source=self.src,
          fill_alpha=0.0,
          line_color='black',
          line_width=1)

      self.r_quad_sel = fig.quad(
          top='top_sel',
          bottom='bottom_sel',
          left='left',
          right='right',
          source=self.src,
          fill_alpha=self.style['selected']['fill_alpha'],
          fill_color=self.style['selected']['color'],
          line_alpha=0.0,
          line_width=0.0)

      self.renderers = [self.r_quad, self.r_quad_sel]

      fig.xaxis.axis_label = self.x
      fig.yaxis.axis_label = self.y_label

      if self.tooltips is not None:
          tooltips = [
              (self.x,"@{left}{1.1111}"),
              ("selected","@{top_sel}{1.1111}"),
              ("unselected","@{top_unsel}{1.1111}"),
              ("total","@{top}{1.1111}")
          ]

          fig.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[self.r_quad]))
      else:
          fig.add_tools(bokeh.models.HoverTool(tooltips=None, renderers=self.renderers))


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
          print("Unknown or unsupported plot type:", _type.name)
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
    sel: a CASM Selection used to make the figure
    score_id: a UUID to use for the 'score' column name in the sel.src ColumnDataSource
    rank_id: a UUID to use for the 'rank' column name in the sel.src ColumnDataSource
    max_score: the maximum score
    min_score: the minimum score
  """
  def __init__(self, data=None, project=None, selection="MASTER", to_query=[],
    scoring_query=None, scoring_module=None, scoring_function=None,
    tooltips=[], index=0, type=None,
    style={}, series_name=None, name=None):
    """
    Arguments:
      data: PlottingData instance
      project: str (optional, default=casm.project.project_path())
      selection: str (optional, default="MASTER")
      query: List[str]
      scoring_query: str
      scoring_module: str
      scoring_function: str
      tooltips: List[str]
      tooltips_exclude: List[str]
      index: int
      type: str
    """
    if data is None:
        self.data = PlottingData()
    self.data = data

    if project is None:
        project = casm.project.project_path()

    if selection is None:
      selection = "MASTER"
    self.sel = self.data.selection(project, selection, proj_kwargs={"verbose":False})

    self.scoring_query = scoring_query

    self.to_query = to_query
    self.scoring_module = scoring_module
    self.scoring_function = scoring_function

    if self.scoring_query is not None:
        if self.scoring_query not in self.to_query:
            self.to_query.append(self.scoring_query)

    self.style = rankplot_style(copy.deepcopy(style))

    self.index = index
    if series_name is None:
        series_name = str(self.index)
    self.series_name = series_name

    self.tooltips = tooltips

    if name is None:
        name = scoring_query
    if name is None:
        name = scoring_module + '.' + scoring_function
    self.name = name

    self.score_id = str(uuid.uuid4())
    self.rank_id = str(uuid.uuid4())

    self._max_score = None
    self._min_score = None


  def plot(self, fig=None, tap_action=None):

    update_scatter_glyphs(self.sel, self.style, self.series_name)

    self.r_stem = fig.segment(
      self.rank_id, 0, self.rank_id, self.score_id, source=self.sel.src,
      line_color='line_color.' + self.series_name,
      line_width='line_width.' + self.series_name,
      line_alpha='line_alpha.' + self.series_name,
      hover_alpha=self.style['hover_alpha'],
      hover_color=self.style['hover_color'])

    self.r_circ = fig.circle(
      self.rank_id,
      self.score_id,
      source=self.sel.src,
      size='radii.' + self.series_name,
      fill_color='color.' + self.series_name,
      fill_alpha='fill_alpha.' + self.series_name,
      line_width=0.,
      line_alpha=0.,
      hover_alpha=self.style['hover_alpha'],
      hover_color=self.style['hover_color'])

    fig.xaxis.axis_label = "Rank"
    fig.yaxis.axis_label = self.name


    self.renderers = [self.r_stem, self.r_circ]

    if tap_action is not None:
        tap_action.add_callback(self.sel, casm.plotting.view_on_tap, [self.r_circ])

    if self.tooltips is not None:
        tooltips = []

        for col in self.tooltips + [self.score_id, self.rank_id]:

            if col == self.score_id:
                tooltips.append(("score","@{" + self.score_id + "}{1.1111}"))
            elif col == self.rank_id:
                tooltips.append(("rank","@{" + self.rank_id + "}{1.1111}"))
            elif col in float_dtypes:
                tooltips.append((col,"@{" + col + "}{1.1111}"))
            else:
                tooltips.append((col,"@{" + col + "}"))

        fig.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[self.r_circ]))
    else:
        fig.add_tools(bokeh.models.HoverTool(tooltips=None, renderers=[self.r_circ]))

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

  def score(self):
    """
    Calculate 'score' and 'rank', but do not change selection. May also change the
    'scoring' function.
    """
    if self.scoring_query is not None:
        try:
            self.sel.data.loc[:,self.score_id] = self.sel.data.loc[:,self.scoring_query]
        except Exception as e:
            print("scoring_query:", self.scoring_query)
            print("columns:", self.sel.data.columns)
            raise e
    else:
        f, filename, description = imp.find_module(self.scoring_module)
        try:
            module = imp.load_module(self.scoring_module, f, filename, description)
            self.sel.data.loc[:,self.score_id] = getattr(module, self.scoring_function)(self.sel)
        finally:
            if f:
                f.close()

    # convert to numeric (NaN)
    self.sel.data.loc[:,self.score_id] = pandas.to_numeric(self.sel.data[self.score_id], errors='coerce')

    # add 'score' to src, as self.score_id
    add_src_data(self.sel, self.score_id, self.sel.data.loc[:,self.score_id], force=True)

    # rank
    self.sel.data.loc[self.sel.data.sort_values([self.score_id, 'configname']).index, self.rank_id] = range(self.sel.data.shape[0])

    # add 'rank' to src as self.rank_id
    add_src_data(self.sel, self.rank_id, self.sel.data.loc[:,self.rank_id], force=True)

    self._min_score = None
    self._max_score = None

  def query(self):
    columns = copy.deepcopy(self.to_query)
    if self.tooltips is not None:
        columns += self.tooltips
    self.sel.query(columns)

    # add project and selection path
    N = self.sel.data.shape[0]
    self.sel.add_data('project', [self.sel.proj.path]*N)
    self.sel.add_data('selection', [self.sel.path]*N)
    self.sel.add_data('configname')

    for col in self.sel.data.columns:
        add_src_data(self.sel, col, self.sel.data.loc[:,col])

    self.score()


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
    self.tools = "crosshair,pan,reset,box_zoom,wheel_zoom,save"
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
    _hplot = bk_layouts.row(self.input['cutoff'].widget, self.select_mode, self.select_action)
    self.layout = bk_layouts.column(self.p, _hplot, self.msg.widget, width=self.p.plot_width)


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
      var data = src.data;
      var loc_data = loc.data;
      data['y'][0] = loc_data['y'][0];
      data['y'][1] = loc_data['y'][0];
      src.data = data;
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
      print("Error applying scoring function")
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
        var data = src.data;
        data['y'][0] = cb_data.geometry.y;
        data['y'][1] = cb_data.geometry.y;
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
        except Exception as e:
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
    set_src_data(self.cutoff_line_src, 'y', [self.input['cutoff'].value, self.input['cutoff'].value], force=True)
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
      with open(self.input_filename, 'wb') as file:
        file.write(six.u(json.dumps(self.fit_input, indent=2)).encode('utf-8'))
    else:
      self.fit_input = json.loads(open(self.input_filename, 'rb').read().decode('utf-8'))

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

    for key, widget in six.iteritems(self.input):
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
      with open(self.input_filename, 'wb') as file:
        file.write(six.u(json.dumps(self.fit_input, indent=2)).encode('utf-8'))
      self.msg.value = "Saved input file: " + self.input_filename
    except Exception as e:
      self.msg.value = str(e)

  def _apply_and_save(self):
    try:
      self._apply()
      self._save()
      self.msg.value = "Applied sample weights and saved input file: " + self.input_filename
    except Exception as e:
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
        set_src_data(self.ref_line_src, 'x', [self.x_min, self.x_max], force=True)
        set_src_data(self.ref_line_src, 'y', [self.input['Eref'].value, self.input['Eref'].value], force=True)
      elif self.select_method.value == "wEmin":
        set_src_data(self.ref_line_src, 'x', [self.x_min, self.x_max], force=True)
        set_src_data(self.ref_line_src, 'y', [self.Ef_min, self.Ef_min], force=True)
      elif self.select_method.value == "wHullDist":
        set_src_data(self.ref_line_src, 'x', [], force=True)
        set_src_data(self.ref_line_src, 'y', [], force=True)
    if self.select_method.value == "wEref":
      set_src_data(self.ref_line_src, 'y', [self.input['Eref'].value, self.input['Eref'].value], force=True)

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
    row = bk_layouts.row(bk_layouts.column(*self._widget_layout()), grid)
    return bk_layouts.column(row, self.msg.widget, width=self.whullplot.p.plot_width*3 + 300)

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
    self.fit_input = json.loads(open(self.input_filename, 'rb').read().decode('utf-8'))
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
    self.eci = pandas.DataFrame(eci, columns=[str(i) for i in range(len(self.hall))])

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

    self.layout = bk_layouts.row(self.cv_rankplot.layout, self.eci_stemplot)


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
      print("tap!")
      if len(sel.src.selected['1d']['indices']) == 1:
        print("update!")
        index = sel.src.selected['1d']['indices'][0]
        set_src_data(sel.eci_src, 'current', sel.eci_src.data[str(index)], force=True)

    # self, sel, scoring, name="Score", mode='set', tooltips=None, on_tap=view_on_tap)
    self.cv_rankplot = RankPlot( self.sel, score, name="CV Score",
      tooltips=tooltips, on_tap=_update_eci_stemplot)

    self.cv_rankplot.p.xaxis.axis_label = "Individual"
    self.cv_rankplot.p.plot_width = 400

  def _eci_stemplot(self, index):

    self.sel.eci_src = bokeh.models.ColumnDataSource(data=self.sel.eci)

    with open(self.proj.dir.basis(self.proj.settings.default_clex), 'rb') as f:
      self.basis = json.loads(f.read().decode('utf-8'))

    set_src_data(sel.eci_src, 'index', range(self.sel.eci.shape[0]), force=True)
    set_src_data(sel.eci_src, 'branch', [j["orbit"][0] for j in self.basis["cluster_functions"] ], force=True)
    set_src_data(sel.eci_src, 'mult', [j["mult"] for j in self.basis["cluster_functions"] ], force=True)
    set_src_data(sel.eci_src, 'max_length', [j["prototype"]["max_length"] for j in self.basis["cluster_functions"] ], force=True)
    set_src_data(sel.eci_src, 'min_length', [j["prototype"]["min_length"] for j in self.basis["cluster_functions"] ], force=True)
    colormap = ["black", "blue", "red", "green", "cyan", "magenta", "yellow", "orange"]
    set_src_data(sel.eci_src, 'color', list(map(lambda x: colormap[x % len(colormap)], self.sel.eci_src.data["branch"])), force=True)
    set_src_data(sel.eci_src, 'current', self.sel.eci_src.data[str(index)], force=True)

    # TOOLS in bokeh plot
    _tools = ["crosshair,pan,reset,box_zoom,wheel_zoom,save"]

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
    self.layout = bk_layouts.row(self.widget, width=200)

  def on_change(self, attrname, old, new):
    for f in self.update:
      f(attrname, old, new)

  @property
  def value(self):
    return self.widget.value

  @value.setter
  def value(self, _value):
    self.widget.value = _value
