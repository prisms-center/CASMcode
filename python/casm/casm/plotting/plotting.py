import bokeh.plotting
import bokeh.models
import numpy as np
import pandas
import uuid
import casm.project

from bokeh.plotting import hplot, vplot
from bokeh.client import push_session
from bokeh.io import curdoc, vform


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


def hull_dist():
  """Equivalent to 'hull_dist(CALCULATED,comp)'"""
  return 'hull_dist(CALCULATED,comp)'

def hull_dist_per_atom():
  """Equivalent to 'hull_dist(CALCULATED,atom_frac)'"""
  return 'hull_dist(CALCULATED,atom_frac)'

def add_src_data(sel, name, data, force=False):
  if sel.src is None:
    sel.src = bokeh.models.ColumnDataSource(data={name:data})
    sel.src_update = []
    def callbacks(attrname, old, new):
      for f in sel.src_update:
        f(attrname, old, new)
    sel.src.on_change('data', callbacks)
  else:
    if name not in sel.src.data.keys() or force == True:
      sel.src.data[name] = data
    
def add_hull_data(sel, hull_tol=1e-8):
  
  sel.add_data(hull_dist())
  add_src_data(sel, hull_dist(), sel.data.loc[:,hull_dist()])
  
  sel.add_data(hull_dist_per_atom())
  
  sel.add_data('on_hull', sel.data.loc[:,hull_dist_per_atom()] < hull_tol, force=True)

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


class ConfigStyle(object):
  def __init__(self):
    ### style
    self.selected_color = "blue"
    self.selected_radii = 10
    self.unselected_color = "green"
    self.unselected_radii = 7
    self.fill_alpha = 0.5
    
    self.cutoff_color = "red"
    self.cutoff_line_width = 1
    self.cutoff_alpha = 1.0

    self.off_hull_line_color = "blue"
    self.on_hull_line_color = "#FF4105"
    self.off_hull_line_width = 0.0
    self.on_hull_line_width = 2.0
    self.on_hull_line_alpha = 0.8
    self.off_hull_line_alpha = 0.0

    self.hover_color = "orange"
    self.hover_alpha = 0.7


def update_glyphs(sel, selected=None, on_hull=None):
  """
  Use the selected and on_hull iterables to update glyph styles based on sel.style.
  """
  
  if selected is None:
    selected = sel.data['selected']
  
  if on_hull is None:
    on_hull = sel.data['on_hull']
  
  style = sel.style
  
  # colors to plot selected (True), and unselected (False) configurations
  colormap = dict({True:style.selected_color, False:style.unselected_color})
  
  # colors to plot on_hull (True), and not on_hull (False) configurations
  linecolormap = dict({True:style.on_hull_line_color, False:style.off_hull_line_color})
  
  # width to plot on_hull (True), and not on_hull (False) configurations
  linewidthmap = dict({True:style.on_hull_line_width, False:style.off_hull_line_width})
  
  # width to plot on_hull (True), and not on_hull (False) configurations
  linealphamap = dict({True:style.on_hull_line_alpha, False:style.off_hull_line_alpha})
  
  # radii to plot selected (True), and unselected (False) configurations
  radiimap = dict({True:style.selected_radii, False:style.unselected_radii})
  
  sel.src.data['color'] = map(lambda x: colormap[x], selected)
  sel.src.data['radii'] = map(lambda x: radiimap[x], selected)
  sel.src.data['line_color'] = map(lambda x: linecolormap[x], on_hull)
  sel.src.data['line_width'] = map(lambda x: linewidthmap[x], on_hull)
  sel.src.data['line_alpha'] = map(lambda x: linealphamap[x], on_hull)
  

class ConvexHullPlot(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
    x_id: a UUID to use for the x-axis column name in the sel.src ColumnDataSource
  """
  
  def __init__(self, sel, x='comp(a)'):
    """
    Arguments:
      sel: A CASM Selection
      x: (str) column name to use as x axis. Default='comp(a)'.
    """
    self.sel = sel
    self.x_id = str(uuid.uuid4())
    self.x = x
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    
    #add data to sel.src if necessary
    self._query()
    
    # will store the plot object in self.p_
    self._plot()
    
    # store plot in layout
    self.layout = self.p
  
  @property
  def p(self):
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'style'):
      self.sel.style = ConfigStyle()
      update_glyphs(self.sel)
    style = self.sel.style
    
    # plot formation energy vs comp, with convex hull states
    self.p_ = bokeh.plotting.Figure(plot_width=800, plot_height=400, tools=self.tools)
    cr = self.p_.circle(self.x_id, 'formation_energy', source=self.sel.src, 
      size='radii', fill_color='color', fill_alpha=style.fill_alpha, 
      line_color='line_color', line_width='line_width', line_alpha='line_alpha',
      hover_alpha=style.hover_alpha, hover_color=style.hover_color, 
      legend="Formation energy")
    hull_line = self.p_.line(self.x_id, 'formation_energy', source=self.calc_hull_src, 
      line_color=style.on_hull_line_color, line_width=1.0, 
      legend="Ground states")
#    self.p_.circle(self.x_id, 'formation_energy',  source=self.calc_hull_src, 
#      size=style.radiimap[True], fill_color=None, line_color=style.hull_color, line_width=1, 
#      legend="Ground states")

    self.p_.legend.location = "bottom_left"
    self.p_.xaxis.axis_label = self.x
    self.p_.yaxis.axis_label = "Formation energy"

    # hover over a point to see 'configname', 'Ef', and 'comp(a)'
    tooltips = [
        ("configname","@configname"), 
        ("Ef","@formation_energy{1.1111}"),
        ("hull_dist", "@{" + hull_dist() + "}{1.1111}")
    ]
    
    i = ord('a')
    while "comp(" + chr(i) + ")" in self.sel.src.data.keys():
      s = "comp(" + chr(i) + ")"
      tooltips.append((s, "@{" + s + "}{1.1111}"))
      i += 1
    
    self.tap_action = ConfigurationTapAction(self.sel)

    self.p_.add_tools(self.tap_action.tool(view_on_tap, [cr]))
    self.p_.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[cr]))
    self.p_.add_tools(bokeh.models.BoxSelectTool(renderers=[cr]))
    self.p_.add_tools(bokeh.models.LassoSelectTool(renderers=[cr]))
  
  def _query(self):
    #add data to sel.data if necessary
    col = ['configname',
           'selected',
           'comp',
           'comp_n',
           'formation_energy']
    
    df = casm.project.query(self.sel.proj, col, self.sel)
    for c in df.columns:
      self.sel.add_data(c, df[c])
      add_src_data(self.sel, c, self.sel.data.loc[:,c])
    
    add_hull_data(self.sel)
    
    # add a column to selection src to use for the x-axis / x-data
    # we use a uuid for the name so we can enable interactively changing it later
    add_src_data(self.sel, self.x_id, self.sel.data.loc[:,self.x])
    
    # get hull data, sorted by 'x' axis value for easy plotting
    df_calc_hull = self.sel.data[self.sel.data['on_hull'] == True].sort_values([hull_dist(), self.x])
    self.calc_hull_src = bokeh.models.ColumnDataSource(data=df_calc_hull)
    self.calc_hull_src.data[self.x_id] = self.calc_hull_src.data[self.x]


class Scatter(object):
  """
  Attributes:
    p: a bokeh Figure containing formation energies and the convex hull
    layout: a bokeh layout element holding p
    sel: a CASM Selection used to make the figure
    x: the value of the x-axis, 'comp(a)' by default
    x_id: a UUID to use for the x-axis column name in the sel.src ColumnDataSource
    y: the value of the y-axis, 'comp(a)' by default
    y_id: a UUID to use for the y-axis column name in the sel.src ColumnDataSource
  """
  
  def __init__(self, sel, x='comp(a)', y='comp(a)', show_gs=True):
    """
    Arguments:
      sel: A CASM Selection
      x: (str) property name to use as x axis, 'comp(a)' by default.
      y: (str) property name to use as y axis, 'comp(a)' by default.
      show_gs: (Boolean) whether to outline ground states
    """
    print "begin Scatter"
    self.sel = sel
    
    self.x = x
    self.sel.add_data('configname')
    add_src_data(self.sel, 'configname', self.sel.data['configname'])
    self.sel.add_data(x)
    # add a column to selection src to use for the x-axis / x-data
    # we use a uuid for the name so we can enable interactively changing it later
    self.x_id = str(uuid.uuid4())
    add_src_data(self.sel, self.x_id, self.sel.data.loc[:,self.x])
    
    self.y = y
    self.sel.add_data(y)
    # add a column to selection src to use for the y-axis / y-data
    # we use a uuid for the name so we can enable interactively changing it later
    self.y_id = str(uuid.uuid4())
    add_src_data(self.sel, self.y_id, self.sel.data.loc[:,self.y])
    
    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    
    # get hull data, sorted by 'x' axis value for easy plotting
    self.show_gs = show_gs
    if self.show_gs:
      add_hull_data(self.sel)
    
    self._plot()
    
    # store plot in layout
    self.layout = self.p
  
  @property
  def p(self):
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'style'):
      self.sel.style = ConfigStyle()
      update_glyphs(self.sel)
    
    style = self.sel.style
    
    # plot formation energy vs comp, with convex hull states
    if True:
      self.p_ = bokeh.plotting.Figure(plot_width=400, plot_height=400, tools=self.tools)
      
      cr = self.p_.circle(self.x_id, self.y_id, source=self.sel.src, 
        size='radii', fill_color='color', fill_alpha=style.fill_alpha,
        line_color='line_color', line_width='line_width', line_alpha='line_alpha',
        hover_alpha=style.hover_alpha, hover_color=style.hover_color)
    
    else:
      # rug plot
      x_max = max(self.sel.src.data[self.x_id])
      x_min = min(self.sel.src.data[self.x_id])
      y_max = max(self.sel.src.data[self.y_id])
      y_min = min(self.sel.src.data[self.y_id])
      x_rug_length = (y_max - y_min)/20.0
      y_rug_length = (x_max - x_min)/20.0
      
      self.p_ = bokeh.plotting.Figure(
        plot_width=400, plot_height=400, tools=self.tools,
        x_range=(x_min-y_rug_length, x_max+2*y_rug_length), y_range=(y_min-x_rug_length, y_max+2*x_rug_length))
      
      cr = self.p_.circle(self.x_id, self.y_id, source=self.sel.src, 
        size='radii', fill_color='color', fill_alpha=style.fill_alpha,
        line_color='line_color', line_width='line_width', line_alpha='line_alpha',
        hover_alpha=style.hover_alpha, hover_color=style.hover_color)
      
      x_rug = self.p_.segment(
        x0=self.x_id, x1=self.x_id,
        y0=y_max + x_rug_length, y1=y_max + 2.0*x_rug_length,
        source=self.sel.src,
        line_color='line_color', line_width=1.0)
      
      y_rug = self.p_.segment(
        x0=x_max + y_rug_length, x1=x_max + 2.0*y_rug_length,
        y0=self.y_id, y1=self.y_id,
        source=self.sel.src,
        line_color='line_color', line_width=1.0)
    
    
    self.p_.xaxis.axis_label = self.x
    self.p_.yaxis.axis_label = self.y

    # hover over a point to see 'configname', 'x', and 'y'
    tooltips = [
        ("configname","@configname"), 
        (self.x,"@{" + self.x_id + "}{1.1111}"),
        (self.y,"@{" + self.y_id + "}{1.1111}")
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
    x_id: a UUID to use for the x-axis column name in the sel.src ColumnDataSource
  """
  
  def __init__(self, sel, x='comp(a)', **kwargs):
    """
    Arguments:
      sel: A CASM Selection
      x: (str) property name to use as x axis, 'comp(a)' by default.
      kwargs: kwargs to pass to numpy.histogram
    """
    print "begin Histogram"
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
    self.sel.src_update.append(self.update)
    
    # store plot in layout
    self.layout = self.p
  
  @property
  def p(self):
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'style'):
      self.sel.style = ConfigStyle()
      update_glyphs(self.sel)
    
    style = self.sel.style
    
    self.p_ = bokeh.plotting.Figure(plot_width=400, plot_height=400, tools=self.tools)
    
    # set data in 'self.src'
    self.update(0, 0, 0)
    
    qd = self.p_.quad(top='top', bottom='bottom', left='left', right='right', source=self.src, 
      fill_alpha=0.0, line_color='black', line_width=1)
    
    qd_sel = self.p_.quad(top='top_sel', bottom='bottom_sel', left='left', right='right', source=self.src, 
      fill_alpha=style.fill_alpha, fill_color=style.selected_color,
      line_alpha=0.0, line_width=0.0)
        
    self.p_.xaxis.axis_label = self.x
    self.p_.yaxis.axis_label = self.y_label

    # hover over a point to see 'x', and 'Freq.'
    tooltips = [
        (self.x,"@{left}{1.1111}"),
        (self.y_label,"@{top}{1.1111}")
    ]
    
    self.p_.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[qd]))
  
  def update(self, attrname, old, new):
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
    x_id: a UUID to use for the x-axis column name in the sel.src ColumnDataSource
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
    print "begin GridPlot"
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
    self.sel.src_update.append(self.update)
    
    # layout
    self.layout = bokeh.models.GridPlot(children=self.p_)
    
    print "finish GridPlot"
  
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
  
  def update(self, attrname, old, new):
    for row in self.casm_p_:
      for casm_p_ in row:
        if hasattr(casm_p_, 'update'):
          casm_p_.update(attrname, old, new)

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
  
  def __init__(self, sel, scoring, session, name="Score", mode='set'):
    """
    Arguments:
      sel: A CASM Selection
      scoring: (function) a scoring function to operate on a row of sel.data
      session: a bokeh push_session, if gui==True
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
    
    # a 'RankSelectCutoff' instance
    self.gui = True
    self.session = session

    # TOOLS in bokeh plot
    self.tools = "crosshair,pan,reset,resize,box_zoom"
    
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
    if self.gui:
      _hplot = hplot(self.input['cutoff'].widget, self.select_mode, self.select_action)
      self.layout = vplot(self.p, _hplot, self.msg.widget, width=self.p.plot_width)
    else:
      self.layout = vplot(self.p)
  
  @property
  def p(self):
    if self.p_ is None:
      self._plot()
    return self.p_
  
  def _plot(self):
    
    # add 'style' to the selection if not already existing
    if not hasattr(self.sel, 'style'):
      self.sel.style = ConfigStyle()
      update_glyphs(self.sel)
    
    style = self.sel.style
    
    if self.gui:
      _tools = [self.mouse, self.tools]
    else:
      _tools = [self.tools]
    
    self.p_ = bokeh.plotting.Figure(plot_width=800, plot_height=400, tools=_tools)
    p_stem = self.p_.segment(self.rank_id, 0, self.rank_id, self.score_id, source=self.sel.src, 
      line_color='color', line_width=1.0,
      hover_alpha=style.hover_alpha, hover_color=style.hover_color)
    
    p_circ = self.p_.circle(self.rank_id, self.score_id, source=self.sel.src, 
      size='radii', fill_color='color', fill_alpha=style.fill_alpha,
      line_color='line_color', line_width='line_width', line_alpha='line_alpha',
      hover_alpha=style.hover_alpha, hover_color=style.hover_color)

    self.p_.xaxis.axis_label = "Rank"
    self.p_.yaxis.axis_label = "Score"

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
    
    
    if self.gui:
      self.tap_action = ConfigurationTapAction(self.sel)
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
      line_width=style.cutoff_line_width, line_alpha=style.cutoff_alpha, color=style.cutoff_color)

    self.p_.add_tools(bokeh.models.HoverTool(tooltips=tooltips, renderers=[p_circ]))
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
      self.df = pandas.DataFrame(self.sel.data.apply(self.scoring,axis='columns'), columns=['score'])
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
    if self.gui:
      
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
      self.select_mode.on_change(
        'value', 
        select_mode_f)
      
      self.select_action = bokeh.models.Select(
        value="Select action", 
        options=["Select action", "Revert", "Apply", "Save", "Close"])
      def select_action_f(attrname, old, new):
        if new == "Revert":
          self.sel.data.loc[:,'to_be_selected'] = self.sel.data.loc[:,'selected']
          update_glyphs(self.sel, selected=self.sel.data['to_be_selected'])
          self.msg.value = "Reverted to original selection"
        elif new == "Apply":
          self.sel.data.loc[:,'selected'] = self.sel.data.loc[:,'to_be_selected']
          update_glyphs(self.sel, selected=self.sel.data['to_be_selected'])
          self.msg.value = "Applied selection"
        elif new == "Save":
          try:
            self.sel.save(force=True)
            self.msg.value = "Saved selection: " + self.sel.path
          except Exception, e:
            self.msg.value = str(e)
        elif new == "Close":
          self.msg.value = "Closing server connection"
          self.session.close()
        self.select_action.value = "Select action"
      self.select_action.on_change('value', select_action_f)
  
  def _post_score_input(self):
    """
    Construct gui widgets
    """
    if self.gui:
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
    self.widget: a bokeh.models.SelectInput containing the cutoff value
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
    
    self.widget = bokeh.models.SelectInput(value=value, title=title, options=options)
    self.update = dict()
    self.widget.on_change('value', self.on_change)
  
  def on_change(self, attrname, old, new):
    if new in self.update.keys():
      for f in self.update[new]:
        f(attrname, old, new)
  
  @property
  def value(self):
    return self.widget.value
  
  @value.setter
  def value(self, _value):
    self.widget.value = _value


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
 

