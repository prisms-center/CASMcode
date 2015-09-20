#ifndef BP_Plot_CC
#define BP_Plot_CC

#include "casm/BP_C++/BP_Plot.hh"

namespace BP {

  ////////////////////////////////////////////////////
  /// BP_Plot_Data Functions

  void BP_Plot_Data::set_ps_data(double *axis_pos, double *axis_size, double *axis_lim) {
    //axis_pos in inches
    //axis_size in inches
    //axis_lim in (data units)

    double dx = axis_lim[1] - axis_lim[0];
    double dy = axis_lim[3] - axis_lim[2];

    if(dx == 0.0) dx = 1.0;
    if(dy == 0.0) dy = 1.0;

    double fx = axis_size[0] * 72.0 / (dx);		// ps pts / data units (x)
    double fy = axis_size[1] * 72.0 / (dy);		// ps pts / data units (x)
    double ps_x_init = axis_pos[0] * 72.0;
    double ps_y_init = axis_pos[1] * 72.0;

    ps_xdata.clear();
    ps_ydata.clear();
    for(unsigned long int i = 0; i < xdata.size(); i++) {
      ps_xdata.add(ps_x_init + (xdata[i] - axis_lim[0])*fx);
      ps_ydata.add(ps_y_init + (ydata[i] - axis_lim[2])*fy);

    }


  };


  void BP_Plot_Data::write(BP_Write &file, double *axis_pos, double *axis_size, double *axis_lim) {
    set_ps_data(axis_pos, axis_size, axis_lim);

    if(line) {
      file << "\n";
      file << "gsave" << std::endl;
      file << "newpath" << std::endl;

      for(int i = 0; i < xdata.size(); i++) {
        if(i == 0) {
          file << ps_xdata[i] << " " << ps_ydata[i] << " moveto" << std::endl;
        }
        else {
          file << ps_xdata[i] << " " << ps_ydata[i] << " lineto" << std::endl;
        }

      }

      //if( linestyle == 1) file << "[4 4] 0 setdash" << std::endl;				//dashed
      //else if( linestyle == 2) file << "[0.5 4] 0 setdash" << std::endl;		//dotted
      //else if( linestyle == 3) file << "[0.5 4 4 4] 0 setdash" << std::endl;	//dash dotted

      if(linestyle == 1) file << "[" << linewidth * 4 << " " << linewidth * 4 << "] 0 setdash" << std::endl;				//dashed
      else if(linestyle == 2) file << "[" << linewidth * 0.5 << " " << linewidth * 4 << "] 0 setdash" << std::endl;		//dotted
      else if(linestyle == 3) file << "[" << linewidth * 0.5 << " " << linewidth * 4 << " " << linewidth * 4 << " " << linewidth * 4 << "] 0 setdash" << std::endl;	//dash dotted
      else if(linestyle == 4) {
        file << "[";
        for(int i = 0; i < custom_linestyle.size(); i++) {
          file << custom_linestyle[i] << " ";
        }
        file << "] 0 setdash" << std::endl;
      }

      file << linecolor << " setrgbcolor" << std::endl;
      file << linewidth << " setlinewidth" << std::endl;
      file << "1 setlinecap" << std::endl;
      file << "stroke" << std::endl;
      file << "grestore" << std::endl;
    }

    //bool	pointface;
    //bool	pointedge;
    //BP_Vec<int>		pointstyle;	//0: circles, 1: square, 2: plus, 3: cross (x), 4: diamond, 5: triangle up, 6: triangle down
    //BP_Vec<double>	pointsize;	// pts (diameter-ish)
    //BP_Vec<BP_RGB>  pointfacecolor;
    //BP_Vec<BP_RGB>  pointedgecolor;
    //BP_Vec<double>  pointedgewidth;


    if(points) {
      // faces:  %color, x, y, r, style
      // edges:  %edgewidth, color, x, y, r, style

      file << "/circ_face {\n\tnewpath\n\t0 360 arc closepath\n\tsetrgbcolor\n\tfill\n} def\n" << std::endl;
      file << "/circ_edge {\n\tnewpath\n\t0 360 arc closepath\n\tsetrgbcolor\n\tsetlinewidth\n\tstroke\n} def\n" << std::endl;
      file << "/square_face {\n\tnewpath\n\t/r exch def\n\tmoveto\n\t-1 r mul -1 r mul rmoveto\n\t2 r mul /d exch def\n\td 0 rlineto\n\t0 d rlineto\n\t-1 d mul 0 rlineto\n\tclosepath\n\tsetrgbcolor\n\tfill\n} def\n" << std::endl;
      file << "/square_edge {\n\tnewpath\n\t/r exch def\n\tmoveto\n\t-1 r mul -1 r mul rmoveto\n\t2 r mul /d exch def\n\td 0 rlineto\n\t0 d rlineto\n\t-1 d mul 0 rlineto\n\tclosepath\n\tsetrgbcolor\n\tsetlinewidth\n\tstroke\n} def\n" << std::endl;


      file << "\n";
      file << "gsave" << std::endl;
      for(int i = 0; i < xdata.size(); i++) {
        if(pointface[i]) {
          file << pointfacecolor[i] << " " << ps_xdata[i] << " " << ps_ydata[i] << " " << pointsize[i] * 0.5;
          switch(pointstyle[i]) {
          case 0:
            file << " circ_face" << std::endl;
            break;
          case 1:
            file << " square_face" << std::endl;
            break;
          case 2:
            file << " plus_face" << std::endl;
            break;
          case 3:
            file << " cross_face" << std::endl;
            break;
          case 4:
            file << " diamond_face" << std::endl;
            break;
          case 5:
            file << " tup_face" << std::endl;
            break;
          case 6:
            file << " tdown_face" << std::endl;
            break;

          }
        }

        if(pointedge[i]) {
          file << pointedgewidth[i] << " " << pointedgecolor[i] << " " << ps_xdata[i] << " " << ps_ydata[i] << " " << pointsize[i] * 0.5;
          switch(pointstyle[i]) {
          case 0:
            file << " circ_edge" << std::endl;
            break;
          case 1:
            file << " square_edge" << std::endl;
            break;
          case 2:
            file << " plus_edge" << std::endl;
            break;
          case 3:
            file << " cross_edge" << std::endl;
            break;
          case 4:
            file << " diamond_edge" << std::endl;
            break;
          case 5:
            file << " tup_edge" << std::endl;
            break;
          case 6:
            file << " tdown_edge" << std::endl;
            break;

          }
        }
      }
      file << "grestore" << std::endl;

    }
  };


  ////////////////////////////////////////////////////
  /// BP_Axes Functions

  void BP_Axes::add_line(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata, std::string name) {
    add();
    last().set_data(xdata, ydata, name);

    // color, style, linewidth
    last().set_lineprops(size() - 1, 0, data_linewidth);

    // track data range
    if(size() == 1) {
      data_lim[0] = min(xdata);
      data_lim[1] = max(xdata);
      data_lim[2] = min(ydata);
      data_lim[3] = max(ydata);
    }
    else {
      double d;
      d = min(xdata);
      if(d < data_lim[0]) data_lim[0] = d;
      d = max(xdata);
      if(d > data_lim[1]) data_lim[1] = d;
      d = min(ydata);
      if(d < data_lim[2]) data_lim[2] = d;
      d = max(ydata);
      if(d > data_lim[3]) data_lim[3] = d;

    }

  };

  void BP_Axes::add_points(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata, std::string name) {
    add();
    last().set_data(xdata, ydata, name);

    // color, style, psize, pointedgesize
    last().set_pointprops(size() - 1, 0, data_pointsize, data_pointedgewidth);

    // track data range
    if(size() == 1) {
      data_lim[0] = min(xdata);
      data_lim[1] = max(xdata);
      data_lim[2] = min(ydata);
      data_lim[3] = max(ydata);
    }
    else {
      double d;
      d = min(xdata);
      if(d < data_lim[0]) data_lim[0] = d;
      d = max(xdata);
      if(d > data_lim[1]) data_lim[1] = d;
      d = min(ydata);
      if(d < data_lim[2]) data_lim[2] = d;
      d = max(ydata);
      if(d > data_lim[3]) data_lim[3] = d;

    }

  };

  void BP_Axes::add_line_wpoints(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata, std::string name) {
    add();
    last().set_data(xdata, ydata, name);

    // color, style, linewidth
    last().set_lineprops(size() - 1, 0, data_linewidth);

    // color, style, psize, pointedgesize
    last().set_pointprops(size() - 1, 0, data_pointsize, data_pointedgewidth);

    // track data range
    if(size() == 1) {
      data_lim[0] = min(xdata);
      data_lim[1] = max(xdata);
      data_lim[2] = min(ydata);
      data_lim[3] = max(ydata);
    }
    else {
      double d;
      d = min(xdata);
      if(d < data_lim[0]) data_lim[0] = d;
      d = max(xdata);
      if(d > data_lim[1]) data_lim[1] = d;
      d = min(ydata);
      if(d < data_lim[2]) data_lim[2] = d;
      d = max(ydata);
      if(d > data_lim[3]) data_lim[3] = d;

    }

  };

  void BP_Axes::set_tick_ps() {
    //axis_pos in inches
    //axis_size in inches
    //axis_lim in (data units)

    //BP_Vec<double> xtick_ps_major;
    //BP_Vec<double> ytick_ps_major;
    //BP_Vec<double> xtick_ps_minor;
    //BP_Vec<double> ytick_ps_minor;

    double fx = axis_size[0] * 72.0 / (axis_lim[1] - axis_lim[0]);		// ps pts / data units (x)
    double fy = axis_size[1] * 72.0 / (axis_lim[3] - axis_lim[2]);		// ps pts / data units (x)
    double ps_x_init = axis_pos[0] * 72.0;
    double ps_y_init = axis_pos[1] * 72.0;

    xtick_ps_major.clear();
    for(unsigned long int i = 0; i < xtick_val_major.size(); i++) {
      xtick_ps_major.add(ps_x_init + (xtick_val_major[i] - axis_lim[0])*fx);

    }

    xtick_ps_minor.clear();
    for(unsigned long int i = 0; i < xtick_val_minor.size(); i++) {
      xtick_ps_minor.add(ps_x_init + (xtick_val_minor[i] - axis_lim[0])*fx);

    }

    ytick_ps_major.clear();
    for(unsigned long int i = 0; i < ytick_val_major.size(); i++) {
      ytick_ps_major.add(ps_y_init + (ytick_val_major[i] - axis_lim[2])*fy);

    }

    ytick_ps_minor.clear();
    for(unsigned long int i = 0; i < ytick_val_minor.size(); i++) {
      ytick_ps_minor.add(ps_y_init + (ytick_val_minor[i] - axis_lim[2])*fy);

    }


  };
}

#endif // BP_Plot_CC

