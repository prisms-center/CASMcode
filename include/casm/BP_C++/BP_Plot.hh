#ifndef BP_Plot_HH
#define BP_Plot_HH

#include "casm/BP_C++/BP_Vec.hh"
#include "casm/BP_C++/BP_Parse.hh"

namespace BP {

  class BP_RGB {
    double val[3];	//range 0-1

  public:

    BP_RGB() {
      val[0] = 1;
      val[1] = 1;
      val[2] = 1;
    };

    BP_RGB(double i0, double i1, double i2) {
      val[0] = i0;
      val[1] = i1;
      val[2] = i2;
    };

    BP_RGB(std::string s) {
      set(s);
    };

    BP_RGB(int i) {
      set(i);
    };

    void set(double i0, double i1, double i2) {
      val[0] = i0;
      val[1] = i1;
      val[2] = i2;

    };

    void set(std::string s) {
      if(s == "k" || s == "black")		set(0, 0, 0);			//0: (k) black
      else if(s == "b" || s == "blue")	set(0, 0, 1);		//1: (b) blue
      else if(s == "r" || s == "red")	set(1, 0, 0);		//2: (r) red
      else if(s == "g" || s == "green")	set(0, 1, 0);		//3: (g) green

      else if(s == "y" || s == "yellow")	 set(1, 1, 0);		//4: (y) yellow
      else if(s == "m" || s == "magenta") set(1, 0, 1);		//5: (m) magenta
      else if(s == "c" || s == "cyan")	 set(0, 1, 1);		//6: (c) cyan

      else if(s == "o" || s == "orange")	set(1, 0.5, 0);		//7: (o) orange
      else if(s == "l" || s == "lime")	set(0.5, 1, 0);		//8: (l) lime

      else if(s == "gray")	set(0.5, 0.5, 0.5);				//9: (gray) gray
      else if(s == "w" || s == "white")	set(1, 1, 1);		//-: (w) white
      else set(0, 0, 0);
    };

    void set(int i) {
      i = i % 10;

      if(i == 0)			set(0, 0, 0);			//0: (k) black
      else if(i == 1)	set(0, 0, 1);		//1: (b) blue
      else if(i == 2)	set(1, 0, 0);		//2: (r) red
      else if(i == 3)	set(0, 1, 0);		//3: (g) green

      else if(i == 4)	set(1, 1, 0);		//4: (y) yellow
      else if(i == 5)	set(1, 0, 1);		//5: (m) magenta
      else if(i == 6)	set(0, 1, 1);		//6: (c) cyan

      else if(i == 7)	set(1, 0.5, 0);		//7: (o) orange
      else if(i == 8)	set(0.5, 1, 0);		//8: (l) lime

      else if(i == 9)	set(0.5, 0.5, 0.5);	//9: (gray) gray

    };

    double &operator[](int i1) {
      return val[i1];
    };

    const double &operator[](int i1) const {
      return val[i1];
    };

    friend std::ostream &operator<<(std::ostream &outstream, const BP_RGB &col) {
      outstream << col[0] << " " << col[1] << " " << col[2] ;
      return outstream;
    };

    void set_color_from_map(double d, std::string map) {
      // give a value from 0 to 1, use to set the color

      if(map == "gray") {
        val[0] = val[1] = val[2] = d;
      }
      else if(map == "jet") {


        if(d < 0.125) {
          val[0] = 0.0;
          val[1] = 0.0;
          val[2] = 0.5 + 4.0 * d;
        }
        else if(d < 0.375) {
          val[0] = 0.0;
          val[1] = -0.5 + 4.0 * d;
          val[2] = 1.0;
        }
        else if(d < 0.625) {
          val[0] = -1.5 + 4.0 * d;
          val[1] = 1.0;
          val[2] = 2.5 - 4.0 * d;
        }
        else if(d < 0.875) {
          val[0] = 1.0;
          val[1] = 3.5 - 4.0 * d;
          val[2] = 0.0;
        }
        else {
          val[0] = 4.5 - 4.0 * d;
          val[1] = 0.0;
          val[2] = 0.0;
        }

        /*if( d < 0.0)
        {
        	val[0] = val[1] = val[2] = 0;
        }
        else if( d < 0.5)
        {
        	val[0] = 0.0;
        	val[1] = d;
        	val[2] = 1.0 - 2.0*d;
        }
        else if( d < 1.0)
        {
        	val[0] = 2.0*d - 1.0;
        	val[1] = 2.0 - 2.0*d;
        	val[2] = 0.0;
        }
        else
        {
        	val[0] = val[1] = val[2] = 0.9;
        }*/
      }
      else {	// unknown map
        val[0] = val[1] = val[2] = 0.0;
      }

    };

  };


  class BP_Text {

    std::string text;
    std::string fontname;
    double fontsize;
    BP_RGB textcolor;

    double	textpos[2];
    int		textpos_mode;	//0: graphic (x,y) (inches), 1: axis (x,y) (inches), 2: axis (x,y) (axis units)
    int		text_align;		//0: left, 1: center, 2: right
    double  text_angle;		// counter-clockwise

  public:

    BP_Text() {
      reset();
    };

    BP_Text(std::string s_text, std::string s_font, double d_fontsize) {
      reset();
      set_text(s_text);
      set_fontname(s_font);
      set_fontsize(d_fontsize);
    };

    void reset() {
      fontname = "Helvetica";
      fontsize = 8;
      text = "";
      textcolor.set(0, 0, 0);

      textpos_mode = 0;	//0: relative to figure, 1: relative to axes, 2: relative to axes (in axes units)
      //textpos[2];
      text_angle = 0;
      text_align = 0;
    };

    void write(BP_Write &file) {
      file << "gsave" << std::endl;
      write_format(file);
      write_text(file);
      file << "grestore" << std::endl;
    };

    void write_format(BP_Write &file) {
      file << fontname << " findfont" << std::endl;
      file << fontsize << " scalefont" << std::endl;
      file << "setfont" << std::endl;
      file << "0 0 0 setrgbcolor" << std::endl;

    };

    void write_text(BP_Write &file) {
      //double	textpos[2];
      //int		textpos_mode;	//0: graphic (x,y) (inches), 1: axis (x,y) (inches), 2: axis (x,y) (axis units)
      //int		text_align;		//0: left, 1: center, 2: right
      //double  text_angle;		// counter-clockwise

      //file << "gsave\n";
      //if( textpos_mode == 0)
      //file <<
      //str = dtos(ytick_val_major[i],5);
      //file << axis_pos[0]*72-7 << " (" << str << ") stringwidth pop sub " << ytick_ps_major[i]-((fontsize-2)*0.5) << " moveto" << std::endl;
      //file << "(" << str << ") show" << std::endl;
    };

    void set_pos(double d_pos0, double d_pos1, int i_pos_mode, int i_align, double d_angle) {
      set_textpos(d_pos0, d_pos1);
      set_textpos_mode(i_pos_mode);
      set_text_align(i_align);
      set_text_angle(d_angle);
    };

    void set_pos(double d_pos0, double d_pos1, int i_pos_mode, std::string s_align, double d_angle) {
      set_textpos(d_pos0, d_pos1);
      set_textpos_mode(i_pos_mode);
      set_text_align(s_align);
      set_text_angle(d_angle);
    };

    void set_fontname(std::string s) {
      fontname = "/" + s;
    };

    void set_fontsize(double i1) {
      fontsize = i1;
    };

    void set_text(std::string s) {
      text = s;
    };

    void set_color(BP_RGB c) {
      textcolor = c;
    };

    void set_color(double i0, double i1, double i2) {
      textcolor.set(i0, i1, i2);
    };

    void set_color(std::string s) {
      textcolor.set(s);
    };

    void set_color(int i1) {
      textcolor.set(i1);
    };

    void set_textpos_mode(int i1) {
      textpos_mode = i1;
    };

    void set_textpos(double d1, double d2) {
      textpos[0] = d1;
      textpos[1] = d2;
    };

    void set_text_angle(double d1) {
      text_angle = d1;
    };

    void set_text_align(int i1) {
      text_align = i1;
    };

    void set_text_align(std::string s) {
      if(s == "l" || s == "left" || s == "L" || s == "LEFT" || s == "Left") text_align = 0;
      else if(s == "c" || s == "center" || s == "C" || s == "CENTER" || s == "Center") text_align = 1;
      else if(s == "r" || s == "right" || s == "R" || s == "RIGHT" || s == "Right") text_align = 2;

    };

  };

  class BP_Plot_Data {
  private:
    std::string  name;

    BP_Vec<double>	xdata;
    BP_Vec<double>	ydata;

    BP_Vec<double>	ps_xdata;	// x pos in ps pts
    BP_Vec<double>	ps_ydata;	// y pos in ps pts

    bool	line;		//
    int		linestyle;	//0: line, 1: dashed, 2: dotted, 3: dash-dotted, 4: custom
    BP_Vec<double> custom_linestyle;
    double	linewidth;	// pts
    BP_RGB	linecolor;

    bool	points;
    BP_Vec<bool>	pointface;
    BP_Vec<bool>	pointedge;
    BP_Vec<int>		pointstyle;	//0: circles, 1: square, 2: plus, 3: x, 4: diamond, 5: triangle up, 6: triangle down
    BP_Vec<double>	pointsize;	// pts (diameter-ish)
    BP_Vec<BP_RGB>  pointfacecolor;
    BP_Vec<BP_RGB>  pointedgecolor;
    BP_Vec<double>	pointedgewidth;	// pts (diameter-ish)

  public:
    BP_Plot_Data() {
      reset();
    };

    void reset() {
      name = "";

      line = 0;
      linestyle = 0;
      linewidth = 0;
      linecolor.set(0);

      points = 0;

      //pointface = 1;
      //pointedge = 1;
      //pointstyle = 0;
      //pointsize = 3;
      //pointfacecolor.set(0);
      //pointedgecolor.set(0);
    };

    void set_ps_data(double *axis_pos, double *axis_size, double *axis_lim);
    void write(BP_Write &file, double *axis_pos, double *axis_size, double *axis_lim);

    void set_data(const BP_Vec<double> &xdat, const BP_Vec<double> &ydat, const std::string &nam) {
      xdata = xdat;
      ydata = ydat;
      name = nam;

    };

    // line properties
    void set_lineprops(int color, int style, double width) {
      line = 1;
      linecolor.set(color);
      linestyle = style;
      linewidth = width;
    };
    void set_lineprops(std::string color, int style, double width) {
      line = 1;
      linecolor.set(color);
      linestyle = style;
      linewidth = width;
    };
    void set_lineprops(BP_RGB color, int style, double width) {
      line = 1;
      linecolor = color;
      linestyle = style;
      linewidth = width;
    };

    void set_line(bool i1) {
      line = i1;
    };
    void set_linestyle(int i1) {
      linestyle = i1;
    };
    void set_custom_linestyle(BP_Vec<double> &i1) {
      custom_linestyle = i1;
      linestyle = 4;
    };
    void set_linewidth(double i1) {
      linewidth = i1;
    };

    void set_linecolor(double i0, double i1, double i2) {
      linecolor.set(i0, i1, i2);
    };
    void set_linecolor(std::string i1) {
      linecolor.set(i1);
    };
    void set_linecolor(int i1) {
      linecolor.set(i1);
    };
    void set_linecolor(BP_RGB i1) {
      linecolor = i1;
    };

    // point properties
    void set_pointprops(int color, int style, double psize, double pewidth) {
      points = 1;
      set_pointedge(1);
      set_pointedgecolor(color);
      set_pointfacecolor(color);
      set_pointface(0);
      set_pointstyle(style);
      set_pointsize(psize);
      set_pointedgewidth(pewidth);
    };
    void set_pointprops(std::string color, int style, double psize, double pewidth) {
      points = 1;
      set_pointedge(1);
      set_pointedgecolor(color);
      set_pointfacecolor(color);
      set_pointface(0);
      set_pointstyle(style);
      set_pointsize(psize);
      set_pointedgewidth(pewidth);
    };
    void set_pointprops(BP_RGB color, int style, double psize, double pewidth) {
      points = 1;
      set_pointedge(1);
      set_pointedgecolor(color);
      set_pointfacecolor(color);
      set_pointface(0);
      set_pointstyle(style);
      set_pointsize(psize);
      set_pointedgewidth(pewidth);
    };

    void set_points(bool i1) {
      points = i1;
    };
    void set_pointface(bool i1) {
      pointface =  BP_Vec<bool>(xdata.size(), i1);
    };
    void set_pointedge(bool i1) {
      pointedge =  BP_Vec<bool>(xdata.size(), i1);
    };
    void set_pointstyle(int i1) {
      pointstyle = BP_Vec<int>(xdata.size(), i1);
    };
    void set_pointsize(double i1) {
      pointsize = BP_Vec<double>(xdata.size(), i1);
    };
    void set_pointedgewidth(double i1) {
      pointedgewidth = BP_Vec<double>(xdata.size(), i1);
    };

    void set_pointfacecolor(double i0, double i1, double i2) {
      set_pointface(1);
      pointfacecolor = BP_Vec<BP_RGB>(xdata.size(), BP_RGB(i0, i1, i2));
    };
    void set_pointfacecolor(std::string i1) {
      set_pointface(1);
      pointfacecolor = BP_Vec<BP_RGB>(xdata.size(), BP_RGB(i1));
    };
    void set_pointfacecolor(int i1) {
      set_pointface(1);
      pointfacecolor = BP_Vec<BP_RGB>(xdata.size(), BP_RGB(i1));
    };
    void set_pointfacecolor(BP_RGB i1) {
      set_pointface(1);
      pointfacecolor = BP_Vec<BP_RGB>(xdata.size(), i1);
    };

    void set_pointedgecolor(double i0, double i1, double i2) {
      set_pointedge(1);
      pointedgecolor = BP_Vec<BP_RGB>(xdata.size(), BP_RGB(i0, i1, i2));
    };
    void set_pointedgecolor(std::string i1) {
      set_pointedge(1);
      pointedgecolor = BP_Vec<BP_RGB>(xdata.size(), BP_RGB(i1));
    };
    void set_pointedgecolor(int i1) {
      set_pointedge(1);
      pointedgecolor = BP_Vec<BP_RGB>(xdata.size(), BP_RGB(i1));
    };
    void set_pointedgecolor(BP_RGB i1) {
      set_pointedge(1);
      pointedgecolor = BP_Vec<BP_RGB>(xdata.size(), i1);
    };

    void set_pointstyle(BP_Vec<int> i1) {
      pointstyle = i1;
    };
    void set_pointsize(BP_Vec<double> i1) {
      pointsize = i1;
    };
    void set_pointedgewidth(BP_Vec<double> i1) {
      pointedgewidth = i1;
    };
    void set_pointfacecolor(BP_Vec<BP_RGB> i1) {
      set_pointface(1);
      pointfacecolor = i1;
    };
    void set_pointedgecolor(BP_Vec<BP_RGB> i1) {
      set_pointedge(1);
      pointfacecolor = i1;
    };

  };

  class BP_Axes : public BP_Vec<BP_Plot_Data> {
  private:
    double	axis_size[2];		//inches: [width][height]
    double	axis_pos[2];		//inches: [x from left edge][y from bottom edge]

    bool	auto_lim[4];		// automatically determine axis_lim [xmin][xmax][ymin][ymax]
    double	axis_lim[4];		//[xmin][xmax][ymin][ymax]
    double  data_lim[4];		//[xmin][xmax][ymin][ymax]

    bool	axis_edge_on[4];	//[left][right][bottom][top]
    bool	axis_tick_on[4];	//[left][right][bottom][top]

    bool auto_tick[2];					//x,y
    BP_Vec<double> xtick_val_major;
    BP_Vec<double> ytick_val_major;
    BP_Vec<double> xtick_val_minor;
    BP_Vec<double> ytick_val_minor;

    double		xtick_pos_major[2];		//[pts outside xaxis][pts inside xaxis]
    double		ytick_pos_major[2];		//[pts outside yaxis][pts inside yaxis]
    double		xtick_pos_minor[2];		//[pts outside xaxis][pts inside xaxis]
    double		ytick_pos_minor[2];		//[pts outside yaxis][pts inside yaxis]

    double		axis_edge_linewidth;
    double		xtick_major_linewidth;
    double		ytick_major_linewidth;
    double		xtick_minor_linewidth;
    double		ytick_minor_linewidth;

    double		target_tick_major_sep;	//inches

    double		data_linewidth;
    double		data_pointsize;
    double		data_pointedgewidth;


    BP_Vec<double> xtick_ps_major;
    BP_Vec<double> ytick_ps_major;
    BP_Vec<double> xtick_ps_minor;
    BP_Vec<double> ytick_ps_minor;


    std::string fontname;					//:"Times-Roman", "Courier", "Helvetica", etc.
    double fontsize;

    BP_RGB	axis_background_col;		//[R][G][B], background
    BP_Vec<BP_RGB>	data_col_list;		// default color order: (k)black, (b)blue, (r)red, (g)green, (c)cyan, (m)magenta, (y)yellow

  public:

    BP_Text		xticklabel;
    BP_Text		yticklabel;
    BP_Text		xaxislabel;
    BP_Text		yaxislabel;
    BP_Text		title;


    void set_tick_ps();

    BP_Axes() {
      reset();
    };

    void set_fontname(std::string s) {
      fontname = "/" + s;

    };

    void set_fontsize(double d) {
      fontsize = d;

    };

    void reset() {
      set_axis_size(4.5, 4.5);	//inches: [width][height]
      set_axis_pos(1.0, 1.0);	//inches: [x from left edge][y from bottom edge]

      set_axis_lim(0, 1, 0, 1);
      set_auto_lim();

      set_axis_edge_on(1, 1, 1, 1);
      set_auto_tick();
      set_axis_tick_on(1, 1, 1, 1);

      set_target_tick_major_sep(0.75);	//inches

      set_xtick_pos_major(0, 5);	//outside, inside
      set_xtick_pos_minor(0, 2);  //outside, inside
      set_ytick_pos_major(0, 5);  //outside, inside
      set_ytick_pos_minor(0, 2);  //outside, inside


      set_axis_edge_linewidth(0.25);
      set_tick_linewidth(0.25);

      set_data_linewidth(1);
      set_data_pointsize(3);
      set_data_pointedgewidth(1);

      set_fontname("Helvetica");
      set_fontsize(8);

      xticklabel.set_text_align("center");
      yticklabel.set_text_align("right");
      xaxislabel.set_pos(3.25, 0.25, 0, "center", 0);
      yaxislabel.set_pos(0.5, 3.25, 0, "center", 90);
      title.set_pos(3.25, 5.55, 0, "center", 0);


      axis_background_col.set("white");
      for(int i = 0; i < 10; i++) {
        data_col_list.add();
        data_col_list[i].set(i);
      }

      xtick_val_major.clear();

    };

    void get_auto_ticks(int axis, double lo, double hi, BP_Vec<double> &tick_val_major, BP_Vec<double> &tick_val_minor) {
      //std::cout << std::endl << std::endl;
      // axis=0: x-axis, 1: y-axis

      BP_Vec<double> j_list;
      j_list.add(1.0);
      j_list.add(2.0);
      j_list.add(5.0);
      j_list.add(10.0);
      j_list.add(20.0);
      j_list.add(50.0);

      //std::cout << "lo: " << lo << " hi: " << hi << std::endl;
      // major tick marks
      {
        int Ntarget = (int)(axis_size[axis] / target_tick_major_sep);
        if(Ntarget == 0) Ntarget = 1;
        double dx = hi - lo;
        if(dx == 0.0) dx = 1.0;
        double stepsize = dx / (1.0 * Ntarget);
        double xstep = 1;
        while(xstep < stepsize) xstep *= 10.0;
        while(stepsize < xstep) xstep /= 10.0;

        //std::cout << "xstep: " << xstep << std::endl;

        int j_index = 0;
        int Ninc = (int) round(dx / xstep);
        //std::cout << "  Ninc: " << Ninc << "  Ntarget: " << Ntarget << "  j_index: " << j_index << std::endl;
        int Nerr1 = std::abs(Ntarget - Ninc);
        j_index += 1;
        Ninc = (int) round(dx / (j_list[j_index] * xstep));
        //std::cout << "  Ninc: " << Ninc << "  Ntarget: " << Ntarget << "  j_index: " << j_index << std::endl;
        int Nerr2 = std::abs(Ntarget - Ninc);
        //std::cout << "Ninc: " << Ninc << "  Ntarget: " << Ntarget << "  Nerr1: " << Nerr1 << "  Nerr2: " << Nerr2 << "  j_index: " << j_index << std::endl;
        while(Nerr2 < Nerr1 && j_index + 1 < j_list.size()) {
          Nerr1 = Nerr2;
          j_index += 1;
          Ninc = (int) round(dx / (j_list[j_index] * xstep));
          Nerr2 = std::abs(Ntarget - Ninc);
          //std::cout << "Ninc: " << Ninc << "  Ntarget: " << Ntarget << "  Nerr1: " << Nerr1 << "  Nerr2: " << Nerr2 << "  j_index: " << j_index << std::endl;

        }
        //std::cout << "final  j_index: " << j_index-1 << std::endl;

        xstep = j_list[j_index - 1] * xstep;

        //std::cout << "xstep: " << xstep << std::endl;


        double xmin;
        if(lo == 0) {
          xmin = 0;
        }
        else if(lo < 0) {
          // mag of min is greatest
          xmin = -1;
          while(xmin < lo) xmin /= 10.0;
          while(lo < xmin) xmin *= 10.0;
          //while( xmin + xstep <= lo) xmin += xstep;
          while(xmin < lo) xmin += xstep;
        }
        else {
          // mag of max is greatest
          xmin = 1;
          while(xmin < lo) xmin *= 10.0;
          //std::cout << "1 xmin: " << xmin << std::endl;
          while(lo < xmin) xmin /= 10.0;
          //std::cout << "2 xmin: " << xmin << std::endl;
          //while( xmin + xstep <= lo) xmin += xstep;
          if(xmin < xstep) xmin = 0;
          while(xmin < lo) xmin += xstep;
          //std::cout << "3 xmin: " << xmin << std::endl;

        }
        //std::cout << "Ntarget: " << Ntarget << "  stepsize: " << stepsize << "  xmin: " << xmin << "  xstep: " << xstep << " j_index: " << j_index << std::endl;


        for(double d = xmin; d < hi + (hi - lo) * (1e-4); d += xstep)
          tick_val_major.add(d);


        //std::cout << "tick_val_major: " << tick_val_major << std::endl;
      }
    };

    void set_lims_and_ticks() {
      // at the end of this function the tick mark positions are in xtick_ps_major, ytick_ps_major, xtick_ps_minor, and ytick_ps_minor

      // set axis limits
      for(int i = 0; i < 4; i++)
        if(auto_lim[i])
          axis_lim[i] = data_lim[i];

      // x-axis ticks
      if(auto_tick[0]) {
        get_auto_ticks(0, axis_lim[0], axis_lim[1], xtick_val_major, xtick_val_minor);
      }

      // y-axis ticks
      if(auto_tick[1]) {
        get_auto_ticks(1, axis_lim[2], axis_lim[3], ytick_val_major, ytick_val_minor);
      }

      // set the tick positions in ps pts
      set_tick_ps();

      //std::cout << "xtick_ps_major: " << xtick_ps_major << std::endl;
      //std::cout << "ytick_ps_major: " << ytick_ps_major << std::endl;


    };

    void write_background(BP_Write &file) {
      // background

      file << "\n";
      file << "gsave" << std::endl;
      file << "newpath" << std::endl;
      file << axis_pos[0] << " inch " << axis_pos[1] << " inch moveto" << std::endl;
      file << axis_size[0] << " inch 0 inch rlineto" << std::endl;
      file << "0 inch " << axis_size[1] << " inch rlineto" << std::endl;
      file << -axis_size[0] << " inch 0 inch rlineto" << std::endl;
      file << "closepath" << std::endl;
      file << axis_background_col << " setrgbcolor" << std::endl;
      file << "fill" << std::endl;
      file << "grestore" << std::endl;

    };

    void write_edges(BP_Write &file) {
      if(axis_edge_on[0]) {	//left
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << axis_pos[0] << " inch " << axis_pos[1] << " inch moveto" << std::endl;
        file << "0 inch " << axis_size[1] << " inch rlineto" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << axis_edge_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }

      if(axis_edge_on[1]) {	//right
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << axis_pos[0] << " " << axis_size[0] << " add inch " << axis_pos[1] << " inch moveto" << std::endl;
        file << "0 inch " << axis_size[1] << " inch rlineto" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << axis_edge_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }

      if(axis_edge_on[2]) {	//bottom
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << axis_pos[0] << " inch " << axis_pos[1] << " inch moveto" << std::endl;
        file << axis_size[0] << " inch 0 inch rlineto" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << axis_edge_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }

      if(axis_edge_on[3]) {	//top
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << axis_pos[0] << " inch " << axis_pos[1] << " " << axis_size[1] << " add inch moveto" << std::endl;
        file << axis_size[0] << " inch 0 inch rlineto" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << axis_edge_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }
    };

    void write_ticks_and_ticklabels(BP_Write &file) {
      // ytick major	// left
      if(axis_edge_on[0]) {
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << ytick_major_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        for(int i = 0; i < ytick_ps_major.size(); i++) {
          file << axis_pos[0] * 72 - ytick_pos_major[0] << " " << ytick_ps_major[i] << " moveto" << std::endl;
          file << ytick_pos_major[0] + ytick_pos_major[1] << " 0 rlineto" << std::endl;

        }
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }

      // ytick major	// right
      if(axis_edge_on[1]) {
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << ytick_major_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        for(int i = 0; i < ytick_ps_major.size(); i++) {
          file << (axis_pos[0] + axis_size[0]) * 72 + ytick_pos_major[0] << " " << ytick_ps_major[i] << " moveto" << std::endl;
          file << -ytick_pos_major[0] - ytick_pos_major[1] << " 0 rlineto" << std::endl;

        }
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }

      // xtick major
      if(axis_edge_on[2]) {
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << xtick_major_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        for(int i = 0; i < xtick_ps_major.size(); i++) {
          file << xtick_ps_major[i] << " " << axis_pos[1] * 72 - xtick_pos_major[0] << " moveto" << std::endl;
          file << "0 " << xtick_pos_major[0] + xtick_pos_major[1] << " rlineto" << std::endl;

        }
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }

      // xtick major
      if(axis_edge_on[3]) {
        file << "\n";
        file << "gsave" << std::endl;
        file << "newpath" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        file << xtick_major_linewidth << " setlinewidth" << std::endl;
        file << "2 setlinecap" << std::endl;
        for(int i = 0; i < xtick_ps_major.size(); i++) {
          file << xtick_ps_major[i] << " " << (axis_pos[1] + axis_size[1]) * 72 + xtick_pos_major[0] << " moveto" << std::endl;
          file << "0 " << -xtick_pos_major[0] - xtick_pos_major[1] << " rlineto" << std::endl;

        }
        file << "stroke" << std::endl;
        file << "grestore" << std::endl;
      }


      // labels
      std::string str;

      // ytick labels // left
      if(axis_edge_on[0]) {
        file << "\n%%ytick labels\n";
        file << "gsave" << std::endl;
        file << fontname << " findfont" << std::endl;
        file << fontsize << " scalefont" << std::endl;
        file << "setfont" << std::endl;
        //file << "newpath" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        for(int i = 0; i < ytick_ps_major.size(); i++) {
          str = dtos(ytick_val_major[i], 5);
          file << axis_pos[0] * 72 - 7 << " (" << str << ") stringwidth pop sub " << ytick_ps_major[i] - ((fontsize - 2) * 0.5) << " moveto" << std::endl;
          file << "(" << str << ") show" << std::endl;

        }
        file << "grestore" << std::endl;
      }

      // xtick labels // bottom
      if(axis_edge_on[2]) {
        file << "\n";
        file << "gsave" << std::endl;
        file << fontname << " findfont" << std::endl;
        file << fontsize << " scalefont" << std::endl;
        file << "setfont" << std::endl;
        //file << "newpath" << std::endl;
        file << "0 0 0 setrgbcolor" << std::endl;
        for(int i = 0; i < xtick_ps_major.size(); i++) {
          str = dtos(xtick_val_major[i], 5);
          file << xtick_ps_major[i] << " (" << str << ") stringwidth pop 0.5 mul sub " << axis_pos[1] * 72 - 5 - fontsize << " moveto" << std::endl;
          file << "(" << str << ") show" << std::endl;

        }
        file << "grestore" << std::endl;
      }
    };

    void write_axis_labels(BP_Write &file) {

    };

    void write_data(BP_Write &file) {
      // data
      for(int i = 0; i < size(); i++)
        (*this)[i].write(file, axis_pos, axis_size, axis_lim);

    };

    void write_legend(BP_Write &file) {

    };

    void write(BP_Write &file) {
      //std::cout << "data_lim: " << data_lim[0] << " " << data_lim[1] << " " << data_lim[2] << " " << data_lim[3] << std::endl;
      set_lims_and_ticks();

      write_background(file);

      write_edges(file);

      write_ticks_and_ticklabels(file);

      write_axis_labels(file);

      write_data(file);

      write_legend(file);


    };

    void add_data(const BP_Plot_Data &data);
    void add_line_wpoints(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata);
    void add_line_wpoints(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata, std::string name);
    void add_line(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata);
    void add_line(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata, std::string name);
    void add_points(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata);
    void add_points(const BP_Vec<double> &xdata,  const BP_Vec<double> &ydata, std::string name);
    using BP_Vec<BP_Plot_Data>::operator[];
    //using BP_Vec<BP_Plot_Data>::operator[] const;
    BP_Plot_Data &operator[](const std::string &name);
    const BP_Plot_Data &operator[](const std::string &name) const;

    void set_auto_tick() {
      auto_tick[0] = 1;
      auto_tick[1] = 1;
    };
    void set_auto_xtick() {
      auto_tick[0] = 1;
    };
    void set_auto_ytick() {
      auto_tick[1] = 1;
    };
    void set_xtick_major(BP_Vec<double> i1) {
      auto_tick[0] = 0;
      xtick_val_major = i1;
    };
    void set_ytick_major(BP_Vec<double> i1) {
      auto_tick[1] = 0;
      ytick_val_major = i1;
    };

    void set_auto_lim() {	// set all to auto
      auto_lim[0] = auto_lim[1] = auto_lim[2] = auto_lim[3] = 1;
    };
    void set_auto_lim(bool i1, bool i2, bool i3, bool i4) {	// any 'true' set to auto, any 'false' leave as is
      if(i1) auto_lim[0] = i1;
      if(i2) auto_lim[1] = i2;
      if(i3) auto_lim[2] = i3;
      if(i4) auto_lim[3] = i4;
    };

    void set_axis_lim(int i, double d) {
      auto_lim[i] = 0;
      axis_lim[i] = d;
    };
    void set_axis_lim(double i0, double i1, double i2, double i3) {
      set_axis_lim(0, i0);
      set_axis_lim(1, i1);
      set_axis_lim(2, i2);
      set_axis_lim(3, i3);
    };
    void set_xaxis_lim(double i0, double i1) {
      set_axis_lim(0, i0);
      set_axis_lim(1, i1);
    };
    void set_yaxis_lim(double i2, double i3) {
      set_axis_lim(2, i2);
      set_axis_lim(3, i3);
    };

    void set_data_linewidth(double i1) {
      data_linewidth = i1;
    };
    void set_data_pointsize(double i1) {
      data_pointsize = i1;
    };
    void set_data_pointedgewidth(double i1) {
      data_pointedgewidth = i1;
    };
    void set_axis_size(double i1, double i2) {
      axis_size[0] = i1;
      axis_size[1] = i2;
    };
    void set_axis_pos(double i1, double i2) {
      axis_pos[0] = i1;
      axis_pos[1] = i2;
    };
    void get_axis_size(double &i1, double &i2) {
      i1 = axis_size[0];
      i2 = axis_size[1];
    };
    void get_axis_pos(double &i1, double &i2) {
      i1 = axis_pos[0];
      i2 = axis_pos[1];
    };

    void set_axis_edge_on(bool i0, bool i1, bool i2, bool i3) {
      axis_edge_on[0] = i0;
      axis_edge_on[1] = i1;
      axis_edge_on[2] = i2;
      axis_edge_on[3] = i3;
    };
    void set_axis_tick_on(bool i0, bool i1, bool i2, bool i3) {
      axis_tick_on[0] = i0;
      axis_tick_on[1] = i1;
      axis_tick_on[2] = i2;
      axis_tick_on[3] = i3;
    };
    void set_target_tick_major_sep(double i1) {
      target_tick_major_sep = i1;
    };

    void set_xtick_pos_major(double i0, double i1) {
      xtick_pos_major[0] = i0;
      xtick_pos_major[1] = i1;
    };
    void set_ytick_pos_major(double i0, double i1) {
      ytick_pos_major[0] = i0;
      ytick_pos_major[1] = i1;
    };
    void set_xtick_pos_minor(double i0, double i1) {
      xtick_pos_minor[0] = i0;
      xtick_pos_minor[1] = i1;
    };
    void set_ytick_pos_minor(double i0, double i1) {
      ytick_pos_minor[0] = i0;
      ytick_pos_minor[1] = i1;
    };
    void set_axis_edge_linewidth(double i0) {
      axis_edge_linewidth = i0;
    };
    void set_tick_linewidth(double i0) {
      xtick_major_linewidth = i0;
      ytick_major_linewidth = i0;
      xtick_minor_linewidth = i0;
      ytick_minor_linewidth = i0;
    };
    void set_tick_major_linewidth(double i0) {
      xtick_major_linewidth = i0;
      ytick_major_linewidth = i0;
    };
    void set_tick_minor_linewidth(double i0) {
      xtick_minor_linewidth = i0;
      ytick_minor_linewidth = i0;
    };


  };


  ///		BP_Plot is a class for creating eps plots
  ///
  ///		BP_Plot is a vector of BP_Axes. It is initialized with one BP_Axes.
  ///		Each BP_Axes is a vector of BP_Plot_Data
  ///
  ///		example:
  ///			BP_Plot plot;					// initialize a BP_Plot object with one BP_Axes
  ///			plot[0].add_line(x,y1,"y=x");	// creates a line with x-values 'x', and y-values 'y'
  ///			plot[0][0].set_linestyle(1);	// makes a dashed line
  ///			plot.write("example");			// writes the file "example.eps"
  ///
  class BP_Plot : public BP_Vec<BP_Axes> {
    //int curr_axes;
    double	graphic_size[2];	//inches: [width][height]
    BP_Vec<BP_Text>	text_list;

  public:

    BP_Plot() {
      reset();
    };

    void reset() {
      graphic_size[0] = 6;
      graphic_size[1] = 6;
      erase();
      add();
      text_list.erase();
      //curr_axes = 0;
    };

    void set_graphic_size(double x, double y) {
      graphic_size[0] = x;
      graphic_size[1] = y;
    };

    void set_graphic_margin(double m) {
      double x, y, tx, ty, maxx, maxy;

      maxx = 0;
      maxy = 0;
      for(int i = 0; i < size(); i++) {
        (*this)[i].get_axis_pos(tx, ty);
        x = tx;
        y = ty;
        (*this)[i].get_axis_size(tx, ty);
        x += tx;
        y += ty;

        if(x > maxx)
          maxx = x;

        if(y > maxy)
          maxy = y;
      }

      graphic_size[0] = maxx + m;
      graphic_size[1] = maxy + m;
    };

    void write(std::string filename) {

      int i;

      std::string epsname = filename + ".eps";

      BP_Write file(epsname);
      file.newfile();


      file << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
      file << "%%BoundingBox: 0 0 " << graphic_size[0] * 72 << " " << graphic_size[1] * 72 << std::endl;
      //file << "%!" << std::endl << std::endl;
      file << "/inch {72 mul} def" << std::endl;


      for(i = 0; i < size(); i++) {
        (*this)[i].write(file);
      }


    };

    void add_xlabel();

    void add_ylabel();

  };

}

#endif // BP_Plot_HH

