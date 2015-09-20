#ifndef ARGS_HH
#define ARGS_HH

#include <iostream>
#include <string>
#include <sstream>

#include "casm/container/Array.hh"

namespace CASM {

  /// A class for reading command line arguments
  ///
  /// Typical Usage:
  /// For command line input such as:
  ///   casmtools -funcname -opt1 val11 val12 -opt2 val21 val22 val23
  /// This class will create a list of the 'val' associated with each 'opt'
  ///
  /// Special named '-opt' are: -PRIM, -CSPECS, -HSPECS, -LCSPECS, -STRUC, -MATRIX
  /// Also special are: -COORD (expects either 'FRAC' or 'CART')
  ///                 : -TOL (expects double)
  /// Other '-opt' are stored in 'keywords'
  ///   and the associated Array of values can be accessed using Args::get_keyvals('key')
  ///
  /// Collect command line arguments in an Array<string> object called '_args'
  /// Parse using: Args( Array<string> _args, cout);
  ///
  /// Then access using:
  ///    const std::string &get_function_name() const ;
  ///    const Array<std::string> &get_prim_filenames() const ;
  ///    const Array<std::string> &get_cspecs_filenames() const ;
  ///    const Array<std::string> &get_hspecs_filenames() const ;
  ///    const Array<std::string> &get_lcspecs_filenames() const ;
  ///    const Array<std::string> &get_structure_filenames() const ;
  ///    const Array<std::string> &get_matrix_filenames() const ;
  ///    const Array<std::string> &get_keyvals(std::string key) const ;
  ///    double get_tol() const ;
  ///    COORD_TYPE get_coord_mode() const ;
  ///
  /// If unexpected input is found, print the input so the user sees what went wrong:
  ///    void print(std::ostream &out) const ;


  class Args {
    bool all_ok;


    // ////////////////////////////////////////
    // Input file names
    Array<std::string> input_file_options;

    Array<std::string> prim_filenames;
    Array<std::string> cspecs_filenames;
    Array<std::string> hspecs_filenames;
    Array<std::string> lcspecs_filenames;
    Array<std::string> structure_filenames;
    Array<std::string> matrix_filenames;

    // collect -keyword value1 value2 ...
    Array<std::string> keywords;
    Array< Array<std::string> > values;
    Array<std::string> null_values;

    // ////////////////////////////////////////
    // Additional settings
    double tol;
    COORD_TYPE coord_mode;

    // ////////////////////////////////////////
    // Function to call (only 1 allowed at present)
    std::string function_name;

    // ////////////////////////////////////////
    // Strings I don't understand
    Array<std::string> err_list;


  public:

    // ////////////////////////////////////////
    // ////////////////////////////////////////
    // Member functions

    Args();
    Args(const Array<std::string> args, std::ostream &out);

    void set(const Array<std::string> args, std::ostream &out);

    void clear();

    // ////////////////////////////////////////
    // Accessors

    bool is_ok() const;

    const std::string &get_function_name() const ;

    const Array<std::string> &get_prim_filenames() const ;
    const Array<std::string> &get_cspecs_filenames() const ;
    const Array<std::string> &get_hspecs_filenames() const ;
    const Array<std::string> &get_lcspecs_filenames() const ;
    const Array<std::string> &get_structure_filenames() const ;
    const Array<std::string> &get_matrix_filenames() const ;
    const Array<std::string> &get_keyvals(std::string key) const ;


    double get_tol() const ;
    COORD_TYPE get_coord_mode() const ;

    // ////////////////////////////////////////
    // Mutators

    void set_function_name(const std::string &in) ;

    void set_prim_filenames(const Array<std::string> &in) ;
    void set_cspecs_filenames(const Array<std::string> &in) ;
    void set_hspecs_filenames(const Array<std::string> &in) ;
    void set_lcspecs_filenames(const Array<std::string> &in) ;
    void set_structure_filenames(const Array<std::string> &in) ;
    void set_matrix_filenames(const Array<std::string> &in) ;
    void set_keyvals(const std::string &key, const Array<std::string> &invals) ;


    Array<std::string> &get_prim_filenames() ;
    Array<std::string> &get_cspecs_filenames() ;
    Array<std::string> &get_hspecs_filenames() ;
    Array<std::string> &get_lcspecs_filenames() ;
    Array<std::string> &get_structure_filenames() ;
    Array<std::string> &get_matrix_filenames() ;

    void set_tol(double in) ;
    void set_coord_mode(COORD_TYPE in) ;

    // ////////////////////////////////////////
    // Print

    void print(std::ostream &out) const ;


  };

};

#endif
