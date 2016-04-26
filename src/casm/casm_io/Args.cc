#include "casm/casm_io/Args.hh"

namespace CASM {

  // ////////////////////////////////////////
  // ////////////////////////////////////////
  // Member functions

  Args::Args() {
    clear();
  };
  Args::Args(const Array<std::string> args, std::ostream &out) {
    clear();
    set(args, out);
  };

  void Args::set(const Array<std::string> args, std::ostream &out) {
    Index i;
    std::string arg;

    // args[0] shold be "casmtools or /path/to/casmtools"
    if(args.size() < 2) {
      all_ok = false;
      return;
    }

    // args[1] shold be function to call (-clust, -hclust, etc.)
    if(args[1][0] != '-') {
      all_ok = false;
      err_list.push_back(args[1]);
    }
    else {
      function_name = args[1];
    }

    // read the rest of the arguments
    for(i = 2; i < args.size(); i++) {
      if(args[i] == "-PRIM") {
        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          prim_filenames.push_back(args[i]);
        }
        if(prim_filenames.size() == 0)
          all_ok = false;
      }
      else if(args[i] == "-CSPECS") {
        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          cspecs_filenames.push_back(args[i]);
        }
        if(cspecs_filenames.size() == 0)
          all_ok = false;
      }
      else if(args[i] == "-HSPECS") {
        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          hspecs_filenames.push_back(args[i]);
        }
        if(hspecs_filenames.size() == 0)
          all_ok = false;
      }
      else if(args[i] == "-LCSPECS") {
        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          lcspecs_filenames.push_back(args[i]);
        }
        if(lcspecs_filenames.size() == 0)
          all_ok = false;
      }
      else if(args[i] == "-STRUC") {
        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          structure_filenames.push_back(args[i]);
        }
        if(structure_filenames.size() == 0)
          all_ok = false;
      }
      else if(args[i] == "-MATRIX") {
        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          matrix_filenames.push_back(args[i]);
        }
        if(matrix_filenames.size() == 0)
          all_ok = false;
      }
      else if(args[i] == "-TOL") {
        if(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          tol = std::stod(args[i]);
        }
        else {
          all_ok = false;
        }
      }
      else if(args[i] == "-COORD") {
        if(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          if(args[i] == "FRAC")
            coord_mode = CASM::FRAC;
          else if(args[i] == "CART")
            coord_mode = CASM::CART;
          else {
            all_ok = false;
            err_list.push_back(args[i]);
          }

        }
        else {
          all_ok = false;
        }
      }
      else if(args[i][0] == '-') { // keywords not pre-specified, "-keyword val0 val1 ..."
        Index index = keywords.find(args[i]);
        bool new_keyword = false;
        if(index == keywords.size()) {
          keywords.push_back(args[i]);
          new_keyword = true;
        }

        while(((i + 1) < args.size()) && (args[i + 1][0] != '-')) {
          i++;
          if(new_keyword) {
            values.push_back(Array<std::string>(1, args[i]));
            new_keyword = false;
          }
          else {
            values[index].push_back(args[i]);
          }
        }

      }
      else {
        all_ok = false;
        err_list.push_back(args[i]);
      }

    }

    if(!all_ok) {
      out << "ERROR. Error reading arguments." << std::endl;
      out << "Expect format: " << "casmtools -functionname -OPTION1 Val -OPTION2 Val" << std::endl;
      out << std::endl;
      out << "INPUT: ";
      for(i = 0; i < args.size(); i++)
        out << "'" << args[i] << "' ";
      out << std::endl;
      out << std::endl;
      print(out);
      exit(1);
    }
  };

  void Args::clear() {
    // set all to default

    prim_filenames.clear();
    cspecs_filenames.clear();
    hspecs_filenames.clear();
    lcspecs_filenames.clear();
    structure_filenames.clear();
    matrix_filenames.clear();

    all_ok = true;

    tol = CASM::TOL;
    coord_mode = CASM::CART;

  };

  // ////////////////////////////////////////
  // Accessors

  bool Args::is_ok() const {
    return all_ok;
  };

  const std::string &Args::get_function_name() const {
    return function_name;
  };

  const Array<std::string> &Args::get_prim_filenames() const {
    return prim_filenames;
  };
  const Array<std::string> &Args::get_cspecs_filenames() const {
    return cspecs_filenames;
  };
  const Array<std::string> &Args::get_hspecs_filenames() const {
    return hspecs_filenames;
  };
  const Array<std::string> &Args::get_lcspecs_filenames() const {
    return lcspecs_filenames;
  };
  const Array<std::string> &Args::get_structure_filenames() const {
    return structure_filenames;
  };
  const Array<std::string> &Args::get_matrix_filenames() const {
    return matrix_filenames;
  };
  const Array<std::string> &Args::get_keyvals(std::string key) const {
    Index index = keywords.find(key);
    if(index == keywords.size())
      return null_values;
    else
      return values[index];
  };


  double Args::get_tol() const {
    return tol;
  };
  COORD_TYPE Args::get_coord_mode() const {
    return coord_mode;
  };

  // ////////////////////////////////////////
  // Mutators

  void Args::set_function_name(const std::string &in) {
    function_name = in;
  };

  void Args::set_prim_filenames(const Array<std::string> &in) {
    prim_filenames = in;
  };
  void Args::set_cspecs_filenames(const Array<std::string> &in) {
    cspecs_filenames = in;
  };
  void Args::set_hspecs_filenames(const Array<std::string> &in) {
    hspecs_filenames = in;
  };
  void Args::set_lcspecs_filenames(const Array<std::string> &in) {
    lcspecs_filenames = in;
  };
  void Args::set_structure_filenames(const Array<std::string> &in) {
    structure_filenames = in;
  };
  void Args::set_matrix_filenames(const Array<std::string> &in) {
    matrix_filenames = in;
  };
  void Args::set_keyvals(const std::string &key, const Array<std::string> &invals) {
    Index index = keywords.find(key);
    if(index == keywords.size()) {
      keywords.push_back(key);
      values.push_back(invals);
    }
    else {
      for(Index i = 0; i < invals.size(); i++)
        values[index].push_back(invals[i]);
    }
  };


  Array<std::string> &Args::get_prim_filenames() {
    return prim_filenames;
  };
  Array<std::string> &Args::get_cspecs_filenames() {
    return cspecs_filenames;
  };
  Array<std::string> &Args::get_hspecs_filenames() {
    return hspecs_filenames;
  };
  Array<std::string> &Args::get_lcspecs_filenames() {
    return lcspecs_filenames;
  };
  Array<std::string> &Args::get_structure_filenames() {
    return structure_filenames;
  };
  Array<std::string> &Args::get_matrix_filenames() {
    return matrix_filenames;
  };

  void Args::set_tol(double in) {
    tol = in;
  };
  void Args::set_coord_mode(COORD_TYPE in) {
    coord_mode = in;
  };

  // ////////////////////////////////////////
  // Print

  void Args::print(std::ostream &out) const {
    out << "# Input file names" << std::endl;
    out << "prim_filenames:      " << prim_filenames << std::endl;
    out << "cspec_filenames:     " << cspecs_filenames << std::endl;
    out << "hspec_filenames:     " << hspecs_filenames << std::endl;
    out << "lcspec_filenames:    " << lcspecs_filenames << std::endl;
    out << "structure_filenames: " << structure_filenames << std::endl;
    out << "matrix_filenames:    " << matrix_filenames << std::endl;
    out << std::endl;
    out << "# Additional settings" << std::endl;
    if(coord_mode == CASM::FRAC)
      out << "coord_mode: " << "FRAC" << std::endl;
    else if(coord_mode == CASM::CART)
      out << "coord_mode: " << "CART" << std::endl;
    out << "tol: " << tol << std::endl;
    out << std::endl;
    out << "# Additional Keywords:" << std::endl;
    for(Index i = 0; i < keywords.size(); i++) {
      out << keywords[i] << ": ";
      if(values.size() >= i)
        for(Index j = 0; j < values[i].size(); j++) {
          out << values[i][j] << " ";
        }
      out << "\n" << std::flush;
    }
    out << std::endl;
    out << "# Function to call" << std::endl;
    out << "function_name: " << function_name << std::endl;
    out << std::endl;
    if(!all_ok) {
      out << "# Strings I don't understand" << std::endl;
      out << "Err: " << err_list << std::endl;
      out << std::endl;
    }
  }

};

