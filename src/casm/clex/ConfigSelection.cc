#include "casm/clex/ConfigSelection.hh"

#include <regex>

#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"

namespace CASM {

  template<>
  ConfigSelection<false>::ConfigSelection(typename ConfigSelection<false>::PrimClexType &_primclex)
    : m_primclex(&_primclex) {
    m_col_headers.clear();
    for(auto it = _primclex.config_begin(); it != _primclex.config_end(); ++it) {
      m_config[it->name()] = it->selected();
    }
  }

  template<>
  ConfigSelection<true>::ConfigSelection(typename ConfigSelection<true>::PrimClexType &_primclex)
    : m_primclex(&_primclex) {
    m_col_headers.clear();
    for(auto it = _primclex.config_cbegin(); it != _primclex.config_cend(); ++it) {
      m_config[it->name()] = it->selected();
    }
  }

  //***********************************************************
  /**
   *   Set 'selected?' to 'true' for configurations that meet 'criteria'.
   *
   *   Used in 'casm select --set'.
   *
   *   'criteria' is a list of operators and arguments
   *      listed in reverse polish notation:
   *
   *   operators:
   *       on: *if this is first, then set matching to 'true'
   *       off: *if this is first, then set matching to 'false'
   *       re: regular expression matching
   *       eq: ==
   *       ne: !=
   *       lt: <
   *       le: <=
   *       gt: >
   *       ge: >=
   *       also: add, sub, div, mult, pow
   *       also: AND, OR, XOR, NOT
   *    arguments:
   *       "scell", "config", "configid", "calculated", "groundstate" (mol name), (string), (double)
   *
   *    example:
   *       casm select --set on
   *         - use to set 'selected?'=true for all configuations
   *         - calls set_selection with criteria = ['on']
   *
   *       casm select --set off scell SCEL3 re
   *         - use to set 'selected?'=false for all volume 3 configuations
   *         - calls set_selection with criteria = ['scell', 'SCEL3', 're']
   *
   *
   *       casm select --set on scell SCEL3 re A le 0.5 AND
   *         - use to set 'selected?'=true for all volume 3 configuations with A concentration less than or equal to 0.5
   *         - calls set_selection with criteria = ['scell', 'SCEL3', 're', 'A', 'lt', '0.501', 'AND']
   *
   */
  //***********************************************************
  bool get_selection(const Array<std::string> &criteria, const Configuration &config, bool init_selection) {

    using namespace ConfigSelection_impl;

    if(criteria.size() == 0) {
      std::cerr << "Error in Supercell::set_selection(const Array<std::string> &criteria, Configuration &config)" << std::endl;
      std::cerr << "  criteria.size() must be > 0." << std::endl;
      exit(1);
    }

    // The first 'criteria' should be whether to select or unselect if the expressions comes out as true ("1")
    bool mk;
    if(criteria[0] == "on")
      mk = true;
    else if(criteria[0] == "off")
      mk = false;
    else {
      std::cerr << "Error in Supercell::set_selection(const Array<std::string> &criteria, Configuration &config)" << std::endl;
      std::cerr << "  criteria[0] must be \"on\" or \"off\", but you gave \"" << criteria[0] << "\"" << std::endl;
      exit(1);
    }

    std::string q;
    std::string A, B;
    Array<std::string> stack;


    // If just 'on' or 'off', then select or unselect all
    if(criteria.size() == 1) {
      return mk;
    }

    // Evaluate the expression by operating on the values in the 'criteria' list one by one
    //   For each value:
    //     check if it is an operator, if true: evaluate it and add result to stack
    //     else check if it is a variable, if it is: substitute in the appropriate value
    //   By the end of the 'criteria' list, the 'stack' should have only one entry which is "1" or "0"
    stack.clear();
    for(Index j = 1; j < criteria.size(); j++) {
      q = criteria[j];
      if(is_operator(q)) {
        if(is_unary(q)) {
          A = stack.back();
          stack.pop_back();
          stack.push_back(operate(q, A));
        }
        else {
          B = stack.back();
          stack.pop_back();
          A = stack.back();
          stack.pop_back();
          stack.push_back(operate(q, A, B));
        }
      }
      else {
        stack.push_back(convert_variable(q, config));
      }
    }

    if(stack.size() != 1) {
      std::cerr << "Error in Supercell::set_selection(const Array<std::string> &criteria, Configuration &config)" << std::endl;
      std::cerr << "  stack.size() != 1, check your criteria." << std::endl;
      std::cerr << "  stack: " << stack << std::endl;
      std::cerr << "  criteria: " << criteria << std::endl;
      exit(1);
    }

    if(stack[0] == "1") {
      return mk;
    }
    else if(stack[0] != "0") {
      std::cerr << "Error in Supercell::set_selection(const Array<std::string> &criteria, Configuration &config)" << std::endl;
      std::cerr << "  stack[0]: " << stack[0] << std::endl;
      exit(1);
    }
    return init_selection;
  }

  namespace ConfigSelection_impl {

    //*******************************************************************************

    bool is_operator(const std::string &q) {
      //std::cout << "is_operator: " << q << std::endl;

      if(q == "NOT")
        return true;
      if(q == "AND")
        return true;
      if(q == "OR")
        return true;
      if(q == "XOR")
        return true;
      if(q == "re")
        return true;
      if(q == "rs")
        return true;
      if(q == "eq")
        return true;
      if(q == "ne")
        return true;
      if(q == "lt")
        return true;
      if(q == "le")
        return true;
      if(q == "gt")
        return true;
      if(q == "ge")
        return true;
      if(q == "add")
        return true;
      if(q == "sub")
        return true;
      if(q == "div")
        return true;
      if(q == "mult")
        return true;
      if(q == "pow")
        return true;
      return false;
    }

    //*******************************************************************************

    std::string operate(const std::string &q, const std::string &A) {
      //std::cout << "operate1: " << q << "  A: " << A << std::endl;

      try {
        if(q == "NOT")
          return A == "0" ? "1" : "0";
        throw std::invalid_argument(std::string("  q: ") + q + " is not recognized");
      }
      catch(...) {
        std::cerr << "Error performing operation: '" << A << " " << q << "'." << std::endl;
        throw;
      }
    }

    //*******************************************************************************

    std::string operate(const std::string &q, const std::string &A, const std::string &B) {
      //std::cout << "operate2: " << q << "  A: " << A << "  B: " << B << std::endl;

      try {
        if(q == "AND")
          return (A == "0" || B == "0") ? "0" : "1";
        if(q == "OR")
          return (A == "0" && B == "0") ? "0" : "1";
        if(q == "XOR")
          return ((A == "0") != (B == "0")) ? "0" : "1";
        if(q == "re") {
          std::regex e(B);
          //std::cout << "A: " << A << "  B: " << B << "  regex_match: " << boost::regex_match(A, e) << std::endl;
          return std::regex_match(A, e) ? "1" : "0";
        }
        if(q == "rs") {
          std::regex e(B);
          //std::cout << "A: " << A << "  B: " << B << "  regex_search: " << boost::regex_search(A, e) << std::endl;
          return std::regex_search(A, e) ? "1" : "0";
        }
        if(q == "eq")
          return (A == B) ? "1" : "0";
        if(q == "ne")
          return (A != B) ? "1" : "0";
        if(q == "lt")
          return (std::stod(A) < std::stod(B)) ? "1" : "0";
        if(q == "le")
          return (std::stod(A) <= std::stod(B)) ? "1" : "0";
        if(q == "gt")
          return (std::stod(A) > std::stod(B)) ? "1" : "0";
        if(q == "ge")
          return (std::stod(A) >= std::stod(B)) ? "1" : "0";
        if(q == "add")
          return std::to_string(std::stod(A) + std::stod(B));
        if(q == "sub")
          return std::to_string(std::stod(A) - std::stod(B));
        if(q == "div")
          return std::to_string(std::stod(A) / std::stod(B));
        if(q == "mult")
          return std::to_string(std::stod(A) * std::stod(B));
        if(q == "pow")
          return std::to_string(pow(std::stod(A), std::stod(B)));
        throw std::invalid_argument(std::string("  q: ") + q + " is not recognized");
      }
      catch(...) {
        std::cerr << "Error performing operation: '" << A << " " << B << " " << q << "'." << std::endl;
        throw;
      }

    }

    //*******************************************************************************

    bool is_unary(const std::string &q) {
      return q == "NOT" ? true : false;
    }

    //***********************************************************
    /**  convert variables
     *
     *   if q is:
     *      "scelname", returns config.get_supercell().get_name()
     *      "configname", returns config.name()
     *      the name of a molecule allowed in the prim,
     *        returns that molecule's 'true_composition'
     *      "calculated", returns:
     *        config.calc_properties().contains(get_primclex().get_curr_property().begin(), get_primclex().get_curr_property().end())
     *   else:
     *      not a variable, returns q
     */
    //***********************************************************
    std::string convert_variable(const std::string &q, const Configuration &config) {

      const Supercell &scel = config.get_supercell();
      const PrimClex &primclex = scel.get_primclex();

      // check for a variable
      if(q == "scelname")
        return scel.get_name();
      if(q == "configname")
        return config.name();
      if(q == "scel_size")
        return std::to_string(scel.volume());
      if(q == "is_groundstate") {
        if(!config.generated_properties().contains("is_groundstate")) {
          return "unknown";
        }
        bool is_groundstate;
        config.generated_properties()["is_groundstate"].get(is_groundstate);
        return is_groundstate ? "1" : "0";
      }
      if(q == "is_calculated") {
        return std::all_of(primclex.get_curr_property().begin(),
                           primclex.get_curr_property().end(),
        [&](const std::string & key) {
          return config.calc_properties().contains(key);
        }) ? "1" : "0";
      }
      if(q == "dist_from_hull") {
        if(!config.generated_properties().contains("dist_from_hull")) {
          //std::cerr << "WARNING: Configuration " << config.get_path() << std::endl
          //          << "         Does not have an initialized 'dist_from_hull' field in its generated properties." << std::endl;
          return "unknown";
        }
        else
          return config.generated_properties()["dist_from_hull"].get<std::string>();
      }
      if(q == "formation_energy") {
        if(!config.delta_properties().contains("relaxed_energy")) {
          //std::cerr << "WARNING: Configuration " << config.get_path() << std::endl
          //          << "         Does not have an initialized 'relaxed_energy' field in its delta properties." << std::endl
          //          << "         Please verify that all desired properties are specified in the project settings." << std::endl;
          return "unknown";
        }
        else
          return config.delta_properties()["relaxed_energy"].get<std::string>();
      }

      // try matching 'clex(clex_name)', 'comp(x)', 'true_comp(x)', 'mol_comp(x)'
      std::smatch sm;

      // scalar cluster expansion property
      auto clex_e = std::regex("clex\\((.*)\\)");
      std::regex_match(q, sm, clex_e);
      if(sm.size()) {
        std::string ss = sm[1];
        Clexulator clexulator = primclex.global_clexulator();
        ECIContainer eci = primclex.global_eci(ss);
        return std::to_string(eci * correlations(config, clexulator));
      }

      // parametric composition
      auto comp_e = std::regex("comp\\((.*)\\)");
      std::regex_match(q, sm, comp_e);
      if(sm.size()) {
        std::string ss = sm[1];
        if(ss.size() > 1) {
          throw std::runtime_error(
            std::string("Error in selecting: '") + q + "'.\n" +
            "  Composition '" + ss + "' is not valid.");
        }
        int index = ((int) ss[0]) - ((int) 'a');
        int Nind = primclex.composition_axes().independent_compositions();
        if(index < 0) {
          throw std::runtime_error(
            std::string("Error in selecting: '") + q + "'.\n" +
            "  Composition '" + ss + "' is not valid.");
        }
        if(index >= Nind) {
          throw std::runtime_error(
            std::string("Error in selecting: '") + q + "'.\n" +
            "  looking for '" + ss[0] + "', with composition index: " + std::to_string(index) +
            "  but, # independent compositions: " + std::to_string(Nind));
        }
        return std::to_string(config.get_param_composition()[index]);
      }

      Array<Molecule> struc_molecule = primclex.get_prim().get_struc_molecule();

      // site fraction i.e. include vacancies in the count
      auto site_frac_e = std::regex("site_frac\\((.*)\\)");
      std::regex_match(q, sm, site_frac_e);
      if(sm.size()) {
        std::string ss = sm[1];
        for(int i = 0; i < struc_molecule.size(); i++)
          if(struc_molecule[i].name == ss) {
            return std::to_string(config.get_true_composition()[i]);
          }
        throw std::runtime_error(
          std::string("Error in selecting: '") + q + "'.\n" +
          "  Attempting to get site fraction, but could not find atom '" + ss + "'");
      }

      // mole fraction i.e. do not include vacancies in the count
      auto atom_frac_e = std::regex("atom_frac\\((.*)\\)");
      std::regex_match(q, sm, atom_frac_e);
      if(sm.size()) {
        std::string ss = sm[1];
        for(int i = 0; i < struc_molecule.size(); i++)
          if(struc_molecule[i].name == ss)
            return std::to_string(config.get_composition()[i]);
        throw std::runtime_error(
          std::string("Error in selecting: '") + q + "'.\n" +
          "  Attempting to get mole fraction, but could not find atom '" + ss + "'");
      }

      // else not a variable:
      return q;
    }

  }

}
