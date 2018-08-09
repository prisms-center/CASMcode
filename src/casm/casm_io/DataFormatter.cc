#include <vector>
#include <string>
#include "casm/external/boost.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  //****************************************************************************************

  std::string::const_iterator end_of_literal(std::string::const_iterator it, std::string::const_iterator end_it) {

    for(++it; it != end_it; ++it) {
      if((*it) == '\'')
        return ++it;
    }
    return it;
  }

  //****************************************************************************************
  /*
   * break 'input' string into a list of format tags and their (option) arguments
   */
  void split_formatter_expression(const std::string &input_expr,
                                  std::vector<std::string> &tag_names,
                                  std::vector<std::string> &sub_exprs) {
    std::string::const_iterator it(input_expr.cbegin()), it_end(input_expr.cend()), t_it1, t_it2;
    while(it != it_end) {
      while(it != it_end && (isspace(*it) || boost::is_any_of(",#")(*it)))
        ++it;
      if(it == it_end)
        break;
      // Identified a formatter tag, save starting iterator
      t_it1 = it;
      if((*it) == '\'') {
        it = end_of_literal(it, it_end);
        t_it2 = it;
        --t_it2;
        if(t_it2 == t_it1 || (*t_it2) != '\'') {
          throw std::runtime_error("Mismatched quotation marks in expression:\n    \"" + input_expr + "\"\n");
        }
      }
      else {
        // find end of formatter tag
        while(it != it_end && !(isspace(*it) || (*it) == ',') && (*it) != '(' && (*it) != ')')
          ++it;

        if(*it == ')')
          throw std::runtime_error("Mismatched parentheses in expression:\n    \"" + input_expr + "\"\n");
      }
      // push_back formatter tag, and push_back an empty string for its optional arguments
      tag_names.push_back(std::string(t_it1, it));
      sub_exprs.push_back(std::string());


      // no argument, we've reached beginning of new tag, so continue
      if(it == it_end || (*it) != '(')
        continue;

      // from here on, we're parsing arguments:
      while(it != it_end && isspace(*(++it)));//skipspace

      if(*it == ',')
        throw std::runtime_error("Stray comma in expression:\n    \"" + input_expr + "\"\n");

      // start of argument
      t_it1 = it;

      // stop at end of string, or if we hit the non-nested closing paren
      Index paren_depth(0);
      while(it != it_end && !((*it) == ')' && paren_depth == 0)) {
        if((*it) == '(')
          ++paren_depth;
        else if((*it) == ')')
          --paren_depth;

        ++it;
      }

      if(it == it_end)
        throw std::runtime_error("Mismatched parentheses in expression:\n    \"" + input_expr + "\"\n");

      t_it2 = it;
      while(isspace(*(--t_it2)));//remove trailing space

      if(*t_it2 == ',')
        throw std::runtime_error("Stray comma in expression:\n    \"" + input_expr + "\"\n");

      sub_exprs.back() = std::string(t_it1, ++t_it2);

      ++it;
    }

  }
}
