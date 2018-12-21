#include <vector>
#include <iostream>
#include <casm/casm_io/stream_io/container.hh>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  //****************************************************************************************

  std::pair<std::vector<long>, std::vector<long> > index_expression_to_bounds(const std::string &_expr) {
    //std::cout << "Parsing index expression: " << _expr << "\n";
    typedef boost::tokenizer<boost::char_separator<char> >
    tokenizer;
    boost::char_separator<char> sep1(","), sep2(" \t", ":");
    tokenizer tok1(_expr, sep1);
    std::vector<std::string> split_expr(tok1.begin(), tok1.end());
    std::vector<long> ind_vec_begin(split_expr.size());
    std::vector<long> ind_vec_end(split_expr.size());
    for(Index i = 0; i < split_expr.size(); i++) {
      tokenizer tok2(split_expr[i], sep2);
      std::vector<std::string> ind_expr(tok2.begin(), tok2.end());
      if(ind_expr.size() == 1) {
        if(ind_expr[0][0] == ':')
          ind_vec_begin[i] = -1;
        else
          ind_vec_begin[i] = boost::lexical_cast<long>(ind_expr[0]);
        ind_vec_end[i] = ind_vec_begin[i];
      }
      else if(ind_expr.size() == 3) {
        ind_vec_begin[i] = boost::lexical_cast<long>(ind_expr[0]);
        ind_vec_end[i] = boost::lexical_cast<long>(ind_expr[2]);
      }
      else
        throw std::runtime_error(std::string("In index_expression_to_bounds(), invalid index expression \"")
                                 + _expr + "\"");
    }
    //std::cout << "Lower bound: " << ind_vec_begin << "\n";
    //std::cout << "upper bound: " << ind_vec_end << "\n";
    return std::make_pair(std::move(ind_vec_begin), std::move(ind_vec_end));
  }
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
