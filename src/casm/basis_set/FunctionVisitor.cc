#include "casm/basis_set/FunctionVisitor.hh"

#include "casm/basis_set/OccupantFunction.hh"

namespace CASM {


  OccFuncLabeler::OccFuncLabeler(const std::string &_template) {
    // parse _template into the Array m_sub_strings
    if(_template.size())
      m_sub_strings.push_back(std::string());
    for(Index i = 0; i < _template.size(); i++) {
      if(_template[i] == '%') {
        if(m_sub_strings.back().size())
          m_sub_strings.push_back(std::string());
        m_sub_strings.back() = _template.substr(i, 2);
        if((++i) + 1 < _template.size())
          m_sub_strings.push_back(std::string());
      }
      else
        m_sub_strings.back().push_back(_template[i]);
    }
    //std::cout << "substring expression: " << m_sub_strings << '\n';
  }

  //************************************************************


  bool OccFuncLabeler::visit(OccupantFunction &host)const {
    m_ss.str(std::string());
    m_ss.clear();

    for(Index i = 0; i < m_sub_strings.size(); i++) {
      if(m_sub_strings[i] == "%n") {
        if(valid_index(host.dof().ID()))
          m_ss << host.dof().ID();
        else
          m_ss << '?';
      }
      else if(m_sub_strings[i] == "%f") {
        if(valid_index(host.occ_func_ind()))
          m_ss << host.occ_func_ind();
        else
          m_ss << '?';
      }
      else if(m_sub_strings[i] == "%b") {
        if(valid_index(host.basis_ind()))
          m_ss << host.basis_ind();
        else
          m_ss << '?';
      }
      else
        m_ss << m_sub_strings[i];
    }
    //std::cout << "Paying a visit. Formula before: " << host.formula() << "\n";
    host.set_formula(m_ss.str());
    //std::cout << "                Formula after:  " << host.formula() << "\n";
    return true;
  }

  //************************************************************

  bool OccFuncBasisIndexer::visit(OccupantFunction &host)const {
    host.set_basis_ind(m_new_index);
    return true;
  };

}

