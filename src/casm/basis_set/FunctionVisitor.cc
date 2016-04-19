#include "casm/basis_set/FunctionVisitor.hh"

#include "casm/basis_set/OccupantFunction.hh"
#include "casm/basis_set/Variable.hh"
#include "casm/basis_set/PolynomialFunction.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/misc/CASM_math.hh"
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

  //*******************************************************************************************


  bool OccFuncLabeler::visit(OccupantFunction &host, BasisSet const *bset_ptr)const {
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

  //*******************************************************************************************

  VariableLabeler::VariableLabeler(const std::string &_template) {
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

  //*******************************************************************************************


  bool VariableLabeler::visit(Variable &host, BasisSet const *bset_ptr)const {

    std::stringstream tformula, ttex;
    Array<int> var_ind;

    for(Index i = 0; i < host.var_compon().size(); i++) {
      if(!almost_zero(host.coeffs()[i])) {
        var_ind.push_back(i);
      }
    }

    if(!var_ind.size()) {
      host.set_formula("0");
      host.set_tex_formula("0");
      return false;
    }

    double var_scale(host.coeffs()[var_ind[0]]);
    if(almost_zero(var_scale + 1)) {
      ttex << '-';
    }
    else if(!almost_zero(var_scale - 1)) {
      ttex << irrational_to_tex_string(var_scale, 10 * host.coeffs().size() * host.coeffs().size());
    }

    if(var_ind.size() > 1) {
      tformula << "(";
      ttex << "(";
    }

    double coeff;
    for(Index i = 0; i < var_ind.size(); i++) {
      const ContinuousDoF &var_compon(host.var_compon()[var_ind[i]]);
      coeff = host.coeffs()[var_ind[i]];

      if(i > 0 && coeff > 0) {
        tformula << '+';
      }
      if(almost_zero(coeff + 1)) {
        tformula << '-';
      }
      if(!almost_zero(std::abs(coeff) - 1)) {
        tformula << coeff;
        tformula << '*';
      }

      if(i > 0 && coeff / var_scale > 0) {
        ttex << '+';
      }
      if(almost_zero(coeff / var_scale + 1)) {
        ttex << '-';
      }
      if(!almost_zero(std::abs(coeff / var_scale) - 1)) {
        ttex << irrational_to_tex_string(coeff / var_scale, 10 * host.coeffs().size() * host.coeffs().size()) << ' ';
      }

      for(Index j = 0; j < m_sub_strings.size(); j++) {
        if(m_sub_strings[j] == "%n") {
          if(valid_index(var_compon.ID())) {
            ttex << var_compon.ID();
            tformula << var_compon.ID();
          }
          else {
            //std::cout << "type_name is " << var_compon.type_name() << ", ID is " << var_compon.ID() << "\n";
            ttex << '?';
            tformula << '?';
          }
        }
        else if(m_sub_strings[j] == "%p") {
          std::string prefix = var_compon.type_name_prefix();
          if(prefix.empty())
            prefix = "?";
          ttex << prefix;
          tformula << prefix;
        }
        else if(m_sub_strings[j] == "%s") {
          std::string suffix = var_compon.type_name_suffix();
          if(suffix.empty())
            suffix = "?";
          ttex << suffix;
          tformula << suffix;
        }
        else {
          ttex << m_sub_strings[j];
          tformula << m_sub_strings[j];
        }
      }
    }
    if(var_ind.size() > 1) {
      tformula << ")";
      ttex << ')';
    }
    host.set_tex_formula(ttex.str());
    host.set_formula(tformula.str());
    return true;
  }

  //*******************************************************************************************

  bool OccFuncBasisIndexer::visit(OccupantFunction &host, BasisSet const *bset_ptr)const {
    host.set_basis_ind(m_new_index);
    return true;
  }

  //*******************************************************************************************
  SubExpressionLabeler::SubExpressionLabeler(const std::string &_bset_name, const std::string &_template) : m_bset_name(_bset_name) {
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
    //std::cout << "substring expression: " << m_sub_strings << " and bset_name is " << m_bset_name << '\n';
  }

  //*******************************************************************************************

  bool SubExpressionLabeler::visit(Variable &host, BasisSet const *bset_ptr)const {
    return _generic_visit(host, bset_ptr);
  }

  //*******************************************************************************************

  bool SubExpressionLabeler::visit(OccupantFunction &host, BasisSet const *bset_ptr)const {
    return _generic_visit(host, bset_ptr);
  }

  //*******************************************************************************************

  bool SubExpressionLabeler::visit(PolynomialFunction &host, BasisSet const *bset_ptr)const {
    return _generic_visit(host, bset_ptr);
  }

  //*******************************************************************************************

  bool SubExpressionLabeler::_generic_visit(Function &host, BasisSet const *bset_ptr) const {

    if(bset_ptr == NULL || (bset_ptr->name()).find(m_bset_name) != 0) {
      //std::cout << "_generic_visit(): bset_ptr is " << bset_ptr << "\n";
      return false;
    }

    m_ss.str(std::string());
    m_ss.clear();

    for(Index i = 0; i < m_sub_strings.size(); i++) {
      if(m_sub_strings[i] == "%n") {
        if(bset_ptr->dof_IDs().size() == 0) {
          m_ss << '?';
        }
        else {
          for(Index i = 0; i < bset_ptr->dof_IDs().size(); i++) {
            m_ss << (bset_ptr->dof_IDs())[i];
            if(i + 1 < bset_ptr->dof_IDs().size())
              m_ss << '_';
          }
        }
      }
      else if(m_sub_strings[i] == "%f") {
        Index f = bset_ptr->find(&host);
        if(f < bset_ptr->size())
          m_ss << f;
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
}

