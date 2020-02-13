#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ConfigDoFValues.hh"

#include "casm/global/definitions.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/basis_set/Adapter.hh"


namespace CASM {

  /*  ConfigDoF::ConfigDoF(Index _N, double _tol) :
    m_occupation(_N, 0),
    m_tol(_tol) {
    }*/

  //*******************************************************************************

  // ConfigDoF::ConfigDoF(const std::vector<int> &_occ, double _tol):
  //  m_occupation(_occ),
  //  m_tol(_tol) {
  // }


  //*******************************************************************************

  void ConfigDoF::clear() {
    m_occupation.values().setZero();
    m_global_dofs.clear();
    m_local_dofs.clear();
  }

  //*******************************************************************************

  void ConfigDoF::swap(ConfigDoF &RHS) {
    std::swap(m_occupation, RHS.m_occupation);
    std::swap(m_local_dofs, RHS.m_local_dofs);
    std::swap(m_global_dofs, RHS.m_global_dofs);
    std::swap(m_tol, RHS.m_tol);
  }


  //*******************************************************************************

  void ConfigDoF::set_local_dof(DoFKey const &_key, Eigen::Ref<const Eigen::MatrixXd> const &_val) {
    auto it = m_local_dofs.find(_key);
    if(it == m_local_dofs.end()) {
      throw std::runtime_error("Attempting to assign local ConfigDoF values of type " + _key + " but no such value type has been allocated.\n");
    }
    assert(_val.cols() == it->second.values().cols() && _val.rows() == it->second.values().rows());
    it->second.values() = _val;
  }
  //*******************************************************************************

  void ConfigDoF::set_global_dof(DoFKey const &_key, Eigen::Ref<const Eigen::VectorXd> const &_val) {
    auto it = m_global_dofs.find(_key);
    if(it == m_global_dofs.end()) {
      throw std::runtime_error("Attempting to assign global ConfigDoF values of type " + _key + " but no such value type has been allocated.\n");
    }
    assert(_val.rows() == it->second.values().rows());
    //std::cout << "SETTING " << _key << " to " << _val << "\n";
    it->second.values() = _val;

  }
  //*******************************************************************************

  // Calculate transformed ConfigDoF from PermuteIterator via
  //   *this = permute_iterator * (*this)
  ConfigDoF &ConfigDoF::apply_sym(PermuteIterator const &it) {
    for(auto &dof : m_global_dofs) {
      dof.second.values() = *(it.global_dof_rep(dof.first).MatrixXd()) * dof.second.values();
    }

    Permutation tperm(it.combined_permute());
    if(occupation().size()) {
      if(it.sym_info().has_aniso_occs()) {
        Index l = 0;
        for(Index b = 0; b < n_sublat(); ++b) {
          for(Index n = 0; n < n_vol(); ++n, ++l) {
            occ(l) = (*(it.occ_rep(b).permutation()))[occ(l)];
          }
        }
      }
      set_occupation(tperm * occupation());
    }

    for(auto &dof : m_local_dofs) {
      LocalContinuousConfigDoFValues tmp = dof.second;

      for(Index b = 0; b < tmp.n_sublat(); ++b)
        tmp.sublat(b) = *(it.local_dof_rep(dof.first, b).MatrixXd()) * dof.second.sublat(b);
      for(Index l = 0; l < size(); ++l) {
        dof.second.site_value(l) = tmp.site_value(tperm[l]);
      }

    }

    return *this;
  }

  //*******************************************************************************

  // Calculate transformed ConfigDoF from PermuteIterator, using only the effect of symmetry on the value at each site
  ConfigDoF &ConfigDoF::apply_sym_no_permute(SymOp const &_op) {
    for(auto &dof : m_global_dofs) {
      dof.second.values() = *(_op.representation(dof.second.info().symrep_ID()).MatrixXd()) * dof.second.values();
    }

    if(occupation().size()) {
      Index l = 0;
      for(Index b = 0; b < n_sublat(); ++b) {
        if(!m_occupation.symrep_IDs()[b].is_identity()) {
          SymPermutation const &permrep(*_op.get_permutation_rep(m_occupation.symrep_IDs()[b]));
          l = b * n_vol();
          for(Index n = 0; n < n_vol(); ++n, ++l) {
            occ(l) = (*permrep.permutation())[occ(l)];
          }
        }
      }
    }

    for(auto &dof : m_local_dofs) {
      LocalContinuousConfigDoFValues tmp = dof.second;

      for(Index b = 0; b < tmp.n_sublat(); ++b)
        dof.second.sublat(b) = *(_op.representation(dof.second.info()[b].symrep_ID()).MatrixXd()) * dof.second.sublat(b);
    }

    return *this;
  }

  //*******************************************************************************

  jsonParser &ConfigDoF::to_json(jsonParser &json) const {
    json = jsonParser::object();
    if(occupation().size())
      json["occ"] = m_occupation;
    if(!m_local_dofs.empty()) {
      json["local_dofs"] = m_local_dofs;
    }
    if(!m_global_dofs.empty()) {
      json["global_dofs"] = m_global_dofs;
    }

    return json;
  }

  //*******************************************************************************
  void ConfigDoF::from_json(const jsonParser &json) { //, Index NB) {
    if(json.contains("occupation")) {
      //For Backwards compatibility
      CASM::from_json(m_occupation, json["occupation"]);
    }
    else if(json.contains("occ")) {
      CASM::from_json(m_occupation, json["occ"]);
    }
    else {
      throw std::runtime_error("JSON serialization of ConfigDoF must contain field \"occ\"\n");
    }


    if(json.contains("local_dofs")) {
      auto end_it = json["local_dofs"].end();
      for(auto it = json["local_dofs"].begin(); it != end_it; ++it)
        CASM::from_json(m_local_dofs[it.name()], *it);
      //.emplace(it.name(), LocalValueType(DoF::traits(it.name()), NB, m_occupation.size() / NB, (*it)["values"].get<Eigen::MatrixXd>()));
    }

    if(json.contains("global_dofs")) {
      auto end_it = json["global_dofs"].end();
      for(auto it = json["global_dofs"].begin(); it != end_it; ++it)
        CASM::from_json(m_global_dofs[it.name()], *it);
      //m_global_dofs.emplace(it.name(), GlobalValueType(DoF::traits(it.name()), NB, m_occupation.size() / NB, (*it)["values"].get<Eigen::VectorXd>()));
    }

  }

  //*******************************************************************************

  ConfigDoF jsonConstructor<ConfigDoF>::from_json(const jsonParser &json, Index NB, std::map<DoFKey, DoFSetInfo> const &global_info, std::map<DoFKey, std::vector<DoFSetInfo> > const &local_info, std::vector<SymGroupRepID> const &_occ_symrep_IDs, double _tol) {

    ConfigDoF result(NB, 0, global_info, local_info,  _occ_symrep_IDs, _tol);
    result.from_json(json);//, NB);

    return result;
  }

  //*******************************************************************************

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json) {
    return value.to_json(json);
  }

  //*******************************************************************************

  void from_json(ConfigDoF &value, const jsonParser &json) { //, Index NB) {
    value.from_json(json);
  }

  //*******************************************************************************

  void swap(ConfigDoF &A, ConfigDoF &B) {
    A.swap(B);
  }

  std::vector<SymGroupRepID> occ_symrep_IDs(xtal::Structure const &_struc) {
    _struc.generate_factor_group();
    return _struc.occupant_symrepIDs();
  }

}

