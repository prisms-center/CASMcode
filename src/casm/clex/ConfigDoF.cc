#include "casm/clex/ConfigDoF.hh"

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/crystallography/Structure.hh"


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
    m_occupation.clear();
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

  void ConfigDoF::set_occupation(const std::vector<int> &new_occupation) {

    if(occupation().size() != new_occupation.size()) {
      std::cerr << "CRITICAL ERROR: In ConfigDoF::set_occupation(), attempting to set occupation to size " << new_occupation.size() << ",\n"
                << "                which does not match initialized size of ConfigDoF -> " << size() << "\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    m_occupation = new_occupation;

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
    it->second.values() = _val;

  }
  //*******************************************************************************

  // Calculate transformed ConfigDoF from PermuteIterator via
  //   *this = permute_iterator * (*this)
  ConfigDoF &ConfigDoF::apply_sym(const PermuteIterator &it) {
    Eigen::Matrix3d fg_cart_op = it.sym_op().matrix();
    for(auto &dof : m_global_dofs) {
      dof.second.values() = *(it.global_dof_rep(dof.first).MatrixXd()) * dof.second.values();
    }

    Permutation tperm(it.combined_permute());
    if(occupation().size()) {
      set_occupation(tperm * occupation());
    }

    for(auto &dof : m_local_dofs) {
      LocalContinuousConfigDoFValues tmp = dof.second;

      for(Index b = 0; b < tmp.n_basis(); ++b)
        tmp.sublat(b) = *(it.local_dof_rep(dof.first, b).MatrixXd()) * dof.second.sublat(b);
      for(Index l = 0; l < size(); ++l) {
        dof.second.site_value(l) = tmp.site_value(tperm[l]);
      }

    }

    return *this;
  }

  //*******************************************************************************

  jsonParser &ConfigDoF::to_json(jsonParser &json) const {
    std::cout << "IN ConfigDoF::to_json()\n";
    json = jsonParser::object();
    if(occupation().size())
      json["occupation"] = occupation();
    if(!m_local_dofs.empty()) {
      json["local_dofs"] = m_local_dofs;
    }
    if(!m_global_dofs.empty()) {
      std::cout << "WRITING GLOBAL DOFS" << std::endl;
      json["global_dofs"] = m_global_dofs;
    }

    return json;
  }

  //*******************************************************************************
  void ConfigDoF::from_json(const jsonParser &json, Index NB) {
    if(!json.contains("occupation")) {
      throw std::runtime_error("JSON serialization of ConfigDoF must contain field \"occupation\"\n");
    }
    CASM::from_json(m_occupation, json["occupation"]);

    if(json.contains("local_dofs")) {
      auto end_it = json["local_dofs"].end();
      for(auto it = json["local_dofs"].begin(); it != end_it; ++it)
        m_local_dofs.emplace(it.name(), LocalValueType(DoF::traits(it.name()), NB, m_occupation.size() / NB, it->get<Eigen::MatrixXd>()));
    }

    if(json.contains("global_dofs")) {
      auto end_it = json["global_dofs"].end();
      for(auto it = json["global_dofs"].begin(); it != end_it; ++it)
        m_global_dofs.emplace(it.name(), GlobalValueType(DoF::traits(it.name()), NB, m_occupation.size() / NB, it->get<Eigen::VectorXd>()));
    }

  }

  //*******************************************************************************

  ConfigDoF jsonConstructor<ConfigDoF>::from_json(const jsonParser &json, Index NB) {

    ConfigDoF result(0, 0, std::map<DoFKey, Index>(), std::map<DoFKey, Index>(), 0.);
    result.from_json(json, NB);

    return result;
  }

  //*******************************************************************************

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json) {
    return value.to_json(json);
  }

  //*******************************************************************************

  void from_json(ConfigDoF &value, const jsonParser &json, Index NB) {
    value.from_json(json, NB);
  }

  //*******************************************************************************

  void swap(ConfigDoF &A, ConfigDoF &B) {
    A.swap(B);
  }


}

