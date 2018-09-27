#include "casm/clex/ConfigDoF.hh"

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/crystallography/Structure.hh"


namespace CASM {

  ConfigDoF::ConfigDoF(Index _N, double _tol) :
    m_occupation(_N, 0),
    m_tol(_tol) {
  }

  //*******************************************************************************

  ConfigDoF::ConfigDoF(const std::vector<int> &_occ, double _tol):
    m_occupation(_occ),
    m_tol(_tol) {
  }


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
  /*
  // Calculate transformed ConfigDoF from PermuteIterator via
  //   *this = permute_iterator * (*this)
  ConfigDoF &ConfigDoF::apply_sym(const PermuteIterator &it) {
    Eigen::Matrix3d fg_cart_op = it.sym_op().matrix();
    //if(has_deformation()) {
      //set_deformation(fg_cart_op * deformation() * fg_cart_op.transpose());
      //}
    Permutation tperm(it.combined_permute());
    if(occupation().size()) {
      set_occupation(tperm * occupation());
    }
    //if(displacement().cols()) {
      //Eigen::MatrixXd new_disp = fg_cart_op * displacement();
      //set_displacement(Eigen::MatrixXd(3, size()));
      //for(Index i = 0; i < size(); i++)
        //disp(i) = new_disp.col(tperm[i]);
        }

    return *this;
  }
  */
  //*******************************************************************************

  jsonParser &ConfigDoF::to_json(jsonParser &json) const {
    json = jsonParser::object();
    if(occupation().size())
      json["occupation"] = occupation();
    if(!m_local_dofs.empty()) {
      json["local_dofs"] = m_local_dofs;
    }
    if(!m_global_dofs.empty()) {
      json["global_dofs"] = m_global_dofs;
    }

    return json;
  }

  //*******************************************************************************
  void ConfigDoF::from_json(const jsonParser &json) {
    if(!json.contains("occupation")) {
      throw std::runtime_error("JSON serialization of ConfigDoF must contain field \"occupation\"\n");
    }
    CASM::from_json(m_occupation, json["occupation"]);

    if(json.contains("local_dofs")) {
      auto end_it = json["local_dofs"].end();
      for(auto it = json["local_dofs"].begin(); it != end_it; ++it)
        m_local_dofs.emplace(it.name(), LocalValueType(DoF::traits(it.name()), it->get<Eigen::MatrixXd>()));
    }

    if(json.contains("global_dofs")) {
      auto end_it = json["global_dofs"].end();
      for(auto it = json["global_dofs"].begin(); it != end_it; ++it)
        m_global_dofs.emplace(it.name(), GlobalValueType(DoF::traits(it.name()), it->get<Eigen::VectorXd>()));
    }

  }

  //*******************************************************************************

  ConfigDoF jsonConstructor<ConfigDoF>::from_json(const jsonParser &json) {

    ConfigDoF result(0);
    result.from_json(json);

    return result;
  }

  //*******************************************************************************

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json) {
    return value.to_json(json);
  }

  //*******************************************************************************

  void from_json(ConfigDoF &value, const jsonParser &json) {
    value.from_json(json);
  }

  //*******************************************************************************

  void swap(ConfigDoF &A, ConfigDoF &B) {
    A.swap(B);
  }


}

