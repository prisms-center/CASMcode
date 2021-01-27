#include "casm/clex/ConfigIOStrucScore.hh"

#include <boost/filesystem.hpp>
#include <functional>

#include "casm/app/AppIO.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/dataformatter/EigenDataStream.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"

namespace CASM {

namespace {
fs::path _calc_properties_path(const Configuration &config) {
  const PrimClex &primclex = config.primclex();
  return primclex.dir().calculated_properties(
      config.name(), primclex.settings().default_clex().calctype);
}
}  // namespace

namespace ConfigIO {

StrucScore::StrucScore()
    : VectorXdAttribute<Configuration>(
          "struc_score",
          "Evaluates the mapping of a configuration onto an arbitrary "
          "primitive structure, specified by its path. Allowed options are [ "
          "'basis_score' (mean-square site displacement) | 'lattice_score' "
          "(lattice deformation metric having units Angstr.^2) | 'total_score' "
          "(w*lattice_score+(1.0-w)*basis_score) ].  The struc_score weighting "
          "parameter 'w' can be provided as an optional decimal parameter from "
          "0.0 to 1.0 (default 0.5). Ex: struc_score(path/to/PRIM, "
          "basis_score, 0.4)"),
      m_strain_weight(0.5){

      };

StrucScore::StrucScore(const StrucScore &RHS)
    : VectorXdAttribute<Configuration>(RHS),
      m_altprim((RHS.m_altprim == nullptr)
                    ? nullptr
                    : new BasicStructure(*(RHS.m_altprim))),
      m_strain_weight(RHS.m_strain_weight),
      m_prim_path(RHS.m_prim_path),
      m_prop_names(RHS.m_prop_names) {
  if (RHS.m_strucmapper)
    m_strucmapper = notstd::make_unique<StrucMapper>(*RHS.m_strucmapper);
}

bool StrucScore::parse_args(const std::string &args) {
  std::vector<std::string> splt_vec;
  double _strain_weight = 0.5;
  bool already_initialized = !m_prim_path.empty();
  int pushed_args = 0;
  boost::split(splt_vec, args, boost::is_any_of(", "),
               boost::token_compress_on);
  if (splt_vec.size() < 2 || splt_vec.size() > 4) {
    throw std::runtime_error(
        "Attempted to initialize format tag " + name() + " with " +
        std::to_string(splt_vec.size()) + " arguments (" + args +
        "). You must provide at least 2 arguments, but no more than 4.\n");
    return false;
  }
  if (!m_prim_path.empty() && fs::path(splt_vec[0]) != m_prim_path)
    return false;
  if (m_prim_path.empty()) {
    m_prim_path = splt_vec[0];
    if (!fs::exists(m_prim_path)) {
      throw std::runtime_error("Attempted to initialize format tag " + name() +
                               " invalid file path '" +
                               fs::absolute(m_prim_path).string() +
                               "'. File does not exist.\n");
    }
  }
  for (Index i = 1; i < splt_vec.size(); ++i) {
    if (splt_vec[i] != "basis_score" && splt_vec[i] != "lattice_score" &&
        splt_vec[i] != "total_score") {
      try {
        _strain_weight = std::stod(splt_vec[i]);
      } catch (...) {
        throw std::runtime_error("Attempted to initialize format tag " +
                                 name() + " with invalid argument '" +
                                 splt_vec[i] +
                                 "'. Valid arguments are [ basis_score | "
                                 "lattice_score | total_score ]\n");
      }
      if (already_initialized &&
          !almost_equal(_strain_weight, m_strain_weight)) {
        for (; pushed_args > 0; pushed_args--) m_prop_names.pop_back();
        return false;
      }
    } else {
      m_prop_names.push_back(splt_vec[i]);
      ++pushed_args;
    }
  }
  return true;
}

//****************************************************************************************

/// \brief If not yet initialized, use the default clexulator from the PrimClex
bool StrucScore::init(const Configuration &_tmplt) const {
  PrimClex const &pclex(_tmplt.primclex());
  m_altprim.reset(
      new BasicStructure(read_prim(m_prim_path, _tmplt.crystallography_tol(),
                                   &(pclex.settings().hamiltonian_modules()))));

  m_strucmapper = notstd::make_unique<StrucMapper>(
      PrimStrucMapCalculator(*m_altprim), m_strain_weight);
  return true;
}

//****************************************************************************************

bool StrucScore::validate(const Configuration &_config) const {
  return fs::exists(_calc_properties_path(_config));
}

//****************************************************************************************

std::vector<std::string> StrucScore::col_header(
    const Configuration &_tmplt) const {
  std::vector<std::string> col;
  for (Index i = 0; i < m_prop_names.size(); i++) {
    std::stringstream t_ss;
    t_ss << "    " << name() << '(' << m_prim_path.string() << ','
         << m_prop_names[i] << ',' << m_strucmapper->lattice_weight() << ')';
    col.push_back(t_ss.str());
  }
  return col;
}

//****************************************************************************************

std::string StrucScore::short_header(const Configuration &_tmplt) const {
  std::stringstream t_ss;
  t_ss << name() << '(' << m_prim_path.string();
  for (Index i = 0; i < m_prop_names.size(); i++)
    t_ss << ',' << m_prop_names[i];
  t_ss << ',' << m_strucmapper->lattice_weight() << ')';
  return t_ss.str();
}

//****************************************************************************************
Eigen::VectorXd StrucScore::evaluate(const Configuration &_config) const {
  SimpleStructure relaxed_struc;

  from_json(relaxed_struc, jsonParser(_calc_properties_path(_config)));

  auto result = m_strucmapper->map_deformed_struc(relaxed_struc);
  Eigen::VectorXd result_vec(m_prop_names.size());
  if (result.empty()) {
    result_vec.setConstant(1e9);
    return result_vec;
  }

  MappingNode const &mapping(*result.begin());

  for (Index i = 0; i < m_prop_names.size(); i++) {
    if (m_prop_names[i] == "basis_score") {
      result_vec[i] = mapping.atomic_node.cost;
    } else if (m_prop_names[i] == "lattice_score") {
      result_vec[i] = mapping.lattice_node.cost;
    } else if (m_prop_names[i] == "total_score") {
      result_vec[i] = mapping.cost;
    }
  }

  return result_vec;
}

}  // namespace ConfigIO
}  // namespace CASM
