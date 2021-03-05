#include "casm/enumerator/DoFSpace.hh"

#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/VectorSpaceSymReport.hh"

// for make_symmetry_adapted_dof_space error_report:
#include "casm/casm_io/Log.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"

namespace CASM {

namespace DoFSpace_impl {

void throw_if_missing_local_dof_requirements(
    DoFKey const &dof_key,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites) {
  if (!transformation_matrix_to_super.has_value() || !sites.has_value()) {
    std::stringstream msg;
    msg << "Error: local DoF '" << dof_key
        << "' require transformation_matrix_to_super and sites" << std::endl;
    throw std::runtime_error(msg.str());
  }
}

/// Print information in case constructing VectorSpaceSymReport fails
void error_report(DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
                  std::vector<PermuteIterator> const &group, bool calc_wedges,
                  std::optional<VectorSpaceSymReport> &symmetry_report) {
  CASM::err_log() << "prim:" << std::endl;
  COORD_TYPE mode = FRAC;
  bool include_va = false;
  jsonParser prim_json;
  write_prim(*dof_space.shared_prim(), prim_json, mode, include_va);
  CASM::err_log() << prim_json << std::endl;

  CASM::err_log() << "transformation_matrix_to_super:" << std::endl;
  CASM::err_log() << sym_info.transformation_matrix_to_super() << std::endl;

  CASM::err_log() << "supercell factor group: " << std::endl;
  SymInfoOptions opt{CART};
  brief_description(log(), sym_info.factor_group(),
                    sym_info.supercell_lattice(), opt);

  jsonParser dof_space_json;
  std::optional<std::string> identifier;
  std::optional<ConfigEnumInput> input_state;
  to_json(dof_space, dof_space_json, identifier, input_state, symmetry_report);
  CASM::err_log() << dof_space_json << std::endl;

  CASM::err_log() << "dof_space.basis().cols(): " << dof_space.basis().cols()
                  << std::endl;
  CASM::err_log() << "dof_space basis: " << std::endl
                  << dof_space.basis() << std::endl;
  if (symmetry_report.has_value()) {
    CASM::err_log() << "symmetry_adapted_dof_subspace.cols(): "
                    << symmetry_report->symmetry_adapted_dof_subspace.cols()
                    << std::endl;
    CASM::err_log() << "symmetry_adapted_dof_subspace: " << std::endl
                    << symmetry_report->symmetry_adapted_dof_subspace
                    << std::endl;
    Index i = 0;
    for (auto const &irrep : symmetry_report->irreps) {
      CASM::err_log() << "irrep[" << i << "].trans_mat: " << std::endl;
      CASM::err_log() << irrep.trans_mat << std::endl;
      i++;
    }
  }
};

/// Print information in case constructing VectorSpaceSymReport fails
void error_report_v2(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<SymRepTools_v2::VectorSpaceSymReport> &symmetry_report) {
  CASM::err_log() << "prim:" << std::endl;
  COORD_TYPE mode = FRAC;
  bool include_va = false;
  jsonParser prim_json;
  write_prim(*dof_space.shared_prim(), prim_json, mode, include_va);
  CASM::err_log() << prim_json << std::endl;

  CASM::err_log() << "transformation_matrix_to_super:" << std::endl;
  CASM::err_log() << sym_info.transformation_matrix_to_super() << std::endl;

  CASM::err_log() << "supercell factor group: " << std::endl;
  SymInfoOptions opt{CART};
  brief_description(log(), sym_info.factor_group(),
                    sym_info.supercell_lattice(), opt);

  jsonParser dof_space_json;
  std::optional<std::string> identifier;
  std::optional<ConfigEnumInput> input_state;
  to_json(dof_space, dof_space_json, identifier, input_state, symmetry_report);
  CASM::err_log() << dof_space_json << std::endl;

  CASM::err_log() << "dof_space.basis().cols(): " << dof_space.basis().cols()
                  << std::endl;
  CASM::err_log() << "dof_space basis: " << std::endl
                  << dof_space.basis() << std::endl;
  if (symmetry_report.has_value()) {
    CASM::err_log() << "symmetry_adapted_dof_subspace.cols(): "
                    << symmetry_report->symmetry_adapted_subspace.cols()
                    << std::endl;
    CASM::err_log() << "symmetry_adapted_dof_subspace: " << std::endl
                    << symmetry_report->symmetry_adapted_subspace << std::endl;
    Index i = 0;
    for (auto const &irrep : symmetry_report->irreps) {
      CASM::err_log() << "irrep[" << i << "].trans_mat: " << std::endl;
      CASM::err_log() << irrep.trans_mat << std::endl;
      i++;
    }
  }
};

}  // namespace DoFSpace_impl

/// DoFSpace constructor
///
/// \param _shared_prim The prim structure
/// \param _dof_key The DoF type
/// \param _transformation_matrix_to_super Specifies the supercell for a local
/// DoF space. Required
///        for local DoF. Ignored for global DoF.
/// \param _sites The sites included in a local DoF space. For local DoF,
/// default value if none
///        given is all sites in the supercell. Ignored for global DoF.
/// \param _basis The DoFSpace basis, as a column vector matrix. May be a
/// subspace (cols <= rows).
///        If no value, will be set to identity matrix of dimension matching
///        result of `get_dof_space_dimension`. If has value, dimension must
///        match result of `get_dof_space_dimension`.
///
DoFSpace::DoFSpace(
    std::shared_ptr<Structure const> const &_shared_prim,
    DoFKey const &_dof_key,
    std::optional<Eigen::Matrix3l> const &_transformation_matrix_to_super,
    std::optional<std::set<Index>> const &_sites,
    std::optional<Eigen::MatrixXd> const &_basis)
    : m_shared_prim(_shared_prim),
      m_dof_key(_dof_key),
      m_transformation_matrix_to_super(_transformation_matrix_to_super),
      m_sites(_sites) {
  if (AnisoValTraits(m_dof_key).global()) {
    m_transformation_matrix_to_super.reset();
    m_sites.reset();
  } else {
    if (!m_transformation_matrix_to_super.has_value()) {
      std::stringstream msg;
      msg << "Error constructing DoFSpace: Local DoF '" << m_dof_key
          << "' requires transformation_matrix_to_super." << std::endl;
      throw std::runtime_error(msg.str());
    }
    if (!m_sites.has_value()) {
      xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
          *transformation_matrix_to_super(), shared_prim()->basis().size());
      m_sites = std::set<Index>();
      for (Index i = 0; i < unitcellcoord_index_converter.total_sites(); ++i) {
        m_sites->insert(i);
      }
    }
  }

  Index dof_space_dimension = get_dof_space_dimension(
      dof_key(), *shared_prim(), transformation_matrix_to_super(), sites());

  if (!_basis.has_value()) {
    m_basis =
        Eigen::MatrixXd::Identity(dof_space_dimension, dof_space_dimension);
  } else {
    m_basis = _basis.value();
  }
  if (m_basis.rows() != dof_space_dimension) {
    std::stringstream msg;
    msg << "Error constructing DoFSpace: # basis rows (" << m_basis.rows()
        << ") != expected dimension (" << dof_space_dimension << ").";
    throw std::runtime_error(msg.str());
  }
  if (m_basis.cols() > m_basis.rows()) {
    std::stringstream msg;
    msg << "Error constructing DoFSpace: # basis columns (" << m_basis.cols()
        << ") > expected dimension (" << m_basis.rows() << ").";
    throw std::runtime_error(msg.str());
  }

  make_dof_space_axis_info(
      dof_key(), *shared_prim(), transformation_matrix_to_super(), sites(),
      m_axis_glossary, m_axis_site_index, m_axis_dof_component);
}

/// Shared prim structure
std::shared_ptr<Structure const> const &DoFSpace::shared_prim() const {
  return m_shared_prim;
}

/// Type of degree of freedom that is under consideration (e.g., "disp",
/// "Hstrain", "occ")
DoFKey const &DoFSpace::dof_key() const { return m_dof_key; }

/// Specifies the supercell for a local DoF space. Has value for local DoF.
std::optional<Eigen::Matrix3l> const &DoFSpace::transformation_matrix_to_super()
    const {
  return m_transformation_matrix_to_super;
}

/// The sites included in a local DoF space. Has value for local DoF.
std::optional<std::set<Index>> const &DoFSpace::sites() const {
  return m_sites;
}

/// True, if local DoF space with all sites in the supercell included
bool DoFSpace::includes_all_sites() const {
  return transformation_matrix_to_super().has_value() && sites().has_value() &&
         (sites()->size() == transformation_matrix_to_super()->determinant() *
                                 shared_prim()->basis().size());
}

/// The DoF space basis, as a column vector matrix. May be a subspace (cols <=
/// rows).
Eigen::MatrixXd const &DoFSpace::basis() const { return m_basis; }

/// The DoF space dimensions (equals to number of rows in basis).
Index DoFSpace::dim() const { return m_basis.rows(); }

/// The DoF subspace dimension (equal to number of columns in basis).
Index DoFSpace::subspace_dim() const { return m_basis.cols(); }

/// Names the DoF corresponding to each dimension (row) of the basis
std::vector<std::string> const &DoFSpace::axis_glossary() const {
  return m_axis_glossary;
}

/// The supercell site_index corresponding to each dimension (row) of the basis.
/// Has value for local DoF.
std::optional<std::vector<Index>> const &DoFSpace::axis_site_index() const {
  return m_axis_site_index;
}

/// The local DoF site DoFSet component index corresponding to each dimension
/// (row) of the basis. Has value for local DoF.
std::optional<std::vector<Index>> const &DoFSpace::axis_dof_component() const {
  return m_axis_dof_component;
}

/// Return true if `dof_space` is valid for `config`
///
/// Checks that:
/// - The prim are equivalent
/// - For local DoF, that the transformation_matrix_to_super are equivalent
bool is_valid_dof_space(Configuration const &config,
                        DoFSpace const &dof_space) {
  if (!AnisoValTraits(dof_space.dof_key()).global()) {
    if (config.supercell().shared_prim() != dof_space.shared_prim()) {
      return false;
    }
    if (config.supercell().sym_info().transformation_matrix_to_super() !=
        dof_space.transformation_matrix_to_super()) {
      return false;
    }
  }
  return true;
}

/// Throw if `!is_valid_dof_space(config, dof_space)`
void throw_if_invalid_dof_space(Configuration const &config,
                                DoFSpace const &dof_space) {
  if (!is_valid_dof_space(config, dof_space)) {
    std::stringstream msg;
    msg << "Error: DoFSpace is not valid for given configuration." << std::endl;
    throw std::runtime_error(msg.str());
  }
}

/// Return `config` DoF value as a coordinate in the DoFSpace basis
///
/// TODO: handle DoF values not in the basis subspace
Eigen::VectorXd get_normal_coordinate(Configuration const &config,
                                      DoFSpace const &dof_space) {
  using namespace DoFSpace_impl;
  throw_if_invalid_dof_space(config, dof_space);

  auto const &dof_key = dof_space.dof_key();
  auto const &basis = dof_space.basis();

  if (AnisoValTraits(dof_key).global()) {
    auto const &dof_values = config.configdof().global_dof(dof_key).values();
    return basis.colPivHouseholderQr().solve(dof_values);
  } else {
    if (dof_key == "occ") {
      auto const &dof_values = config.configdof().occupation().cast<double>();
      return basis.colPivHouseholderQr().solve(dof_values);
    } else {
      Eigen::VectorXd vector_values = Eigen::VectorXd::Zero(dof_space.dim());
      Eigen::MatrixXd const &matrix_values =
          config.configdof().local_dof(dof_key).values();

      auto const &axis_dof_component = dof_space.axis_dof_component().value();
      auto const &axis_site_index = dof_space.axis_site_index().value();

      for (Index i = 0; i < dof_space.dim(); ++i) {
        vector_values[i] =
            matrix_values(axis_dof_component[i], axis_site_index[i]);
      }
      return basis.colPivHouseholderQr().solve(vector_values);
    }
  }
}

/// Set `config` DoF value from a coordinate in the DoFSpace basis
void set_dof_value(Configuration &config, DoFSpace const &dof_space,
                   Eigen::VectorXd const &normal_coordinate) {
  using namespace DoFSpace_impl;
  throw_if_invalid_dof_space(config, dof_space);

  if (normal_coordinate.size() != dof_space.subspace_dim()) {
    std::stringstream msg;
    msg << "Error in set_dof_value: normal coordinate size ("
        << normal_coordinate.size() << ") != # basis axes ("
        << dof_space.subspace_dim() << ").";
    throw std::runtime_error(msg.str());
  }

  auto const &dof_key = dof_space.dof_key();
  auto const &basis = dof_space.basis();

  if (AnisoValTraits(dof_key).global()) {
    config.configdof().set_global_dof(dof_key, basis * normal_coordinate);
  } else {
    if (dof_key == "occ") {
      std::stringstream msg;
      msg << "Error: set_dof_value is not supported for occupation."
          << std::endl;
      throw std::runtime_error(msg.str());
    }

    auto &local_dof = config.configdof().local_dof(dof_key);
    Eigen::VectorXd vector_values = basis * normal_coordinate;
    Eigen::MatrixXd matrix_values = local_dof.values();

    auto const &axis_dof_component = dof_space.axis_dof_component().value();
    auto const &axis_site_index = dof_space.axis_site_index().value();

    for (Index i = 0; i < dof_space.dim(); ++i) {
      matrix_values(axis_dof_component[i], axis_site_index[i]) =
          vector_values[i];
    }

    local_dof.set_values(matrix_values);
  }
}

Index get_dof_space_dimension(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites) {
  if (AnisoValTraits(dof_key).global()) {
    return prim.global_dof(dof_key).dim();
  } else {
    using namespace DoFSpace_impl;
    throw_if_missing_local_dof_requirements(
        dof_key, transformation_matrix_to_super, sites);

    xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
        *transformation_matrix_to_super, prim.basis().size());
    Index dof_space_dimension = 0;
    for (Index site_index : *sites) {
      Index sublattice_index =
          unitcellcoord_index_converter(site_index).sublattice();
      xtal::Site const &site = prim.basis()[sublattice_index];
      if (dof_key == "occ") {
        dof_space_dimension += site.occupant_dof().size();
      } else if (site.has_dof(dof_key)) {
        dof_space_dimension += site.dof(dof_key).dim();
      }
    }
    return dof_space_dimension;
  }
}

/// The axis_glossary gives names to an un-symmetrized coordinate system
std::vector<std::string> make_axis_glossary(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites) {
  std::vector<std::string> axis_glossary;

  if (AnisoValTraits(dof_key).global()) {
    // Global DoF, axis_glossary comes straight from the DoF
    axis_glossary = component_descriptions(prim.global_dof(dof_key));
  } else {
    using namespace DoFSpace_impl;
    throw_if_missing_local_dof_requirements(
        dof_key, transformation_matrix_to_super, sites);

    xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
        *transformation_matrix_to_super, prim.basis().size());
    // Generate full axis_glossary for all active sites of the config_region
    for (Index site_index : *sites) {
      Index sublattice_index =
          unitcellcoord_index_converter(site_index).sublattice();
      xtal::Site const &site = prim.basis()[sublattice_index];
      if (dof_key == "occ") {
        for (auto const &molecule : site.occupant_dof()) {
          axis_glossary.push_back("occ[" + std::to_string(site_index + 1) +
                                  "][" + molecule.name() + "]");
        }
      } else if (site.has_dof(dof_key)) {
        std::vector<std::string> tdescs =
            component_descriptions(site.dof(dof_key));
        for (std::string const &desc : tdescs) {
          axis_glossary.push_back(desc + "[" + std::to_string(site_index + 1) +
                                  "]");
        }
      }
      if (!site.has_dof(dof_key)) continue;
    }
  }
  return axis_glossary;
}

/// The axis_glossary gives names to an un-symmetrized coordinate system
void make_dof_space_axis_info(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites,
    std::vector<std::string> &axis_glossary,
    std::optional<std::vector<Index>> &axis_site_index,
    std::optional<std::vector<Index>> &axis_dof_component) {
  if (AnisoValTraits(dof_key).global()) {
    axis_glossary.clear();
    axis_site_index = std::nullopt;
    axis_dof_component = std::nullopt;

    // Global DoF, axis_glossary comes straight from the DoF
    axis_glossary = component_descriptions(prim.global_dof(dof_key));
  } else {
    using namespace DoFSpace_impl;
    throw_if_missing_local_dof_requirements(
        dof_key, transformation_matrix_to_super, sites);

    axis_glossary.clear();
    axis_site_index = std::vector<Index>();
    axis_dof_component = std::vector<Index>();

    xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter(
        *transformation_matrix_to_super, prim.basis().size());
    // Generate full axis_glossary for all active sites of the config_region
    for (Index site_index : *sites) {
      Index sublattice_index =
          unitcellcoord_index_converter(site_index).sublattice();
      xtal::Site const &site = prim.basis()[sublattice_index];
      if (dof_key == "occ") {
        Index i = 0;
        for (auto const &molecule : site.occupant_dof()) {
          axis_glossary.push_back("occ[" + std::to_string(site_index + 1) +
                                  "][" + molecule.name() + "]");
          axis_dof_component->push_back(i);
          axis_site_index->push_back(site_index);
          ++i;
        }
      } else if (site.has_dof(dof_key)) {
        std::vector<std::string> tdescs =
            component_descriptions(site.dof(dof_key));
        Index i = 0;
        for (std::string const &desc : tdescs) {
          axis_glossary.push_back(desc + "[" + std::to_string(site_index + 1) +
                                  "]");
          axis_dof_component->push_back(i);
          axis_site_index->push_back(site_index);
          ++i;
        }
      }
      if (!site.has_dof(dof_key)) continue;
    }
  }
}

/// Make DoFSpace using `dof_key`, transformation matrix and sites from
/// `input_state`, and `basis`
DoFSpace make_dof_space(DoFKey dof_key, ConfigEnumInput const &input_state,
                        std::optional<Eigen::MatrixXd> const &basis) {
  auto const &supercell = input_state.configuration().supercell();

  return DoFSpace(supercell.shared_prim(), dof_key,
                  supercell.sym_info().transformation_matrix_to_super(),
                  input_state.sites(), basis);
}

/// Make a SymGroupRep for a DoFSpace
///
/// Usage:
/// \code
/// DoFSpace dof_space = ... ;
/// SupercellSymInfo const &sym_info = ... ;
/// std::vector<PermuteIterator> invariant_group = ... ;
/// MasterSymGroup symrep_master_group;
/// SymGroupRepID id;
/// SymGroupRep const &rep = make_dof_space_symrep(
///    dof_space, sym_info, invariant_group, symrep_master_group, id);
/// for (SymOp const &op : group) {
///   Eigen::MatrixXd matrix_rep = *(rep.MatrixXd(op));
///   // transformed_dof_values = matrix_rep * dof_values_in_dof_space_basis;
// }
/// \endcode
///
SymGroupRep const &make_dof_space_symrep(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group,
    MasterSymGroup &symrep_master_group, SymGroupRepID &id) {
  MasterSymGroup &g = symrep_master_group;
  DoFKey dof_key = dof_space.dof_key();
  xtal::BasicStructure const &prim_struc = dof_space.shared_prim()->structure();

  // Populate temporary objects for two cases
  // CASE 1: DoF is global (use prim_struc's list of global DoFs, rather than
  // relying on val_traits.is_global()
  if (prim_struc.global_dofs().count(dof_key)) {
    // Global DoF, use point group only
    SymGroup pointgroup = make_point_group(group, sym_info.supercell_lattice());
    g = make_master_sym_group(pointgroup, sym_info.supercell_lattice());

    id = g.allocate_representation();
    SymGroupRep const &rep = *(sym_info.global_dof_symrep(dof_key).rep_ptr());
    for (Index i = 0; i < pointgroup.size(); ++i) {
      Index fg_ix = pointgroup[i].index();
      g[i].set_rep(id, *rep[fg_ix]);
    }

  }
  // CASE 2: DoF is local
  else {
    if (!dof_space.sites().has_value()) {
      throw std::runtime_error(
          "Error in make_dof_space_symrep: Local DoF, but no sites");
    }
    auto group_and_ID = make_collective_dof_symrep(dof_space.sites().value(),
                                                   sym_info, dof_key, group);
    g = group_and_ID.first;
    g.is_temporary_of(group_and_ID.first);
    id = group_and_ID.second;
  }

  return g.representation(id);
}

/// Make VectorSpaceSymReport
///
/// \param dof_space DoFSpace to make VectorSpaceSymReport for
/// \param sym_info Supercell symmetry info
/// \param group Group used for vector space symmetry report
/// \param calc_wedges If true, calculate the irreducible wedges for the vector
/// space. This may take a long time.
VectorSpaceSymReport vector_space_sym_report(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges) {
  // We need a temporary mastersymgroup to manage the symmetry representation
  // for the DoF
  MasterSymGroup g;
  SymGroupRepID id;
  SymGroupRep const &rep =
      make_dof_space_symrep(dof_space, sym_info, group, g, id);

  // Generate report, based on constructed inputs
  VectorSpaceSymReport result =
      vector_space_sym_report(rep, g, dof_space.basis(), calc_wedges);
  result.axis_glossary = dof_space.axis_glossary();

  return result;
}

/// Make VectorSpaceSymReport
///
/// \param dof_space DoFSpace to make VectorSpaceSymReport for
/// \param sym_info Supercell symmetry info
/// \param group Group used for vector space symmetry report
/// \param calc_wedges If true, calculate the irreducible wedges for the vector
/// space. This may take a long time.
SymRepTools_v2::VectorSpaceSymReport vector_space_sym_report_v2(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges) {
  // We need a temporary mastersymgroup to manage the symmetry representation
  // for the DoF
  MasterSymGroup g;
  SymGroupRepID id;
  SymGroupRep const &rep =
      make_dof_space_symrep(dof_space, sym_info, group, g, id);

  // Generate report, based on constructed inputs
  SymRepTools_v2::VectorSpaceSymReport result =
      vector_space_sym_report_v2(rep, g, dof_space.basis(), calc_wedges);
  result.axis_glossary = dof_space.axis_glossary();

  return result;
}

/// Make DoFSpace with symmetry adapated basis
DoFSpace make_symmetry_adapted_dof_space(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<VectorSpaceSymReport> &symmetry_report) {
  using namespace DoFSpace_impl;

  try {
    symmetry_report =
        vector_space_sym_report(dof_space, sym_info, group, calc_wedges);
  } catch (std::exception &e) {
    error_report(dof_space, sym_info, group, calc_wedges, symmetry_report);
    CASM::err_log() << "Error constructing vector space symmetry report: "
                    << e.what() << std::endl;
    throw e;
  }
  // check for error occuring for "disp"
  if (symmetry_report->symmetry_adapted_dof_subspace.cols() <
      dof_space.basis().cols()) {
    error_report(dof_space, sym_info, group, calc_wedges, symmetry_report);
    std::stringstream msg;
    msg << "Error in make_symmetry_adapted_dof_space: "
        << "symmetry_adapted_dof_subspace.cols() < dof_space.basis().cols()";
    throw make_symmetry_adapted_dof_space_error(msg.str());
  }

  return DoFSpace(dof_space.shared_prim(), dof_space.dof_key(),
                  sym_info.transformation_matrix_to_super(), dof_space.sites(),
                  symmetry_report->symmetry_adapted_dof_subspace);
}

/// Make DoFSpace with symmetry adapated basis
DoFSpace make_symmetry_adapted_dof_space_v2(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<SymRepTools_v2::VectorSpaceSymReport> &symmetry_report) {
  using namespace DoFSpace_impl;

  try {
    symmetry_report =
        vector_space_sym_report_v2(dof_space, sym_info, group, calc_wedges);
  } catch (std::exception &e) {
    error_report_v2(dof_space, sym_info, group, calc_wedges, symmetry_report);
    CASM::err_log() << "Error constructing vector space symmetry report: "
                    << e.what() << std::endl;
    throw e;
  }
  // check for error occuring for "disp"
  if (symmetry_report->symmetry_adapted_subspace.cols() <
      dof_space.basis().cols()) {
    error_report_v2(dof_space, sym_info, group, calc_wedges, symmetry_report);
    std::stringstream msg;
    msg << "Error in make_symmetry_adapted_dof_space_v2: "
        << "symmetry_adapted_subspace.cols() < dof_space.basis().cols()";
    throw make_symmetry_adapted_dof_space_error(msg.str());
  }

  return DoFSpace(dof_space.shared_prim(), dof_space.dof_key(),
                  sym_info.transformation_matrix_to_super(), dof_space.sites(),
                  symmetry_report->symmetry_adapted_subspace);
}

/// Removes the homogeneous mode space from the local continuous DoFSpace basis.
///
/// \param dof_space DoF space to remove the homogeneous mode space from.
/// Must be a DoF space for a local continuous DoF and include all sites in the
/// supercell (`dof_space.includes_all_sites==true`), else will throw.
///
/// \returns A copy of dof_space with basis modified to remove homogeneous
/// modes.
DoFSpace exclude_homogeneous_mode_space(DoFSpace const &dof_space) {
  if (AnisoValTraits(dof_space.dof_key()).global() ||
      dof_space.dof_key() == "occ" || !dof_space.includes_all_sites()) {
    std::stringstream msg;
    msg << "Error in exclude_homogeneous_mode_space: Must be a DoF space for a "
           "local continuous degrees of freedom that includes all sites in the "
           "supercell.";
    throw std::runtime_error(msg.str());
  }

  Eigen::MatrixXd null_space =
      make_homogeneous_mode_space(dof_space).transpose().fullPivLu().kernel();

  return DoFSpace{dof_space.shared_prim(), dof_space.dof_key(),
                  dof_space.transformation_matrix_to_super(), dof_space.sites(),
                  null_space};
}

/// Make the homogeneous mode space of a local DoFSpace
///
/// \param dof_space DoF space to find the homogeneous mode space of. Should be
/// a DoF space for a local continuous DoF and include all sites in the
/// supercell.
///
/// \returns The column vector space of allowed homogeneous modes (i.e. allowed
/// rigid translations)
///
/// The prim DoF basis may not be either equal to the standard DoF basis, or
/// the same on all sublattices. The homoegeneous mode space then is limited
/// to the common part of dof_info.basis() for each site in the prim.
///
/// Ex: most typical, all displacements allowed on all sites:
///     prim_dof_info[0].basis(): [[dx], [dy], [dz]],
///     prim_dof_info[1].basis(): [[dx], [dy], [dz]]
///     Common basis: Rigid translations in [[dx], [dy], [dz]]
///
/// Ex: 2d displacements:
///     prim_dof_info[0].basis(): [[dx], [dy]],
///     prim_dof_info[1].basis(): [[dx], [dy]]
///     Common basis: Rigid translations are only allowed in [[dx], [dy]]
///
/// Ex: 2d displacements of differing orientation, 1d common basis:
///     prim_dof_info[0].basis(): [[dx, dy], [dz]],
///     prim_dof_info[1].basis(): [[-dx, dy], [dz]]
///     Common basis: Rigid translations are only allowed in [[dz]]
///
/// Ex: some fixed sites:
///     prim_dof_info[0].basis(): <not allowed>,
///     prim_dof_info[1].basis(): [[dx], [dy], [dz]]
///     Common basis: No rigid translations allowed
///
/// The common basis, in the standard DoF basis, is:
///     common_standard_basis = nullspace( Prod_i P(i) - I ),
/// where projector P(i), is:
///     P(i) = prim_dof_info[i].basis() * prim_dof_info[i].inv_basis()
///
/// If the common_standard_basis is not empty, the homogeneous mode space is
/// the column vector space constructed by transforming the
/// common_standard_basis back into the basis for each site in the DoFSpace:
///
///     [[ sites_dof_info[0].inv_basis() * common_standard_basis ],
///      [ sites_dof_info[1].inv_basis() * common_standard_basis ],
///      ...,
///      [ sites_dof_info[n_sites-1].inv_basis() * common_standard_basis ]]
///
/// Note that the lines above represent blocks equal to the dimension of the
/// DoF basis on each site.
Eigen::MatrixXd make_homogeneous_mode_space(DoFSpace const &dof_space) {
  if (AnisoValTraits(dof_space.dof_key()).global() ||
      dof_space.dof_key() == "occ" || !dof_space.includes_all_sites()) {
    std::stringstream msg;
    msg << "Error in make_homogeneous_mode_space: Must be a DoF space for a "
           "local continuous degrees of freedom that includes all sites in the "
           "supercell.";
    throw std::runtime_error(msg.str());
  }

  auto const &dof_key = dof_space.dof_key();
  auto const &prim = *dof_space.shared_prim();
  auto const &T = *dof_space.transformation_matrix_to_super();
  auto const &sites = *dof_space.sites();

  /// DoFSetInfo for each sublattice
  std::vector<DoFSetInfo> prim_dof_info = local_dof_info(prim)[dof_key];

  /// DoFSetInfo for each site in the DoFSpace with 'dof_key'
  std::vector<DoFSetInfo> sites_dof_info;
  // b: sublattice index
  // l: linear index in supercell
  // bijk: UnitCellCoord, integral site coordinates
  xtal::UnitCellCoordIndexConverter l_to_bijk(T, prim.basis().size());
  for (Index l : sites) {
    Index b = l_to_bijk(l).sublattice();
    xtal::Site const &site = prim.basis()[b];
    if (site.has_dof(dof_key)) {
      sites_dof_info.push_back(prim_dof_info[b]);
    }
  }

  // standard_dof_values = dof_info.basis() * prim_dof_values
  // dof_info.inv_basis() * standard_dof_values = prim_dof_values

  // find the common basis, in the standard DoF basis, among all sites:
  int standard_basis_dim = prim_dof_info[0].basis().rows();
  Eigen::MatrixXd I =
      Eigen::MatrixXd::Identity(standard_basis_dim, standard_basis_dim);
  Eigen::MatrixXd prod = I;
  for (auto const &sublat_dof : prim_dof_info) {
    prod = sublat_dof.basis() * sublat_dof.inv_basis() * prod;
  }
  // common_standard_basis is nullspace of (prod - I):
  Eigen::MatrixXd common_standard_basis =
      (prod - I).transpose().fullPivLu().kernel();

  if (common_standard_basis.isZero(TOL)) {
    return Eigen::MatrixXd::Zero(dof_space.basis().rows(), 0);
  }

  // construct homogeneous_mode_space by transforming common_standard_basis into
  // the site basis values, and combining for each site
  Eigen::MatrixXd homogeneous_mode_space{dof_space.basis().rows(),
                                         common_standard_basis.cols()};

  Index row = 0;  // block starting row
  Index col = 0;  // block starting column
  for (auto const &site_dof : sites_dof_info) {
    auto const &values = site_dof.inv_basis() * common_standard_basis;
    int n_rows = values.rows();
    int n_cols = values.cols();
    homogeneous_mode_space.block(row, col, n_rows, n_cols) = values;
    row += site_dof.dim();
  }

  return homogeneous_mode_space;
}

VectorSpaceMixingInfo::VectorSpaceMixingInfo(
    Eigen::MatrixXd const &column_vector_space, Eigen::MatrixXd const &subspace,
    double tol) {
  // Make a projection operator out of homogeneous mode space and project each
  // of the basis vectors onto it If they have a partial projection (not full or
  // zero) => translational modes are mixed between irreps
  Eigen::MatrixXd proj_operator = subspace * subspace.transpose();
  for (Index i = 0; i < column_vector_space.cols(); ++i) {
    Eigen::VectorXd col_projection = proj_operator * column_vector_space.col(i);
    if (col_projection.isZero(tol)) {
      axes_not_in_subspace.push_back(i);
    } else if (CASM::almost_equal(col_projection.normalized(),
                                  column_vector_space.col(i).normalized(),
                                  tol)) {
      axes_in_subspace.push_back(i);
    }

    else {
      axes_mixed_with_subspace.push_back(i);
    }
  }

  if (axes_mixed_with_subspace.size() == 0) {
    are_axes_mixed_with_subspace = false;
  }
}
}  // namespace CASM
