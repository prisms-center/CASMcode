#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/DoFSpace_impl.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"

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

/// Make DoFSpace with symmetry adapated basis
DoFSpace make_symmetry_adapted_dof_space(
    DoFSpace const &dof_space, ConfigEnumInput const &input_state,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<VectorSpaceSymReport> &symmetry_report) {
  symmetry_report = vector_space_sym_report(
      dof_space, input_state, group.begin(), group.end(), calc_wedges);

  // check for error occuring for "disp"
  if (symmetry_report->symmetry_adapted_dof_subspace.cols() <
      dof_space.basis().cols()) {
    std::stringstream msg;
    msg << "Error in make_symmetry_adapted_dof_space: "
        << "symmetry_adapted_dof_subspace.cols() < dof_space.basis().cols()";
    throw make_symmetry_adapted_dof_space_error(msg.str());
  }

  return make_dof_space(dof_space.dof_key(), input_state,
                        symmetry_report->symmetry_adapted_dof_subspace);
}

}  // namespace CASM
