#ifndef CASM_ProjectBuilder
#define CASM_ProjectBuilder

#include <set>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
}

class ProjectSettings;
class Structure;

/// Build a CASM project directory tree and populate with standard files
///
/// Notes:
/// - The project is built at `project_settings.root_dir()`
/// - Throws if `!project_settings.has_dir()`
/// - Throws if a CASM project already exists at that location (subprojects are
/// OK)
void build_project(ProjectSettings const &project_settings,
                   Structure const &prim);

/// Construct default project settings for a given prim structure, without
/// setting a root dir
///
/// Notes:
/// - Uses `default_nlist_weight_matrix` and `default_nlist_sublat_indices`
/// - Create default clex using `default_configuration_clex`
/// - Sets required_properties for Configuration to "energy"
ProjectSettings make_default_project_settings(xtal::BasicStructure const &prim,
                                              std::string project_name);

/// Construct default project settings for a given prim structure, including a
/// root dir
///
/// Notes:
/// - Uses `default_nlist_weight_matrix` and `default_nlist_sublat_indices`
/// - Create default clex using `default_configuration_clex`
/// - Sets required_properties for Configuration to "energy"
/// - Sets root dir in project settings
ProjectSettings make_default_project_settings(xtal::BasicStructure const &prim,
                                              std::string project_name,
                                              fs::path root_dir);

/// Make default weight matrix for approximately spherical neighborhood in
/// Cartesian coordinates
///
/// Equivalent to:
/// \code
/// PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10,
/// tol()); \endcode
Eigen::Matrix3l default_nlist_weight_matrix(xtal::BasicStructure const &prim,
                                            double tol);

/// Make default list of sublattice indices that will be included in the
/// neighbor list
///
/// Includes sublattices with either:
/// - 2 or more occupants (Site::occupant_dof().size() >= 2)
/// - Continuous degrees of freedom (Site::dof_size() > 0)
std::set<int> default_nlist_sublat_indices(xtal::BasicStructure const &prim);

}  // namespace CASM

#endif
