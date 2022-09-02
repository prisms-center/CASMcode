#ifndef CASM_DoFSpace
#define CASM_DoFSpace

#include <optional>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

struct VectorSpaceSymReport;

namespace SymRepTools_v2 {
struct VectorSpaceSymReport;
}

/// DoFSpace
///
/// The DoFSpace class is used to specify a particular degree of freedom space.
/// DoFSpace have:
///   - shared_prim (std::shared_ptr<Structure const>): The prim structure
///   - dof_key (DoFKey): A string indicating which DoF type (e.g., "disp",
///   "Hstrain", "occ")
///   - transformation_matrix_to_super (std::optional<Eigen::Matrix3l>):
///   Specifies the supercell for a local DoF space. Has value for local DoF.
///   - sites (std::optional<std::set<Index>>): The sites included in a local
///   DoF space. Has value for local DoF.
///   - basis (Eigen::MatrixXd): Allows specifying a subspace of the space
///   determined from dof_key, and for local DoF, transformation_matrix_to_super
///   and sites. The rows of `basis` correspond to prim DoF basis axes, so the
///   following relations apply:
///
///       standard_dof_values = dof_set.basis() * prim_dof_values
///       prim_dof_values = dof_space.basis() * normal_coordinate
///
///    Examples:
///
///       For dof_key=="disp", and sites.value().size()==4 (any
///       transformation_matrix_to_super):
///         The default DoF space has dimension 12 corresponding to (x,y,z) for
///         each of the 4 sites. Then, basis is a matrix with 12 rows and 12 or
///         fewer columns.
///
///       For dof_key=="occ", and sites().value().size()==4 (any
///       transformation_matrix_to_super), with the particular sites selected
///       having allowed occupations ["A", "B"], ["A", "B", "C"], ["A, B"], and
///       ["A", "B", "C"] :
///         The default DoF space has dimension 10 = 2 + 3 + 2 + 3.
///         Then, basis is a matrix with 10 rows and 10 or fewer columns.
///
///       For dof_key=="GLstrain":
///         The default DoF space has dimension 6 for the six independent
///         GLstrain components. Then, basis is a matrix with 6 rows and 6 or
///         fewer columns.
//
///     By default, basis is the full DoF space (identity matrix with dimensions
///     matching the full space for the particular combination of config_region
///     and dof_key).
class DoFSpace {
 public:
  DoFSpace(std::shared_ptr<Structure const> const &_shared_prim,
           DoFKey const &_dof_key,
           std::optional<Eigen::Matrix3l> const
               &_transformation_matrix_to_super = std::nullopt,
           std::optional<std::set<Index>> const &_sites = std::nullopt,
           std::optional<Eigen::MatrixXd> const &_basis = std::nullopt);

  /// Shared prim structure
  std::shared_ptr<Structure const> const &shared_prim() const;

  /// Type of degree of freedom that is under consideration (e.g., "disp",
  /// "Hstrain", "occ")
  DoFKey const &dof_key() const;

  /// Specifies the supercell for a local DoF space. Has value for local DoF.
  std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super() const;

  /// The sites included in a local DoF space. Has value for local DoF.
  std::optional<std::set<Index>> const &sites() const;

  /// True, if local DoF space with all sites in the supercell included
  bool includes_all_sites() const;

  /// The DoF space basis, as a column vector matrix. May be a subspace (cols <=
  /// rows).
  Eigen::MatrixXd const &basis() const;

  /// The pseudo-inverse of the DoFSpace basis.
  Eigen::MatrixXd const &basis_inv() const;

  /// Solve for normal coordinates resulting in specified DoF values
  Eigen::MatrixXd solve(Eigen::MatrixXd const &prim_dof_values) const;

  /// The DoF space dimension (equal to number of rows in basis).
  Index dim() const;

  /// The DoF subspace dimension (equal to number of columns in basis).
  Index subspace_dim() const;

  /// Names the DoF corresponding to each dimension (row) of the basis
  std::vector<std::string> const &axis_glossary() const;

  /// \brief The supercell site_index corresponding to each dimension (row) of
  /// the basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the column index in
  /// `ConfigDoF.local_dof(dof_key).values()` corresponding to each row in the
  /// DoFSpace basis.
  /// - For "occ" DoF this gives the site index in `ConfigDoF.occupation()`
  /// corresponding to each row in the DoFSpace basis.
  std::optional<std::vector<Index>> const &axis_site_index() const;

  /// \brief The local DoF site DoFSet component index corresponding to each
  /// dimension (row) of the basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the row index in
  /// `ConfigDoF.local_dof(dof_key).values()` corresponding to each row in the
  /// DoFSpace basis.
  /// - For "occ" DoF this gives the index into `Site.occupant_dof()`, which is
  /// the value of `ConfigDoF.occupation()[site_index]`.
  std::optional<std::vector<Index>> const &axis_dof_component() const;

  /// \brief Gives the index of basis row corresponding to a given
  /// axis_site_index and axis_dof_component.
  Index basis_row_index(Index site_index, Index dof_component) const;

  /// \brief Gives the index of basis row corresponding to a given
  /// axis_site_index and axis_dof_component.
  std::optional<std::vector<std::vector<Index>>> const &basis_row_index() const;

 private:
  /// Shared prim structure
  std::shared_ptr<Structure const> m_shared_prim;

  /// Type of degree of freedom that is under consideration (e.g., "disp",
  /// "Hstrain", "occ")
  DoFKey m_dof_key;

  /// Specifies the supercell for a local DoF space. Required for local DoF.
  std::optional<Eigen::Matrix3l> m_transformation_matrix_to_super;

  /// The sites included in a local DoF space. Required for local DoF.
  std::optional<std::set<Index>> m_sites;

  /// The DoF space basis, as a column vector matrix. May be a subspace (cols <=
  /// rows).
  Eigen::MatrixXd m_basis;

  /// The pseudo-inverse of the DoFSpace basis.
  Eigen::MatrixXd m_basis_inv;

  /// QR factorization of basis
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> m_qr;

  /// Names the DoF corresponding to each dimension (row) of the basis
  std::vector<std::string> m_axis_glossary;

  /// The supercell site_index corresponding to each dimension (row) of the
  /// basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the column index in
  /// `ConfigDoF.local_dof(dof_key).values()` corresponding to each row in the
  /// DoFSpace basis.
  /// - For "occ" DoF this gives the site index in `ConfigDoF.occupation()`
  /// corresponding to each row in the DoFSpace basis.
  std::optional<std::vector<Index>> m_axis_site_index;

  /// The local DoF site DoFSet component index corresponding to each dimension
  /// (row) of the basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the row index in
  /// `ConfigDoF.local_dof(dof_key).values()` corresponding to each row in the
  /// DoFSpace basis.
  /// - For "occ" DoF this gives the index into `Site.occupant_dof()`, which is
  /// the value of `ConfigDoF.occupation()[site_index]`.
  std::optional<std::vector<Index>> m_axis_dof_component;

  /// Use the site_index and dof_component to lookup the corresponding basis
  /// row index. Has value for local DoF.
  ///
  /// Usage:
  /// Index basis_row_index = (*m_basis_row_index)[site_index][dof_component]
  std::optional<std::vector<std::vector<Index>>> m_basis_row_index;
};

/// Return true if `dof_space` is valid for `config`
///
/// Checks that:
/// - The prim are equivalent
/// - For local DoF, that the transformation_matrix_to_super are equivalent
bool is_valid_dof_space(Configuration const &config, DoFSpace const &dof_space);

/// Throw if `!is_valid_dof_space(config, dof_space)`
void throw_if_invalid_dof_space(Configuration const &config,
                                DoFSpace const &dof_space);

/// Perform conversions between DoFSpace site indices to config site indices
struct DoFSpaceIndexConverter {
  DoFSpaceIndexConverter(Configuration const &config,
                         DoFSpace const &dof_space);

  std::shared_ptr<Structure const> shared_prim;
  xtal::UnitCellCoordIndexConverter const &config_index_converter;
  xtal::UnitCellCoordIndexConverter dof_space_index_converter;

  /// Perform conversion from Coordinate to DoFSpace site index
  Index dof_space_site_index(xtal::Coordinate const &coord,
                             double tol = TOL) const;

  /// Perform conversion from DoFSpace site index to config site index
  Index dof_space_site_index(Index config_site_index) const;

  /// Perform conversion from Coordinate to config site index
  Index config_site_index(xtal::Coordinate const &coord,
                          double tol = TOL) const;

  /// Perform conversion from DoFSpace site index to config site index
  Index config_site_index(Index dof_space_site_index) const;

  /// \brief Perform conversion from DoFSpace site index to config site index,
  /// with additional translation within config
  Index config_site_index(Index dof_space_site_index,
                          UnitCell const &translation) const;
};

/// Return `config` DoF value as a unrolled coordinate in the prim basis
Eigen::VectorXd get_dof_value(Configuration const &config,
                              DoFSpace const &dof_space);

/// Return `config` DoF value as a coordinate in the DoFSpace basis
Eigen::VectorXd get_normal_coordinate(Configuration const &config,
                                      DoFSpace const &dof_space);

/// \brief Return DoF value in the DoFSpace basis of the subset of the
/// configuration located at a particular coordinate
Eigen::VectorXd get_normal_coordinate_at(
    Configuration const &config, DoFSpace const &dof_space,
    DoFSpaceIndexConverter const &index_converter,
    UnitCell integral_lattice_coordinate);

/// Return mean DoF value of configuration in the DoFSpace basis
Eigen::VectorXd get_mean_normal_coordinate(Configuration const &config,
                                           DoFSpace const &dof_space);

/// Set `config` DoF value from a coordinate in the DoFSpace basis
void set_dof_value(Configuration &config, DoFSpace const &dof_space,
                   Eigen::VectorXd const &normal_coordinate);

/// Return dimension of DoFSpace
Index get_dof_space_dimension(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super =
        std::nullopt,
    std::optional<std::set<Index>> const &sites = std::nullopt);

/// Make DoFSpace axis glossary
std::vector<std::string> make_dof_space_axis_glossary(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super =
        std::nullopt,
    std::optional<std::set<Index>> const &sites = std::nullopt);

/// Make DoFSpace axis glossary, axis site index, and axis dof component
void make_dof_space_axis_info(
    DoFKey dof_key, xtal::BasicStructure const &prim,
    std::optional<Eigen::Matrix3l> const &transformation_matrix_to_super,
    std::optional<std::set<Index>> const &sites,
    std::vector<std::string> &axis_glossary,
    std::optional<std::vector<Index>> &axis_site_index,
    std::optional<std::vector<Index>> &axis_dof_component,
    std::optional<std::vector<std::vector<Index>>> &basis_row_index);

/// Make DoFSpace using `dof_key`, transformation matrix and sites from
/// `input_state`, and `basis`
DoFSpace make_dof_space(
    DoFKey dof_key, ConfigEnumInput const &input_state,
    std::optional<Eigen::MatrixXd> const &basis = std::nullopt);

/// Make a SymGroupRep for a DoFSpace
SymGroupRep const &make_dof_space_symrep(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group,
    MasterSymGroup &symrep_master_group, SymGroupRepID &id);

/// Make VectorSpaceSymReport
VectorSpaceSymReport vector_space_sym_report(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges = false);

/// Make VectorSpaceSymReport
SymRepTools_v2::VectorSpaceSymReport vector_space_sym_report_v2(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges = false);

/// Make DoFSpace with symmetry adapated basis
DoFSpace make_symmetry_adapted_dof_space(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<VectorSpaceSymReport> &symmetry_report);

/// Make DoFSpace with symmetry adapated basis
DoFSpace make_symmetry_adapted_dof_space_v2(
    DoFSpace const &dof_space, SupercellSymInfo const &sym_info,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<SymRepTools_v2::VectorSpaceSymReport> &symmetry_report);

class make_symmetry_adapted_dof_space_error : public std::runtime_error {
 public:
  make_symmetry_adapted_dof_space_error(std::string _what)
      : std::runtime_error(_what) {}
  virtual ~make_symmetry_adapted_dof_space_error() {}
};

/// Removes the default occupation modes from the DoFSpace basis
DoFSpace exclude_default_occ_modes(DoFSpace const &dof_space);

/// Removes the homogeneous mode space from the DoFSpace basis
DoFSpace exclude_homogeneous_mode_space(DoFSpace const &dof_space);

/// Make the homogeneous mode space of a local DoFSpace
Eigen::MatrixXd make_homogeneous_mode_space(DoFSpace const &dof_space);

/// \brief A struct which gives out the mixing information of given
/// column_vector_space and a subspace
///
/// For example, consider a column_vector_space of 3 dimensions.
/// Let it be [q1, q2, q3] where q1, q2 and q3 are columns of the said vector
/// space Let's say q1 - has a zero projection onto the given subspace, q2 - has
/// partial projection onto the subspace, q3 - full projection onto the subspace
/// Then axes_not_in_subspace will be {0} - where 0 is the index of the column
/// of the column_vector_space which has zero projection onto the given subspace
/// Similarly, axes_in_subspace will be {2} & axes_mixed_with_subspace will be
/// {1} If any of the columns of the given column_vector_space has a partial
/// projection onto the subspace, are_axes_mixed_with_subspace will be true;
/// otherwise it will be false
struct VectorSpaceMixingInfo {
  VectorSpaceMixingInfo(Eigen::MatrixXd const &column_vector_space,
                        Eigen::MatrixXd const &subspace, double tol);
  std::vector<Index> axes_not_in_subspace;
  std::vector<Index> axes_in_subspace;
  std::vector<Index> axes_mixed_with_subspace;
  bool are_axes_mixed_with_subspace = true;
};

}  // namespace CASM

#endif
