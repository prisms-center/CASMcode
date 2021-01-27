#ifndef CASM_DoFSpace
#define CASM_DoFSpace

#include "casm/crystallography/DoFDecl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

struct VectorSpaceSymReport;

/// DoFSpace
///
/// The DoFSpace class is used to specify a particular degree of freedom space.
/// DoFSpace have:
///   - shared_prim (std::shared_ptr<Structure const>): The prim structure
///   - dof_key (DoFKey): A string indicating which DoF type (e.g., "disp",
///   "Hstrain", "occ")
///   - transformation_matrix_to_super (std::optional<Eigen::Matrix3l>):
///   Specifies the supercell
///     for a local DoF space. Has value for local DoF.
///   - sites (std::optional<std::set<Index>>): The sites included in a local
///   DoF space. Has value
///     for local DoF.
///   - basis (Eigen::MatrixXd): Allows specifying a subspace of the space
///   determined from
///     dof_key, and for local DoF, transformation_matrix_to_super and sites.
///     Examples:
///
///       For dof_key=="disp", and sites.value().size()==4 (any
///       transformation_matrix_to_super):
///         The default DoF space has dimension 12 corresponding to (x,y,z) for
///         each of the 4 sites. Then, basis is a matrix with 12 rows and 12 or
///         fewer columns.
///
///       For dof_key=="occ", and sites().value().size()==4 (any
///       transformation_matrix_to_super), with the particular sites selected
///       having allowed occupations ["A", "B"],
///       ["A", "B", "C"], ["A, B"], and ["A", "B", "C"] :
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

  /// The DoF space basis, as a column vector matrix. May be a subspace (cols <=
  /// rows).
  Eigen::MatrixXd const &basis() const;

  /// The DoF space dimension (equal to number of rows in basis).
  Index dim() const;

  /// The DoF subspace dimension (equal to number of columns in basis).
  Index subspace_dim() const;

  /// The inverse of the DoF space basis.
  Eigen::MatrixXd const &basis_inverse() const;

  /// Names the DoF corresponding to each dimension (row) of the basis
  std::vector<std::string> const &axis_glossary() const;

  /// The supercell site_index corresponding to each dimension (row) of the
  /// basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the column index in
  /// `ConfigDoF.local_dof(dof_key).values()`
  ///   corresponding to each row in the DoFSpace basis.
  /// - For "occ" DoF this gives the site index in `ConfigDoF.occupation()`
  /// corresponding to
  ///   each row in the DoFSpace basis.
  std::optional<std::vector<Index>> const &axis_site_index() const;

  /// The local DoF site DoFSet component index corresponding to each dimension
  /// (row) of the basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the row index in
  /// `ConfigDoF.local_dof(dof_key).values()`
  ///   corresponding to each row in the DoFSpace basis.
  /// - For "occ" DoF this gives the index into `Site.occupant_dof()`, which is
  /// the value of
  ///   of `ConfigDoF.occupation()[site_index]`.
  std::optional<std::vector<Index>> const &axis_dof_component() const;

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

  /// The pseudoinverse of the DoF space basis.
  Eigen::MatrixXd m_basis_inverse;

  /// Names the DoF corresponding to each dimension (row) of the basis
  std::vector<std::string> m_axis_glossary;

  /// The supercell site_index corresponding to each dimension (row) of the
  /// basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the column index in
  /// `ConfigDoF.local_dof(dof_key).values()`
  ///   corresponding to each row in the DoFSpace basis.
  /// - For "occ" DoF this gives the site index in `ConfigDoF.occupation()`
  /// corresponding to
  ///   each row in the DoFSpace basis.
  std::optional<std::vector<Index>> m_axis_site_index;

  /// The local DoF site DoFSet component index corresponding to each dimension
  /// (row) of the basis. Has value for local DoF.
  ///
  /// Note:
  /// - For continuous DoF this gives the row index in
  /// `ConfigDoF.local_dof(dof_key).values()`
  ///   corresponding to each row in the DoFSpace basis.
  /// - For "occ" DoF this gives the index into `Site.occupant_dof()`, which is
  /// the value of
  ///   of `ConfigDoF.occupation()[site_index]`.
  std::optional<std::vector<Index>> m_axis_dof_component;
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

/// Return `config` DoF value as a coordinate in the DoFSpace basis
Eigen::VectorXd get_normal_coordinate(Configuration const &config,
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
    std::optional<std::vector<Index>> &axis_dof_component);

/// Make DoFSpace using `dof_key`, transformation matrix and sites from
/// `input_state`, and `basis`
DoFSpace make_dof_space(
    DoFKey dof_key, ConfigEnumInput const &input_state,
    std::optional<Eigen::MatrixXd> const &basis = std::nullopt);

/// Make VectorSpaceSymReport
template <typename PermuteIteratorIt>
VectorSpaceSymReport vector_space_sym_report(DoFSpace const &dof_space,
                                             ConfigEnumInput const &input_state,
                                             PermuteIteratorIt group_begin,
                                             PermuteIteratorIt group_end,
                                             bool calc_wedges = false);

/// Make DoFSpace with symmetry adapated basis
DoFSpace make_symmetry_adapted_dof_space(
    DoFSpace const &dof_space, ConfigEnumInput const &input_state,
    std::vector<PermuteIterator> const &group, bool calc_wedges,
    std::optional<VectorSpaceSymReport> &symmetry_report);

class make_symmetry_adapted_dof_space_error : public std::runtime_error {
 public:
  make_symmetry_adapted_dof_space_error(std::string _what)
      : std::runtime_error(_what) {}
  virtual ~make_symmetry_adapted_dof_space_error() {}
};

}  // namespace CASM

#endif
