#ifndef CASM_DoFSet
#define CASM_DoFSet

#include <string>
#include <unordered_set>
#include <vector>

#include "casm/basis_set/DoF.hh"

namespace CASM {

class AnisoValTraits;
class DoFSet;

/// Data structure storing the SymGroupRepID and coordinates axes for continuous
/// vector DoF
///
/// Note:
/// - This is the minimal data required for performing most crystallographic
/// symmetry analysis
struct DoFSetInfo {
  /// Construct DoFSetInfo with SymGroupRepID and DoFSet coordinate axes matrix
  /// (basis)
  ///
  /// \param _id SymGroupRepID allows lookup of the DoFSet symmetry
  /// representation that transforms
  ///        DoF values.
  /// \param _basis DoFSet coordinate axes. Columns are vectors in the standard
  /// vector space that
  ///        specify the alignment of DoFSet components.
  ///
  /// A particular DoFSet column-vector value, 'v', can be represented in the
  /// standard coordinate space by:
  ///
  ///     v_standard = _basis * v
  ///
  /// Notes:
  /// - "_basis" is not required to be full rank. This allows DoF values to be
  /// restricted to
  ///   to a subspace of the standard values.
  /// - Caches the pseudoinverse (which is used for symmetry analysis, among
  /// other things)
  ///
  DoFSetInfo(SymGroupRepID _id, Eigen::Ref<const Eigen::MatrixXd> const &_basis)
      : m_symrep_ID(_id), m_basis(_basis) {
    if (basis().cols() > 0 && basis().rows() > 0)
      m_inv_basis = basis()
                        .transpose()
                        .colPivHouseholderQr()
                        .solve(Eigen::MatrixXd::Identity(dim(), dim()))
                        .transpose();
  }

  /// Returns SymGroupRepID of associated DoFSet
  SymGroupRepID const &symrep_ID() const { return m_symrep_ID; }

  /// Sets SymGroupRepID
  void set_symrep_ID(SymGroupRepID _id) { m_symrep_ID = _id; }

  /// Const reference to the DoFSet coordinate axes
  ///
  /// A particular DoFSet column-vector value, 'v', can be represented in the
  /// standard coordinate space by:
  ///
  ///     v_standard = basis * v
  ///
  /// Note:
  /// - The basis is not required to be full rank (i.e. basis.cols() <
  /// basis().rows() is valid).
  ///   This allows DoF values to be restricted to to a subspace of the standard
  ///   values.
  Eigen::MatrixXd const &basis() const { return m_basis; }

  /// Set the DoFSet coordinate axes
  ///
  /// A particular DoFSet column-vector value, 'v', can be represented in the
  /// standard coordinate space by:
  ///
  ///     v_standard = basis * v
  ///
  /// Note:
  /// - The basis is not required to be full rank (i.e. basis.cols() <
  /// basis().rows() is valid).
  ///   This allows DoF values to be restricted to to a subspace of the standard
  ///   values.
  /// - Caches the pseudoinverse (which is used for symmetry analysis, among
  /// other things)
  Eigen::MatrixXd set_basis(Eigen::Ref<const Eigen::MatrixXd> const &_basis) {
    m_basis = _basis;
    if (basis().cols() > 0 && basis().rows() > 0)
      m_inv_basis = basis()
                        .transpose()
                        .colPivHouseholderQr()
                        .solve(Eigen::MatrixXd::Identity(dim(), dim()))
                        .transpose();
    return m_basis;
  }

  /// The pseudoinverse of DoFSet coordinate axes
  ///
  /// The pseudoinverse is useful for transforming from standard basis to this
  /// DoFSet's basis. A standard coordinate column-vector value, 'v_standard',
  /// can be represented in this DoFSet's space by:
  ///
  ///     v = inv_basis * v_standard
  ///
  Eigen::MatrixXd const &inv_basis() const { return m_inv_basis; }

  /// Dimension of the DoFSet, equivalent to basis().cols()
  Index dim() const { return basis().cols(); }

 private:
  SymGroupRepID m_symrep_ID;
  Eigen::MatrixXd m_basis;
  Eigen::MatrixXd m_inv_basis;
};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// DoFSet: A set of degrees of freedom (DoF)
///
/// DoFSet specifies all identifying information for a vector of continuous
/// independent variables (may be site or global variables). For example, a
/// DoFSet may be used to represent allowed displacements on one site, or
/// allowed lattice strains.
///
/// A DoFSet has:
///     - AnisoValTraits, which provides the DoF type name, a standard
///     coordinate system (the
///       "standard basis"), and specifies how values transform under
///       application of symmetry.
///     - a "DoF basis", a set of named basis vectors which are denoted relative
///     to the standard
///       basis, allowing the user to specify the DoFSet components, name them,
///       and restrict DoF values to a particular subspace
///         - The basis is stored as a matrix and components of the basis are
///         also represented as
///           a `std::vector<ContinuousDoF>`. The `ContinuousDoF` are used for
///           constructing basis functions.
///     - the SymGroupRepID of the DoFSet, which is a key for retrieving the
///     SymGroupRep that
///       encodes how the DoFSet transforms with symmetry from the prim
///       Structure factor group
///     - a list of site occupants for which the DoF does not apply
///     ("excluded_occupants",
///       std::unordered_set<std::string>). As an example, this could be used if
///       some allowed site occupant molecules have magnetic spin, but other
///       allowed occupants do not.
///
/// Examples of standard basis specified by AnisoValTraits:
/// - "disp" -> (dx, dy, dz) -> displacement components relative to fixed
/// laboratory frame
/// - "strain" -> (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy)
/// -> tensor elements
///
class DoFSet {
 public:
  using BasicTraits = AnisoValTraits;

  using TypeFunc = std::function<notstd::cloneable_ptr<AnisoValTraits>()>;

  using Container = std::vector<ContinuousDoF>;

  using const_iterator = std::vector<ContinuousDoF>::const_iterator;

  /// Construct a DoFSet having standard dimension and coordinate basis for
  /// given type.
  ///
  /// Example usage:
  /// \code
  /// DoFSet my_disp_dof = DoFSet::make_default("disp");
  /// \endcode
  static DoFSet make_default(AnisoValTraits const &_type);

  /// General case DoFSet constructor
  ///
  /// \param _type AnisoValTraits specifying the type of DoF.
  /// \param components Name of each component of the DoFSet. (i.e. "dx", "dy",
  /// "dz" for "disp" DoF) \param _info DoFSetInfo containing the SymGroupRepID
  /// and DoFSet coordinate axes matrix (basis) \param excluded_occupants A list
  /// of site occupants for which the DoF does not apply
  ///
  /// Note: Use the static function `DoFSet::make_default` to construct a DoFSet
  /// having standard dimension and coordinate basis.
  ///
  DoFSet(AnisoValTraits const &_type, std::vector<std::string> components,
         DoFSetInfo _info,
         std::unordered_set<std::string> excluded_occupants = {});

  /// Returns the number of components in this DoFSet
  Index size() const { return m_components.size(); }

  /// Returns the DoF type name
  ///
  /// Note: Will be a standardized DoF type name (e.g., "disp", "magspin",
  /// "GLstrain") from the DoFSet's AnisoValTraits object.
  std::string const &type_name() const { return m_traits.name(); }

  /// Returns the AnisoValTraits object for the DoF type of this DoFSet
  AnisoValTraits const &traits() const { return m_traits; }

  /// Set identifying index (ID) of components in this DoFSet
  ///
  /// Notes:
  /// - ID is context dependent but is used to different distinct DoFSets that
  /// are otherwise
  ///   equivalent (i.e., have same type_name(), components, etc).
  /// - For example, "disp" DoFSet on two sites in the primitive cell may have
  /// the same basis and
  ///   be indistringuishable except for being associated with different sites.
  ///   Their components are given different IDs in order to distinguish the
  ///   DoFSets from each other.
  void set_ID(Index _ID) {
    for (auto &c : m_components) c.set_ID(_ID);
  }

  /// Locks IDs of components in this DoFSet so they can no longer be updated
  void lock_IDs();

  /// Sets the IDs of each component of the DoFSet to be their index
  /// within the set. That is, components of the set get labeled [0, 1, ...,
  /// n-1]
  void set_sequential_IDs();

  /// Returns reference to DoFSetInfo (which contains SymGroupRepID and basis
  /// vectors of the DoFSet coordinate system)
  DoFSetInfo const &info() const { return m_info; }

  /// Return i'th component of DoFSet
  ///
  /// DoFSet components are represented as ContinuousDoF, which, in addition to
  /// a type_name also has a variable name.
  ContinuousDoF const &operator[](Index i) const { return m_components[i]; }

  /// Iterator pointing to the first ContinuousDoF component
  const_iterator begin() const { return m_components.cbegin(); }

  /// Iterator pointing one past last ContinuousDoF component
  const_iterator end() const { return m_components.cend(); }

  /// Const iterator to first ContinuousDoF component
  const_iterator cbegin() const { return m_components.cbegin(); }

  /// Const iterator pointing one past last ContinuousDoF component
  const_iterator cend() const { return m_components.cend(); }

  /// Returns true if DoFSet is inactive (e.g., takes zero values) when
  /// specified occupant is present
  bool is_excluded_occ(std::string const &_occ_name) const {
    return m_excluded_occs.count(_occ_name);
  }

  /// Dimension of the DoFSet, equivalent to basis().cols()
  Index dim() const { return m_info.dim(); }

  /// Const reference to the DoFSet coordinate axes
  ///
  /// A particular DoFSet column-vector value, 'v', can be represented in the
  /// standard coordinate space by:
  ///
  ///     v_standard = basis * v
  ///
  /// Note:
  /// - The basis is not required to be full rank (i.e. basis.cols() <
  /// basis().rows() is valid).
  ///   This allows DoF values to be restricted to to a subspace of the standard
  ///   values.
  Eigen::MatrixXd const &basis() const { return m_info.basis(); }

  /// SymGroupRepID of this DoFSet
  ///
  /// Note:
  /// - SymGroupRepID is the key for retrieving the SymGroupRep that encodes how
  /// the DoFSet
  ///   transforms with symmetry from the prim Structure factor group
  SymGroupRepID const &symrep_ID() const { return m_info.symrep_ID(); }

  /// Returns true if "rhs" has identical components and basis to this DoFSet
  bool identical(DoFSet const &rhs) const;

  /// Equivalent to `m_basis=trans_mat*m_basis`. Invalidates SymGroupRepID.
  void transform_basis(Eigen::Ref<const Eigen::MatrixXd> const &trans_mat);

  /// Update component IDs
  ///
  /// Notes:
  /// - This is a convenience function used when DoFSet IDs are being set en
  /// masse
  /// - Each component with ID in `before_IDs[i]` will get changed to the
  /// corresponding ID in
  ///   `after_IDs[i]`.
  /// - If any component ID is updated, returns true; otherwise returns false.
  bool update_IDs(const std::vector<Index> &before_IDs,
                  const std::vector<Index> &after_IDs);

 private:
  /// AnisoValTraits. Describes the type of DoF, and can convert Cartesian
  /// symmetry representations into the appropriate representation
  BasicTraits m_traits;

  /// ContinuousDoF components, one for each column in `m_info->basis()`
  std::vector<ContinuousDoF> m_components;

  /// DoFSetInfo containing the SymGroupRepID and DoFSet coordinate axes matrix
  /// (basis)
  DoFSetInfo m_info;

  /// A list of site occupants for which the DoF does not apply
  std::unordered_set<std::string> m_excluded_occs;
};

//********************************************************************

inline bool operator==(DoFSet const &A, DoFSet const &B) {
  return A.identical(B);
}

//********************************************************************

inline bool operator!=(DoFSet const &A, DoFSet const &B) {
  return !A.identical(B);
}

}  // namespace CASM

#endif
