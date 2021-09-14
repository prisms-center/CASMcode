#ifndef CASM_XTAL_DOFSET
#define CASM_XTAL_DOFSET

#include <map>
#include <unordered_set>
#include <vector>

#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/external/Eigen/Core"

namespace CASM {
namespace xtal {
struct SymOp;

/**
 * DoFSet specifies all identifying information for a vector of continuous
 * independent variables (Degrees of Freedom / DoFs) DoFSets are associated with
 * a specific DoF 'type', which has a predefined 'standard' coordinate system
 * ex: displacement -> 3-vector (x,y,z) -> displacement components (relative to
 * fixed laboratory frame) strain -> 6-vector (e_xx, e_yy, e_zz, sqrt(2)*e_yz,
 * sqrt(2)*e_xz, sqrt(2)*e_xy) -> tensor elements DoFSets have a typename, which
 * specifies the type, and a set of basis vectors, which are denoted relative to
 * the DoF type's standard axes. This allows the DoFSet components to be
 * specified by the user, including the ability to only allow DoF values within
 * a subspace of the standard values. DoFSet records  the DoF typename, the
 * names of the vector components, and the axes of the vector components
 * (relative to a set of standard axes)
 */
class DoFSet {
 public:
  using BasicTraits = AnisoValTraits;

  DoFSet(const BasicTraits &init_traits,
         const std::vector<std::string> &init_component_names,
         const Eigen::MatrixXd &init_basis)
      : m_traits(init_traits),
        m_component_names(init_component_names),
        m_basis(init_basis) {
    assert(m_component_names.size() == this->dim());
    assert(m_basis.cols() == this->dim());
    assert(m_basis.rows() ==
           this->traits().dim());  // TODO: This makes sense I think?
    if (basis().cols() > 0 && basis().rows() > 0) {
      m_inv_basis = basis()
                        .transpose()
                        .colPivHouseholderQr()
                        .solve(Eigen::MatrixXd::Identity(dim(), dim()))
                        .transpose();
    }
  }

  DoFSet(const BasicTraits &init_traits)
      : DoFSet(
            init_traits, init_traits.standard_var_names(),
            Eigen::MatrixXd::Identity(init_traits.dim(), init_traits.dim())) {}

  /// \brief Returns type_name of DoFSet, which should be a standardized DoF
  /// type (e.g., "disp", "magspin", "GLstrain")
  const std::string &type_name() const { return m_traits.name(); }

  /// Returns the names of each of the component axes
  const std::vector<std::string> &component_names() const {
    return this->m_component_names;
  }

  /// \brief Returns traits object for the DoF type of this DoFSet
  BasicTraits const &traits() const { return m_traits; }

  /// Returns the number of dimension of the DoF, corresponding to the number of
  /// axes in the vector space
  Index dim() const { return this->basis().cols(); }

  /// Matrix that relates DoFSet variables to a conventional coordiante system
  Eigen::MatrixXd const &basis() const { return m_basis; }

  /// The pseudoinverse of DoFSet coordinate axes
  ///
  /// The pseudoinverse is useful for transforming from standard basis to this
  /// DoFSet's basis. A standard coordinate column-vector value, 'v_standard',
  /// can be represented in this DoFSet's space by:
  ///
  ///     v = inv_basis * v_standard
  ///
  Eigen::MatrixXd const &inv_basis() const { return m_inv_basis; }

 private:
  /// AnisoValTraits. Describes the type of DoF, and can convert Cartesian
  /// symmetry representations into the appropriate representation
  BasicTraits m_traits;

  /// Names for each axis of the basis, for example "dx", "dy", "dz" for
  /// displacement
  std::vector<std::string> m_component_names;

  /// The basis defines the space of the DoF, which should be a linear
  /// combination of the AnisoValTraits conventional coordinates. For example,
  /// you may want to define displacements that only happen along a particular
  /// direction
  Eigen::MatrixXd m_basis;

  /// The pseudoinverse of DoFSet coordinate axes
  Eigen::MatrixXd m_inv_basis;
};

/// Returns descriptive names of the components in a DoFSet, using
/// AnisoValTraits::variable_descriptors()
std::vector<std::string> component_descriptions(DoFSet const &dofset);

/**
 * Identical to xtal::DoFSet, but also keeps track of a list of molecule names
 * that the DoFSet does not apply to. For example, don't apply displacements to
 * a vacancy.
 */

class SiteDoFSet : public DoFSet {
 public:
  SiteDoFSet(const DoFSet &init_dofset,
             const std::unordered_set<std::string> &init_exclude_occs)
      : DoFSet(init_dofset), m_excluded_occs(init_exclude_occs) {}

  SiteDoFSet(const DoFSet &init_dofset) : DoFSet(init_dofset) {}

  SiteDoFSet(const BasicTraits &init_traits,
             const std::vector<std::string> &init_component_names,
             const Eigen::MatrixXd &init_basis,
             const std::unordered_set<std::string> &init_exclude_occs)
      : SiteDoFSet(DoFSet(init_traits, init_component_names, init_basis),
                   init_exclude_occs) {}

  SiteDoFSet(const BasicTraits &init_traits,
             const std::unordered_set<std::string> &init_exclude_occs = {})
      : SiteDoFSet(DoFSet(init_traits), init_exclude_occs) {}

  /// \brief Returns true if DoFSet is inactive (e.g., takes zero values) when
  /// specified occupant is present
  bool is_excluded_occ(std::string const &_occ_name) const {
    return m_excluded_occs.count(_occ_name);
  }

  /// Return all occupants that the DoFSet should not be applied to
  const std::unordered_set<std::string> &excluded_occupants() const {
    return this->m_excluded_occs;
  }

 private:
  std::unordered_set<std::string> m_excluded_occs;
};

template <typename DoFSetType>
std::map<std::string, DoFSetType> make_dofset_map(
    std::vector<DoFSetType> const &dofset_vec) {
  std::map<std::string, DoFSetType> dofset_map;
  for (auto const &dofset : dofset_vec) {
    dofset_map.emplace(dofset.type_name(), dofset);
  }
  return dofset_map;
}

/**
 * Comparator class for checking equivalence of two DoFSet values.
 * Evaluate by constructing object with one of the values, and then pass
 * the other DoFSet to the oveloaded operator().
 *
 * DoFSets are considered equivalent if:
 * - The traits names match
 * - Basis vectors span the same space
 */

struct DoFSetIsEquivalent_f {
 public:
  DoFSetIsEquivalent_f(const DoFSet &reference_value, double tol)
      : m_reference_dofset(reference_value), m_tol(tol) {}

  /// Returns true if the passed is equivalent to the stored value that *this
  /// was constructed with
  bool operator()(const DoFSet &other_value) const;

 private:
  /// Values passed to operator() will be compared against this
  DoFSet m_reference_dofset;

  /// Tolerance value for making comparisons
  double m_tol;

  /// Returns true if the traits match. Only compares the names
  bool _traits_match(const DoFSet &other_value) const;

  /// Returns true if the axis span the same space. For example the basis would
  /// be considered equivalent if they are the same but have been rotated
  bool _basis_spans_same_space(const DoFSet &other_value) const;
};

/**
 * Comparator class for checking equivalence of two DoFSet values.
 * Behaves exactly like DoFSetIsEquivalent, but also checks for excluded
 * occupants. Evaluate by constructing object with one of the values, and then
 * pass the other DoFSet to the oveloaded operator().
 *
 * DoFSets are considered equivalent if:
 * - The traits names match
 * - They have the same names for each of the axes
 * - Excluded occupants match
 */

struct SiteDoFSetIsEquivalent_f : private DoFSetIsEquivalent_f {
  SiteDoFSetIsEquivalent_f(const SiteDoFSet &reference_value, double tol)
      : DoFSetIsEquivalent_f(reference_value, tol),
        m_reference_excluded_occs(reference_value.excluded_occupants()) {}

  /// Returns true if the passed is equivalent to the stored value that *this
  /// was constructed with
  bool operator()(const SiteDoFSet &other_value) const {
    return DoFSetIsEquivalent_f::operator()(other_value) &&
           this->_excluded_occupants_match(other_value);
  }

 private:
  mutable Eigen::MatrixXd m_U;

  std::unordered_set<std::string> m_reference_excluded_occs;

  bool _excluded_occupants_match(const SiteDoFSet &other_value) const {
    return this->m_reference_excluded_occs == other_value.excluded_occupants();
  }
};

/// Create the symmtery representation for going from one basis to another
Eigen::MatrixXd dofset_transformation_matrix(const Eigen::MatrixXd &from_basis,
                                             const Eigen::MatrixXd &to_basis,
                                             double tol);

}  // namespace xtal
}  // namespace CASM

namespace CASM {
// This is how it's currently organized for xtal::Lattice, but maybe we want
// something else  (see SymTools.hh)
namespace sym {
/// \brief Copy and apply SymOp to a DoFSet
xtal::DoFSet copy_apply(const xtal::SymOp &op, const xtal::DoFSet &_dof);
}  // namespace sym
}  // namespace CASM

#endif
