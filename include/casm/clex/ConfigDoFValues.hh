#ifndef CASM_ConfigDoFValues
#define CASM_ConfigDoFValues

#include "casm/basis_set/DoFSet.hh"
#include "casm/crystallography/AnisoValTraits.hh"
namespace CASM {

class ConfigDoFValues {
 public:
  ConfigDoFValues(DoFKey const &_type_name, Index _n_sublat, Index _n_vol);

  std::string const &type_name() const;

  Index n_vol() const;

  Index n_sublat() const;

 private:
  DoFKey m_type_name;
  Index m_n_sublat;
  Index m_n_vol;
};

/// Access a vector with local discrete (occupation) DoF values
///
/// Notes:
/// - Values are stored in a vector of size=n_vol()*n_sublat()
/// - The value at site_index l indicates the type of Molecule occupying the
///   site, via `prim.basis()[b].occupant_dof()[this->occ(l)]`, where the
///   sublattice, `b = l / this->n_vol()`.
/// - This class does not hold the values directly, but includes a reference to
///   the plain data structure `clexulator::ConfigDoFValues` which does.
///
/// See ConfigDoF documentation for more details on how DoF values are stored.
///
/// The vector size is fixed at construction and attempts to change it, other
/// than via a complete copy, will cause an exception to be thrown.
///
class LocalDiscreteConfigDoFValues : public ConfigDoFValues {
 public:
  typedef Eigen::VectorXi ValueType;
  typedef Eigen::VectorXi &Reference;
  typedef Eigen::VectorXi const &ConstReference;

  typedef typename ValueType::Scalar SiteValueType;
  typedef int &SiteReference;
  typedef const int &ConstSiteReference;

  typedef ValueType SublatValueType;
  typedef typename ValueType::SegmentReturnType SublatReference;
  typedef typename ValueType::ConstSegmentReturnType ConstSublatReference;

  LocalDiscreteConfigDoFValues(DoFKey const &_type_name, Index _n_sublat,
                               Index _n_vol,
                               std::vector<SymGroupRepID> const &_symrep_IDs,
                               ValueType &_values);

  LocalDiscreteConfigDoFValues(LocalDiscreteConfigDoFValues const &other) =
      delete;

  LocalDiscreteConfigDoFValues &operator=(
      LocalDiscreteConfigDoFValues const &other) = delete;

  /// Reference occupation value on site i
  int &occ(Index i) {  // keep inline
    return m_vals(i);
  }

  /// Const reference occupation value on site i
  int const &occ(Index i) const {  // keep inline
    return m_vals(i);
  }

  /// Set occupation values (values are indices into Site::occupant_dof())
  void set_values(Eigen::Ref<ValueType const> const &_values);

  /// Set occupation values to zero
  void setZero();

  /// Const access occupation values (values are indices into
  /// Site::occupant_dof())
  Eigen::VectorXi const &values() const {  // keep inline
    return m_vals;
  }

  /// Access vector block of values for all sites on one sublattice
  SublatReference sublat(Index b);

  /// Const access vector block of values for all sites on one sublattice
  ConstSublatReference sublat(Index b) const;

  /// Provides the symmetry representations for transforming `values` (i.e. due
  /// to molecule orientation, not for permuting sites)
  std::vector<SymGroupRepID> const &symrep_IDs() const;

 private:
  void _throw_if_invalid_size(Eigen::Ref<ValueType const> const &_values) const;

  ValueType &m_vals;
  std::vector<SymGroupRepID> m_symrep_IDs;
};

/// Access a matrix with local continuous DoF values
///
/// Notes:
/// - Values are stored in a matrix of size rows=dim(), cols=n_vol()*n_sublat(),
///   - The value of dim() is the maximum, over all sublattices in the prim, of
///     the site basis dimension for this type of DoF. In case the site basis
///     dimension varies from site to site (for example, displacements are
///     restricted to a 2d plane on some sites, but not others), the columns
///     corresponding to sublattices with a basis dimension less than dim()
///     have a tail of zeros which should not be modified.
/// - Each column represents a site DoF value in the prim DoF basis
/// - The prim DoF basis for each sublattice can be accessed by
///   `this->info()[b]`, where `b` is the sublattice index, `b = column_index /
///   this->n_vol()`
/// - This class does not hold the values directly, but includes a reference to
///   the plain data structure `clexulator::ConfigDoFValues` which does.
///
/// See ConfigDoF documentation for more details on how DoF values are stored.
///
/// The matrix size is fixed at construction and attempts to change it, other
/// than via a complete copy, will cause an exception to be thrown.
///
class LocalContinuousConfigDoFValues : public ConfigDoFValues {
 public:
  typedef Eigen::MatrixXd ValueType;
  typedef Eigen::MatrixXd &Reference;
  typedef Eigen::MatrixXd const &ConstReference;

  typedef Eigen::VectorXd SiteValueType;
  typedef typename ValueType::ColXpr SiteReference;
  typedef const typename ValueType::ConstColXpr ConstSiteReference;

  typedef ValueType SublatValueType;
  typedef typename Eigen::Block<ValueType> SublatReference;
  typedef const typename Eigen::Block<const ValueType> ConstSublatReference;

  /// local continuous DoF values matrix has #rows == max( DoFSetInfo::dim() )
  static Index matrix_dim(std::vector<DoFSetInfo> const &_info);

  LocalContinuousConfigDoFValues(DoFKey const &_type_name, Index _n_sublat,
                                 Index _n_vol,
                                 std::vector<DoFSetInfo> const &_info,
                                 ValueType &_values);

  LocalContinuousConfigDoFValues(LocalContinuousConfigDoFValues const &other) =
      delete;

  LocalContinuousConfigDoFValues &operator=(
      LocalContinuousConfigDoFValues const &other) = delete;

  /// maximum DoF vector representation size (max of DoFSetInfo::dim())
  Index dim() const;

  /// Access site DoF values (prim DoF basis, matrix representing all sites)
  void set_values(Eigen::Ref<const ValueType> const &_values);

  /// Set DoF values to zero
  void setZero();

  /// Const access DoF values (prim DoF basis, matrix representing all sites)
  Eigen::MatrixXd const &values() const;

  /// Set local DoF values from standard DoF values
  void from_standard_values(
      Eigen::Ref<const Eigen::MatrixXd> const &_standard_values);

  /// Get local DoF values as standard DoF values
  Eigen::MatrixXd standard_values() const;

  /// Access site DoF value vector
  SiteReference site_value(Index l);

  /// Const access site DoF value vector
  ConstSiteReference site_value(Index l) const;

  /// Access matrix block of values for all sites on one sublattice
  SublatReference sublat(Index b);

  /// Const access matrix block of values for all sites on one sublattice
  ConstSublatReference sublat(Index b) const;

  /// DoFSetInfo provides the basis and symmetry representations for `values`
  std::vector<DoFSetInfo> const &info() const;

 private:
  void _throw_if_invalid_size(Eigen::Ref<ValueType const> const &_values) const;

  Index m_dim;
  std::vector<DoFSetInfo> m_info;
  ValueType &m_vals;
};

// TODO: It might be confusing that ValueType, Reference, and ConstReference are
// in terms of Eigen::VectorXd, but dim(), from_standard_values(), and
// standard_values() are written using Eigen::MatrixXd. Since Eigen::VectorXd is
// a 1 column Eigen::MatrixXd it could be all written in terms of
// Eigen::VectorXd or Eigen::MatrixXd. Not sure if it is more useful to express
// in terms of Eigen::MatrixXd and match the LocalContinuousConfigDoFValues
// interface or use Eigen::VectorXd consistently to express that the value is 1
// dimensional.

/// Access a vector with global continuous DoF values
///
/// Dof values are stored as values in the prim DoF basis. The conversion to
/// the standard DoF basis values is:
///
///     this->standard_values() = info().basis() * this->values()
///
/// See ConfigDoF documentation for more details on how DoF values are stored.
///
/// The matrix size is fixed at construction and attempts to change it, other
/// than via a complete copy, will cause an exception to be thrown.
///
/// Notes:
/// - This class does not hold the values directly, but includes a reference to
///   the plain data structure `clexulator::ConfigDoFValues` which does.
///
class GlobalContinuousConfigDoFValues : public ConfigDoFValues {
 public:
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::VectorXd &Reference;
  typedef Eigen::VectorXd const &ConstReference;

  typedef typename ValueType::Scalar SiteValueType;
  typedef int &SiteReference;
  typedef const int &ConstSiteReference;

  GlobalContinuousConfigDoFValues(DoFKey const &_name, Index _n_sublat,
                                  Index _n_vol, DoFSetInfo const &_info,
                                  ValueType &_values);

  GlobalContinuousConfigDoFValues(
      GlobalContinuousConfigDoFValues const &other) = delete;

  GlobalContinuousConfigDoFValues &operator=(
      GlobalContinuousConfigDoFValues const &other) = delete;

  /// Global DoF vector representation dimension
  Index dim() const;

  /// Set global DoF values
  void set_values(Eigen::Ref<const Eigen::MatrixXd> const &_values);

  /// Set DoF values to zero
  void setZero();

  /// Const access global DoF values
  Eigen::VectorXd const &values() const {  // keep inline
    return m_vals;
  }

  /// Set global DoF values from standard DoF values
  void from_standard_values(
      Eigen::Ref<const Eigen::MatrixXd> const &_standard_values);

  /// Get global DoF values as standard DoF values
  Eigen::MatrixXd standard_values() const;

  /// DoFSetInfo provides the basis and symmetry representations for `values`
  DoFSetInfo const &info() const;

 private:
  void _throw_if_invalid_size(Eigen::Ref<ValueType const> const &_values) const;

  DoFSetInfo m_info;
  ValueType &m_vals;
};

/// Access site DoF value vector
///
/// Note:
/// - This is the vector of DoF values associated with a single site, in the
///   prim DoF basis
/// - If the prim DoF basis dimension for this site is less than dim(), this
///   includes a tail of zeros (for rows >= this->info()[b].dim(), where b is
///   the sublattice index for site_index l). The tail of zeros should not be
///   modified.
inline  // keep inline
    LocalContinuousConfigDoFValues::SiteReference
    LocalContinuousConfigDoFValues::site_value(Index l) {
  return m_vals.col(l);
}

/// Const access site DoF value vector
///
/// Note:
/// - This is the vector of DoF values associated with a single site, in the
///   prim DoF basis
/// - If the prim DoF basis dimension for this site is less than dim(), this
///   includes a tail of zeros (for rows >= this->info()[b].dim(), where b is
///   the sublattice index for site_index l). The tail of zeros should not be
///   modified.
inline  // keep inline
    LocalContinuousConfigDoFValues::ConstSiteReference
    LocalContinuousConfigDoFValues::site_value(Index l) const {
  return const_cast<ValueType const &>(m_vals).col(l);
}

/// Set global DoF values
///
/// \throws std::runtime_error ("Invalid size in
/// GlobalContinuousConfigDoFValues...") if size is not valid
inline  // keep inline
    void
    GlobalContinuousConfigDoFValues::set_values(
        Eigen::Ref<const Eigen::MatrixXd> const &_values) {
  _throw_if_invalid_size(_values);
  m_vals = _values;
}

inline  // keep inline
    void
    GlobalContinuousConfigDoFValues::_throw_if_invalid_size(
        Eigen::Ref<ValueType const> const &_values) const {
  if (_values.size() != m_vals.size()) {
    std::stringstream msg;
    msg << "Invalid size in LocalContinuousConfigDoFValues: "
        << "Expected size=" << m_info.basis().cols()
        << ", received size=" << _values.size();
    throw std::runtime_error(msg.str());
  }
}

}  // namespace CASM

#endif
