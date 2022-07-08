#include "casm/clex/ConfigDoFValues.hh"

#include "casm/clexulator/ConfigDoFValuesTools.hh"

namespace CASM {

ConfigDoFValuesBase::ConfigDoFValuesBase(DoFKey const &_type_name,
                                         Index _n_sublat, Index _n_vol)
    : m_type_name(_type_name), m_n_sublat(_n_sublat), m_n_vol(_n_vol) {}

std::string const &ConfigDoFValuesBase::type_name() const {
  return m_type_name;
}

Index ConfigDoFValuesBase::n_vol() const { return m_n_vol; }

Index ConfigDoFValuesBase::n_sublat() const { return m_n_sublat; }

LocalDiscreteConfigDoFValues::LocalDiscreteConfigDoFValues(
    DoFKey const &_type_name, Index _n_sublat, Index _n_vol,
    std::vector<SymGroupRepID> const &_symrep_IDs, ValueType &_values)
    : ConfigDoFValuesBase(_type_name, _n_sublat, _n_vol),
      m_vals(_values),
      m_symrep_IDs(_symrep_IDs) {}

/// Set occupation values (values are indices into Site::occupant_dof())
///
/// \throws std::runtime_error ("Invalid size in
/// LocalDiscreteConfigDoFValues...") if size does not match _n_sublat *
/// _n_vol
void LocalDiscreteConfigDoFValues::set_values(
    Eigen::Ref<ValueType const> const &_values) {
  _throw_if_invalid_size(_values);
  m_vals = _values;
}

/// Set occupation values to zero
void LocalDiscreteConfigDoFValues::setZero() { m_vals.setZero(); }

/// Access vector block of values for all sites on one sublattice
LocalDiscreteConfigDoFValues::SublatReference
LocalDiscreteConfigDoFValues::sublat(Index b) {
  return m_vals.segment(b * n_vol(), n_vol());
}

/// Const access vector block of values for all sites on one sublattice
LocalDiscreteConfigDoFValues::ConstSublatReference
LocalDiscreteConfigDoFValues::sublat(Index b) const {
  return const_cast<ValueType const &>(m_vals).segment(b * n_vol(), n_vol());
}

/// Provides the symmetry representations for transforming `values` (i.e. due
/// to molecule orientation, not for permuting sites)
std::vector<SymGroupRepID> const &LocalDiscreteConfigDoFValues::symrep_IDs()
    const {
  return m_symrep_IDs;
}

void LocalDiscreteConfigDoFValues::_throw_if_invalid_size(
    Eigen::Ref<ValueType const> const &_values) const {
  if (_values.size() != n_vol() * n_sublat()) {
    std::stringstream msg;
    msg << "Invalid size in LocalDiscreteConfigDoFValues: "
        << "Expected size=" << n_vol() * n_sublat()
        << ", received size=" << _values.size();
    throw std::runtime_error(msg.str());
  }
}

/// local continuous DoF values matrix has #rows == max( DoFSetInfo::dim() )
Index LocalContinuousConfigDoFValues::matrix_dim(
    std::vector<DoFSetInfo> const &_info) {
  return clexulator::max_dim(_info);
}

LocalContinuousConfigDoFValues::LocalContinuousConfigDoFValues(
    DoFKey const &_type_name, Index _n_sublat, Index _n_vol,
    std::vector<DoFSetInfo> const &_info, ValueType &_values)
    : ConfigDoFValuesBase(_type_name, _n_sublat, _n_vol),
      m_dim(matrix_dim(_info)),
      m_info(_info),
      m_vals(_values) {
  _throw_if_invalid_size(m_vals);
}

/// maximum DoF vector representation size (max of DoFSetInfo::dim())
Index LocalContinuousConfigDoFValues::dim() const { return m_dim; }

/// Access site DoF values (prim DoF basis, matrix representing all sites)
///
/// Notes:
/// - Matrix of size:
///   - rows=dim() (the max, over all sublattices, prim DoF basis dimension)
///   - cols=n_vol()*n_sublat()
/// - Each column represents a site DoF value in the prim DoF basis
/// - The prim DoF basis can be accessed by `this->info()[b].basis()`, where
///   `b` is the sublattice index, `b = column_index / this->n_vol()`
/// - If the prim DoF basis dimension (this->info()[b].basis().cols()) is less
///   than the standard DoF basis dimension, the matrix includes blocks of
///   zeros for the corresponding sublattice.
///
/// \throws std::runtime_error ("Invalid size in
/// LocalContinuousConfigDoFValues...") if size is not valid
void LocalContinuousConfigDoFValues::set_values(
    Eigen::Ref<const ValueType> const &_values) {
  _throw_if_invalid_size(_values);
  m_vals = _values;
}

/// Set DoF values to zero
void LocalContinuousConfigDoFValues::setZero() { m_vals.setZero(); }

/// Const access DoF values (prim DoF basis, matrix representing all sites)
///
/// Notes:
/// - Matrix of size rows=dim(), cols=n_vol()*n_sublat(),
///   - In case the site basis dimension varies from site to site (for
///     example, displacements are restricted to a 2d plane on some sites,
///     but not others), the value of dim() is the maximum site basis
///     dimension for this DoF. The columns corresponding to sublattices with
///     a smaller basis dimension than dim() have a tail of zeros which should
///     not be modified.
/// - Each column represents a site DoF value in the prim DoF basis
/// - The prim DoF basis can be accessed by `this->info()[b]`, where `b` is
///   the sublattice index, `b = column_index / this->n_vol()`
///
Eigen::MatrixXd const &LocalContinuousConfigDoFValues::values() const {
  return m_vals;
}

/// Set local DoF values from standard DoF values
///
/// Notes:
/// - Standard DoF values are those expressed according to the standard DoF
///   basis, i.e. coordinates whose values corrspond to
///   AnisoValTraits::standard_var_names().
void LocalContinuousConfigDoFValues::from_standard_values(
    Eigen::Ref<const Eigen::MatrixXd> const &_standard_values) {
  if (_standard_values.rows() != info().front().basis().rows() ||
      _standard_values.cols() != n_vol() * n_sublat()) {
    std::stringstream msg;
    msg << "Invalid standard values input size in "
           "LocalContinuousConfigDoFValues: "
        << "Expected rows=" << info().front().basis().rows()
        << ", received rows=" << _standard_values.rows()
        << ", expected cols=" << n_vol() * n_sublat()
        << ", received cols=" << _standard_values.cols();
    throw std::runtime_error(msg.str());
  }
  for (Index b = 0; b < n_sublat(); ++b){
    if (info()[b].symrep_ID().empty()){
        continue;
    }
    sublat(b).topRows(info()[b].dim()) =
        info()[b].inv_basis() *
        _standard_values.block(0, b * n_vol(), m_vals.rows(), n_vol());
  }
}

/// Get local DoF values as standard DoF values
///
/// Notes:
/// - Standard DoF values are those expressed according to the standard DoF
///   basis, i.e. coordinates whose values corrspond to
///   AnisoValTraits::standard_var_names().
Eigen::MatrixXd LocalContinuousConfigDoFValues::standard_values() const {
  Index rows = m_info.front().basis().rows();
  Eigen::MatrixXd result(rows, m_vals.cols());
  for (Index b = 0; b < n_sublat(); ++b) {
    result.block(0, b * n_vol(), rows, n_vol()) =
        info()[b].basis() * sublat(b).topRows(info()[b].dim());
  }
  return result;
}

/// Access matrix block of values for all sites on one sublattice
LocalContinuousConfigDoFValues::SublatReference
LocalContinuousConfigDoFValues::sublat(Index b) {
  return m_vals.block(0, b * n_vol(), m_vals.rows(), n_vol());
}

/// Const access matrix block of values for all sites on one sublattice
LocalContinuousConfigDoFValues::ConstSublatReference
LocalContinuousConfigDoFValues::sublat(Index b) const {
  return const_cast<ValueType const &>(m_vals).block(0, b * n_vol(),
                                                     m_vals.rows(), n_vol());
}

/// DoFSetInfo provides the basis and symmetry representations for `values`
std::vector<DoFSetInfo> const &LocalContinuousConfigDoFValues::info() const {
  return m_info;
}

void LocalContinuousConfigDoFValues::_throw_if_invalid_size(
    Eigen::Ref<ValueType const> const &_values) const {
  if (_values.rows() != dim() || _values.cols() != n_vol() * n_sublat()) {
    std::stringstream msg;
    msg << "Invalid size in LocalContinuousConfigDoFValues: "
        << "Expected rows=" << dim() << ", received rows=" << _values.rows()
        << ", expected cols=" << n_vol() * n_sublat()
        << ", received cols=" << _values.cols();
    throw std::runtime_error(msg.str());
  }
}

GlobalContinuousConfigDoFValues::GlobalContinuousConfigDoFValues(
    DoFKey const &_name, Index _n_sublat, Index _n_vol, DoFSetInfo const &_info,
    ValueType &_values)
    : ConfigDoFValuesBase(_name, _n_sublat, _n_vol),
      m_info(_info),
      m_vals(_values) {}

/// Global DoF vector representation dimension
Index GlobalContinuousConfigDoFValues::dim() const { return m_vals.rows(); }

/// Set DoF values to zero
void GlobalContinuousConfigDoFValues::setZero() { m_vals.setZero(); }

/// Set global DoF values from standard DoF values
///
/// Notes:
/// - Standard DoF values are those expressed according to the standard DoF
///   basis, i.e. coordinates whose values corrspond to
///   AnisoValTraits::standard_var_names().
void GlobalContinuousConfigDoFValues::from_standard_values(
    Eigen::Ref<const Eigen::MatrixXd> const &_standard_values) {
  if (_standard_values.size() != m_info.basis().rows()) {
    std::stringstream msg;
    msg << "Invalid standard values input size in "
           "GlobalContinuousConfigDoFValues: "
        << "Expected size=" << m_info.basis().rows()
        << ", received size=" << _standard_values.size();
    throw std::runtime_error(msg.str());
  }
  m_vals = info().inv_basis() * _standard_values;
}

/// Get global DoF values as standard DoF values
///
/// Notes:
/// - Standard DoF values are those expressed according to the standard DoF
/// basis, i.e.
///   coordinates whose values corrspond to
///   AnisoValTraits::standard_var_names().
Eigen::MatrixXd GlobalContinuousConfigDoFValues::standard_values() const {
  return m_info.basis() * m_vals;
}

/// DoFSetInfo provides the basis and symmetry representations for `values`
DoFSetInfo const &GlobalContinuousConfigDoFValues::info() const {
  return m_info;
}

}  // namespace CASM
