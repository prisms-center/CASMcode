#include "casm/crystallography/SuperlatticeEnumerator.hh"

#include "casm/crystallography/HermiteCounter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

//*******************************************************************************************************************//
// ScelEnumProps

//*******************************************************************************************************************//
// SuperlatticeIterator

SuperlatticeIterator::SuperlatticeIterator(
    const SuperlatticeEnumerator &enumerator, int volume, int dims)
    : m_enum(&enumerator),
      m_current(notstd::make_cloneable<HermiteCounter>(volume, dims)),
      m_super_updated(false),
      m_matrix_updated(false) {
  if (enumerator.begin_volume() > enumerator.end_volume()) {
    throw std::runtime_error(
        "The beginning volume of the SuperlatticeEnumerator cannot be greater "
        "than the end volume!");
  }

  if (dims < 1) {
    throw std::runtime_error(
        "Dimensions to count over must be greater than 0!");
  }

  _advance_if_invalid();
}

SuperlatticeIterator &SuperlatticeIterator::operator=(
    const SuperlatticeIterator &B) {
  m_enum = B.m_enum;
  m_current = B.m_current;
  m_super_updated = false;
  m_matrix_updated = false;

  m_canon_hist = B.m_canon_hist;
  return *this;
}

bool SuperlatticeIterator::operator==(const SuperlatticeIterator &B) const {
  return (m_enum == B.m_enum) && ((volume() >= m_enum->end_volume() &&
                                   B.volume() >= m_enum->end_volume()) ||
                                  (matrix() - B.matrix()).isZero());
}

bool SuperlatticeIterator::operator!=(const SuperlatticeIterator &B) const {
  return !(*this == B);
}

typename SuperlatticeIterator::reference SuperlatticeIterator::operator*()
    const {
  if (!m_super_updated) {
    m_super = make_superlattice(m_enum->unit(), matrix());
    m_super_updated = true;
  }
  return m_super;
}

typename SuperlatticeIterator::pointer SuperlatticeIterator::operator->()
    const {
  if (!m_super_updated) {
    m_super = make_superlattice(m_enum->unit(), matrix());
    m_super_updated = true;
  }
  return &m_super;
}

HermiteCounter::value_type SuperlatticeIterator::volume() const {
  return m_current->determinant();
}

const SuperlatticeEnumerator &SuperlatticeIterator::enumerator() const {
  return *m_enum;
}

// prefix
SuperlatticeIterator &SuperlatticeIterator::operator++() {
  _increment();
  return *this;
}

void SuperlatticeIterator::_increment() {
  // check if this == end
  if (m_current->determinant() >= m_enum->end_volume()) {
    return;
  }
  m_canon_hist.push_back(matrix());
  _advance_one();
  _advance_if_invalid();
}

/// \brief Advance m_current by one, updating flags and history
void SuperlatticeIterator::_advance_one() {
  HermiteCounter::value_type last_determinant = m_current->determinant();
  ++(*m_current);
  if (last_determinant != m_current->determinant()) {
    m_canon_hist.clear();
  }
  m_super_updated = false;
  m_matrix_updated = false;
}

/// \brief Advance m_current if it is invalid, updating flags and history
void SuperlatticeIterator::_advance_if_invalid() {
  while (!_current_is_valid_and_unique()) {
    if (m_current->determinant() >= m_enum->end_volume()) {
      break;
    }
    _advance_one();
  }
}

Eigen::Matrix3i const &SuperlatticeIterator::matrix() const {
  if (!m_matrix_updated) {
    Eigen::Matrix3i expanded = HermiteCounter_impl::_expand_dims(
        m_current->current(), m_enum->gen_mat());
    m_matrix = canonical_hnf(expanded, m_enum->point_group(), m_enum->unit());
    m_matrix_updated = true;
  }
  return m_matrix;
}

/// \brief Check m_current for uniqueness and diagonal_only, fixed_shape
/// validity
bool SuperlatticeIterator::_current_is_valid_and_unique() const {
  Eigen::MatrixXi H = m_current->current();

  // check diagonal_only criteria
  if (m_enum->diagonal_only() && !H.isDiagonal()) {
    return false;
  }

  // check fixed_shape criteria
  if (m_enum->fixed_shape()) {
    int dim = m_enum->dimension();
    bool fixed_shape = true;
    for (int d = 1; d < dim; ++d) {
      fixed_shape = fixed_shape && H(d, d) == H(0, 0);
    }
    if (!fixed_shape) {
      return false;
    }
  }

  // check canonical hnf criteria for uniqueness
  if (std::find(m_canon_hist.begin(), m_canon_hist.end(), matrix()) !=
      m_canon_hist.end()) {
    return false;
  }

  return true;
}

//*******************************************************************************************************************//
// SuperlatticeEnumerator

SuperlatticeEnumerator::SuperlatticeEnumerator(const Lattice &unit,
                                               const SymOpVector &point_grp,
                                               const ScelEnumProps &enum_props)
    : m_unit(unit),
      m_point_group(point_grp),
      m_begin_volume(enum_props.begin_volume()),
      m_end_volume(enum_props.end_volume()),
      m_gen_mat(enum_props.generating_matrix()),
      m_dims(enum_props.dims()),
      m_diagonal_only(enum_props.diagonal_only()),
      m_fixed_shape(enum_props.fixed_shape()) {
  if (m_gen_mat.determinant() < 1) {
    throw std::runtime_error(
        "The transformation matrix to expand into a 3x3 matrix must have a "
        "positive determinant!");
  }
}

const Lattice &SuperlatticeEnumerator::unit() const { return m_unit; }

const SymOpVector &SuperlatticeEnumerator::point_group() const {
  return m_point_group;
}

const Eigen::Matrix3i &SuperlatticeEnumerator::gen_mat() const {
  return m_gen_mat;
}

int SuperlatticeEnumerator::dimension() const { return m_dims; }

bool SuperlatticeEnumerator::diagonal_only() const { return m_diagonal_only; }

bool SuperlatticeEnumerator::fixed_shape() const { return m_fixed_shape; }

typename SuperlatticeEnumerator::size_type
SuperlatticeEnumerator::begin_volume() const {
  return m_begin_volume;
}

typename SuperlatticeEnumerator::size_type SuperlatticeEnumerator::end_volume()
    const {
  return m_end_volume;
}

typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::begin()
    const {
  return const_iterator(*this, m_begin_volume, dimension());
}

typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::end()
    const {
  return const_iterator(*this, m_end_volume, dimension());
}

typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::cbegin()
    const {
  return const_iterator(*this, m_begin_volume, dimension());
}

typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::cend()
    const {
  return const_iterator(*this, m_end_volume, dimension());
}

typename SuperlatticeEnumerator::const_iterator
SuperlatticeEnumerator::citerator(size_type volume) const {
  return SuperlatticeIterator(*this, volume, dimension());
}

//*******************************************************************************************************************//
// Functions

Eigen::Matrix3i enforce_min_volume(const Lattice &unit,
                                   const Eigen::Matrix3i &T,
                                   const SymOpVector &point_grp, Index volume,
                                   bool fix_shape) {
  if (fix_shape) {
    auto init_vol = T.determinant();
    Index m = 1;
    while (m * m * m * init_vol < volume) {
      ++m;
    }

    return m * Eigen::Matrix3i::Identity();
  } else {
    auto init_vol = T.determinant();
    Index m = 1;
    while (m * init_vol < volume) {
      ++m;
    }

    auto compactness = [](const Lattice &lat) {
      Eigen::Matrix3d L = lat.lat_column_mat();
      return (L.transpose() * L).trace();
    };

    auto compare = [&](const Lattice &A, const Lattice &B) {
      return compactness(A) < compactness(B);
    };

    SuperlatticeEnumerator scel(unit, point_grp, ScelEnumProps(m, m + 1));
    auto best_it = std::min_element(scel.begin(), scel.end(), compare);
    return best_it.matrix();
  }
}

Eigen::Matrix3i canonical_hnf(const Eigen::Matrix3i &T,
                              const SymOpVector &effective_pg,
                              const Lattice &ref_lattice) {
  Eigen::Matrix3d lat = ref_lattice.lat_column_mat();

  // get T in hermite normal form
  // H is the canonical form of the initial T matrix
  const Eigen::Matrix3i H = hermite_normal_form(T).first;

  // H_best will be the most canonical version and is returned
  Eigen::Matrix3i H_best;
  H_best = H;

  for (const auto &op : effective_pg) {
    Eigen::Matrix3i transformed =
        iround(lat.inverse() * get_matrix(op) * lat) * H;
    Eigen::Matrix3i H_transformed = hermite_normal_form(transformed).first;

    // If you fall in here then transformed was greater than H
    if (HermiteCounter_impl::_canonical_compare(H_best, H_transformed) == 1) {
      H_best = H_transformed;
    }
  }

  return H_best;
  // return std::make_pair<Eigen::Matrix3i, Eigen::MatrixXd>(H_best,
  // effective_pg[i_canon].matrix());
}

}  // namespace xtal
}  // namespace CASM
