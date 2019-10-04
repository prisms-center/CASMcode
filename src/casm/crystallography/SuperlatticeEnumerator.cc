#include "casm/crystallography/SuperlatticeEnumerator.hh"

#include "casm/crystallography/HermiteCounter.hh"
#include "casm/external/Eigen/Dense"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  //*******************************************************************************************************************//
  //ScelEnumProps

  //*******************************************************************************************************************//
  //SuperlatticeIterator

  SuperlatticeIterator::SuperlatticeIterator(const SuperlatticeEnumerator &enumerator,
                                             int volume,
                                             int dims):
    m_super_updated(false),
    m_enum(&enumerator),
    m_current(notstd::make_cloneable<HermiteCounter>(volume, dims)) {
    if(enumerator.begin_volume() > enumerator.end_volume()) {
      throw std::runtime_error("The beginning volume of the SuperlatticeEnumerator cannot be greater than the end volume!");
    }

    if(dims < 1) {
      throw std::runtime_error("Dimensions to count over must be greater than 0!");
    }
  }

  SuperlatticeIterator &SuperlatticeIterator::operator=(const SuperlatticeIterator &B) {
    m_enum = B.m_enum;
    m_current = B.m_current;
    m_super_updated = false;

    m_canon_hist = B.m_canon_hist;
    return *this;
  }

  bool SuperlatticeIterator::operator==(const SuperlatticeIterator &B) const {
    return (m_enum == B.m_enum) && (matrix() - B.matrix()).isZero();
  }

  bool SuperlatticeIterator::operator!=(const SuperlatticeIterator &B) const {
    return !(*this == B);
  }

  typename SuperlatticeIterator::reference SuperlatticeIterator::operator*() const {
    if(!m_super_updated) {
      m_super = make_supercell(m_enum->unit(), matrix());
      m_super_updated = true;
    }
    return m_super;
  }

  typename SuperlatticeIterator::pointer SuperlatticeIterator::operator->() const {
    if(!m_super_updated) {
      m_super = make_supercell(m_enum->unit(), matrix());
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
    m_canon_hist.push_back(matrix());
    HermiteCounter::value_type last_determinant = m_current->determinant();
    ++(*m_current);

    if(last_determinant != m_current->determinant()) {
      m_canon_hist.clear();
    }

    while(std::find(m_canon_hist.begin(), m_canon_hist.end(), matrix()) != m_canon_hist.end()) {
      ++(*m_current);
    }

    m_super_updated = false;
  }

  Eigen::Matrix3i SuperlatticeIterator::matrix() const {
    Eigen::Matrix3i expanded = HermiteCounter_impl::_expand_dims((*m_current)(), m_enum->gen_mat());
    return canonical_hnf(expanded, m_enum->point_group(), m_enum->unit());
    /* return canonical_hnf(expanded, m_enum->point_group(), m_enum->lattice()); */
  }

  //*******************************************************************************************************************//
  //SuperlatticeEnumerator

  SuperlatticeEnumerator::SuperlatticeEnumerator(const Lattice &unit,
                                                 const SymGroupType &point_grp,
                                                 const ScelEnumProps &enum_props) :
    m_unit(unit),
    m_point_group(point_grp),
    m_begin_volume(enum_props.begin_volume()),
    m_end_volume(enum_props.end_volume()),
    m_gen_mat(enum_props.generating_matrix()),
    m_dims(enum_props.dims()) {

    if(m_gen_mat.determinant() < 1) {
      throw std::runtime_error("The transformation matrix to expand into a 3x3 matrix must have a positive determinant!");
    }
  }

  const Lattice &SuperlatticeEnumerator::unit() const {
    return m_unit;
  }

  const SuperlatticeEnumerator::SymGroupType &SuperlatticeEnumerator::point_group() const {
    return m_point_group;
  }

  const Eigen::Matrix3i &SuperlatticeEnumerator::gen_mat() const {
    return m_gen_mat;
  }

  int SuperlatticeEnumerator::dimension() const {
    return m_dims;
  }

  typename SuperlatticeEnumerator::size_type SuperlatticeEnumerator::begin_volume() const {
    return m_begin_volume;
  }

  typename SuperlatticeEnumerator::size_type SuperlatticeEnumerator::end_volume() const {
    return m_end_volume;
  }

  typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::begin() const {
    return const_iterator(*this, m_begin_volume, dimension());
  }

  typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::end() const {
    return const_iterator(*this, m_end_volume, dimension());
  }

  typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::cbegin() const {
    return const_iterator(*this, m_begin_volume, dimension());
  }

  typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::cend() const {
    return const_iterator(*this, m_end_volume, dimension());
  }

  typename SuperlatticeEnumerator::const_iterator SuperlatticeEnumerator::citerator(size_type volume) const {
    return SuperlatticeIterator(*this, volume, dimension());
  }

  //*******************************************************************************************************************//
  //Functions

  Eigen::Matrix3i make_unit_cell(const PrimClex &primclex, const jsonParser &input) {

    // read generating matrix (unit cell)
    Eigen::Matrix3i generating_matrix;
    if(input.is_null() || !input.contains("unit_cell")) {
      generating_matrix = Eigen::Matrix3i::Identity();
    }
    else if(input["unit_cell"].is_array()) {
      from_json(generating_matrix, input["unit_cell"]);
    }
    else if(input["unit_cell"].is_string()) {
      generating_matrix = primclex.db<Supercell>().find(input["unit_cell"].get<std::string>())->transf_mat();
    }
    else {
      throw std::invalid_argument(
        "Error reading unit cell from JSON input: 'unit_cell' must be a 3x3 integer matrix or supercell name");
    }
    return generating_matrix;
  }

  ScelEnumProps make_scel_enum_props(const PrimClex &primclex, const jsonParser &input) {

    // read volume range
    ScelEnumProps::size_type min_vol;
    ScelEnumProps::size_type max_vol;

    if(!input.contains("min")) {
      min_vol = 1;
    }
    else {
      from_json(min_vol, input["min"]);
    }
    if(!(min_vol > 0)) {
      throw std::invalid_argument(
        "Error in ScelEnumProps JSON input: 'min_volume' must be >0");
    }

    // read "max" scel size, or by default use largest existing supercell
    ScelEnumProps::size_type max_scel_size = 1;
    for(const auto &scel : primclex.db<Supercell>()) {
      if(scel.volume() > max_scel_size) {
        max_scel_size = scel.volume();
      }
    }

    input.get_else(max_vol, "max", max_scel_size);
    if(!(max_vol >= min_vol)) {
      throw std::invalid_argument(
        "Error in ScelEnumProps JSON input: 'max' must be greater than or equal to 'min'");
    }

    // read generating matrix (unit cell)
    Eigen::Matrix3i generating_matrix = make_unit_cell(primclex, input);

    std::string dirs;
    input.get_else<std::string>(dirs, "dirs", "abc");

    return ScelEnumProps(min_vol, max_vol + 1, dirs, generating_matrix);
  }

  jsonParser &to_json(const ScelEnumProps &props, jsonParser &json) {
    json.put_obj();
    json["min"] = props.begin_volume();
    json["max"] = props.end_volume() - 1;
    json["dirs"] = props.dirs();
    json["unit_cell"] = props.generating_matrix();
    return json;
  }

  ScelEnumProps jsonConstructor<ScelEnumProps>::from_json(const jsonParser &json, const PrimClex &primclex) {
    return make_scel_enum_props(primclex, json);
  }

  Eigen::Matrix3i enforce_min_volume(
    const Lattice &unit,
    const Eigen::Matrix3i &T,
    const SuperlatticeEnumerator::SymGroupType &point_grp,
    Index volume,
    bool fix_shape) {

    if(fix_shape) {
      auto init_vol = T.determinant();
      Index m = 1;
      while(m * m * m * init_vol < volume) {
        ++m;
      }

      return m * Eigen::Matrix3i::Identity();
    }
    else {
      auto init_vol = T.determinant();
      Index m = 1;
      while(m * init_vol < volume) {
        ++m;
      }

      auto compactness = [](const Lattice & lat) {
        Eigen::Matrix3d L = lat.lat_column_mat();
        return (L.transpose() * L).trace();
      };

      auto compare = [&](const Lattice & A, const Lattice & B) {
        return compactness(A) < compactness(B);
      };

      SuperlatticeEnumerator scel(unit, point_grp, ScelEnumProps(m, m + 1));
      auto best_it = std::min_element(scel.begin(), scel.end(), compare);
      return best_it.matrix();
    }

  }

  Eigen::Matrix3i canonical_hnf(const Eigen::Matrix3i &T, const SuperlatticeEnumerator::SymGroupType &effective_pg, const Lattice &ref_lattice) {
    Eigen::Matrix3d lat = ref_lattice.lat_column_mat();

    //get T in hermite normal form
    //H is the canonical form of the initial T matrix
    const Eigen::Matrix3i H = hermite_normal_form(T).first;

    //H_best will be the most canonical version and is returned
    Eigen::Matrix3i H_best;
    H_best = H;

    for(const auto &op : effective_pg) {
      Eigen::Matrix3i transformed = iround(lat.inverse() * Adapter::get_matrix(op) * lat) * H;
      Eigen::Matrix3i H_transformed = hermite_normal_form(transformed).first;

      //If you fall in here then transformed was greater than H
      if(HermiteCounter_impl::_canonical_compare(H_best, H_transformed) == 1) {
        H_best = H_transformed;
      }
    }

    return H_best;
    //return std::make_pair<Eigen::Matrix3i, Eigen::MatrixXd>(H_best, effective_pg[i_canon].matrix());
  }

  //*******************************************************************************************************************//



}
