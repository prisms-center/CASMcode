#include "casm/crystallography/Molecule.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  void from_json(AtomSpecies &_species, jsonParser const &json) {
    _species = AtomSpecies(json.get<std::string>());
  }

  jsonParser &to_json(AtomSpecies const &_species, jsonParser &json) {
    return json = _species.name();
  }


  //****************************************************
  //
  //****************************************************

  bool AtomPosition::identical(AtomPosition const &RHS, double _tol) const {
    return species() == RHS.species() && almost_equal(cart(), RHS.cart(), _tol);
  }

  //****************************************************
  //
  //****************************************************

  void AtomPosition::print(std::ostream &stream,
                           Eigen::Ref<const Eigen::Vector3d> const &trans,
                           Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
                           int spaces,
                           bool print_sd_flags /* = false */) const {
    for(int i = 0; i < spaces; i++) {
      stream << ' ';
    }

    stream << (cart2frac * (cart() + trans)).transpose();

    if(print_sd_flags) {
      for(int i = 0; i < 3; i++) {
        if(m_sd_flag[i]) stream << "  T";
        else stream << "  F";
      }
    }

    stream << "   " << name();

    return;
  }

  //****************************************************
  //
  //****************************************************

  AtomPosition &AtomPosition::apply_sym(const SymOp &op) {
    m_position = op.matrix() * m_position;
    return *this;
  }

  //****************************************************
  //
  //****************************************************

  jsonParser &to_json(AtomPosition const &apos, jsonParser &json, Eigen::Matrix3d const &c2f_mat) {
    json.put_obj();
    to_json_array(c2f_mat * apos.cart(), json["coordinate"]);
    json["species"] = apos.species();
    json["SD_flag"] = apos.sd_flag();
    return json;
  }

  //****************************************************
  //
  //****************************************************

  void from_json(AtomPosition &apos, const jsonParser &json, Eigen::Matrix3d const &f2c_mat) {
    std::string _name;
    Eigen::Vector3d _pos(0., 0., 0.);
    AtomPosition::sd_type _SD_flag({false, false, false});
    if(json.is_obj()) {
      _name = json["species"].get<std::string>();
      if(json.contains("coordinate"))
        _pos = f2c_mat * json["coordinate"].get<Eigen::Vector3d>();
      if(json.contains("SD_flag"))
        _SD_flag = json["SD_flag"].get<typename AtomPosition::sd_type>();
    }
    else if(json.is_string()) {
      _name = json.get<std::string>();
    }
    else
      throw std::runtime_error("Invalid JSON input encountered. Unable to parse AtomPosition object.\n");
    apos = AtomPosition(_pos, AtomSpecies(_name), _SD_flag);
    return;
  }

  /// \brief Read from json [b, i, j, k], using 'unit' for AtomPosition::unit()
  AtomPosition jsonConstructor<AtomPosition>::from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat) {
    std::string _name;
    Eigen::Vector3d _pos(0., 0., 0.);
    AtomPosition::sd_type _SD_flag({false, false, false});
    if(json.is_obj()) {
      _name = json["species"].get<std::string>();
      if(json.contains("coordinate"))
        _pos = f2c_mat * json["coordinate"].get<Eigen::Vector3d>();
      if(json.contains("SD_flag"))
        _SD_flag = json["SD_flag"].get<typename AtomPosition::sd_type>();
    }
    else if(json.is_string()) {
      _name = json.get<std::string>();
    }
    else
      throw std::runtime_error("Invalid JSON input encountered. Unable to parse AtomPosition object.\n");
    return AtomPosition(_pos, AtomSpecies(_name), _SD_flag);
  }

  bool Molecule::is_vacancy() const {
    //return m_atoms.empty();
    return CASM::is_vacancy(m_atoms[0].name());
  }

  //****************************************************
  // Applies symmetry operation to a Molecule
  //****************************************************

  Molecule &Molecule::apply_sym(SymOp const &op) {
    for(Index i = 0; i < size(); i++) {
      m_atoms[i].apply_sym(op);
    }
    for(auto it = m_attribute_map.begin(); it != m_attribute_map.end(); ++it)
      (it->second).apply_sym(op);
    return *this;
  }

  //****************************************************
  //
  //****************************************************

  bool Molecule::identical(Molecule const &RHS, double _tol) const {
    // compare number of atoms
    if(size() != RHS.size())
      return false;

    // compare number of attributes
    if(m_attribute_map.size() != RHS.m_attribute_map.size())
      return false;

    // compare atoms, irrespective of order
    for(Index i = 0; i < RHS.size(); i++) {
      Index j = 0;
      for(j = 0; j < size(); j++) {
        if(atom(i).identical(RHS.atom(i), _tol))
          break;
      }
      if(j == size())
        return false;
    }

    // compare attributes
    auto it(m_attribute_map.cbegin()), end_it(m_attribute_map.cend());
    for(; it != end_it; ++it) {
      auto it_RHS = RHS.m_attribute_map.find(it->first);
      if(it_RHS == RHS.m_attribute_map.cend() || !(it->second).identical(it_RHS->second, _tol))
        return false;
    }
    return true;
  }

  //****************************************************
  //
  //****************************************************

  bool Molecule::contains(std::string const &_name) const {
    for(Index i = 0; i < size(); i++)
      if(atom(i).name() == _name)
        return true;
    return false;
  }

  //****************************************************
  //
  //****************************************************

  void Molecule::print(std::ostream &stream,
                       Eigen::Ref<const Eigen::Vector3d> const &trans,
                       Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
                       int spaces,
                       char delim,
                       bool print_sd_flags /* = false */) const {
    for(Index i = 0; i < size(); i++) {
      atom(i).print(stream, trans, cart2frac, spaces, print_sd_flags);
      stream << delim;
    }
    return;
  }

  //****************************************************
  //   Write Molecule to json.
  //****************************************************

  jsonParser &Molecule::to_json(jsonParser &json, Eigen::Matrix3d const &f2c_mat) const {
    json.put_obj();
    CASM::to_json(m_atoms, json["atoms"], f2c_mat);
    json["attributes"] = m_attribute_map;
    json["name"] = name();
    return json;
  }

  //****************************************************
  //
  //    Read Molecule from json.
  //
  //****************************************************

  void Molecule::from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat) {
    m_atoms.clear();
    m_attribute_map.clear();
    if(json.contains("atoms")) {
      CASM::from_json(m_atoms, json["atoms"], f2c_mat);
    }

    if(json.contains("attributes")) {
      auto it = json["attributes"].cbegin(), end_it = json["attributes"].cend();
      for(; it != end_it; ++it) {
        auto result_pair = m_attribute_map.emplace(it.name(), it.name());
        (result_pair.first->second).from_json(*it);
      }
    }
  }

  //****************************************************
  /// \brief Return an atomic Molecule with specified name and Lattice
  //****************************************************

  Molecule Molecule::make_atom(std::string const &atom_name, AtomPosition::sd_type const &_sd_flags /*=AtomPosition::sd_type{false,false,false}*/) {
    //if(CASM::is_vacancy(atom_name))
    //  return make_vacancy();
    return Molecule(atom_name, {AtomPosition(0., 0., 0., atom_name)});
  }

  //****************************************************
  //
  //****************************************************

  Molecule Molecule::make_vacancy() {
    //return Molecule("Va", {});
    return make_atom("Va");
  }

  //****************************************************
  //
  //****************************************************

  jsonParser &to_json(const Molecule &mol, jsonParser &json, Eigen::Matrix3d const &c2f_mat) {
    return mol.to_json(json, c2f_mat);
  }

  //****************************************************
  //
  //****************************************************

  void from_json(Molecule &mol, const jsonParser &json, Eigen::Matrix3d const &f2c_mat) {
    mol.from_json(json, f2c_mat);
  }

  Molecule jsonConstructor<Molecule>::from_json(const jsonParser &json, Eigen::Matrix3d const &f2c_mat) {
    return json.get<Molecule>(f2c_mat);
  }


};
