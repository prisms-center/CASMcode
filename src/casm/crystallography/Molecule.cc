#include "casm/crystallography/Molecule.hh"

#include "casm/symmetry/SymOp.hh"

#include "casm/casm_io/json_io/container.hh"

namespace CASM {


  jsonParser &to_json(const Specie &specie, jsonParser &json) {
    json.put_obj();
    json["name"] = specie.name;
    json["mass"] = specie.mass;
    json["magmom"] = specie.magmom;
    json["U"] = specie.U;
    json["J"] = specie.J;
    return json;
  }

  void from_json(Specie &specie, const jsonParser &json) {
    try {
      from_json(specie.name, json["name"]);
      from_json(specie.mass, json["mass"]);
      from_json(specie.magmom, json["magmom"]);
      from_json(specie.U, json["U"]);
      from_json(specie.J, json["J"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

  //****************************************************

  //****************************************************
  /**
   *
   */
  //****************************************************

  bool AtomPosition::operator==(const AtomPosition &RHS) const {
    return specie == RHS.specie && Coordinate::operator==(RHS);


  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  void AtomPosition::print(std::ostream &stream, const Coordinate &trans, int spaces, bool SD_is_on) const {
    for(int i = 0; i < spaces; i++) {
      stream << ' ';
    }

    (Coordinate(*this) + trans).print(stream);

    if(SD_is_on) {
      for(int i = 0; i < 3; i++) {
        if(SD_flag[i]) stream << "  T";
        else stream << "  F";
      }
    }

    stream << "   " << specie.name; // << '\n';

    return;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  AtomPosition &AtomPosition::apply_sym(const SymOp &op) {
    Coordinate::apply_sym(op);
    return *this;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  AtomPosition &AtomPosition::apply_sym_no_trans(const SymOp &op) {
    Coordinate::apply_sym_no_trans(op);
    return *this;
  }


  jsonParser &to_json(const AtomPosition &apos, jsonParser &json) {
    json.put_obj();
    json["coordinate"].put<Coordinate>(apos);
    json["specie"] = apos.specie;
    json["SD_flag"] = apos.SD_flag;
    return json;
  }

  // Lattice must be set separately
  void from_json(AtomPosition &apos, const jsonParser &json) {
    try {
      Coordinate &coord_ref = apos;
      from_json(coord_ref, json["coordinate"]);
      from_json(apos.specie, json["specie"]);
      from_json(apos.SD_flag, json["SD_flag"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }



  //****************************************************
  /**
   * Applies symmetry operation to a Molecule
   *
   * @param SymOp &op
   */
  //****************************************************

  Molecule &Molecule::apply_sym(const SymOp &op) {
    // TODO: Is this the right way to apply_sym?

    center.apply_sym(op);
    for(Index i = 0; i < size(); i++) {
      //using apply_sym from Atomposition
      //since AtomPosition inherits from Coordinate
      at(i).apply_sym(op);
    }
    return *this;
  }

  //****************************************************
  /**
   * Applies symmetry operation to a Molecule
   * without adding translation
   * @param SymOp &op
   */
  //****************************************************

  Molecule &Molecule::apply_sym_no_trans(const SymOp &op) {
    // TODO: Is this the right way to apply_sym?

    center.apply_sym_no_trans(op);

    for(Index i = 0; i < size(); i++) {
      //using apply_sym from Coordinate
      //since AtomPosition inherits from Coordinate
      at(i).apply_sym_no_trans(op);
    }
    return *this;
  }


  //****************************************************
  /**
   *
   */
  //****************************************************

  void Molecule::set_lattice(const Lattice &new_lat, COORD_TYPE invariant_mode) {
    m_home = &new_lat;
    center.set_lattice(new_lat, invariant_mode);

    for(Index i = 0; i < size(); i++) {
      at(i).set_lattice(new_lat, invariant_mode);
    }
  };

  //****************************************************
  /**
   *
   */
  //****************************************************
  //TODO:
  bool Molecule::operator==(const Molecule &RHS) const {

    Index i, j;
    if(center != RHS.center)
      return false;

    for(i = 0; i < RHS.size(); i++) {

      for(j = 0; j < size(); j++) {
        if(RHS.at(i) == at(j))
          break;
      }

      if(j == size())
        return false;
    }
    return true;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************
  bool Molecule::contains(const std::string &name) const {
    for(Index i = 0; i < size(); i++)
      if(at(i).specie.name == name)
        return true;
    return false;
  }

  //****************************************************
  /**
   *
   */
  //****************************************************

  void Molecule::print(std::ostream &stream, const Coordinate &trans, int spaces, char delim, bool SD_is_on) const {
    for(Index i = 0; i < size(); i++) {
      at(i).print(stream, trans, spaces, SD_is_on);
      stream << delim;
    }
    return;
  }

  //****************************************************
  /**
   *   Write Molecule to json. Does not write lattice.
   */
  //****************************************************

  jsonParser &Molecule::to_json(jsonParser &json) const {
    json.put_obj();
    const Array<AtomPosition> &apos_array_ref = (*this);
    json["atomposition"] = apos_array_ref;
    if(point_group.size() > 0) {
      json["point_group"] = point_group;
    }
    json["center"] = center;
    json["name"] = name;
    return json;
  }

  //****************************************************
  /**
   *    Read Molecule from json. Assumes current lattice.
   */
  //****************************************************

  void Molecule::from_json(const jsonParser &json) {

    try {

      AtomPosition apos(*home());
      for(int i = 0; i < json["atomposition"].size(); i++) {
        CASM::from_json(apos, json["atomposition"][i]);
        push_back(apos);
      }

      point_group.clear();
      if(json.contains("point_group")) {
        for(int i = 0; i < json["point_group"].size(); i++) {
          point_group.push_back(json["point_group"][i].get<SymOp>());
        }
      }
      CASM::from_json(center, json["center"]);
      CASM::from_json(name, json["name"]);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }


  /// \brief Return an atomic Molecule with specified name and Lattice
  Molecule make_atom(std::string atom_name, const Lattice &lat) {
    Molecule atom(lat);
    atom.name = atom_name;
    atom.push_back(AtomPosition(0.0, 0.0, 0.0, atom.name, lat, CART));
    return atom;
  }

  /// \brief Return an vacancy Molecule with specified Lattice
  Molecule make_vacancy(const Lattice &lat) {
    return make_atom("Va", lat);
  }

  /// \brief Return true if Molecule name matches 'name', including Va name checks
  ///
  /// Equivalent to:
  /// \code
  /// ( (is_vacancy(mol.name) && is_vacancy(name)) || (mol.name == name) )
  /// \endcode
  ///
  bool is_molecule_name(const Molecule &mol, std::string name) {
    return (is_vacancy(mol.name) && is_vacancy(name)) || (mol.name == name);
  }

  jsonParser &to_json(const Molecule &mol, jsonParser &json) {
    return mol.to_json(json);
  }

  // Lattice must be set separately
  void from_json(Molecule &mol, const jsonParser &json) {
    try {
      mol.from_json(json);
    }
    catch(...) {
      /// re-throw exceptions
      throw;
    }
  }

};
