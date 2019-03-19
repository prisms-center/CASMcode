#include "casm/kinetics/OccupationTransformation.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/Configuration.hh"
#include "casm/basis_set/DoF.hh"

#include "casm/kinetics/DoFTransformation_impl.hh"

namespace CASM {


  namespace Kinetics {

    template class DoFTransformation<CRTPBase<Kinetics::OccupationTransformation> >;

    // OccupationTransformation

    OccupationTransformation::OccupationTransformation(
      const UnitCellCoord &_uccoord,
      Index _from_value,
      Index _to_value) :
      from_value(_from_value),
      to_value(_to_value),
      uccoord(_uccoord) {}

    const UnitCellCoord::UnitType &OccupationTransformation::prim() const {
      return uccoord.unit();
    }

    const UnitCellCoord OccupationTransformation::coord() const {
      return uccoord;
    }

    const Index OccupationTransformation::from_occ() const {
      return from_value;
    }

    const Index OccupationTransformation::to_occ() const {
      return to_value;
    }

    const Molecule &OccupationTransformation::from_mol() const {
      return this->uccoord.sublat_site().site_occupant()[from_value];
    }

    const Molecule &OccupationTransformation::to_mol() const {
      return this->uccoord.sublat_site().site_occupant()[to_value];
    }

    bool OccupationTransformation::operator<(const OccupationTransformation &B) const {
      return _tuple() < B._tuple();
    }

    OccupationTransformation &OccupationTransformation::operator+=(UnitCell frac) {
      uccoord += frac;
      return *this;
    }

    OccupationTransformation &OccupationTransformation::apply_sym(const SymOp &op) {
      uccoord.apply_sym(op);
      return *this;
    }

    Configuration &OccupationTransformation::apply_to(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), to_value);
      return config;
    }

    Configuration &OccupationTransformation::apply_reverse_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), from_value);
      return config;
    }

    void OccupationTransformation::reverse() {
      using std::swap;
      swap(from_value, to_value);
    }

    std::tuple<UnitCellCoord, Index, Index> OccupationTransformation::_tuple() const {
      return std::make_tuple(uccoord, from_value, to_value);
    }

    /// \brief Print OccupationTransformation to stream, using default Printer<Kinetics::OccupationTransformation>
    std::ostream &operator<<(std::ostream &sout, const OccupationTransformation &trans) {
      Printer<Kinetics::OccupationTransformation> printer;
      Log out(sout);
      printer.print(trans, out);
      return sout;
    }
  }

  std::map<std::string, Index> empty_species_count(const UnitCellCoord::UnitType &prim) {
    auto species = struc_species(prim);
    std::map<std::string, Index> _species_count;
    for(const std::string &s : species) {
      _species_count[s] = 0;
    }
    return _species_count;
  }

  template<typename OccTransfIt>
  std::map<std::string, Index> from_species_count(OccTransfIt begin, OccTransfIt end) {
    if(begin == end) {
      return std::map<std::string, Index>();
    }
    std::map<std::string, Index> _species_count = empty_species_count(begin->prim());
    for(; begin != end; ++begin) {
      const OccupationTransformation &t = *begin;
      const Molecule &mol = t.uccoord.sublat_site().site_occupant()[t.from_value];
      for(const AtomPosition &species_pos : mol.atoms()) {
        _species_count[species_pos.name()]++;
      }
    }
    return _species_count;
  }
  typedef std::vector<OccupationTransformation>::iterator OccTransfVecIt;
  typedef std::vector<OccupationTransformation>::const_iterator OccTransfVecConstIt;
  template std::map<std::string, Index> from_species_count<OccTransfVecIt>(
    OccTransfVecIt begin,
    OccTransfVecIt end);
  template std::map<std::string, Index> from_species_count<OccTransfVecConstIt>(
    OccTransfVecConstIt begin,
    OccTransfVecConstIt end);

  template<typename OccTransfIt>
  std::map<std::string, Index> to_species_count(OccTransfIt begin, OccTransfIt end) {
    if(begin == end) {
      return std::map<std::string, Index>();
    }
    std::map<std::string, Index> _species_count = empty_species_count(begin->prim());
    for(; begin != end; ++begin) {
      const OccupationTransformation &t = *begin;
      const Molecule &mol = t.uccoord.sublat_site().site_occupant()[t.to_value];
      for(const AtomPosition &species_pos : mol.atoms()) {
        _species_count[species_pos.name()]++;
      }
    }
    return _species_count;
  }
  template std::map<std::string, Index> to_species_count<OccTransfVecIt>(
    OccTransfVecIt begin,
    OccTransfVecIt end);
  template std::map<std::string, Index> to_species_count<OccTransfVecConstIt>(
    OccTransfVecConstIt begin,
    OccTransfVecConstIt end);

  /// \brief Write OccupationTransformation to JSON object
  jsonParser &to_json(const Kinetics::OccupationTransformation &trans, jsonParser &json) {
    json.put_obj();
    json["uccoord"] = trans.uccoord;
    json["from_value"] = trans.from_value;
    json["to_value"] = trans.to_value;
    return json;
  }

  Kinetics::OccupationTransformation jsonConstructor<Kinetics::OccupationTransformation>::from_json(const jsonParser &json, const UnitCellCoord::UnitType &prim) {
    return Kinetics::OccupationTransformation(
             jsonConstructor<UnitCellCoord>::from_json(json["uccoord"], prim),
             json["from_value"].get<Index>(),
             json["to_value"].get<Index>());
  }

  void from_json(Kinetics::OccupationTransformation &fill_value, const jsonParser &read_json) {
    fill_value.from_value = read_json["from_value"].get<Index>();
    fill_value.to_value = read_json["to_value"].get<Index>();
  }

  const std::string Printer<Kinetics::OccupationTransformation>::element_name = "OccupationTransformation";

  void Printer<Kinetics::OccupationTransformation>::print(const Kinetics::OccupationTransformation &trans, Log &out) {
    if(!out.print()) {
      return;
    }
    COORD_MODE printer_mode(this->opt.coord_type);

    out << out.indent_str();
    out << trans.uccoord << " : ";
    out << trans.from_value << " (" << trans.from_mol().name() << ")";
    out << "  ->  ";
    out << trans.to_value << " (" << trans.to_mol().name() << ")";
    if(this->opt.delim)
      out << this->opt.delim;
    out << std::flush;
  }

}
