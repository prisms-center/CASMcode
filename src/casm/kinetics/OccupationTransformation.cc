#include "casm/kinetics/OccupationTransformation.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  namespace Kinetics {

    // OccupationTransformation

    OccupationTransformation::OccupationTransformation(
      const UnitCellCoord &uccoord,
      Index from_value,
      Index to_value) :
      SiteDoFTransformation(uccoord),
      from_value(from_value),
      to_value(to_value) {}

    OccupationTransformation &OccupationTransformation::apply_sym(const SymOp &op) {
      static_cast<SiteDoFTransformation &>(*this).apply_sym(op);
      return *this;
    }

    std::unique_ptr<OccupationTransformation> OccupationTransformation::clone() const {
      return std::unique_ptr<OccupationTransformation>(this->_clone());
    }

    Configuration &OccupationTransformation::apply_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), to_value);
      return config;
    }

    Configuration &OccupationTransformation::apply_reverse_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), from_value);
      return config;
    }

    void OccupationTransformation::reverse_impl() {
      using std::swap;
      swap(from_value, to_value);
    }

    OccupationTransformation *OccupationTransformation::_clone() const {
      return new OccupationTransformation(*this);
    }

    /// \brief Print OccupationTransformation to stream, using default Printer<Kinetics::OccupationTransformation>
    std::ostream &operator<<(std::ostream &sout, const OccupationTransformation &trans) {
      Printer<Kinetics::OccupationTransformation> printer;
      printer.print(trans, sout);
      return sout;
    }
  }

  /// \brief Write OccupationTransformation to JSON object
  jsonParser &to_json(const Kinetics::OccupationTransformation &trans, jsonParser &json) {
    json.put_obj();
    json["uccoord"] = trans.uccoord;
    json["from_value"] = trans.from_value;
    json["to_value"] = trans.to_value;
    return json;
  }

  Kinetics::OccupationTransformation jsonConstructor<Kinetics::OccupationTransformation>::from_json(const jsonParser &json, const Structure &prim) {
    return Kinetics::OccupationTransformation(
             json["uccoord"].get<UnitCellCoord>(prim),
             json["from_value"].get<Index>(),
             json["to_value"].get<Index>());
  }

  void from_json(Kinetics::OccupationTransformation &fill_value, const jsonParser &read_json) {
    fill_value.from_value = read_json["from_value"].get<Index>();
    fill_value.to_value = read_json["to_value"].get<Index>();
  }

  const std::string Printer<Kinetics::OccupationTransformation>::element_name = "OccupationTransformation";

  void Printer<Kinetics::OccupationTransformation>::print(const Kinetics::OccupationTransformation &trans, std::ostream &out) {
    COORD_MODE printer_mode(mode);

    out << indent() << indent() << indent();
    out << trans.uccoord << " : ";
    out << trans.from_value << " (" << trans.from_mol().name << ")";
    out << "  ->  ";
    out << trans.to_value << " (" << trans.to_mol().name << ")";
    if(delim)
      out << delim;
    out << std::flush;
  }

}

