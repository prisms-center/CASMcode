#include "casm/app/ClexDescription.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/misc/algorithm.hh"

namespace CASM {

void ClexDescription::print(std::ostream &sout, bool is_default,
                            int indent) const {
  std::string in(' ', indent);
  sout << in << name;
  if (is_default) {
    sout << "*";
  }
  sout << ": \n";
  sout << in << std::setw(16) << "property: " << property << "\n";
  sout << in << std::setw(16) << "calctype: " << calctype << "\n";
  sout << in << std::setw(16) << "ref: " << ref << "\n";
  sout << in << std::setw(16) << "bset: " << bset << "\n";
  sout << in << std::setw(16) << "eci: " << eci << "\n";
  sout << "\n";
}

ClexDescription default_configuration_clex() {
  std::string name = "formation_energy";
  std::string property = "formation_energy";
  std::string calctype = "default";
  std::string ref = "default";
  std::string bset = "default";
  std::string eci = "default";
  return ClexDescription{name, property, calctype, ref, bset, eci};
}

/// \brief Compare using name strings: A.name < B.name
bool operator<(const ClexDescription &A, const ClexDescription &B) {
  return A.name < B.name;
}

jsonParser &to_json(const ClexDescription &desc, jsonParser &json) {
  json.put_obj();
  json["name"] = desc.name;
  json["property"] = desc.property;
  json["calctype"] = desc.calctype;
  json["ref"] = desc.ref;
  json["bset"] = desc.bset;
  json["eci"] = desc.eci;
  return json;
}

void from_json(ClexDescription &desc, const jsonParser &json) {
  from_json(desc.name, json["name"]);
  from_json(desc.property, json["property"]);
  from_json(desc.calctype, json["calctype"]);
  from_json(desc.ref, json["ref"]);
  from_json(desc.bset, json["bset"]);
  from_json(desc.eci, json["eci"]);
}

bool new_dir(const DirectoryStructure &dir, ClexDescription const &desc) {
  bool result{true};
  result &= dir.new_bset_dir(desc.bset);
  result &= dir.new_calc_settings_dir(desc.calctype);
  result &= dir.new_ref_dir(desc.calctype, desc.ref);
  result &= dir.new_eci_dir(desc.property, desc.calctype, desc.ref, desc.bset,
                            desc.eci);
  return result;
}

bool clex_exists(const DirectoryStructure &dir, const ClexDescription &desc) {
  return contains(dir.all_calctype(), desc.calctype) &&
         contains(dir.all_ref(desc.calctype), desc.ref) &&
         contains(dir.all_bset(), desc.bset) &&
         contains(
             dir.all_eci(desc.property, desc.calctype, desc.ref, desc.bset),
             desc.eci);
}

}  // namespace CASM
