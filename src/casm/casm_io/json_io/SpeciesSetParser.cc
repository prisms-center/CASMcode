#include "casm/casm_io/json_io/SpeciesSetParser_impl.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  const std::string traits<ALLOWED_SPECIES_TYPES>::name = "allowed_species_types";

  const std::multimap<ALLOWED_SPECIES_TYPES, std::vector<std::string> > traits<ALLOWED_SPECIES_TYPES>::strval = {
    {ALLOWED_SPECIES_TYPES::ALL, {"ALL", "All", "all"} },
    {ALLOWED_SPECIES_TYPES::ATOM, {"ATOM", "Atom", "atom"} },
    {ALLOWED_SPECIES_TYPES::MOLECULE, {"MOLECULE", "Molecule", "molecule", "MOL", "Mol", "mol"} }
  };

  ENUM_IO_DEF(ALLOWED_SPECIES_TYPES)

  std::string msg(ALLOWED_SPECIES_TYPES types) {
    if(types == ALLOWED_SPECIES_TYPES::ALL) {
      return std::string("atom or molecule");
    }
    else if(types == ALLOWED_SPECIES_TYPES::ATOM) {
      return std::string("atom");
    }
    else if(types == ALLOWED_SPECIES_TYPES::MOLECULE) {
      return std::string("molecule");
    }
    throw std::invalid_argument("Error: Unknown ALLOWED_SPECIES_TYPES value");
  }


  SpeciesSetParser::SpeciesSetParser(
    const PrimClex &_primclex,
    ALLOWED_SPECIES_TYPES _allowed_species_types,
    std::string _option_name,
    jsonParser &_input,
    fs::path _path,
    bool _required) :
    KwargsParser(_input, _path, _required),
    m_primclex(_primclex),
    m_allowed_species_types(_allowed_species_types),
    m_option_name(_option_name) {
    if(exists()) {
      require_valid_species();
    }
  }

  std::set<std::string> SpeciesSetParser::values() const {
    if(exists()) {
      return self.get<std::set<std::string>>();
    }
    else {
      return std::set<std::string>();
    }
  }

  void SpeciesSetParser::require_valid_species() {
    auto ptr = this->optional_at<std::set<std::string>>(fs::path());
    if(ptr) {
      std::set<std::string> allowed_names;

      auto copy = [&](const std::vector<std::string> &names) {
        std::copy(names.begin(), names.end(), std::inserter(allowed_names, allowed_names.end()));
      };
      if(m_allowed_species_types == ALLOWED_SPECIES_TYPES::ATOM || m_allowed_species_types == ALLOWED_SPECIES_TYPES::ALL) {
        copy(prim().struc_species_name());
      }
      if(m_allowed_species_types == ALLOWED_SPECIES_TYPES::MOLECULE || m_allowed_species_types == ALLOWED_SPECIES_TYPES::ALL) {
        copy(prim().struc_molecule_name());
      }
      jsonParser json;

      for(std::string name : *ptr) {
        if(allowed_names.find(name) == allowed_names.end()) {
          error.insert(std::string("Error: '") + name + "' is not a recognized "
                       + msg(m_allowed_species_types) + " name.");
        }
      }
    }
  }

}
