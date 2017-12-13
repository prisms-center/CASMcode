#ifndef CASM_jsonIO_SpeciesSetParser
#define CASM_jsonIO_SpeciesSetParser

#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex_impl.hh"

#include "casm/casm_io/EnumIO.hh"
#include "casm/casm_io/InputParser_impl.hh"
#include "casm/misc/CRTPBase.hh"
#include "casm/clex/HasPrimClex.hh"

namespace CASM {

  enum class ALLOWED_SPECIES_TYPES {
    ALL, ATOM, MOLECULE
  };

  ENUM_IO_DECL(ALLOWED_SPECIES_TYPES)
  ENUM_TRAITS(ALLOWED_SPECIES_TYPES)

  std::string msg(ALLOWED_SPECIES_TYPES types);

  /// \brief Check that input species (atom or molecule) are recognized for given primclex
  class SpeciesSetParser : public KwargsParser, public HasPrimClex<CRTPBase<SpeciesSetParser>> {

  public:

    static const std::string require_all_help;
    static const std::string exclude_all_help;

    SpeciesSetParser(
      const PrimClex &_primclex,
      ALLOWED_SPECIES_TYPES _allowed_species_types,
      std::string _option_name,
      jsonParser &_input,
      fs::path _path,
      bool _required);

    std::set<std::string> values() const;

    const PrimClex &primclex() const {
      return m_primclex;
    }

  private:
    void require_valid_species();

    const PrimClex &m_primclex;
    ALLOWED_SPECIES_TYPES m_allowed_species_types;
    std::string m_option_name;
  };
}

#endif
