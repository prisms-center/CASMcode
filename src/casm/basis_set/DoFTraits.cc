#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/basis_set/StrainDoFTraits.hh"
#include "casm/basis_set/DisplacementDoFTraits.hh"
#include "casm/basis_set/MagSpinDoFTraits.hh"

namespace CASM {

  template<>
  DoFType::TraitsDictionary make_parsing_dictionary<DoF::BasicTraits>() {
    DoF::register_traits(DoFType::occupation());
    DoF::register_traits(DoFType::displacement());
    DoF::register_traits(DoFType::magspin());
    DoF::register_traits(DoFType::EAstrain());
    DoF::register_traits(DoFType::Hstrain());
    DoF::register_traits(DoFType::GLstrain());
    DoFType::TraitsDictionary dict;

    dict.insert(
      DoFType::occupation(),
      DoFType::displacement(),
      DoFType::magspin(),
      DoFType::EAstrain(),
      DoFType::Hstrain(),
      DoFType::GLstrain());
    return dict;
  }
}
