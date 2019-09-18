#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
  AnisoValTraits AnisoValTraits::disp() {
    return AnisoValTraits("disp",
    {"x", "y", "z"},
    LOCAL,
    CartesianSymRepBuilder(),
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::force() {
    return AnisoValTraits("force",
    {"fx", "fy", "fz"},
    LOCAL,
    CartesianSymRepBuilder(),
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::strain(std::string const &_prefix) {
    return AnisoValTraits(_prefix + "strain",
    {"e_1", "e_2", "e_3", "e_4", "e_5", "e_6"},
    GLOBAL,
    Rank2TensorSymRepBuilder(),
    {},
    {"atomize", "disp"});
  }

  AnisoValTraits AnisoValTraits::magspin() {
    return AnisoValTraits("magspin",
    {"sx", "sy", "sz"},
    LOCAL | UNIT_LENGTH,
    AngularMomentumSymRepBuilder(),
    {},
    {"atomize"});
  }

  AnisoValTraits AnisoValTraits::magmom() {
    return AnisoValTraits("magmom",
    {"mx", "my", "mz"},
    LOCAL | UNIT_LENGTH,
    TimeReversalSymRepBuilder(),
    {},
    {"atomize"});
  }

}
