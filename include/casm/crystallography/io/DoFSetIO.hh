#ifndef DOFSETIO_HH
#define DOFSETIO_HH

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {
namespace xtal {
class SiteDoFSet;
class DoFSet;
}  // namespace xtal

class AnisoValTraits;

template <>
struct jsonConstructor<xtal::SiteDoFSet> {
  static xtal::SiteDoFSet from_json(const jsonParser &json,
                                    AnisoValTraits const &traits);
};

template <>
struct jsonConstructor<xtal::DoFSet> {
  static xtal::DoFSet from_json(const jsonParser &json,
                                AnisoValTraits const &traits);
};

jsonParser &to_json(xtal::SiteDoFSet const &_dof, jsonParser &json);
jsonParser &to_json(xtal::DoFSet const &_dof, jsonParser &json);
}  // namespace CASM

#endif
