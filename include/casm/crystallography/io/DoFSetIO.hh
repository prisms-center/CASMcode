#ifndef DOFSETIO_HH
#define DOFSETIO_HH

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {
  namespace xtal {
    class SiteDoFSet;
    class DoFSet;
  }

  template <>
  xtal::SiteDoFSet from_json<xtal::SiteDoFSet>(const jsonParser &json);
  template <>
  xtal::DoFSet from_json<xtal::DoFSet>(const jsonParser &json);

  jsonParser &to_json(xtal::SiteDoFSet const &_dof, jsonParser &json);
  jsonParser &to_json(xtal::DoFSet const &_dof, jsonParser &json);
}

#endif
