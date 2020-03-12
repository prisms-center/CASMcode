#ifndef SUPERLATTICEENUMERATORIO_HH
#define SUPERLATTICEENUMERATORIO_HH

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {
  namespace xtal {
    class ScelEnumProps;
  }

  jsonParser &to_json(const xtal::ScelEnumProps &props, jsonParser &json);
}

#endif
