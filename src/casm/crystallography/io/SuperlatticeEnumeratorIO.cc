#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/io/SuperlatticeEnumeratorIO.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/container/json_io.hh"

namespace CASM {
  jsonParser &to_json(const xtal::ScelEnumProps &props, jsonParser &json) {
    json.put_obj();
    json["min"] = props.begin_volume();
    json["max"] = props.end_volume() - 1;
    json["dirs"] = props.dirs();
    json["unit_cell"] = props.generating_matrix();
    return json;
  }
}
