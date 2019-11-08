#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

  /// \brief Output internal state to JSON
  jsonParser const &SymGroupRepID::from_json(jsonParser const &json) {
    json.get_else(m_group_index, "group_index", Index(-1));
    json.get_else(m_rep_index, "group_index", Index(-1));
    return json;
  }

  jsonParser &to_json(SymGroupRepID const &_id, jsonParser &json) {
    json["group_index"] = _id.group_index();
    json["rep_index"] = _id.rep_index();
    return json;
  }

  jsonParser const &from_json(SymGroupRepID &_id, jsonParser const &json) {
    return _id.from_json(json);
  }
}
