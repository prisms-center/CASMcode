#include "casm/casm_io/json_io/global.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  /// Gets string value as COORD_TYPE
  ///
  /// If first letter is:
  /// - 'c', 'C', 'k', or 'K' -> CART
  /// - 'd', 'D', 'f', or 'F' -> FRAC
  /// - else: COORD_DEFAULT
  ///
  void from_json(COORD_TYPE &value, const jsonParser &json) {
    std::string c_case = from_json<std::string>(json);

    char c = c_case[0];
    c = tolower(c);

    if((c == 'c') || (c == 'k')) {
      value = CART;
    }
    else if((c == 'd') || (c == 'f')) {
      value = FRAC;
    }
    else {
      value = COORD_DEFAULT;
    }

  }


  /// Puts COORD_TYPE value (as string)
  jsonParser &to_json(const COORD_TYPE &value, jsonParser &json) {
    if(value == FRAC) {
      return to_json("FRAC", json);
    }
    else if(value == CART) {
      return to_json("CART", json);
    }
    else {
      return to_json("COORD_DEFAULT", json);
    }
  }

  /// Gets string value as PERIODICITY_TYPE
  void from_json(PERIODICITY_TYPE &value, const jsonParser &json) {
    std::string periodicity = from_json<std::string>(json);

    if(periodicity == "PERIODIC") {
      value = PERIODIC;
    }
    else if(periodicity == "LOCAL") {
      value = LOCAL;
    }
    else if(periodicity == "PERIODICITY_DEFAULT") {
      value = PERIODICITY_DEFAULT;
    }
    else {
      std::stringstream ss;
      ss << "Could not recognize PERIODICITY_TYPE from json: '" << periodicity << "'";
      throw std::runtime_error(ss.str());
    }
  }


  /// Puts PERIODICITY_TYPE value (as string)
  jsonParser &to_json(const PERIODICITY_TYPE &value, jsonParser &json) {
    if(value == PERIODIC) {
      return to_json("PERIODIC", json);
    }
    else if(value == LOCAL) {
      return to_json("LOCAL", json);
    }
    else {
      return to_json("PERIODICITY_DEFAULT", json);
    }
  }

  /// Gets string value as CELL_TYPE
  void from_json(CELL_TYPE &value, const jsonParser &json) {
    std::string type = from_json<std::string>(json);

    if(type == "PRIM") {
      value = PRIM;
    }
    else if(type == "SCEL") {
      value = SCEL;
    }
    else {
      std::stringstream ss;
      ss << "Could not recognize CELL_TYPE from json: '" << type << "'";
      throw std::runtime_error(ss.str());
    }
  }


  /// Puts CELL_TYPE value (as string)
  jsonParser &to_json(const CELL_TYPE &value, jsonParser &json) {
    if(value == PRIM) {
      return to_json("PRIM", json);
    }
    else {
      return to_json("SCEL", json);
    }

  }

}

