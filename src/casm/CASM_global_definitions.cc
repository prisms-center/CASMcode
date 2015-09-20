#ifndef CASM_GLOBAL_DEFINTIONS_CC
#define CASM_GLOBAL_DEFINTIONS_CC

//#include "casm/CASM_classes.hh"

#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {
  bool valid_index(Index i) {
    return 0 <= i;
  };

  std::istream &operator>>(std::istream &sin, COORD_TYPE &coord) {
    std::string s;
    sin >> s;
    if(s == "FRAC" || s == "0") {
      coord = FRAC;
    }
    else if(s == "CART" || s == "1") {
      coord = CART;
    }
    else if(s == "COORD_DEFAULT" || s == "2") {
      coord = COORD_DEFAULT;
    }
    return sin;
  }

  /// Gets string value as COORD_TYPE
  void from_json(COORD_TYPE &value, const jsonParser &json) {
    try {
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
    catch(...) {
      /// re-throw exceptions
      throw;
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
    try {
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
        std::cout << "Could not recognize PERIODICITY_TYPE from json: '" << periodicity << "'" << std::endl;
        exit(1);
      }

    }
    catch(...) {
      /// re-throw exceptions
      throw;
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
    try {
      std::string type = from_json<std::string>(json);

      if(type == "PRIM") {
        value = PRIM;
      }
      else if(type == "SCEL") {
        value = SCEL;
      }
      else {
        std::cout << "Could not recognize CELL_TYPE from json: '" << type << "'" << std::endl;
        exit(1);
      }

    }
    catch(...) {
      /// re-throw exceptions
      throw;
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

  jsonParser &to_json(const MTRand &twister, jsonParser &json) {
    std::cerr << "UNIMPLEMENTED ROUTINE to_json for MTRand" << std::endl;
    std::cerr << "Pssst, look in CASM_global_definitions.cc" << std::endl;
    exit(1);
  }


  void print_splash(std::ostream &out) {

    out << "       .::::::::.        .:::::.         .::::::.          .:.       .:.     \n"
        << "     .:::::::::::.      .:::::::.      .::::::::::.       .:::.     .:::.    \n"
        << "    .:::'     ':::.    .:::' ':::.    .:::'   '::::      .:::::.   .:::::.   \n"
        << "    ::::       ::::   .:::'   ':::.   ::::     '::'     .:::::::. .:::::::.  \n"
        << "    ::::              ::::     ::::   '::::.            ::::'':::.:::''::::  \n"
        << "    ::::              ::::     ::::    '::::::.         ::::  ':::::'  ::::  \n"
        << "    ::::              :::::::::::::      ''::::::.      ::::   ':::'   ::::  \n"
        << "    ::::              :::::::::::::         '::::::     ::::    ':'    ::::  \n"
        << "    ::::      .::::   ::::     ::::            :::::    ::::           ::::  \n"
        << "    ':::.    .::::'   ::::     ::::   .::.     .::::    ::::           ::::  \n"
        << "     '::::::::::'     ::::     ::::   :::::...:::::'    ::::           ::::  \n"
        << "        ':::::'       '::'     '::'   ':::::::::::'     '::'           '::'  \n";
  }


};

#endif

