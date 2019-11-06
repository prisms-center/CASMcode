#ifndef CASM_support_io_json_enum
#define CASM_support_io_json_enum

#include "casm/casm_io/traits/enum.hh"

namespace CASM {

  class jsonParser;

#define ENUM_JSON_IO_DECL(ENUM) \
  jsonParser &to_json(const ENUM &val, jsonParser &json); \
  \
  void from_json(ENUM& val, const jsonParser& json); \

#define ENUM_JSON_IO_DEF(ENUM) \
  jsonParser &to_json(const ENUM &val, jsonParser &json) { \
    return to_json(to_string<ENUM>(val), json); \
  } \
  \
  void from_json(ENUM& val, const jsonParser& json) { \
    val = from_string<ENUM>(json.get<std::string>()); \
  } \

#define ENUM_JSON_IO_INLINE(ENUM) \
  inline jsonParser &to_json(const ENUM &val, jsonParser &json) { \
    return to_json(to_string<ENUM>(val), json); \
  } \
  \
  inline void from_json(ENUM& val, const jsonParser& json) { \
    val = from_string<ENUM>(json.get<std::string>()); \
  } \


}

#endif
