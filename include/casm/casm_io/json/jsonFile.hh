#ifndef CASM_jsonFile
#define CASM_jsonFile

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

/// \brief Slightly more compact construction of jsonParser from file
///
/// \code
/// // enables using:
/// jsonFile json_a("path/to/file.json");
///
/// // instead of:
/// jsonParser json_b(fs::path("path/to/file.json"));
/// \endcode
///
class jsonFile : public jsonParser {
 public:
  /// \brief Convert fs::path from t, then construct jsonParser
  template <typename T>
  jsonFile(const T &t) : jsonParser(fs::path(t)) {}
};

}  // namespace CASM

#endif
