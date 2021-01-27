#ifndef CASM_FormattedDataFile
#define CASM_FormattedDataFile

#include <boost/filesystem.hpp>
#include <iostream>

#include "casm/casm_io/dataformatter/DataFormatterDecl.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/definitions.hh"

namespace CASM {

struct FormattedDataFileOptions {
  FormattedDataFileOptions(fs::path _file_path = "", bool _json_output = false,
                           bool _json_arrays = false, bool _compress = false);

  fs::path file_path;  // location to write
  bool json_output;    // if true output as json, else csv
  bool json_arrays;    // if true, and json_output==true, output json as arrays
  bool compress;       // if true, write compressed .gz file
};

template <typename DataObject>
class FormattedDataFile {
 public:
  FormattedDataFile(FormattedDataFileOptions const &options);

  FormattedDataFile(fs::path file_path, bool json_output, bool json_arrays,
                    bool compress);

  ~FormattedDataFile();

  void operator()(DataFormatter<DataObject> const &formatter,
                  DataObject const &object);

 private:
  bool m_initialized;
  jsonParser m_json;
  bool m_json_output;
  bool m_json_arrays;
  std::unique_ptr<std::ostream> m_ostream;
};

}  // namespace CASM

#endif
