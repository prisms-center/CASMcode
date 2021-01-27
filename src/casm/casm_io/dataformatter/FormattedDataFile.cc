#include "casm/casm_io/dataformatter/FormattedDataFile.hh"

namespace CASM {

FormattedDataFileOptions::FormattedDataFileOptions(fs::path _file_path,
                                                   bool _json_output,
                                                   bool _json_arrays,
                                                   bool _compress)
    : file_path(_file_path),
      json_output(_json_output),
      json_arrays(_json_arrays),
      compress(_compress) {}

}  // namespace CASM
