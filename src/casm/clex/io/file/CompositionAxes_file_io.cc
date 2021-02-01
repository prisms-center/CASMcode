#include "casm/clex/io/file/CompositionAxes_file_io.hh"

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/CompositionAxes.hh"
#include "casm/clex/io/json/CompositionAxes_json_io.hh"

namespace CASM {

CompositionAxes read_composition_axes(fs::path _filename) {
  try {
    jsonParser json{_filename};
    return json.get<CompositionAxes>();

  } catch (...) {
    std::cerr << "Error reading composition axes from " << _filename
              << std::endl;
    throw;
  }
}

/// \brief Write CompositionAxes to file
void write_composition_axes(fs::path _filename,
                            CompositionAxes const& composition_axes) {
  SafeOfstream outfile;
  outfile.open(_filename);
  jsonParser json;
  to_json(composition_axes, json);
  json.print(outfile.ofstream());
  outfile.close();
}

}  // namespace CASM
