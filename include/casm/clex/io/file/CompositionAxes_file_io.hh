#ifndef CASM_clex_CompositionAxes_file_io
#define CASM_clex_CompositionAxes_file_io

#include "casm/global/definitions.hh"

namespace CASM {

struct CompositionAxes;

CompositionAxes read_composition_axes(fs::path _filename);

/// \brief Write CompositionAxes to file
void write_composition_axes(fs::path _filename,
                            CompositionAxes const& composition_axes);

}  // namespace CASM

#endif
