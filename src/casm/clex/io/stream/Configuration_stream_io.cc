#include "casm/clex/io/stream/Configuration_stream_io.hh"

#include "casm/clex/Configuration.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/crystallography/io/VaspIO.hh"

namespace CASM {

void print_poscar(Configuration const &configuration, std::ostream &sout) {
  VaspIO::PrintPOSCAR p{make_simple_structure(configuration),
                        configuration.name()};
  p.sort();
  p.print(sout);
}

}  // namespace CASM
