#ifndef CASM_clex_Configuration_stream_io
#define CASM_clex_Configuration_stream_io

#include <iostream>

namespace CASM {

  /// Print Configuration as a VASP POSCAR
  ///
  /// Note: This makes use of VaspIO::PrintPOSCAR, which can be used directly for additional
  ///       formatting options.
  void print_poscar(Configuration const &configuration, std::ostream &sout);

}

#endif
