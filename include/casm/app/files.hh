#ifndef CASM_files_HH
#define CASM_files_HH

#include <iostream>
#include "casm/casm_io/Log.hh"

namespace CASM {

  class PrimClex;
  
  int files_command(
    int argc,
    char *argv[],
    PrimClex *_primclex = nullptr,
    Log &log = default_log(),
    std::ostream &serr = std::cerr);

}

#endif
