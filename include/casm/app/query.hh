#ifndef QUERY_HH
#define QUERY_HH

#include <iostream>
#include "casm/casm_io/Log.hh"

namespace CASM {

  class PrimClex;
  
  int query_command(
    int argc,
    char *argv[],
    PrimClex *_primclex = nullptr,
    Log &log = default_log(),
    std::ostream &serr = std::cerr);

}

#endif
