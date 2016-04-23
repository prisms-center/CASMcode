#ifndef QUERY_HH
#define QUERY_HH

#include <iostream>

namespace CASM {

  class PrimClex;

  int query_command(
    int argc,
    char *argv[],
    PrimClex *_primclex = nullptr,
    std::ostream &sout = std::cout,
    std::ostream &serr = std::cerr);

}

#endif
