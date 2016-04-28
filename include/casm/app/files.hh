#ifndef CASM_files_HH
#define CASM_files_HH

#include <iostream>

namespace CASM {

  class PrimClex;

  int files_command(
    int argc,
    char *argv[],
    PrimClex *_primclex = nullptr,
    std::ostream &sout = std::cout,
    std::ostream &serr = std::cerr);

}

#endif
