#include <iostream>
#include "casm/app/casm_functions.hh"
#include "casm/casm_io/Log.hh"

using namespace CASM;

// ///////////////////////////////////////
// ccasm main:

int main(int argc, char *argv[]) {
  try {
    PrimClex *_primclex = nullptr;
    CommandArgs args(argc, argv, _primclex, fs::path());

    return casm_api(args);
  }
  catch(std::exception const &e) {
    log() << "Uncaught exception: \n" << e.what();
  }
  return 1;
}
