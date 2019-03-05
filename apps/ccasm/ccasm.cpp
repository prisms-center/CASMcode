#include <iostream>
#include "casm/app/casm_functions.hh"

using namespace CASM;

// ///////////////////////////////////////
// ccasm main:

int main(int argc, char *argv[]) {
  try {
    PrimClex *_primclex = nullptr;
    CommandArgs args(argc, argv, _primclex, fs::path(), default_log(), default_err_log());

    return casm_api(args);
  }
  catch(std::exception const &e) {
    default_log() << "Uncaught exception: \n" << e.what();
  }
  return 1;
}
