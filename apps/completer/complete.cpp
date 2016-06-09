#include <iostream>
#include <iterator>
#include "casm/completer/Complete.hh"
#include "casm/completer/Handlers.hh"
#include "completer_functions.hh"

using namespace CASM;

namespace CASM {
  namespace Completer {
    Engine build_casm_engine() {
      Completer::Engine casm_engine;

      QueryOption dumbquery;
      casm_engine.push_back(Option(dumbquery.tag(), dumbquery.desc()));

      MonteOption dumbmonte;
      casm_engine.push_back(Option(dumbmonte.tag(), dumbmonte.desc()));

      RunOption dumbrun;
      casm_engine.push_back(Option(dumbrun.tag(), dumbrun.desc()));

      BsetOption dumbbset;
      casm_engine.push_back(Option(dumbbset.tag(), dumbbset.desc()));

      CompositionOption dumbcomposition;
      casm_engine.push_back(Option(dumbcomposition.tag(), dumbcomposition.desc()));

      return casm_engine;
    }
  }
}

int main(int argc, char *argv[]) {


  CASM::PrimClex *_primclex = nullptr;

  Completer::Engine casm_engine = Completer::build_casm_engine();


  if(argc == 1) {
    std::cout << casm_engine.probe_options();
  }

  if(argc == 2) {
    std::cout << casm_engine.probe_suboptions(argv[1]);
  }

  if(argc == 3) {
    std::cout << casm_engine.probe_arguments(argv[1], argv[2]);
  }

  return 0;
}
