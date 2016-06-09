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

      //Generate a fake values to pass as arguments to the option descriptions
      std::string dumb_str;
      fs::path dumb_pth;
      std::vector<std::string> dumb_str_vec;
      bool dumb_bool;
      CASM::Index dumb_idx;

      //*******************************************************//

      QueryOption dumbquery("query");
      Completer::Option query_opts("query", dumbquery.desc());
      casm_engine.push_back(query_opts);

      //*******************************************************//

      MonteOption dumbmonte("monte");
      Completer::Option monte_opts("monte", dumbmonte.desc());
      casm_engine.push_back(monte_opts);

      //*******************************************************//

      RunOption dumbrun("run");
      Completer::Option run_opts("run", dumbrun.desc());
      casm_engine.push_back(run_opts);

      //*******************************************************//

      BsetOption dumbbset("bset");
      Completer::Option bset_opts("bset", dumbbset.desc());
      casm_engine.push_back(bset_opts);

      //*******************************************************//

      CompositionOption dumbcomposition("composition");
      Completer::Option composition_opts("composition", dumbcomposition.desc());
      casm_engine.push_back(composition_opts);

      //*******************************************************//

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
