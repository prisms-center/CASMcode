#include <iostream>
#include <iterator>
#include "casm/completer/complete.hh"
#include "casm/completer/handlers.hh"
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

      po::options_description dumbquery("dumb query");
      Completer::add_query_options(dumbquery, dumb_str, dumb_pth, dumb_str_vec, dumb_str_vec, dumb_str_vec, dumb_bool, dumb_bool, dumb_bool, dumb_bool);
      Completer::Option query_opts("query", dumbquery);
      casm_engine.push_back(query_opts);

      //*******************************************************//

      MonteOption dumbmonte("monte");
      Completer::Option monte_opts("monte", dumbmonte.desc());
      casm_engine.push_back(monte_opts);

      //*******************************************************//

      po::options_description dumbrun("dumb run");
      Completer::add_run_options(dumbrun, dumb_str, dumb_str);
      Completer::Option run_opts("run", dumbrun);
      casm_engine.push_back(run_opts);

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
