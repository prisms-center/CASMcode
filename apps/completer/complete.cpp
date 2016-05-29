#include <iostream>
#include <iterator>
#include "casm/completer/complete.hh"
#include "completer_functions.hh"

int main(int argc, char *argv[]) {

  Completer::Engine casm_engine;

  CASM::PrimClex *_primclex = nullptr;

  //Generate a fake values to pass as arguments to the option descriptions
  std::string dumb_str;
  fs::path dumb_pth;
  std::vector<std::string> dumb_str_vec;
  bool dumb_bool;
  CASM::Index dumb_idx;

  po::options_description dumbquery("dumb query");
  Completer::add_query_options(dumbquery, dumb_str, dumb_pth, dumb_str_vec, dumb_str_vec, dumb_str_vec, dumb_bool, dumb_bool, dumb_bool, dumb_bool);
  Completer::Option query_opts("query", dumbquery);
  casm_engine.push_back(query_opts);

  po::options_description dumbmonte("dumb monte");
  Completer::add_monte_options(dumbmonte, dumb_pth, dumb_str, dumb_idx);
  Completer::Option monte_opts("monte", dumbmonte);
  casm_engine.push_back(monte_opts);



  std::cout << casm_engine.probe_options() << std::endl;
  std::cout << casm_engine.probe_suboptions("query") << std::endl;;
  std::cout << casm_engine.probe_suboptions("monte") << std::endl;;

  return 0;
}
