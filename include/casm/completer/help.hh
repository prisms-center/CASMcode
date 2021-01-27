#ifndef CASM_completer_help
#define CASM_completer_help

#include <string>

namespace CASM {

std::string dry_run_help();
std::string coordinate_mode_help();
std::string indent_space_help();
std::string orbit_print_mode_help();
std::string prec_help(std::string what, int default_prec);
std::string verbosity_help();
}  // namespace CASM

#endif
