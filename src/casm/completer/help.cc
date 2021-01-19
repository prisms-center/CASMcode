#include "casm/completer/help.hh"

#include <sstream>

namespace CASM {

std::string dry_run_help() {
  return "  dry_run: bool (optional, default=false)\n"
         "    Perform dry run.\n\n";
}

std::string indent_space_help() {
  return "  indent_space: string (optional, default=6)\n"
         "    Number of spaces to indent for pretty-printing.\n\n";
}

std::string coordinate_mode_help() {
  return "  coordinate_mode: string (optional, default=FRAC)\n"
         "    Coordinate mode (FRAC, CART, INTEGRAL) for printing orbits.\n\n";
}

std::string orbit_print_mode_help() {
  return "  orbit_print_mode: string (optional, default=\"PROTO\")\n"
         "    Mode (FULL, PROTO) to select printing full orbits or just orbit "
         "prototypes.\n\n";
}

std::string prec_help(std::string what, int default_prec) {
  std::stringstream ss;
  ss << "  prec: int (optional, default=" << default_prec << ")\n"
     << "    Precision for printing " << what << "\n\n";
  return ss.str();
}

std::string verbosity_help() {
  return "  verbosity: string or int (optional, default=\"standard\")\n"
         "    Verbosity of output. Options are 'none', 'quiet', 'standard', "
         "'verbose', 'debug',\n"
         "    or an integer 0-100 (0: 'none', 100: 'debug').\n\n";
}

}  // namespace CASM
