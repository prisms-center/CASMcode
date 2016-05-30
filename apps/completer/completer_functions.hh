#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  namespace Completer {
    void add_query_options
    (CASM::po::options_description &desc,
     std::string &selection_str,
     CASM::fs::path &out_path,
     std::vector<std::string> &columns,
     std::vector<std::string> &help_opt_vec,
     std::vector<std::string> &new_alias,
     bool &json_flag,
     bool &no_header,
     bool &verbatim_flag,
     bool &gz_flag);

    void add_monte_options
    (CASM::po::options_description &desc,
     CASM::fs::path &settings_path,
     std::string &verbosity_str,
     CASM::Index &condition_index);

  }
}
