#include "casm/app/enum/enumerate_configurations.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  void print_options(Log &log, ConfigEnumOptions const &options) {
    log << std::boolalpha;
    log.indent() << "primitive_only: " << options.primitive_only << std::endl;
    log.indent() << "filter: " << static_cast<bool>(options.filter) << std::endl;
    if(options.filter) {
      log.indent() << "filter expression: " << options.filter_expression << std::endl;
    }
    log.indent() << "verbosity: " << options.verbosity << std::endl;
    log.indent() << "dry_run: " << options.dry_run << std::endl;
    log.indent() << "output_configurations: " << options.output_configurations << std::endl;
    if(options.output_configurations) {
      auto const &output_options = options.output_options;
      log.increase_indent();
      log.indent() << "path: " << output_options.file_path << std::endl;
      log.indent() << "json: " << output_options.json_output << std::endl;
      log.indent() << "json_array: " << output_options.json_arrays << std::endl;
      log.indent() << "compress: " << output_options.compress << std::endl;
      log.indent() << "output_filtered_configurations: " << options.output_filtered_configurations << std::endl;
      log.decrease_indent();
    }
    log << std::noboolalpha;
  }

}
