#include "casm/app/LogRuntimeLibrary.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  void print_runtime_lib_options_help(std::ostream &sout) {
    sout << "Error compiling clexulator. To fix: \n";
    sout << "  - Check compiler error messages.\n";
    sout << "  - Check compiler options with 'casm settings -l'\n";
    sout << "    - Update compiler options with 'casm settings --set-compile-options '...options...'\n";
    sout << "    - Make sure the casm headers can be found by including '-I/path/to/casm'\n";
  };

  /// Make shared_ptr<RuntimeLibrary>, logging progress and errors
  std::shared_ptr<RuntimeLibrary> log_make_shared_runtime_lib(
    std::string filename_base,
    std::string compile_options,
    std::string so_options,
    std::string compile_msg) {

    log().compiling<Log::standard>(filename_base + ".cc");
    log().begin_lap();
    log() << compile_msg << std::endl;
    try {
      std::shared_ptr<RuntimeLibrary> result = std::make_shared<RuntimeLibrary>(
                                                 filename_base,
                                                 compile_options,
                                                 so_options);
      log() << "compile time: " << log().lap_time() << " (s)\n" << std::endl;
      return result;
    }
    catch(runtime_lib_compile_error &e) {
      e.print(err_log());
      print_runtime_lib_options_help(log());
      throw;
    }
    catch(runtime_lib_shared_error &e) {
      e.print(err_log());
      print_runtime_lib_options_help(log());
      throw;
    }
    catch(std::exception &e) {
      print_runtime_lib_options_help(log());
      throw;
    }
  }
}
