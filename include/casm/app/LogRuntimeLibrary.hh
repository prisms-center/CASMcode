#ifndef LogRuntimeLibrary_HH
#define LogRuntimeLibrary_HH

#include <memory>
#include <string>

namespace CASM {

class RuntimeLibrary;

std::shared_ptr<RuntimeLibrary> log_make_shared_runtime_lib(
    std::string filename_base, std::string compile_options,
    std::string so_options, std::string compile_msg);
}  // namespace CASM

#endif
