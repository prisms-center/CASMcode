#ifndef AUTOTOOLS_HH
#define AUTOTOOLS_HH

#include <string>

namespace autotools {
std::string abs_srcdir();
std::string abs_top_builddir();
/// Get the full path to the compiled ccasm executable
std::string abs_ccasm_path();
/// Get the full path to the directory containing libcasm, libccasm, etc at
/// build time
std::string abs_libdir();
/// Get the full path to the include directory of the *repository*.
/// This is the directory that gets installed on the system upon `make install`
std::string abs_includedir();
}  // namespace autotools

#endif
