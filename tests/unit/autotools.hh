#ifndef AUTOTOOLS_HH
#define AUTOTOOLS_HH

#include <string>

namespace autotools {
  std::string abs_srcdir();
  std::string abs_top_builddir();
  ///Get the full path to the compiled ccasm executable
  std::string abs_ccasm_path();
}

#endif
