#include "autotools.hh"

namespace autotools {

#ifndef ABS_SRCDIR
#define ABS_SRCDIR "BAD_ABS_SRCDIR_FLAG"
#endif

#ifndef ABS_TOP_BUILDDIR
#define ABS_TOP_BUILDDIR "BAD_ABS_TOP_BUILDDIR_FLAG"
#endif

std::string abs_srcdir() {
  static const std::string &srcdir = ABS_SRCDIR;
  return srcdir;
}

std::string abs_top_builddir() {
  static const std::string &top_builddir = ABS_TOP_BUILDDIR;
  return top_builddir;
}

std::string abs_ccasm_path() { return abs_top_builddir() + "/ccasm"; }

std::string abs_libdir() { return abs_top_builddir() + "/.libs"; }

std::string abs_includedir() { return abs_srcdir() + "/include"; }
}  // namespace autotools
