#include "autotools.hh"

namespace autotools {

#ifndef SRCDIR
#define SRCDIR "BAD_FLAG"
#endif

  std::string srcdir_name() {
    static const std::string &srcdir = SRCDIR;
    return srcdir;
  }
}
