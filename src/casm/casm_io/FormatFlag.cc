#include "casm/casm_io/FormatFlag.hh"

namespace CASM {
  int FormatFlag::iword_index() {
    static const int _index = std::ios_base::xalloc();
    return _index;
  }

  std::ostream &operator<<(std::ostream &_stream, const FormatFlag &flag) {
    return flag.set(_stream);
  }

}
