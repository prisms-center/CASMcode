#ifndef CASM_support_enum_stream_io
#define CASM_support_enum_stream_io

#include "casm/casm_io/enum/io_traits.hh"

namespace CASM {

#define ENUM_IO_DECL(ENUM)                                       \
  std::ostream &operator<<(std::ostream &sout, const ENUM &val); \
                                                                 \
  std::istream &operator>>(std::istream &sin, ENUM &val);

#define ENUM_IO_DEF(ENUM)                                         \
  std::ostream &operator<<(std::ostream &sout, const ENUM &val) { \
    sout << to_string<ENUM>(val);                                 \
    return sout;                                                  \
  }                                                               \
                                                                  \
  std::istream &operator>>(std::istream &sin, ENUM &val) {        \
    std::string s;                                                \
    sin >> s;                                                     \
    val = from_string<ENUM>(s);                                   \
    return sin;                                                   \
  }

#define ENUM_IO_INLINE(ENUM)                                             \
  inline std::ostream &operator<<(std::ostream &sout, const ENUM &val) { \
    sout << to_string<ENUM>(val);                                        \
    return sout;                                                         \
  }                                                                      \
                                                                         \
  inline std::istream &operator>>(std::istream &sin, ENUM &val) {        \
    std::string s;                                                       \
    sin >> s;                                                            \
    val = from_string<ENUM>(s);                                          \
    return sin;                                                          \
  }

}  // namespace CASM

#endif
