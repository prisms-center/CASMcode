#ifndef CASM_support_container_stream_io
#define CASM_support_container_stream_io

#include <iostream>
#include <sstream>
#include <vector>

namespace CASM {
template <class T>
std::istream &operator>>(std::istream &_in, std::vector<T> &vec) {
  std::string line;
  std::getline(_in, line, '\n');
  std::stringstream tss(line);
  T tval;
  while (_in) {
    _in >> tval;
    vec.push_back(tval);
  }

  return _in;
}
}  // namespace CASM

namespace std {
template <class T>
ostream &operator<<(ostream &out, const vector<T> &vec) {
  if (vec.size() == 0) out << "[empty]  ";
  for (auto it = vec.cbegin(); it != vec.cend(); ++it) {
    out << *it << "  ";
  }
  return out;
}
}  // namespace std

#endif
