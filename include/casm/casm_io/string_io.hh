#ifndef CASM_string_io
#define CASM_string_io

#include <map>
#include <set>
#include <string>
#include <vector>

namespace CASM {

std::string to_string(const std::string &s, std::string begin = "\"",
                      std::string end = "\"") {
  return begin + s + end;
}

template <typename T1, typename T2>
std::string to_string(const std::pair<T1, T2> &value) {
  using namespace std;
  return to_string(value.first) + ":" + to_string(value.second);
}

template <typename Iterator>
std::string container_to_string(Iterator begin_it, Iterator end_it,
                                std::string begin = "[", std::string end = "]",
                                std::string delim = ", ") {
  std::stringstream out;
  using namespace std;
  out << begin;
  while (begin_it != end_it) {
    out << to_string(*begin_it);
    ++begin_it;
    if (begin_it != end_it) {
      out << delim;
    }
  }
  out << end;
  return out.str();
}

template <class T>
std::string to_string(const std::vector<T> &container, std::string begin = "[",
                      std::string end = "]", std::string delim = ", ") {
  return container_to_string(container.begin(), container.end(), begin, end,
                             delim);
}

template <class T, class Compare>
std::string to_string(const std::set<T, Compare> &container,
                      std::string begin = "[", std::string end = "]",
                      std::string delim = ", ") {
  return container_to_string(container.begin(), container.end(), begin, end,
                             delim);
}

template <class Key, class T, class Compare>
std::string to_string(const std::map<T, Key, Compare> &container,
                      std::string begin = "{", std::string end = "}",
                      std::string delim = ", ") {
  return container_to_string(container.begin(), container.end(), begin, end,
                             delim);
}

}  // namespace CASM

#endif
