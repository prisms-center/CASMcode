#ifndef CASM_algorithm
#define CASM_algorithm

#include <algorithm>

namespace CASM {

  /// \brief Equivalent to std::find(begin, end, value), but with custom comparison
  template<typename Iterator, typename T, typename BinaryCompare>
  Iterator find(Iterator begin, Iterator end, const T &value, BinaryCompare q) {
    return std::find_if(begin, end, [&value, q](const T & other)->bool {return q(value, other);});
  }

  /// \brief Equivalent to std::distance(begin, std::find(begin, end, value))
  template<typename Iterator, typename T>
  Index find_index(Iterator begin, Iterator end, const T &value) {
    return std::distance(begin, std::find(begin, end, value));
  }

  /// \brief Equivalent to std::distance(begin, find(begin, end, value, q))
  template<typename Iterator, typename T, typename BinaryCompare>
  Index find_index(Iterator begin, Iterator end, const T &value, BinaryCompare q) {
    return std::distance(begin, find(begin, end, value, q));
  }

  /// \brief Equivalent to std::distance(begin, std::find_if(begin, end, p))
  template<typename Iterator, typename UnaryPredicate>
  Index find_index_if(Iterator begin, Iterator end, UnaryPredicate p) {
    return std::distance(begin, std::find_if(begin, end, p));
  }

  /// \brief Equivalent to std::distance(begin, std::find_if_not(begin, end, q))
  template<typename Iterator, typename UnaryPredicate>
  Index find_index_if_not(Iterator begin, Iterator end, UnaryPredicate q) {
    return std::distance(begin, std::find_if_not(begin, end, q));
  }


  /// \brief Equivalent to std::distance(container.begin(), find(container.begin(), container.end(), value,q))
  template<typename Container, typename T, typename BinaryCompare>
  Index find_index(const Container &container, const T &value, BinaryCompare q) {
    return std::distance(container.begin(), find(container.begin(), container.end(), value, q));
  }

  /// \brief Equivalent to std::distance(container.begin(), std::find(container.begin(), container.end(), value))
  template<typename Container, typename T>
  Index find_index(const Container &container, const T &value) {
    return std::distance(container.begin(), std::find(container.begin(), container.end(), value));
  }

  /// \brief Equivalent to std::distance(container.begin(), std::find_if(container.begin(), container.end(), p))
  template<typename Container, typename UnaryPredicate>
  Index find_index_if(const Container &container, UnaryPredicate p) {
    return std::distance(container.begin(), std::find_if(container.begin(), container.end(), p));
  }

  /// \brief Equivalent to std::distance(container.begin(), std::find_if_not(container.begin(), container.end(), p))
  template<typename Container, typename UnaryPredicate>
  Index find_index_if_not(const Container &container, UnaryPredicate q) {
    return std::distance(container.begin(), std::find_if_not(container.begin(), container.end(), q));
  }


  /// \brief Equivalent to container.end() != std::find(container.begin(), container.end(), value)
  template<typename Container, typename T>
  bool contains(const Container &container, const T &value) {
    return container.end() != std::find(container.begin(), container.end(), value);
  }

  /// \brief Equivalent to container.end() != find(container.begin(), container.end(), value, q)
  template<typename Container, typename T, typename BinaryCompare>
  bool contains(const Container &container, const T &value, BinaryCompare q) {
    return container.end() != find(container.begin(), container.end(), value, q);
  }

  /// \brief Equivalent to container.end() != std::find_if(container.begin(), container.end(), p)
  template<typename Container, typename UnaryPredicate>
  bool contains_if(const Container &container, UnaryPredicate p) {
    return container.end() != std::find_if(container.begin(), container.end(), p);
  }

  /// \brief Equivalent to container.end() != std::find_if_not(container.begin(), container.end(), q)
  template<typename Container, typename UnaryPredicate>
  bool contains_if_not(const Container &container, UnaryPredicate q) {
    return container.end() != std::find_if_not(container.begin(), container.end(), q);
  }

}

#endif
