#ifndef CASM_TMP_HH
#define CASM_TMP_HH

#include <cassert>
#include <iterator>
#include <string>
#include <tuple>
#include <type_traits>

#include "casm/misc/CRTPBase.hh"

namespace CASM {

template <typename T>
struct traits;

template <bool B, class T = void>
using enable_if_t = typename std::enable_if<B, T>::type;

namespace CASM_TMP {

// ---------------------

/// \brief Enables calling a function on each argument in a parameter pack
template <typename... Args>
void ignore_returnvalues(Args &&...) {}

// ---------------------

/// \brief Alias for void, to help SFINAE work
template <typename... Ts>
struct make_void {
  typedef void type;
};

/// \brief Alias for void, to help SFINAE work
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;

/// \brief Base type inherits from std::false_type if T is not iterator
///
/// - T is considered an iterator if it is incrementable, dereferenceable, and
/// comparable
template <typename T, typename = void>
struct is_iterator : std::false_type {};

/// \brief Specialized case inherits from std::true_type if T is an iterator
///
/// - T is considered an iterator if it is incrementable, dereferenceable, and
/// comparable
template <typename T>
struct is_iterator<
    T, void_t<decltype(++std::declval<T &>()), decltype(*std::declval<T &>()),
              decltype(std::declval<T &>() == std::declval<T &>())>>
    : std::true_type {};

/// \brief Template alias to enable a template function via SFINAE for any
/// iterator
///
/// - example function argument: 'typename
/// CASM_TMP::enable_if_iterator<IteratorType>::type* = nullptr'
template <typename T>
using enable_if_iterator = std::enable_if<is_iterator<T>::type::value, void>;

/// \brief Template alias to enable a template function via SFINAE for an
/// iterator with value_type V
///
/// - example function argument: 'typename
/// CASM_TMP::enable_if_iterator_of<IteratorType, ValueType>::type* = nullptr'
template <typename T, typename V>
using enable_if_iterator_of = std::enable_if<
    is_iterator<T>::type::value &&
        std::is_same<typename std::iterator_traits<T>::value_type,
                     V>::type::value,
    void>;

template <bool IsConst, typename T>
using ConstSwitch = typename std::conditional<IsConst, const T, T>::type;

// ---------------------

/// \brief Unary transformation that behaves as Identity (i.e. transform(arg) ==
/// arg is true)
template <typename T>
struct UnaryIdentity {
  T operator()(T const &arg) const { return arg; }
};

/// \brief N-nary function that behaves as a constant (i.e.
/// transform(arg1,arg2,...) == constant is true)
template <typename OutputType>
struct ConstantFunctor {
  ConstantFunctor(OutputType const &_const) : m_const(_const) {}

  template <typename... Args>
  OutputType operator()(Args const &... args) const {
    return m_const;
  }

 private:
  OutputType m_const;
};

// ---------------------

/// \brief Helper Functor for Counter container access using operator[]
template <typename Container,
          typename _value_type = typename Container::value_type,
          typename _size_type = typename Container::size_type>
struct BracketAccess {
  typedef _value_type value_type;
  typedef _size_type size_type;

  BracketAccess() {}
  BracketAccess &operator=(const BracketAccess &) { return *this; };

  /// \brief Identical to 'return container[index];'
  value_type &operator()(Container &container, size_type index) const {
    return container[index];
  }

  /// \brief Identical to 'return container[index];'
  const value_type &operator()(const Container &container,
                               size_type index) const {
    return container[index];
  }

  /// \brief Identical to 'return container[i][j];'
  value_type &operator()(Container &container, size_type i, size_type j) const {
    return container[i][j];
  }

  /// \brief Identical to 'return container[i][j];'
  const value_type &operator()(const Container &container, size_type i,
                               size_type j) const {
    return container[i][j];
  }

  /// \brief Identical to 'return container[index];'
  static value_type &at(Container &container, size_type index) {
    return container[index];
  }

  /// \brief Identical to 'return container[index];'
  static const value_type &at(const Container &container, size_type index) {
    return container[index];
  }

  /// \brief Identical to 'return container[i][j];'
  static value_type &at(Container &container, size_type i, size_type j) {
    return container[i][j];
  }

  /// \brief Identical to 'return container[i][j];'
  static const value_type &at(const Container &container, size_type i,
                              size_type j) {
    return container[i][j];
  }
};

/// \brief Helper Functor for Counter container access using operator()
template <typename Container,
          typename _value_type = typename Container::value_type,
          typename _size_type = typename Container::size_type>
struct ParenthesesAccess {
  typedef _value_type value_type;
  typedef _size_type size_type;

  ParenthesesAccess() {}
  ParenthesesAccess &operator=(const ParenthesesAccess &) { return *this; };

  /// \brief Identical to 'return container(index);'
  value_type &operator()(Container &container, size_type index) const {
    return container(index);
  }

  /// \brief Identical to 'return container(index);'
  const value_type &operator()(const Container &container,
                               size_type index) const {
    return container(index);
  }

  /// \brief Identical to 'return container(index);'
  value_type &operator()(Container &container, size_type i, size_type j) const {
    return container(i, j);
  }

  /// \brief Identical to 'return container(index);'
  const value_type &operator()(const Container &container, size_type i,
                               size_type j) const {
    return container(i, j);
  }

  /// \brief Identical to 'return container(index);'
  static value_type &at(Container &container, size_type index) {
    return container(index);
  }

  /// \brief Identical to 'return container(index);'
  static const value_type &at(const Container &container, size_type index) {
    return container(index);
  }

  /// \brief Identical to 'return container(i,j);'
  static value_type &at(Container &container, size_type i, size_type j) {
    return container(i, j);
  }

  /// \brief Identical to 'return container(i,j);'
  static const value_type &at(const Container &container, size_type i,
                              size_type j) {
    return container(i, j);
  }
};

// ----------------------------------------

template <int I>
struct for_each_type_impl {
  template <typename TupleType, typename F>
  static void eval(F f) {
    typedef typename std::tuple_element<I - 1, TupleType>::type ElementType;
    f.template eval<ElementType>();
    for_each_type_impl<I - 1>::template eval<TupleType, F>(f);
  }
};

template <>
struct for_each_type_impl<1> {
  template <typename TupleType, typename F>
  static void eval(F f) {
    typedef typename std::tuple_element<0, TupleType>::type ElementType;
    return f.template eval<ElementType>();
  }
};

template <typename TupleType, typename F>
void for_each_type(F f = F()) {
  return for_each_type_impl<std::tuple_size<TupleType>::value>::template eval<
      TupleType, F>(f);
}

// ----------------------------------------

template <int I>
struct for_type_impl {
  template <typename TupleType, typename F>
  static void eval(std::string name, F f) {
    typedef typename std::tuple_element<I - 1, TupleType>::type ElementType;
    if (traits<ElementType>::name == name) {
      f.template eval<ElementType>();
    }
    for_type_impl<I - 1>::template eval<TupleType, F>(name, f);
  }
};

template <>
struct for_type_impl<1> {
  template <typename TupleType, typename F>
  static void eval(std::string name, F f) {
    typedef typename std::tuple_element<0, TupleType>::type ElementType;
    if (traits<ElementType>::name == name) {
      f.template eval<ElementType>();
    }
  }
};

template <typename TupleType, typename F>
void for_type(std::string name, F f = F()) {
  return for_type_impl<std::tuple_size<TupleType>::value>::template eval<
      TupleType, F>(name, f);
}

// ----------------------------------------

template <int I>
struct for_type_short_impl {
  template <typename TupleType, typename F>
  static void eval(std::string short_name, F f) {
    typedef typename std::tuple_element<I - 1, TupleType>::type ElementType;
    if (traits<ElementType>::short_name == short_name) {
      f.template eval<ElementType>();
    }
    for_type_short_impl<I - 1>::template eval<TupleType, F>(short_name, f);
  }
};

template <>
struct for_type_short_impl<1> {
  template <typename TupleType, typename F>
  static void eval(std::string short_name, F f) {
    typedef typename std::tuple_element<0, TupleType>::type ElementType;
    if (traits<ElementType>::short_name == short_name) {
      f.template eval<ElementType>();
    }
  }
};

template <typename TupleType, typename F>
void for_type_short(std::string short_name, F f = F()) {
  return for_type_short_impl<std::tuple_size<TupleType>::value>::template eval<
      TupleType, F>(short_name, f);
}

// ----------------------------------------

template <int I, bool IsLast>
struct for_each_element_impl {
  template <typename TupleType, typename F>
  static void eval(F f) {
    f.template eval<I>();
    for_each_element_impl<std::tuple_size<TupleType>::value == I + 1,
                          I + 1>::template eval<TupleType, F>(f);
  }
};

template <int I>
struct for_each_element_impl<I, true> {
  template <typename TupleType, typename F>
  static void eval(F f) {
    return f.template eval<I>();
  }
};

/// \brief Call f.eval<I>(), for int I in range [0,
/// std::tuple_size<TupleType>::value)
template <typename TupleType, typename F>
void for_each_element(F f = F()) {
  return for_each_element_impl<0, std::tuple_size<TupleType>::value ==
                                      1>::template eval<TupleType, F>(f);
}

// -------------------------------------------

template <typename T, typename Tuple>
struct has_type;

template <typename T>
struct has_type<T, std::tuple<>> : std::false_type {};

template <typename T, typename U, typename... Ts>
struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

template <typename T, typename... Ts>
struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

}  // namespace CASM_TMP

}  // namespace CASM

#endif
