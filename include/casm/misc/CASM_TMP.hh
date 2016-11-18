#ifndef CASM_TMP_HH
#define CASM_TMP_HH

#include <type_traits>
#include <cassert>

namespace CASM {
  namespace CASM_TMP {

    template <typename T>
    struct traits {};

    // ---------------------

    /// \brief Enables calling a function on each argument in a parameter pack
    template<typename...Args>
    void ignore_returnvalues(Args &&...) {}

    // ---------------------

    /// \brief Alias for void, to help SFINAE work
    template<typename... Ts>
    struct make_void {
      typedef void type;
    };

    /// \brief Alias for void, to help SFINAE work
    template<typename... Ts>
    using void_t = typename make_void<Ts...>::type;

    /// \brief Base type inherits from std::false_type if T is not iterator
    ///
    /// - T is considered an iterator if it is incrementable, dereferenceable, and comparable
    template <typename T, typename = void>
    struct is_iterator : std::false_type { };

    /// \brief Specialized case inherits from std::true_type if T is an iterator
    ///
    /// - T is considered an iterator if it is incrementable, dereferenceable, and comparable
    template <typename T>
    struct is_iterator < T,
           void_t < decltype(++std::declval<T &>()),
           decltype(*std::declval<T &>()),
           decltype(std::declval<T &>() == std::declval<T &>()) > >
       : std::true_type { };

    /// \brief Template alias to enable a template function via SFINAE for any iterator
    ///
    /// - example function argument: 'typename CASM_TMP::enable_if_iterator<IteratorType>::type* = nullptr'
    template<typename T>
    using enable_if_iterator = std::enable_if<is_iterator<T>::type::value, void>;

    /// \brief Template alias to enable a template function via SFINAE for an iterator with value_type V
    ///
    /// - example function argument: 'typename CASM_TMP::enable_if_iterator_of<IteratorType, ValueType>::type* = nullptr'
    template<typename T, typename V>
    using enable_if_iterator_of =
      std::enable_if < is_iterator<T>::type::value &&
      std::is_same<typename std::iterator_traits<T>::value_type, V>::type::value, void >;


    /// \brief Helper Functor for Counter container access using operator[]
    template < typename Container,
               typename _value_type = typename Container::value_type,
               typename _size_type = typename Container::size_type >
    struct BracketAccess {

      typedef _value_type value_type;
      typedef _size_type size_type;

      BracketAccess() {}
      BracketAccess &operator=(const BracketAccess &) {
        return *this;
      };

      /// \brief Identical to 'return container[index];'
      value_type &operator()(Container &container, size_type index) const {
        return container[index];
      }

      /// \brief Identical to 'return container[index];'
      const value_type &operator()(const Container &container, size_type index) const {
        return container[index];
      }

      /// \brief Identical to 'return container[i][j];'
      value_type &operator()(Container &container, size_type i, size_type j) const {
        return container[i][j];
      }

      /// \brief Identical to 'return container[i][j];'
      const value_type &operator()(const Container &container, size_type i, size_type j) const {
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
      static const value_type &at(const Container &container, size_type i, size_type j) {
        return container[i][j];
      }
    };

    /// \brief Helper Functor for Counter container access using operator()
    template < typename Container,
               typename _value_type = typename Container::value_type,
               typename _size_type = typename Container::size_type >
    struct ParenthesesAccess {

      typedef _value_type value_type;
      typedef _size_type size_type;

      ParenthesesAccess() {}
      ParenthesesAccess &operator=(const ParenthesesAccess &) {
        return *this;
      };


      /// \brief Identical to 'return container(index);'
      value_type &operator()(Container &container, size_type index) const {
        return container(index);
      }

      /// \brief Identical to 'return container(index);'
      const value_type &operator()(const Container &container, size_type index) const {
        return container(index);
      }

      /// \brief Identical to 'return container(index);'
      value_type &operator()(Container &container, size_type i, size_type j) const {
        return container(i, j);
      }

      /// \brief Identical to 'return container(index);'
      const value_type &operator()(const Container &container, size_type i, size_type j) const {
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
      static const value_type &at(const Container &container, size_type i, size_type j) {
        return container(i, j);
      }

    };


  }

}

#endif
