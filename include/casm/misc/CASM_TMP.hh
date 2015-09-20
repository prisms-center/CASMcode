#include <type_traits>
#include <cassert>

#include "casm/CASM_global_definitions.hh"

#ifndef CASM_TMP_HH
#define CASM_TMP_HH

namespace CASM {
  namespace CASM_TMP {
    template <typename T>
    struct traits {};

    // Definitions for IfIntegralTol
    template <typename tol_type, bool IsIntegral>
    struct IfIntegralTol;

    template <typename tol_type>
    struct IfIntegralTol<tol_type, true> {
      IfIntegralTol() {};
      IfIntegralTol(tol_type) {};
      tol_type tol() const {
        return 0;
      }
    };

    template <typename tol_type>
    struct IfIntegralTol<tol_type, false> {
      IfIntegralTol(tol_type _tol) : m_tol(_tol) {};
      tol_type tol() {
        return m_tol;
      }
    private:
      tol_type m_tol;
    };

    template<typename T>
    using TypedTol = IfIntegralTol<T, std::is_integral<T>::value >;
    // End of IfIntegralTol

    template<bool IsConst, typename T>
    using ConstSwitch = typename std::conditional<IsConst, const T, T>::type;


    // Definitions for MuchLessThan
    template<typename value_type>
    struct IntegralLessThan {
      IntegralLessThan() {};
      IntegralLessThan(value_type) {};
      bool operator()(const value_type &A, const value_type &B) const {
        return A < B;
      };
    };

    template<typename value_type>
    struct FloatingPointLessThan {
      FloatingPointLessThan(value_type _tol = TOL) : m_tol(_tol) {};
      bool operator()(const value_type &A, const value_type &B) const {
        return A + m_tol < B;
      };
    private:
      value_type m_tol;
    };

    template<typename T>
    using MuchLessThan = typename std::conditional<boost::is_integral<T>::value, IntegralLessThan<T>, FloatingPointLessThan<T> >::type;
    // End of MuchLessThan


    /// \brief Helper Functor for Counter container access using operator[]
    template < typename Container,
               typename value_type = typename Container::value_type,
               typename size_type = typename Container::size_type >
    struct BracketAccess {

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
               typename value_type = typename Container::value_type,
               typename size_type = typename Container::size_type >
    struct ParenthesesAccess {

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
