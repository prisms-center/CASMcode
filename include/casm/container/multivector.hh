#ifndef CASM_multivector
#define CASM_multivector

#include <vector>

namespace CASM {

  namespace multivector_impl {
    template<typename T, size_t N>
    struct multivector_tmp;

    template<typename T>
    struct multivector_tmp<T, 0> {
      using type = T;
    };
    template<typename T, size_t N>
    struct multivector_tmp {
      using type = std::vector < typename multivector_tmp < T, N - 1 >::type >;
    };
  }

  /// \brief Shortcut for multidimensional vector (std::vector< std::vector< ...)
  ///
  /// - multivector<Type>::X<D> is a D-dimensional container of std::vector<T>
  template<typename T>
  struct multivector {
    template<size_t N>
    using X = typename multivector_impl::multivector_tmp<T, N>::type;
  };

}

#endif
