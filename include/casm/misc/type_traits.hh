#ifndef CASM_type_traits
#define CASM_type_traits

namespace notstd {
  /// \brief Alias for void, to help SFINAE work
  template<typename... Ts>
  struct make_void {
    typedef void type;
  };

  /// \brief Alias for void, to help SFINAE work
  template<typename... Ts>
  using void_t = typename make_void<Ts...>::type;
}

#endif
