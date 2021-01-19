#ifndef CASM_type_traits
#define CASM_type_traits

#include <type_traits>
#include <utility>

namespace notstd {
/// \brief Alias for void, to help SFINAE work
template <typename... Ts>
struct make_void {
  typedef void type;
};

/// \brief Alias for void, to help SFINAE work
template <typename... Ts>
using void_t = typename make_void<Ts...>::type;

/// \brief Deduces return type and argument types for non-generic
/// labmdas and functors.
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType (ClassType::*)(Args...) const> {
  enum { arity = sizeof...(Args) };

  typedef ReturnType result_type;

  template <std::size_t i>
  struct arg {
    typedef typename std::tuple_element<i, std::tuple<Args...>>::type type;
  };
};

template <typename F>
/* using first_argument_type = std::decay_t<typename
   function_traits<F>::template arg<0>::type>; */
using first_argument_type = typename std::decay<
    typename function_traits<F>::template arg<0>::type>::type;
template <typename F>
/* using second_argument_type = std::decay_t<typename
   function_traits<F>::template arg<1>::type>; */
using second_argument_type = typename std::decay<
    typename function_traits<F>::template arg<1>::type>::type;
template <typename F>
/* using return_type = std::decay_t<typename function_traits<F>::result_type>;
 */
using return_type =
    typename std::decay<typename function_traits<F>::result_type>::type;

}  // namespace notstd

#endif
