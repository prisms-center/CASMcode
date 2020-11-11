#ifndef CASM_GLOBAL_DEFINTIONS_HH
#define CASM_GLOBAL_DEFINTIONS_HH

#include <utility>
#include <casm/misc/type_traits.hh>

namespace boost {
  namespace filesystem {
    class path;
  }
  namespace program_options {
  }
}

#ifdef EIGEN_DEFAULT_DENSE_INDEX_TYPE
#define INDEX_TYPE EIGEN_DEFAULT_DENSE_INDEX_TYPE
#else
#define INDEX_TYPE long
#endif

namespace CASM {

  namespace fs = boost::filesystem;
  namespace po = boost::program_options;

  typedef unsigned int uint;
  typedef unsigned long int ulint;
  typedef long int lint;

  //tolerance
  const double TOL = 0.00001;

  //Boltzmann Constant
  const double KB = 8.6173423E-05; //eV/K

  //Planck's Constant
  const double PLANCK = 4.135667516E-15; //eV-s

  ///For long integer indexing:
  typedef INDEX_TYPE Index;
  bool valid_index(Index i);

  template<typename T> struct traits;

  //************************************************************

  /// For any binary functor class that accepts two arguments, such as
  /// an equals or compare functor, create a unary functor that
  /// stores a reference to the LHS value, and a reference to the
  /// binary functor, which then only needs to be provided the RHS value.
  /// This is useful for cases where one wants to use a binary functor
  /// in the context of something like std::find_if
  ///
  /// Example:
  /// UnaryCompare_f<BinaryCompare_f> equals_target_element(target_element,args,for,binary,comapre);
  /// std::find_if(v.begin(),v.end(),equals_target_element);
  template <typename BinaryCompare>
  class UnaryCompare_f {
  public:
    using argument_type = notstd::first_argument_type<BinaryCompare>;

    template <typename... CompareArgs>
    UnaryCompare_f(const argument_type &lhs, const CompareArgs &... args) : m_lhs(lhs), m_compare_method(args...) {
    }

    bool operator()(const argument_type &rhs) {
      return m_compare_method(m_lhs, rhs);
    }

  private:
    //TODO: Is having this as a reference too scary? We could just make it a copy;
    const argument_type &m_lhs;
    const BinaryCompare m_compare_method;
  };


  //************************************************************

  /// Apply a transformation, in place
  /// - Default is equivalent to \code obj.apply_sym(f, args...) \endcode
  template<typename Object, typename Transform, typename...Args>
  Object &apply(const Transform &f, Object &obj, Args &&...args) {
    obj.apply_sym(f, std::forward<Args>(args)...);
    return obj;
  }

  /// Copy and apply a transformation
  /// - We can also include a specialization for cloneable objects, via SFINAE
  /// - Default is equivalent to \code Object result(obj); apply(f, obj, args...) \endcode
  template<typename Object, typename Transform, typename...Args>
  Object copy_apply(const Transform &f, Object obj, Args &&...args) {
    return CASM::apply(f, obj, std::forward<Args>(args)...);
  }

}

#endif
