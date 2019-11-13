#ifndef CASM_GLOBAL_DEFINTIONS_HH
#define CASM_GLOBAL_DEFINTIONS_HH

#include <utility>

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
    return apply(f, obj, std::forward<Args>(args)...);
  }

}

#endif
