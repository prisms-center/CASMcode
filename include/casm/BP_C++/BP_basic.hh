#ifndef BP_basic_HH
#define BP_basic_HH

#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
///  Useful functions

namespace BP {

  extern const double PI;
  //double PI = atan(1)*4;				///< CONSTANT, defined in BP_basic.cc

  void BP_pause();

  template<class T> inline T sqr(T a) {
    return a * a;
  };

  bool is_even(int i);

  int int_pow(int x, int p);

  unsigned long int ulint_pow(unsigned long int x, unsigned long int p);
}

#endif // BP_basic_HH
