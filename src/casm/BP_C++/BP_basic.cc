#ifndef BP_basic_CC
#define BP_basic_CC

#include "casm/BP_C++/BP_basic.hh"

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
///  Useful functions

namespace BP {

  const double PI = atan(1) * 4;				///< CONSTANT

  void BP_pause() {
    int tmp;
    std::cout << "pause(): " << std::endl;
    std::cin >> tmp;
    std::cout << std::endl;
  };

  bool is_even(int i) {
    if(i % 2 == 0) return true;
    else return false;
  };

  int int_pow(int x, int p) {
    if(p == 0) return 1;
    if(p == 1) return x;
    int result = x;
    for(int i = 1; i < p; i++)
      result *= x;
    return result;
  }

  unsigned long int ulint_pow(unsigned long int x, unsigned long int p) {

    if(p == 0) return 1;
    if(p == 1) return x;
    unsigned long int result = x;
    for(unsigned long int i = 1; i < p; i++)
      result *= x;
    return result;
  }

}

#endif // BP_basic_CC
