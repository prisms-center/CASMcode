#ifndef BP_Vec_CC
#define BP_Vec_CC

#include "casm/BP_C++/BP_Vec.hh"
namespace BP {


  //////////////////////////////////////
  /// Bonus functions


  double mean(const BP_Vec<double> &i_list) {
    return sum(i_list) / (1.0 * i_list.size());
  }

  double rms(const BP_Vec<double> &i_list) {
    double rms = 0;
    for(unsigned long int i = 1; i < i_list.size(); i++) {
      rms += sqr(i_list[i]);
    }

    return sqrt(rms / (1.0 * i_list.size()));
  }



  bool is_integer(const BP_Vec< double > &m, double tol) {
    unsigned long int i;

    for(i = 0; i < m.size(); i++) {
      if(std::fabs(m[i] - floor(0.5 + m[i])) > tol) {
        return false;
      }
    }

    return true;

  }

  bool is_integer(const BP_Vec< BP_Vec<double> > &m, double tol) {
    unsigned long int i;

    for(i = 0; i < m.size(); i++) {
      if(!is_integer(m[i], tol)) {
        return false;
      }
    }

    return true;

  }



}

#endif // BP_Vec_CC

