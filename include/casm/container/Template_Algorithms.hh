#ifndef TEMPLATE_ALGORITHMS_HH
#define TEMPLATE_ALGORITHMS_HH

#include <cmath>

namespace CASM {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

template <class Container, class T>
class GramSchmidt {
 public:
  static void orthogonalize(Array<Container> &rows) {
    Index i, j;
    T tcoeff;
    for (i = 0; i < rows.size(); i++) {
      for (j = 0; j < i; j++) rows[i] -= dot(rows[i], rows[j]) * rows[j];

      tcoeff = norm(rows[i]);
      if (TOL < tcoeff)
        rows[i] /= tcoeff;

      else {
        rows.remove(i);
        i--;
      }
    }
    return;
  }
};

};  // namespace CASM
#endif
