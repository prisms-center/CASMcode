#include "casm/symmetry/SymPermutation.hh"

#include "casm/container/Permutation.hh"
#include "casm/container/io/PermutationIO.hh"

namespace CASM {

  double SymPermutation::character() const {
    int n_fix(0);
    for(Index i = 0; i < m_permute.size(); i++) {
      n_fix += int(i == m_permute[i]);
    }
    return n_fix;
  }

  //*******************************************************************************************

  //Makes permutation matrix from m_permute
  //Permute matrix P reorders elements of a column vector v via
  // v' = P*v
  //Therefore, we loop over m_permute and assign ones to the
  //appropriate row/column

  void SymPermutation::_calc_mat() const {
    m_mat.resize(m_permute.size(), m_permute.size());
    m_mat.setZero();
    for(Index i = 0; i < m_permute.size(); i++)
      m_mat(i, m_permute[i]) = 1;
    return;
  }
}
