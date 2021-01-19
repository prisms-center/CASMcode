#include "casm/misc/CASM_Array_math.hh"

namespace CASM {

// returns 'i' if 'input' is equivalent to 'unique[i]', w.r.t. permutation of
// the equivalent elements of 'input'. equivalent elements are specified by
// 'ind_equiv' if 'input' specifies a new combination of integers, unique.size()
// is returned
Index which_unique_combination(const Array<Index> &input,
                               const Array<Array<Index> > &unique,
                               const Array<Array<Index> > &ind_equiv) {
  Index tval, tcount;
  Index i, j, k, l;
  for (i = 0; i < unique.size(); i++) {
    // Is input equivalent to unique[i] w.r.t. ind_equiv?
    // Loop over groups of equivalent indices
    if (unique[i].size() != input.size()) continue;
    for (j = 0; j < ind_equiv.size(); j++) {
      // Loop over indices that are equivalent
      for (k = 0; k < ind_equiv[j].size(); k++) {
        tval = input[ind_equiv[j][k]];
        tcount = 0;
        for (l = 0; l < ind_equiv[j].size(); l++)
          tcount += int(unique[i][ind_equiv[j][l]] == tval) -
                    int(input[ind_equiv[j][l]] == tval);

        if (tcount != 0) break;
      }
      if (k != ind_equiv[j].size()) break;
    }
    if (j == ind_equiv.size()) return i;
  }
  return i;
}

//****************************************

Index which_unique_combination(const Array<Index> &input,
                               const Array<Array<Index> > &unique) {
  Index tval, tcount;
  Index i, j, l;
  for (i = 0; i < unique.size(); i++) {
    // Is input equivalent to unique[i] w.r.t. ind_equiv?
    // Loop over groups of equivalent indices
    if (unique[i].size() != input.size()) continue;
    for (j = 0; j < input.size(); j++) {
      tval = input[j];
      tcount = 0;
      for (l = 0; l < input.size(); l++)
        tcount += int(unique[i][l] == tval) - int(input[l] == tval);

      if (tcount != 0) break;
    }

    if (j == input.size()) return i;
  }
  return i;
}

//************************************************************
/// Find least common multiple
int lcm(const Array<int> &series) {
  if (!series.size()) return 0;
  int lcm_val(series[0]);
  for (Index i = 1; i < series.size(); i++) lcm_val = lcm(lcm_val, series[i]);

  return lcm_val;
}

//************************************************************

ReturnArray<Array<int> > get_prime_factors(int target) {
  Array<Array<int> > factors_array;
  Array<int> factor_list;

  if (target <= 1) {
    std::cerr << "WARNING in global/definitions::get_prime_factors"
              << std::endl;
    std::cerr << "You're asking for prime factors of " << target
              << ". Returning empty array." << std::endl
              << std::endl;
    return factors_array;
  }
  int factor = 2;

  while (target != 1) {
    while (target % factor == 0) {
      factor_list.push_back(factor);
      target = target / factor;
    }

    factor++;

    if (factor_list.size() > 0) {
      factors_array.push_back(factor_list);
    }
    factor_list.clear();
  }

  return factors_array;
}

}  // namespace CASM
