#include "casm/clex/ConfigDoFValues.hh"

namespace CASM {

  void LocalContinuousConfigDoFValues::from_standard_values(Eigen::Ref<const Eigen::MatrixXd> const &_values) {
    resize_vol(_values.cols() / n_sublat());
    for(Index b = 0; b < n_sublat(); ++b)
      sublat(b).topRows(info()[b].dim()) = info()[b].inv_basis() * _values.block(0, b * n_vol(), dim(), n_vol());
  }

  void GlobalContinuousConfigDoFValues::from_standard_values(Eigen::Ref<const Eigen::MatrixXd> const &_values) {
    values() = info().inv_basis() * _values;

  }

}
