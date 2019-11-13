#include "casm/clex/ConfigDoFValues.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/container/json_io.hh"

namespace CASM {

  void LocalContinuousConfigDoFValues::from_standard_values(Eigen::Ref<const Eigen::MatrixXd> const &_values) {
    resize_vol(_values.cols() / n_sublat());
    for(Index b = 0; b < n_sublat(); ++b)
      sublat(b).topRows(info()[b].dim()) = info()[b].inv_basis() * _values.block(0, b * n_vol(), dim(), n_vol());
  }

  jsonParser &to_json(LocalContinuousConfigDoFValues const &_values, jsonParser &_json) {
    to_json(_values.standard_values().transpose(), _json["values"]);
    return _json;
  }

  void from_json(LocalContinuousConfigDoFValues &_values, jsonParser const &_json) {
    Eigen::MatrixXd tval = _json["values"].get<Eigen::MatrixXd>().transpose();
    _values.from_standard_values(tval);
  }

  jsonParser &to_json(LocalDiscreteConfigDoFValues const &_values, jsonParser &_json) {
    to_json_array(_values.values(), _json);
    return _json;
  }

  void from_json(LocalDiscreteConfigDoFValues &_values, jsonParser const &_json) {
    _values.resize_vol(_json.size() / _values.n_sublat());
    _values.values() = _json.get<LocalDiscreteConfigDoFValues::ValueType>();
  }

  jsonParser &to_json(GlobalContinuousConfigDoFValues const &_values, jsonParser &_json) {
    to_json_array(_values.standard_values(), _json["values"]);
    return _json;
  }

  void GlobalContinuousConfigDoFValues::from_standard_values(Eigen::Ref<const Eigen::MatrixXd> const &_values) {
    values() = info().inv_basis() * _values;

  }

  void from_json(GlobalContinuousConfigDoFValues &_values, jsonParser const &_json) {
    _values.from_standard_values(_json["values"].get<Eigen::VectorXd>());
  }



}
