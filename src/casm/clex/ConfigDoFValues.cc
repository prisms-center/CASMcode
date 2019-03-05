#include "casm/clex/ConfigDoFValues.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  jsonParser &to_json(LocalContinuousConfigDoFValues const &_values, jsonParser &_json) {
    Eigen::MatrixXd tval;

    //Size tval to the size of the native DoF vector space
    tval.setZero(DoF::traits(_values.type_name()).dim(), _values.values().cols());
    for(Index b = 0; b < _values.n_sublat(); ++b)
      tval.block(0, b * _values.n_vol(), tval.rows(), _values.n_vol()) = _values.info()[b].basis() * _values.sublat(b).topRows(_values.info()[b].dim());

    to_json(tval.transpose(), _json["values"]);
    return _json;
  }

  void from_json(LocalContinuousConfigDoFValues &_values, jsonParser const &_json) {
    Eigen::MatrixXd tval = _json["values"].get<Eigen::MatrixXd>().transpose();
    _values.resize_vol(tval.cols() / _values.n_sublat());
    for(Index b = 0; b < _values.n_sublat(); ++b)
      _values.sublat(b).topRows(_values.info()[b].dim()) = _values.info()[b].inv_basis() * tval.block(0, b * _values.n_vol(), _values.dim(), _values.n_vol());
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
    Eigen::VectorXd tval;
    //Size tval to the size of the native DoF vector space
    tval.setZero(DoF::traits(_values.type_name()).dim());
    tval = _values.info().basis() * _values.values();

    to_json_array(tval, _json["values"]);
    return _json;
  }

  void from_json(GlobalContinuousConfigDoFValues &_values, jsonParser const &_json) {
    _values.values() = _values.info().inv_basis() * _json["values"].get<Eigen::VectorXd>();
  }



}
