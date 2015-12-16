#include <functional>
#include <boost/algorithm/string.hpp>
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  namespace ConfigIO {
    bool RelaxationStrain::parse_args(const std::string &args) {

      std::string tmetric_name = "GL";
      std::string index_expr;
      auto it = args.cbegin(), it_end = args.cend();
      while(it != it_end) {
        if(std::isdigit(*it) || (*it) == ':') {
          auto it2 = it;
          while(it2 != it_end && (std::isdigit(*it2) || isspace(*it2) || (*it2) == ':'))
            ++it2;
          index_expr = std::string(it, it2);
          it = it2;
        }
        if(it != it_end && std::isalpha(*it)) {
          auto it2 = it;
          while(it2 != it_end && std::isalpha(*it2))
            ++it2;
          tmetric_name = std::string(it, it2);
          it = it2;
        }
        while(it != it_end && !std::isalnum(*it) && *it != ':')
          ++it;
      }
      if(m_metric_name.size() > 0 && tmetric_name != m_metric_name) {
        return false;
      }
      m_metric_name = tmetric_name;
      if(index_expr.size() > 0) {
        _parse_index_expression(index_expr);
      }

      return true;
    }

    //****************************************************************************************

    void RelaxationStrain::init(const Configuration &_tmplt) const {
      if(m_metric_name.size() == 0)
        m_metric_name = "GL";
      m_straincalc.set_mode(m_metric_name);
      if(_index_rules().size() > 0)
        return;

      for(Index i = 0; i < 6; i++)
        _add_rule(std::vector<Index>({i}));

    }
    //****************************************************************************************

    bool RelaxationStrain::validate(const Configuration &_config) const {
      return _config.calc_properties().contains("relaxation_deformation");
    }

    //****************************************************************************************

    std::string RelaxationStrain::long_header(const Configuration &_tmplt) const {
      std::stringstream t_ss;
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      Index s = max(8 - int(name().size()), 0);
      for(; it != end_it; ++it)
        t_ss << "    " << name() << '(' << m_metric_name << ',' << (*it)[0] << ')';
      return t_ss.str();
    }


    //****************************************************************************************

    std::string RelaxationStrain::short_header(const Configuration &_tmplt) const {
      return name() + "(" + m_metric_name + ")";
    }

    //****************************************************************************************
    Eigen::VectorXd RelaxationStrain::evaluate(const Configuration &_config) const {
      return m_straincalc.unrolled_strain_metric(_config.calc_properties()["relaxation_deformation"].get<Eigen::Matrix3d>());
    }
/*
    // ****************************************************************************************
    void RelaxationStrain::inject(const Configuration &_config, DataStream &_stream, Index) const {
      if(!validate(_config)) {
        _stream << DataStream::failbit << std::vector<double>(_index_rules().size(), NAN);
      }
      else {
        Eigen::VectorXd strain_vec = _evaluate(_config);
        for(auto it = _index_rules().cbegin(); it != _index_rules().cend(); ++it)
          _stream << strain_vec[(*it)[0]];
      }
    }

    // ****************************************************************************************

    void RelaxationStrain::print(const Configuration &_config, std::ostream &_stream, Index) const {
      _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
      _stream.precision(8);

      if(!validate(_config)) {
        for(auto it = _index_rules().cbegin(); it != _index_rules().cend(); ++it)
          _stream << "     unknown";
      }

      else {
        Eigen::VectorXd strain_vec = _evaluate(_config);

        for(auto it = _index_rules().cbegin(); it != _index_rules().cend(); ++it)
          _stream << "     " << strain_vec[(*it)[0]];
      }

    }

    // ****************************************************************************************

    jsonParser &RelaxationStrain::to_json(const Configuration &_config, jsonParser &json)const {
      if(validate(_config))
        json = _evaluate(_config);

      return json;
    }
*/
  }

}
