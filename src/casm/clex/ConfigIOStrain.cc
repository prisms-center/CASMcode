#include <functional>
#include "casm/external/boost.hh"
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

    std::vector<std::string> RelaxationStrain::col_header(const Configuration &_tmplt) const {
      std::vector<std::string> col;
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      // Index s = max(8 - int(name().size()), 0);
      for(; it != end_it; ++it) {
        std::stringstream t_ss;
        t_ss << "    " << name() << '(' << m_metric_name << ',' << (*it)[0] << ')';
        col.push_back(t_ss.str());
      }
      return col;
    }


    //****************************************************************************************

    std::string RelaxationStrain::short_header(const Configuration &_tmplt) const {
      return name() + "(" + m_metric_name + ")";
    }

    //****************************************************************************************
    Eigen::VectorXd RelaxationStrain::evaluate(const Configuration &_config) const {
      return m_straincalc.unrolled_strain_metric(_config.calc_properties()["relaxation_deformation"].get<Eigen::Matrix3d>());
    }

    //****************************************************************************************

    bool DoFStrain::parse_args(const std::string &args) {

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

    void DoFStrain::init(const Configuration &_tmplt) const {
      if(m_metric_name.size() == 0)
        m_metric_name = "GL";
      m_straincalc.set_mode(m_metric_name);
      if(_index_rules().size() > 0)
        return;

      for(Index i = 0; i < 6; i++)
        _add_rule(std::vector<Index>({i}));

    }

    //****************************************************************************************

    std::vector<std::string> DoFStrain::col_header(const Configuration &_tmplt) const {
      std::vector<std::string> col;
      auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
      // Index s = max(8 - int(name().size()), 0);
      for(; it != end_it; ++it) {
        std::stringstream t_ss;
        t_ss << name() << '(' << m_metric_name << ',' << (*it)[0] << ')';
        col.push_back(t_ss.str());
      }
      return col;
    }


    //****************************************************************************************

    std::string DoFStrain::short_header(const Configuration &_tmplt) const {
      return name() + "(" + m_metric_name + ")";
    }

    //****************************************************************************************
    Eigen::VectorXd DoFStrain::evaluate(const Configuration &_config) const {
      return m_straincalc.unrolled_strain_metric(_config.deformation());
    }

  }

}
