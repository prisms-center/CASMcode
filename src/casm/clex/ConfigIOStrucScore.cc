#include <functional>
#include "casm/external/boost.hh"
#include "casm/casm_io/EigenDataStream.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/jsonStruc.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  namespace ConfigIO {
    bool StrucScore::parse_args(const std::string &args) {
      std::vector<std::string> splt_vec;
      double _lattice_weight(0.5);
      bool already_initialized = !m_prim_path.empty();
      int pushed_args = 0;
      boost::split(splt_vec, args, boost::is_any_of(", "), boost::token_compress_on);
      if(splt_vec.size() < 2 || splt_vec.size() > 4) {
        throw std::runtime_error("Attempted to initialize format tag " + name()
                                 + " with " + std::to_string(splt_vec.size()) + " arguments ("
                                 + args + "). You must provide at least 2 arguments, but no more than 4.\n");
        return false;
      }
      if(!m_prim_path.empty() && fs::path(splt_vec[0]) != m_prim_path)
        return false;
      if(m_prim_path.empty()) {
        m_prim_path = splt_vec[0];
        if(!fs::exists(m_prim_path)) {
          throw std::runtime_error("Attempted to initialize format tag " + name()
                                   + " invalid file path '" + fs::absolute(m_prim_path).string() + "'. File does not exist.\n");
        }

        m_altprimclex.reset(new PrimClex(Structure(m_prim_path)));
      }
      for(Index i = 1; i < splt_vec.size(); ++i) {
        if(splt_vec[i] != "basis_score" && splt_vec[i] != "lattice_score" && splt_vec[i] != "total_score") {
          try {
            _lattice_weight = std::stod(splt_vec[i]);
          }
          catch(...) {
            throw std::runtime_error("Attempted to initialize format tag " + name()
                                     + " with invalid argument '" + splt_vec[i] + "'. Valid arguments are [ basis_score | lattice_score | total_score ]\n");
          }
          if(already_initialized && !almost_equal(_lattice_weight, m_configmapper.lattice_weight())) {
            for(; pushed_args > 0; pushed_args--)
              m_prop_names.pop_back();
            return false;
          }
        }
        else {
          m_prop_names.push_back(splt_vec[i]);
          ++pushed_args;
        }
      }
      m_configmapper = ConfigMapper(*m_altprimclex, _lattice_weight);
      return true;
    }

    //****************************************************************************************

    bool StrucScore::validate(const Configuration &_config) const {
      return fs::exists(_config.calc_properties_path());
    }

    //****************************************************************************************

    std::vector<std::string> StrucScore::col_header(const Configuration &_tmplt) const {
      std::vector<std::string> col;
      for(Index i = 0; i < m_prop_names.size(); i++) {
        std::stringstream t_ss;
        t_ss << "    " << name() << '(' << m_prim_path.string() << ',' << m_prop_names[i] << ',' << m_configmapper.lattice_weight() << ')';
        col.push_back(t_ss.str());
      }
      return col;
    }


    //****************************************************************************************

    std::string StrucScore::short_header(const Configuration &_tmplt) const {
      std::stringstream t_ss;
      t_ss << name() << '(' << m_prim_path.string();
      for(Index i = 0; i < m_prop_names.size(); i++)
        t_ss   << ',' << m_prop_names[i];
      t_ss << ',' << m_configmapper.lattice_weight() << ')';
      return t_ss.str();
    }

    //****************************************************************************************
    Eigen::VectorXd StrucScore::evaluate(const Configuration &_config)const {
      std::vector<double> result_vec;

      BasicStructure<Site> relaxed_struc;
      ConfigDoF mapped_configdof;
      Lattice mapped_lat;

      auto lambda = [&](const std::vector<double> &result_vec) {
        Eigen::VectorXd res = Eigen::VectorXd::Zero(result_vec.size());
        for(Index i = 0; i < res.size(); ++i) {
          res(i) = result_vec[i];
        }
        return res;
      };

      from_json(simple_json(relaxed_struc, "relaxed_"), jsonParser(_config.calc_properties_path()));

      if(!m_configmapper.struc_to_configdof(relaxed_struc, mapped_configdof, mapped_lat)) {
        for(Index i = 0; i < m_prop_names.size(); i++) {
          result_vec.push_back(1e9);
        }

        return lambda(result_vec);
      }
      for(Index i = 0; i < m_prop_names.size(); i++) {
        if(m_prop_names[i] == "basis_score")
          result_vec.push_back(ConfigMapping::basis_cost(mapped_configdof, relaxed_struc.basis.size()));
        else if(m_prop_names[i] == "lattice_score")
          result_vec.push_back(ConfigMapping::strain_cost(relaxed_struc.lattice(), mapped_configdof, relaxed_struc.basis.size()));
        else if(m_prop_names[i] == "total_score") {
          double sc = ConfigMapping::strain_cost(relaxed_struc.lattice(), mapped_configdof, relaxed_struc.basis.size());

          double bc = ConfigMapping::basis_cost(mapped_configdof, relaxed_struc.basis.size());

          double w = m_configmapper.lattice_weight();
          result_vec.push_back(w * sc + (1.0 - w)*bc);
        }
      }

      return lambda(result_vec);
    }
    /*
        // ****************************************************************************************
        void StrucScore::inject(const Configuration &_config, DataStream &_stream, Index) const {
          if(!validate(_config)) {
            _stream << DataStream::failbit << std::vector<double>(m_prop_names.size(), NAN);
          }
          else {
            _stream << _evaluate(_config);
          }
        }

        // ****************************************************************************************

        void StrucScore::print(const Configuration &_config, std::ostream &_stream, Index) const {
          _stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
          _stream.precision(8);

          if(!validate(_config)) {
            for(auto it = m_prop_names.cbegin(); it != m_prop_names.cend(); ++it)
              _stream << "     unknown";
          }

          else {
            std::vector<double> result_vec(_evaluate(_config));
            for(auto it = result_vec.cbegin(); it != result_vec.cend(); ++it)
              _stream << "     " << *it;
          }

        }

        // ****************************************************************************************

        jsonParser &StrucScore::to_json(const Configuration &_config, jsonParser &json)const {
          if(validate(_config))
            json = _evaluate(_config);

          return json;
        }
    */
  }
}

