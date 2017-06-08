#include "casm/clex/Calculable.hh"

#include <algorithm>
#include <boost/filesystem.hpp>
#include "casm/clex/PrimClex.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/database/DatabaseTypeDefs.hh"

namespace CASM {

  template<typename Derived>
  const jsonParser &Calculable<Derived>::calc_properties() const {
    return m_calc_properties;
  }

  template<typename Derived>
  void Calculable<Derived>::set_calc_properties(const jsonParser &json) {
    m_calc_properties = json;
  }

  template<typename Derived>
  const jsonParser &Calculable<Derived>::source() const {
    return m_source;
  }

  template<typename Derived>
  void Calculable<Derived>::set_source(const jsonParser &source) {
    if(source.is_null() || source.size() == 0) {
      m_source.put_array();
    }
    else if(!source.is_array()) {
      m_source.put_array();
      m_source.push_back(source);
    }
    else {
      m_source = source;
    }
  }

  template<typename Derived>
  void Calculable<Derived>::push_back_source(const jsonParser &source) {

    if(source.is_null() || source.size() == 0) {
      return;
    }
    if(!source.is_array()) {

      // check if the new source is already listed, if it is do nothing
      for(int i = 0; i < m_source.size(); i++) {
        if(m_source[i] == source)
          return;
      }

      // else, add the new source
      m_source.push_back(source);
    }
    else {

      // check all new sources, if already listed skip, if the any of the new sources is already listed, if it is do nothing

      for(int s = 0; s < source.size(); s++) {

        bool found = false;

        for(int i = 0; i < m_source.size(); i++) {
          if(m_source[i] == source[s]) {
            found = true;
            break;
          }
        }

        if(!found) {
          // else, add the new source
          m_source.push_back(source[s]);
        }
      }
    }
  }

  /// Call in Derived any time DoF may be modified
  template<typename Derived>
  void Calculable<Derived>::_modify_dof() {
    this->clear_name();
    m_calc_properties.put_null();

    if(cache_updated()) {
      cache_clear();
    }

    m_source.put_null();
  }

  /// \brief Return true if all required properties have been been calculated for
  /// the configuration
  template<typename ConfigType>
  bool is_calculated(const ConfigType &config) {
    const auto &props = config.primclex().settings().template properties<ConfigType>();
    return is_calculated(config.calc_properties(), props);
  }

  template<typename ConfigType>
  void reset_properties(ConfigType &config) {
    config.set_calc_properties(jsonParser());
  }

  /// \brief Status of calculation
  template<typename ConfigType>
  std::string calc_status(const ConfigType &config) {
    fs::path p = calc_status_path(config);
    if(fs::exists(p)) {
      jsonParser json(p);
      if(json.contains("status"))
        return json["status"].get<std::string>();
    }
    return("not_submitted");
  }

  // \brief Reason for calculation failure.
  template<typename ConfigType>
  std::string failure_type(const ConfigType &config) {
    fs::path p = calc_status_path(config);
    if(fs::exists(p)) {
      jsonParser json(p);
      if(json.contains("failure_type"))
        return json["failure_type"].get<std::string>();
    }
    return("none");
  }

  template<typename ConfigType>
  bool has_calc_status(const ConfigType &config) {
    return !calc_status(config).empty();
  }

  template<typename ConfigType>
  bool has_failure_type(const ConfigType &config) {
    return !failure_type(config).empty();
  }

  template<typename ConfigType>
  fs::path calc_properties_path(const ConfigType &config) {
    return calc_properties_path(config.primclex(), config.name());
  }

  template<typename ConfigType>
  fs::path pos_path(const ConfigType &config) {
    return pos_path(config.primclex(), config.name());
  }

  template<typename ConfigType>
  fs::path calc_status_path(const ConfigType &config) {
    return calc_status_path(config.primclex(), config.name());
  }

  /// \brief Read properties.calc.json from training_data
  ///
  /// \returns tuple of:
  /// - 0: JSON with calculated properties (or empty object)
  /// - 1: bool indicating there is any data
  /// - 2: bool indicating complete data
  ///
  /// JSON includes:
  /// - file contents verbatim
  /// -  "data_timestamp" with the last write time of the file
  ///
  template<typename ConfigType>
  std::tuple<jsonParser, bool, bool> read_calc_properties(const ConfigType &config) {
    return read_calc_properties<ConfigType>(
             config.primclex(),
             calc_properties_path(config.primclex(), config.name()));
  }

  /// \brief Read properties.calc.json from file
  ///
  /// \returns tuple of:
  /// - 0: JSON with calculated properties (or empty object)
  /// - 1: bool indicating there is any data
  /// - 2: bool indicating complete data
  ///
  /// JSON includes:
  /// - file contents verbatim
  /// - "data_timestamp" with the last write time of the file
  ///
  template<typename ConfigType>
  std::tuple<jsonParser, bool, bool> read_calc_properties(const PrimClex &primclex, const fs::path &filepath) {
    if(!fs::exists(filepath)) {
      return std::make_tuple(jsonParser(), false, false);
    }
    jsonParser props(filepath);
    if(!props.is_obj()) {
      primclex.err_log() << "error parsing: " << filepath << std::endl;
      primclex.err_log() << "not a valid properties.calc.json for Configuration: not a JSON object" << std::endl;
      return std::make_tuple(jsonParser(), false, false);
    }
    props["data_timestamp"] = fs::last_write_time(filepath);

    const auto &prop_vec = primclex.settings().properties<ConfigType>();
    bool is_calc = is_calculated(jsonParser(filepath), prop_vec);
    return std::make_tuple(props, true, is_calc);
  }

  /// \brief Return true if all required properties are included in the JSON
  bool is_calculated(
    const jsonParser &calc_properties,
    const std::vector<std::string> &required_properties) {

    return std::all_of(required_properties.begin(),
                       required_properties.end(),
    [&](const std::string & key) {
      return calc_properties.contains(key);
    });
  }

  fs::path calc_properties_path(const PrimClex &primclex, const std::string &configname) {
    return primclex.dir().calculated_properties(configname, primclex.settings().default_clex().calctype);
  }

  fs::path pos_path(const PrimClex &primclex, const std::string &configname) {
    return primclex.dir().POS(configname);
  }

  fs::path calc_status_path(const PrimClex &primclex, const std::string &configname) {
    return primclex.dir().calc_status(configname, primclex.settings().default_clex().calctype);
  }


#define INST_ConfigType(r, data, type) \
template bool is_calculated(const type &config); \
template void reset_properties(type &config); \
template std::string calc_status(const type &_config); \
template std::string failure_type(const type &config); \
template bool has_calc_status(const type &config); \
template bool has_failure_type(const type &config); \
template fs::path calc_properties_path(const type &config); \
template fs::path pos_path(const type &config); \
template fs::path calc_status_path(const type &config); \
template std::tuple<jsonParser, bool, bool> read_calc_properties<type>(const type &config); \
template std::tuple<jsonParser, bool, bool> read_calc_properties<type>(const PrimClex &primclex, const fs::path &filepath); \
template class Calculable<type>;

  BOOST_PP_SEQ_FOR_EACH(INST_ConfigType, _, CASM_DB_CONFIG_TYPES)

}
