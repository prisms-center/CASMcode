#include "casm/clex/Calculable.hh"

#include <algorithm>
#include <boost/filesystem.hpp>

#include "casm/app/ProjectSettings.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/PropertiesDatabase.hh"

namespace CASM {

/// \brief Return MappedProperties for requested calctype
///
/// - If nothing calculated, returns empty MappedProperties object
template <typename _Base>
MappedProperties const &Calculable<_Base>::calc_properties(
    std::string calctype) const {
  if (calctype.empty() && derived().has_primclex()) {
    calctype = derived().primclex().settings().default_clex().calctype;
  }
  auto it = m_calc_properties_map.find(calctype);
  if (it == m_calc_properties_map.end()) {
    _refresh_calc_properties(calctype);
    it = m_calc_properties_map.find(calctype);
  }

  return it->second;
}

template <typename _Base>
void Calculable<_Base>::set_calc_properties(const MappedProperties &_prop,
                                            std::string calctype) {
  if (calctype.empty()) {
    if (!derived().has_primclex()) {
      m_calc_properties_map.clear();
      return;
    }
    calctype = derived().primclex().settings().default_clex().calctype;
  }
  m_calc_properties_map[calctype] = _prop;
}

template <typename _Base>
void Calculable<_Base>::refresh_calc_properties(std::string calctype) {
  static_cast<const Calculable<_Base> *>(this)->_refresh_calc_properties(
      calctype);
}

template <typename _Base>
const jsonParser &Calculable<_Base>::source() const {
  return m_source;
}

template <typename _Base>
void Calculable<_Base>::set_source(const jsonParser &source) {
  if (source.is_null() || source.size() == 0) {
    m_source.put_array();
  } else if (!source.is_array()) {
    m_source.put_array();
    m_source.push_back(source);
  } else {
    m_source = source;
  }
}

template <typename _Base>
void Calculable<_Base>::push_back_source(const jsonParser &source) {
  if (source.is_null() || source.size() == 0) {
    return;
  }
  if (!source.is_array()) {
    // check if the new source is already listed, if it is do nothing
    for (int i = 0; i < m_source.size(); i++) {
      if (m_source[i] == source) return;
    }

    // else, add the new source
    m_source.push_back(source);
  } else {
    // check all new sources, if already listed skip, if the any of the new
    // sources is already listed, if it is do nothing

    for (int s = 0; s < source.size(); s++) {
      bool found = false;

      for (int i = 0; i < m_source.size(); i++) {
        if (m_source[i] == source[s]) {
          found = true;
          break;
        }
      }

      if (!found) {
        // else, add the new source
        m_source.push_back(source[s]);
      }
    }
  }
}

/// Call in _Base any time DoF may be modified
template <typename _Base>
void Calculable<_Base>::_modify_dof() {
  this->clear_name();
  m_calc_properties_map.clear();

  if (cache_updated()) {
    cache_clear();
  }

  m_source.put_null();
}

template <typename _Base>
void Calculable<_Base>::_refresh_calc_properties(std::string calctype) const {
  if (!derived().has_primclex()) {
    return;
  }
  const PrimClex &primclex = derived().primclex();
  const auto &db = primclex.const_db_props<MostDerived>(calctype);
  if (calctype.empty()) {
    calctype = primclex.settings().default_clex().calctype;
  }
  auto it = db.find_via_to(this->name());
  if (it != db.end()) {
    m_calc_properties_map[calctype] = *it;
  } else {
    m_calc_properties_map[calctype] = MappedProperties();
  }
}

/// \brief Return true if all required properties have been been calculated for
/// the configuration in the calctype specified (current calctype if none
/// specified)
template <typename ConfigType>
bool is_calculated(const ConfigType &config, std::string calctype) {
  auto const &settings = config.primclex().settings();
  if (calctype == "") {
    calctype = settings.default_clex().calctype;
  }
  const auto &props =
      settings.required_properties(traits<ConfigType>::name, calctype);
  return is_calculated(config.calc_properties(calctype), props);
}
template <typename ConfigType>
void reset_properties(ConfigType &config) {
  config.set_calc_properties(MappedProperties(), "");
}

/// \brief Status of calculation
template <typename ConfigType>
std::string calc_status(const ConfigType &config, std::string calctype) {
  if (calctype == "") {
    calctype = config.primclex().settings().default_clex().calctype;
  }
  fs::path p = calc_status_path(config, calctype);
  if (fs::exists(p)) {
    jsonParser json(p);
    if (json.contains("status")) return json["status"].get<std::string>();
  }
  return ("not_submitted");
}

// \brief Reason for calculation failure.
template <typename ConfigType>
std::string failure_type(const ConfigType &config, std::string calctype) {
  if (calctype == "") {
    calctype = config.primclex().settings().default_clex().calctype;
  }
  fs::path p = calc_status_path(config, calctype);
  if (fs::exists(p)) {
    jsonParser json(p);
    if (json.contains("failure_type"))
      return json["failure_type"].get<std::string>();
  }
  return ("none");
}

template <typename ConfigType>
bool has_calc_status(const ConfigType &config, std::string calctype) {
  if (calctype == "") {
    calctype = config.primclex().settings().default_clex().calctype;
  }
  return !calc_status(config, calctype).empty();
}

template <typename ConfigType>
bool has_failure_type(const ConfigType &config, std::string calctype) {
  if (calctype == "") {
    calctype = config.primclex().settings().default_clex().calctype;
  }
  return !failure_type(config, calctype).empty();
}

template <typename ConfigType>
std::string calc_properties_path(const ConfigType &config,
                                 std::string calctype) {
  if (calctype == "") {
    calctype = config.primclex().settings().default_clex().calctype;
  }
  return calc_properties_path(config.primclex(), config.name(), calctype);
}

template <typename ConfigType>
std::string pos_path(const ConfigType &config) {
  return pos_path(config.primclex(), config.name());
}

template <typename ConfigType>
std::string calc_status_path(const ConfigType &config, std::string calctype) {
  if (calctype == "") {
    calctype = config.primclex().settings().default_clex().calctype;
  }
  return calc_status_path(config.primclex(), config.name(), calctype);
}

/// \brief Return true if all required properties are included in the JSON
bool is_calculated(const MappedProperties &calc_properties,
                   const std::vector<std::string> &required_properties) {
  return std::all_of(required_properties.begin(), required_properties.end(),
                     [&](const std::string &key) {
                       return (calc_properties.global.count(key) ||
                               calc_properties.site.count(key));
                     });
}

std::string calc_properties_path(const PrimClex &primclex,
                                 const std::string &configname,
                                 std::string calctype) {
  if (calctype == "") {
    calctype = primclex.settings().default_clex().calctype;
  }
  return primclex.dir().calculated_properties(configname, calctype).string();
}

std::string pos_path(const PrimClex &primclex, const std::string &configname) {
  return primclex.dir().POS(configname).string();
}

std::string calc_status_path(const PrimClex &primclex,
                             const std::string &configname,
                             std::string calctype) {
  if (calctype == "") {
    calctype = primclex.settings().default_clex().calctype;
  }

  if (configname[0] == 'd') {
    jsonParser calcjson;
    std::vector<std::string> name;
    boost::split(name, configname, boost::is_any_of("/"),
                 boost::token_compress_on);
    if (fs::exists(primclex.dir().configuration_calc_settings_dir(configname,
                                                                  calctype) /
                   "calc.json")) {
      calcjson.read(
          primclex.dir().configuration_calc_dir(configname, calctype) /
          "calc.json");
    } else if (fs::exists(primclex.dir().supercell_calc_settings_dir(
                              name[0] + name[1] + name[2], calctype) /
                          "calc.json")) {
      calcjson.read(primclex.dir().supercell_calc_settings_dir(
                        name[0] + name[1] + name[2], calctype) /
                    "calc.json");
    } else {
      calcjson.read(primclex.dir().calc_settings_dir(calctype) / "calc.json");
    }
    int n_images;
    calcjson.get_else<int>(n_images, "n_images", 0);
    fs::path tmp = primclex.dir().calc_status(configname, calctype);
    return (tmp.parent_path() / ("/N_images_" + std::to_string(n_images)) /
            tmp.filename())
        .string();
  }
  return primclex.dir().calc_status(configname, calctype).string();
}

#define INST_ConfigType(r, data, type)                                         \
  template bool is_calculated(const type &config, std::string calctype);       \
  template void reset_properties(type &config);                                \
  template std::string calc_status(const type &_config, std::string calctype); \
  template std::string failure_type(const type &config, std::string calctype); \
  template bool has_calc_status(const type &config, std::string calctype);     \
  template bool has_failure_type(const type &config, std::string calctype);    \
  template std::string calc_properties_path(const type &config,                \
                                            std::string calctype);             \
  template std::string pos_path(const type &config);                           \
  template std::string calc_status_path(const type &config,                    \
                                        std::string calctype);                 \
  template class Calculable<CRTPBase<type>>;

BOOST_PP_SEQ_FOR_EACH(INST_ConfigType, _, CASM_DB_CONFIG_TYPES)

}  // namespace CASM
