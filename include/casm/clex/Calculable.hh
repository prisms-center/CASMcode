#ifndef CASM_Calculable
#define CASM_Calculable

#include <string>
#include <vector>

#include "casm/clex/MappedProperties.hh"
#include "casm/database/Cache.hh"
#include "casm/database/Named.hh"

namespace CASM {

class jsonParser;

///
/// - name and calculated properties should be invalidated whenever the
///   ConfigType DoF are modified. Can do this by calling _modify_dof(). This
///   means all DoF should be accessed via functions.
/// - cache should only be used for DoF-dependent properties, not
///   calctype-dependent properties
///
template <typename _Base>
class Calculable : public DB::Cache, public DB::Indexed<_Base> {
 public:
  typedef typename DB::Indexed<_Base> Base;
  typedef typename Base::MostDerived MostDerived;
  using Base::derived;

  /// \brief Return MappedProperties for requested calctype
  MappedProperties const &calc_properties(std::string calctype = "") const;

  void set_calc_properties(const MappedProperties &_prop,
                           std::string calctype = "");

  /// \brief grabs properties from the indicated calctype and adds info to
  /// calc_properties_map
  void refresh_calc_properties(std::string calctype = "");

  const jsonParser &source() const;

  void set_source(const jsonParser &source);

  void push_back_source(const jsonParser &source);

  std::map<std::string, MappedProperties> calc_properties_map() const {
    return m_calc_properties_map;
  }

 protected:
  /// Call in MostDerived any time DoF may be modified
  void _modify_dof();

  /// \brief grabs properties from the indicated calctype and adds info to
  /// calc_properties_map
  void _refresh_calc_properties(std::string calctype = "") const;

 private:
  mutable std::map<std::string, MappedProperties> m_calc_properties_map;
  jsonParser m_source;
};

/// \brief Return true if all required properties have been been calculated for
/// the configuration
template <typename ConfigType>
bool is_calculated(const ConfigType &config, std::string calctype = "");

template <typename ConfigType>
void reset_properties(ConfigType &config);

/// \brief Status of calculation
/// A calculation that has been run successfully will be marked 'complete'.
/// Mapped or imported configurations may have 'is_calculated = 1' without
/// 'calc_status = complete'.
template <typename ConfigType>
std::string calc_status(const ConfigType &_config, std::string calctype = "");

// \brief Reason for calculation failure.
template <typename ConfigType>
std::string failure_type(const ConfigType &config, std::string calctype = "");

template <typename ConfigType>
bool has_calc_status(const ConfigType &config, std::string calctype = "");

template <typename ConfigType>
bool has_failure_type(const ConfigType &config, std::string calctype = "");

template <typename ConfigType>
std::string calc_properties_path(const ConfigType &config,
                                 std::string calctype = "");

template <typename ConfigType>
std::string pos_path(const ConfigType &config, std::string calctype = "");

template <typename ConfigType>
std::string calc_status_path(const ConfigType &config,
                             std::string calctype = "");

/// \brief Return true if all required properties are included in the JSON
bool is_calculated(const MappedProperties &calc_properties,
                   const std::vector<std::string> &required_properties);

std::string calc_properties_path(const PrimClex &primclex,
                                 const std::string &configname,
                                 std::string calctype = "");

std::string pos_path(const PrimClex &primclex, const std::string &configname);

std::string calc_status_path(const PrimClex &primclex,
                             const std::string &configname,
                             std::string calctype = "");

}  // namespace CASM

#endif
