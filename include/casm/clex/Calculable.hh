#ifndef CASM_Calculable
#define CASM_Calculable

#include <vector>
#include <string>

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
  template<typename _Base>
  class Calculable : public DB::Cache, public DB::Indexed<_Base> {

  public:

    typedef typename DB::Indexed<_Base> Base;
    using Base::derived;

    const jsonParser calc_properties(std::string calctype = "") const;

    void set_calc_properties(const jsonParser &json, std::string calctype);

    /// \brief grabs properties from the indicated calctype and adds info to calc_properties_map
    void refresh_calc_properties(std::string calctype);

    const jsonParser &source() const;

    void set_source(const jsonParser &source);

    void push_back_source(const jsonParser &source);

  protected:

    /// Call in MostDerived any time DoF may be modified
    void _modify_dof();

  private:

    std::map<std::string, jsonParser> m_calc_properties_map;
    jsonParser m_source;
  };

  /// \brief Return true if all required properties have been been calculated for
  /// the configuration
  template<typename ConfigType>
  bool is_calculated(const ConfigType &config);

  template<typename ConfigType>
  void reset_properties(ConfigType &config);

  /// \brief Status of calculation
  template<typename ConfigType>
  std::string calc_status(const ConfigType &_config);

  // \brief Reason for calculation failure.
  template<typename ConfigType>
  std::string failure_type(const ConfigType &config);

  template<typename ConfigType>
  bool has_calc_status(const ConfigType &config);

  template<typename ConfigType>
  bool has_failure_type(const ConfigType &config);

  template<typename ConfigType>
  fs::path calc_properties_path(const ConfigType &config);

  template<typename ConfigType>
  fs::path pos_path(const ConfigType &config);

  template<typename ConfigType>
  fs::path calc_status_path(const ConfigType &config);

  /// \brief Read properties.calc.json from training_data
  template<typename ConfigType>
  std::tuple<jsonParser, bool, bool> read_calc_properties(const ConfigType &config);

  /// \brief Read properties.calc.json from file
  template<typename ConfigType>
  std::tuple<jsonParser, bool, bool> read_calc_properties(const PrimClex &primclex, const fs::path &filepath);


  /// \brief Return true if all required properties are included in the JSON
  bool is_calculated(
    const jsonParser &calc_properties,
    const std::vector<std::string> &required_properties);

  fs::path calc_properties_path(const PrimClex &primclex, const std::string &configname);

  fs::path pos_path(const PrimClex &primclex, const std::string &configname);

  fs::path calc_status_path(const PrimClex &primclex, const std::string &configname);

}

#endif
