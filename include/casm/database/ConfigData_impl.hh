#ifndef CASM_DB_ConfigData_impl
#define CASM_DB_ConfigData_impl

#include "casm/database/ConfigData.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/ProjectSettings.hh"

namespace CASM {
  namespace DB {

    // --- class ConfigData ---
    template<typename ConfigType>
    ConfigData::ConfigData(const PrimClex &_primclex, TypeTag<ConfigType>) :
      m_primclex(_primclex) {
      m_db_props_func = [&]()->PropertiesDatabase& {
        return primclex().template db_props<ConfigType>(
          primclex().settings().default_clex().calctype);
      };
    }


    template<typename ConfigType>
    Database<ConfigType> &ConfigData::db_config() const {
      return primclex().template db<ConfigType>();
    }

    //template<typename ConfigType>
    //PropertiesDatabase &ConfigData<ConfigType>::db_props() const {
    //return primclex().template db_props<ConfigType>(
    //         primclex().settings().default_clex().calctype);
    //}

  }
}

#endif
