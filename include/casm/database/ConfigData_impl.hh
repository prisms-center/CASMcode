#ifndef CASM_DB_ConfigData_impl
#define CASM_DB_ConfigData_impl

#include "casm/database/ConfigData.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {
  namespace DB {

    // --- class ConfigData ---

    template<typename _ConfigType>
    Database<_ConfigType> &ConfigData<_ConfigType>::db_config() const {
      return primclex().template db<ConfigType>();
    }

    template<typename _ConfigType>
    PropertiesDatabase &ConfigData<_ConfigType>::db_props() const {
      return primclex().template db_props<ConfigType>(
               primclex().settings().default_clex().calctype);
    }

  }
}

#endif
