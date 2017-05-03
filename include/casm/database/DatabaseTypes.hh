#ifndef CASM_DatabaseTypes
#define CASM_DatabaseTypes

#include <set>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_TMP.hh"

namespace CASM {

  // add Database types:
  // 1) here
  // 2) in DataObjectTypeTuple & ConfigTypeTuple if appropriate
  // 3) in DatabaseDefs.hh & ConfigDatabaseDefs.hh if appropriate
  // 4) in DatabaseTypeTraits.hh & ConfigTypeTraits.hh if appropriate
  // 5) in DatabaseTypeDefs.hh & ConfigTypeDefs.hh if appropriate

  class Supercell;
  class Configuration;

  namespace DB {

    // List of all types stored in a Database
    typedef std::tuple<Supercell, Configuration> DataObjectTypeTuple;

    // List of all configtypes which have calculations stored in a PropertiesDatabase
    typedef std::tuple<Configuration> ConfigTypeTuple;


    // -- List of Database DataObject and ConfigType --

    template<typename F>
    void for_each_type(F f) {
      return CASM_TMP::for_each_type<DataObjectTypeTuple, F>(f);
    }

    template<typename F>
    void for_type(std::string name, F f) {
      return CASM_TMP::for_type<DataObjectTypeTuple, F>(name, f);
    }

    /// std::set of all QueryTraits<DataObject>::name
    const std::set<std::string> &types();

    /// std::set of all QueryTraits<DataObject>::short_name
    const std::set<std::string> &types_short();



    template<typename F>
    void for_each_config_type(F f) {
      return CASM_TMP::for_each_type<ConfigTypeTuple, F>(f);
    }

    template<typename F>
    void for_config_type(std::string name, F f) {
      return CASM_TMP::for_type<ConfigTypeTuple, F>(name, f);
    }

    /// std::set of all QueryTraits<ConfigType>::name
    const std::set<std::string> &config_types();

    /// std::set of all QueryTraits<ConfigType>::short_name
    const std::set<std::string> &config_types_short();

    /// Total number of configs of all types in a supercell
    Index config_count(std::string scelname);

  }
}

#endif
