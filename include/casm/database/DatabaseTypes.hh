#ifndef CASM_DatabaseTypes
#define CASM_DatabaseTypes

#include <set>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_TMP.hh"
#include "casm/database/DatabaseTypeTraits.hh"

// 1) add Database types here
#define CASM_DB_NONCONFIG_TYPES (Supercell) (Kinetics::PrimPeriodicDiffTransOrbit)
#define CASM_DB_CONFIG_TYPES (Configuration)

// 2) add #include locations as appropriate to:
//    DatabaseTypeDefs.hh, DatabaseTypeTraits.hh,
//    ConfigTypeDefs.hh, ConfigTypeTraits.hh

// --- the rest works based on ^ -----------------------------------------------

#define CASM_DB_TYPES CASM_DB_NONCONFIG_TYPES CASM_DB_CONFIG_TYPES

namespace CASM {

  class PrimClex;

  namespace DB {

    // List of all types stored in a Database
    typedef std::tuple<BOOST_PP_SEQ_ENUM(CASM_DB_TYPES)> DataObjectTypeTuple;

    // List of all configtypes which have calculations stored in a PropertiesDatabase
    typedef std::tuple<BOOST_PP_SEQ_ENUM(CASM_DB_CONFIG_TYPES)> ConfigTypeTuple;


    // -- SFINAE helper --

    template<typename T>
    using IfConfigType = std::enable_if<CASM_TMP::has_type<T, ConfigTypeTuple>::value, T>;

    template<typename T>
    using IfNotConfigType = std::enable_if < !CASM_TMP::has_type<T, ConfigTypeTuple>::value, T >;

    // -- List of Database DataObject and ConfigType --

    template<typename F>
    void for_each_type(F f) {
      return CASM_TMP::for_each_type<DataObjectTypeTuple, F>(f);
    }

    template<typename F>
    void for_type(std::string name, F f) {
      return CASM_TMP::for_type<DataObjectTypeTuple, F>(name, f);
    }

    template<typename F>
    void for_type_short(std::string short_name, F f) {
      return CASM_TMP::for_type_short<DataObjectTypeTuple, F>(short_name, f);
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

    template<typename F>
    void for_config_type_short(std::string short_name, F f) {
      return CASM_TMP::for_type_short<ConfigTypeTuple, F>(short_name, f);
    }

    /// std::set of all QueryTraits<ConfigType>::name
    const std::set<std::string> &config_types();

    /// std::set of all QueryTraits<ConfigType>::short_name
    const std::set<std::string> &config_types_short();

    /// Total number of configs of all types in a supercell
    Index config_count(std::string scelname, const PrimClex &primclex);

    /// Total number of configs of a specific type in a supercell
    Index config_count(std::string configtype, std::string scelname, const PrimClex &primclex);

    /// Total number of calculated configs of all types in a supercell
    Index config_calculated_count(std::string scelname, const PrimClex &primclex);

    /// Total number of calculated configs of a specific type in a supercell
    Index config_calculated_count(std::string configtype, std::string scelname, const PrimClex &primclex);

    /// Total number of configs w/ data or files of all types in a supercell
    Index config_data_count(std::string scelname, const PrimClex &primclex);

    /// Total number of configs w/ data or files of a specific type in a supercell
    Index config_data_count(std::string configtype, std::string scelname, const PrimClex &primclex);

  }
}

#endif
