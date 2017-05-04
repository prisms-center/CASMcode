#include "casm/database/DatabaseTypes.hh"
//#include "casm/clex/PrimClex.hh"
//#include "casm/database/DatabaseTypeTraits.hh"
//#include "casm/database/DatabaseDefs.hh"
//#include "casm/clex/ConfigurationTraits.hh"

namespace CASM {

  template<>
  struct traits<Configuration> {
    static const std::string name;
    static const std::string short_name;
  };

  template<>
  struct traits<Supercell> {
    static const std::string name;
    static const std::string short_name;
  };

  namespace DB {
    namespace {

      // --- helpers -------------------------------------

      struct insert_name {
        insert_name(std::set<std::string> &_s) : s(_s) {}
        std::set<std::string> &s;
        template<typename T> void eval() {
          s.insert(traits<T>::name);
        }
      };

      struct insert_short_name {
        insert_short_name(std::set<std::string> &_s) : s(_s) {}
        std::set<std::string> &s;
        template<typename T> void eval() {
          s.insert(traits<T>::short_name);
        }
      };

      /*
      struct ConfigCountImpl {
        ConfigCountImpl(std::string _scelname, const PrimClex &_primclex) :
          count(0), scelname(_scelname), primclex(_primclex) {}

        template<typename T>
        void eval() {
          count += primclex.db<T>().scel_range(scelname).size();
        }

        Index count;
        std::string scelname;
        const PrimClex &primclex;
      };
      */
    }

    const std::set<std::string> &types() {
      static std::set<std::string> _types;
      if(!_types.size()) {
        for_each_type(insert_name(_types));
      }
      return _types;
    };

    const std::set<std::string> &types_short() {
      static std::set<std::string> _types_short;
      if(!_types_short.size()) {
        for_each_type(insert_short_name(_types_short));
      }
      return _types_short;
    };


    const std::set<std::string> &config_types() {
      static std::set<std::string> _config_types;
      if(!_config_types.size()) {
        for_each_config_type(insert_name(_config_types));
      }
      return _config_types;
    };

    const std::set<std::string> &config_types_short() {
      static std::set<std::string> _config_types_short;
      if(!_config_types_short.size()) {
        for_each_config_type(insert_short_name(_config_types_short));
      }
      return _config_types_short;
    };

    /*
        /// Total number of configs of all types in a supercell
        Index config_count(std::string scelname, const PrimClex &primclex) {
          ConfigCountImpl f(scelname, primclex);
          if(primclex.db<Supercell>().count(scelname)) {
            for_each_config_type(f);
          }
          return f.count;
        }
    */
  }
}
