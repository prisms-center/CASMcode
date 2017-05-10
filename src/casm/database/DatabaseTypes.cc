#include <boost/range/iterator_range.hpp>
#include "casm/clex/PrimClex.hh"
#include "casm/database/Import.hh"
#include "casm/database/DatabaseTypeTraits.hh"
#include "casm/database/DatabaseDefs.hh"
#include "casm/database/DatabaseTypes.hh"

namespace CASM {

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

      struct ConfigCountImplBase {
        ConfigCountImplBase(std::string _scelname, const PrimClex &_primclex) :
          count(0), scelname(_scelname), primclex(_primclex) {}

        Index count;
        std::string scelname;
        const PrimClex &primclex;
      };

      struct ConfigCountImpl : public ConfigCountImplBase {
        using ConfigCountImplBase::ConfigCountImplBase;

        template<typename T>
        void eval() {
          count += primclex.db<T>().scel_range_size(scelname);
        }
      };

      struct ConfigCalculatedCountImpl : public ConfigCountImplBase {
        using ConfigCountImplBase::ConfigCountImplBase;

        template<typename T>
        void eval() {
          for(const auto &config : primclex.db<T>().scel_range(scelname)) {
            if(is_calculated(config)) {
              count += 1;
            }
          }
        }
      };

      struct ConfigDataCountImpl : public ConfigCountImplBase {
        using ConfigCountImplBase::ConfigCountImplBase;

        template<typename T>
        void eval() {
          ConfigData<T> data(primclex, null_log());
          auto it = primclex.db<T>().scel_range(scelname).begin();
          auto end = primclex.db<T>().scel_range(scelname).end();
          for(; it != end; ++it) {
            if(data.has_existing_data_or_files(it.name())) {
              count += 1;
            }
          }
        }
      };

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

    /// Total number of configs of all types in a supercell
    Index config_count(std::string scelname, const PrimClex &primclex) {
      ConfigCountImpl f(scelname, primclex);
      if(primclex.db<Supercell>().count(scelname)) {
        for_each_config_type(f);
      }
      return f.count;
    }

    /// Total number of configs of a specific type in a supercell
    Index config_count(std::string configtype, std::string scelname, const PrimClex &primclex) {
      ConfigCountImpl f(scelname, primclex);
      if(primclex.db<Supercell>().count(scelname)) {
        for_config_type_short(configtype, f);
      }
      return f.count;
    }

    /// Total number of calculated configs of all types in a supercell
    Index config_calculated_count(std::string scelname, const PrimClex &primclex) {
      ConfigCalculatedCountImpl f(scelname, primclex);
      if(primclex.db<Supercell>().count(scelname)) {
        for_each_config_type(f);
      }
      return f.count;
    }

    /// Total number of calculated configs of a specific type in a supercell
    Index config_calculated_count(std::string configtype, std::string scelname, const PrimClex &primclex) {
      ConfigCalculatedCountImpl f(scelname, primclex);
      if(primclex.db<Supercell>().count(scelname)) {
        for_config_type_short(configtype, f);
      }
      return f.count;
    }

    /// Total number of configs w/ data or files of all types in a supercell
    Index config_data_count(std::string scelname, const PrimClex &primclex) {
      ConfigDataCountImpl f(scelname, primclex);
      if(primclex.db<Supercell>().count(scelname)) {
        for_each_config_type(f);
      }
      return f.count;
    }

    /// Total number of configs w/ data or files of a specific type in a supercell
    Index config_data_count(std::string configtype, std::string scelname, const PrimClex &primclex) {
      ConfigDataCountImpl f(scelname, primclex);
      if(primclex.db<Supercell>().count(scelname)) {
        for_config_type_short(configtype, f);
      }
      return f.count;
    }

  }
}
