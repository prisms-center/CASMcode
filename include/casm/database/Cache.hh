#ifndef CASM_Cache
#define CASM_Cache

#include <string>
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  namespace DB {

    /// \brief Store data in JSON
    class Cache {

    public:

      Cache() :
        m_cache_updated(false) {}

      /// Insert data in the configuration cache, which will be saved in the database
      ///
      /// - For data that depends on DoF only, not calculated properties or
      ///   composition axes
      /// - Adding to cache is modeled as const, but a flag is set so the data
      ///   will be written by 'to_json'
      /// - Data only written if 'name' does not already exist in cache
      /// - Sets 'cache_updated()' to true
      template<typename T>
      void cache_insert(std::string name, const T &data) const {
        if(!m_cache.contains(name)) {
          m_cache[name] = data;
          m_cache_updated = true;
        }
      }

      /// Access the configuration cache, which will be saved in the database
      const jsonParser &cache() const {
        return m_cache;
      }

      /// Check if cache updated
      bool cache_updated() const {
        return m_cache_updated;
      }

      /// Clear the cache
      ///
      /// - Sets 'cache_updated()' to true
      void cache_clear() const {
        m_cache.put_null();
        m_cache_updated = true;
      }


    protected:

      /// Access the configuration cache, which will be saved in the database
      ///
      /// - Does not change 'cache_updated()' value
      jsonParser &cache() {
        return m_cache;
      }


    private:

      jsonParser m_cache;
      bool m_cache_updated;
    };
  }

}

#endif
