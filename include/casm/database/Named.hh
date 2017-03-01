#ifndef CASM_Named
#define CASM_Named

#include <string>
#include "casm/database/Database.hh"

namespace CASM {

  namespace DB {

    /// Derived should implement:
    /// - std::string generate_name() const
    /// - const PrimClex& primclex() const
    ///
    template<typename Derived>
    class Named {

    public:

      Named() :
        m_name("") {}

      std::string name() const {
        if(m_name.empty()) {
          m_name = derived().generate_name();
        }
        return m_name;
      }

      /// \brief Return "alias" if object stored in database and alias exists,
      ///        return empty string otherwise
      std::string alias() const {
        return derived().primclex().template db<Derived>().alias(name());
      }

      /// \brief Unset "name", if object is modified
      void clear_name() {
        m_name = "";
      }


    private:

      friend ValDatabase<Derived>;

      const Derived &derived() const {
        return *static_cast<const Derived *>(this);
      }

      mutable std::string m_name;

    };

    /// Similar to 'Named', but includes an incrementing 'id' string
    ///
    /// Setting id should be done through Database<Derived> implementations of
    /// insert or emplace.
    ///
    template<typename Derived>
    class Indexed : public Named<Derived> {

    public:

      Indexed() :
        m_id("none") {}

      std::string id() const {
        return m_id;
      }

      /// \brief Unset "id" and "name", if object is modified
      void clear_name() {
        m_id = "none";
        Named<Derived>::clear_name();
      }

    protected:

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations of
      /// insert or emplace.
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(Index _id) {
        m_id = std::to_string(_id);
      }

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations of
      /// insert or emplace.
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(std::string _id) {
        m_id = _id;
      }


    private:

      friend ValDatabase<Derived>;

      std::string m_id;
    };

  }
}

#endif
