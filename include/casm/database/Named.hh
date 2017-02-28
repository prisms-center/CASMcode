#ifndef CASM_Named
#define CASM_Named

#include <string>
#include "casm/database/Database.hh"

namespace CASM {

  namespace DB {

    /// Derived should implement:
    /// - std::string generate_name() const
    ///
    /// Setting alias should be done through ValDatabase<Derived>::set_alias in
    /// order to prevent multiple objects getting the same alias
    ///
    template<typename Derived>
    class Named {

    public:

      Named() :
        m_name(""),
        m_alias("") {}

      std::string name() const {
        if(m_name.empty()) {
          m_name = derived().generate_name();
        }
        return m_name;
      }

      /// \brief User-specified alternative to 'name'
      std::string alias() const {
        return m_alias;
      }


      /// \brief Return alias, if exists, else name
      std::string alias_or_name() const {
        if(alias().empty()) {
          return name();
        }
        return alias();
      }

      /// \brief Unset "name" and "alias", if object is modified
      void clear_name() {
        m_name = "";
        m_alias.clear();
      }

    protected:

      /// \brief Set alias
      ///
      /// Setting alias should be done through ValDatabase<Derived>::set_alias in
      /// order to prevent multiple objects getting the same alias
      ///
      /// - protected, to allow reading Derived from database and setting alias
      void set_alias(const std::string &_alias) {
        m_alias = _alias;
      }


    private:

      friend ValDatabase<Derived>;

      const Derived &derived() const {
        return *static_cast<const Derived *>(this);
      }

      mutable std::string m_name;
      std::string m_alias;

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

      /// \brief Unset "id", "name", and "alias", if object is modified
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
