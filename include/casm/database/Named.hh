#ifndef CASM_Named
#define CASM_Named

#include <string>
#include "casm/database/Database.hh"

namespace CASM {
  namespace DB {

    /// CRTP Mixin for 'named' database objects with no id
    ///
    /// - MostDerived should implement:
    ///   - std::string MostDerived::generate_name_impl() const
    ///   - const PrimClex& MostDerived::primclex_impl() const
    /// - Use one of Indexed or Named
    ///
    template<typename Base>
    class Named : public Base {

    public:

      typedef typename Base::MostDerived MostDerived;
      using Base::derived;

      Named() :
        m_name("") {}

      std::string name() const {
        if(m_name.empty()) {
          m_name = derived().generate_name_impl();
        }
        return m_name;
      }

      /// \brief Return "alias" if object stored in database and alias exists,
      ///        return empty string otherwise
      std::string alias() const {
        return derived().primclex().template db<MostDerived>().alias(name());
      }

    protected:

      /// \brief Unset "name", if object is modified
      void clear_name() const {
        m_name = "";
      }

    private:

      mutable std::string m_name;

    };

    /// Similar to 'Named', but includes an incrementing 'id' string
    ///
    /// - Setting id should be done through Database<Derived> implementations of
    /// insert or emplace.
    /// - Use one of Indexed or Named
    ///
    template<typename _Base>
    class Indexed : public Named<_Base> {

    public:

      typedef Named<_Base> Base;
      typedef typename Base::MostDerived MostDerived;
      using Base::derived;

      Indexed() :
        m_id("none") {}

      std::string id() const {
        return m_id;
      }

    protected:

      /// \brief Unset "id" and "name", if object is modified
      void clear_name() const {
        m_id = "none";
        Named<_Base>::clear_name();
      }

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations of
      /// insert or emplace.
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(Index _id) const {
        set_id(std::to_string(_id));
      }

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations of
      /// insert or emplace.
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(std::string _id) const {
        Named<_Base>::clear_name();
        m_id = _id;
      }


    private:

      friend ValDatabase<MostDerived>;

      mutable std::string m_id;
    };

  }
}

#endif
