#ifndef CASM_Named_impl
#define CASM_Named_impl

#include "casm/database/Named.hh"
#include "casm/database/Database_impl.hh"

namespace CASM {
  namespace DB {

    template<typename Base>
    Named<Base>::Named() :
      m_name("") {}

    template<typename Base>
    std::string Named<Base>::name() const {
      if(m_name.empty()) {
        regenerate_name();
      }
      return m_name;
    }

    /// \brief Return "alias" if object stored in database and alias exists,
    ///        return empty string otherwise
    template<typename Base>
    std::string Named<Base>::alias() const {
      return derived().primclex().template db<MostDerived>().alias(name());
    }

    /// \brief Unset "name", if object is modified
    template<typename Base>
    void Named<Base>::clear_name() const {
      m_name = "";
    }

    /// \brief Regenerate "name"
    template<typename Base>
    void Named<Base>::regenerate_name() const {
      m_name = derived().generate_name_impl();
    }

    /// \brief Set "name", explicity modified
    template<typename Base>
    void Named<Base>::set_name(std::string _name) const {
      m_name = _name;
    }


    template<typename _Base>
    Indexed<_Base>::Indexed() :
      m_id("none") {}

    template<typename _Base>
    std::string Indexed<_Base>::id() const {
      return m_id;
    }

    // If 'id' is already known, just return configname
    template<typename _Base>
    std::string Indexed<_Base>::name() const {
      if(id() == "none") {
        this->regenerate_name();
      }
      return Named<_Base>::name();
    }

    /// \brief Unset "id" and "name", if object is modified
    template<typename _Base>
    void Indexed<_Base>::clear_name() const {
      m_id = "none";
      Named<_Base>::clear_name();
    }

    /// \brief Set id
    ///
    /// Setting id should be done through Database<Derived> implementations
    ///
    /// - protected, to allow reading Derived from database and setting id
    template<typename _Base>
    void Indexed<_Base>::set_id(Index _id) const {
      set_id(std::to_string(_id));
    }

    /// \brief Set id
    ///
    /// Setting id should be done through Database<Derived> implementations
    ///
    /// - protected, to allow reading Derived from database and setting id
    template<typename _Base>
    void Indexed<_Base>::set_id(std::string _id) const {
      Named<_Base>::clear_name();
      m_id = _id;
    }

  }
}

#endif
