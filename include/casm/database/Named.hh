#ifndef CASM_Named
#define CASM_Named

#include <string>
#include "casm/database/Database.hh"

namespace CASM {
  namespace DB {

    template<typename Derived> class Indexed;

    /// Derived should implement:
    /// - std::string _generate_name() const
    /// - const PrimClex& primclex() const
    ///
    template<typename Derived>
    class Named {

    public:

      Named() :
        m_primclex(nullptr),
        m_name("") {}

      Named(const PrimClex &_primclex) :
        m_primclex(&_primclex),
        m_name("") {}

      std::string name() const {
        if(m_name.empty()) {
          m_name = derived()._generate_name();
        }
        return m_name;
      }

      /// \brief Return "alias" if object stored in database and alias exists,
      ///        return empty string otherwise
      std::string alias() const {
        return primclex().template db<Derived>().alias(name());
      }

      /// \brief Unset "name", if object is modified
      void clear_name() {
        m_name = "";
      }

      /// \brief Get PrimClex
      const PrimClex &primclex() const {
        if(!m_primclex) {
          throw std::runtime_error("Error in Orbit::primclex(): PrimClex not valid");
        }
        return *m_primclex;
      }

    private:

      friend ValDatabase<Derived>;
      friend Indexed<Derived>;

      const Derived &derived() const {
        return *static_cast<const Derived *>(this);
      }

      /// \brief Add PrimClex pointer to objects constructed from Database
      void set_primclex(const PrimClex &_primclex) const {
        m_primclex = &_primclex;
      }

      mutable std::string m_name;
      mutable const PrimClex *m_primclex;

    };

    /// Similar to 'Named', but includes an incrementing 'id' string
    ///
    /// Setting id should be done through Database<Derived> implementations of
    /// insert or emplace.
    ///
    template<typename Derived>
    class Indexed : public Named<Derived> {

    public:
      Indexed(const PrimClex &_primclex) :
        Named<Derived>::Named(_primclex),
        m_id("none") {}

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
      void set_id(Index _id) const {
        m_id = std::to_string(_id);
        this->m_name = "";
      }

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations of
      /// insert or emplace.
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(std::string _id) const {
        m_id = _id;
        this->m_name = "";
      }


    private:

      friend ValDatabase<Derived>;

      mutable std::string m_id;
    };

  }
}

#endif
