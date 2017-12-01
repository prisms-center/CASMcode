#ifndef CASM_Named
#define CASM_Named

#include <string>
#include "casm/CASM_global_definitions.hh"

namespace CASM {
  class PrimClex;

  namespace DB {
    template<typename T> class ValDatabase;

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

      Named() ;

      std::string name() const;

      /// \brief Return "alias" if object stored in database and alias exists,
      ///        return empty string otherwise
      std::string alias() const;

    protected:

      /// \brief Unset "name", if object is modified
      void clear_name() const;

      /// \brief Regenerate "name"
      void regenerate_name() const;

      /// \brief Set "name", explicity
      void set_name(std::string _name) const;

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

      Indexed();

      // ID for object stored in database. Typically an integer as a string, or 'none'.
      std::string id() const;

      // If 'id' != 'none', just return configname; else regenerate name
      std::string name() const;

    protected:

      /// \brief Unset "id" and "name", if object is modified
      void clear_name() const;

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(Index _id) const;

      /// \brief Set id
      ///
      /// Setting id should be done through Database<Derived> implementations
      ///
      /// - protected, to allow reading Derived from database and setting id
      void set_id(std::string _id) const;


    private:

      friend ValDatabase<MostDerived>;

      mutable std::string m_id;
    };

  }
}

#endif
