#ifndef CASM_PropertiesDatabase
#define CASM_PropertiesDatabase

#include <string>
#include <set>
#include <boost/iterator/iterator_facade.hpp>
#include "casm/database/Database.hh"
#include "casm/database/MappedProperties.hh"

namespace CASM {
  namespace DB {

    class PropertiesDatabaseIterator;

    /// Dereferencing PropertiesDatabaseIteratorBase only provides const references, whether
    /// the underlying resource is persistent (jsonPropertiesDatabase) or temporary
    /// (other database types). Changing database entries must be done via copy then
    /// insert, or update.
    ///
    /// PropertiesDatabaseIteratorBase should always be dereferenceable (except
    /// end or when default constructed), though the reference may be invalidated
    /// when a second PropertiesDatabaseIteratorBase is dereferenced, dereferencing
    /// the first again should be valid (though it may require re-allocation):
    ///
    /// \code
    /// DerivedDatabaseIterator A = ... get from DerivedDatabase ...;
    /// const MappedProperties& A_ref1 = *A;  // ok
    /// A_ref1.method();  // ok
    /// A->method(); // ok
    ///
    /// DerivedDatabaseIterator B = ... get from DerivedDatabase ...;
    /// const MappedProperties& B_ref = *B;  // ok, but now A_ref1 may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    ///
    /// const MappedProperties& A_ref2 = *A;  // ok, but now B_ref may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A_ref2.method();  // ok
    /// B_ref.method(); // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    /// \endcode
    ///
    /// Derived classes must implement private methods:
    /// - std::string name() const
    /// - void increment()
    /// - reference dereference() const
    /// - PropertiesDatabaseIteratorBase *_clone() const
    ///
    class PropertiesDatabaseIteratorBase {

      friend PropertiesDatabaseIterator;

    public:

      typedef MappedProperties value_type;
      typedef const value_type &reference;

      PropertiesDatabaseIteratorBase() {}
      virtual ~PropertiesDatabaseIteratorBase() {}

      std::unique_ptr<PropertiesDatabaseIteratorBase> clone() const {
        return std::unique_ptr<PropertiesDatabaseIteratorBase>(this->_clone());
      }

    private:

      virtual void increment() = 0;

      virtual reference dereference() const = 0;

      virtual bool equal(const PropertiesDatabaseIteratorBase &other) const = 0;

      virtual PropertiesDatabaseIteratorBase *_clone() const = 0;

      virtual long distance_to(const PropertiesDatabaseIteratorBase &other) const = 0;

    };


    /// \brief Wrapper class for specializations PropertiesDatabaseIteratorBase
    ///
    /// - Gives all specialized Database<ValueType> the same iterator type
    ///
    /// Dereferencing PropertiesDatabaseIterator only provides const references, whether
    /// the underlying resource is persistent (jsonPropertiesDatabase) or temporary
    /// (other types). Changing database entries must be done via copy then
    /// insert, or update.
    ///
    /// DatabaseIterator should always be dereferenceable (except end or when
    /// default constructed), though the reference may be invalidated when a
    /// second DatabaseIterator is dereferenced, dereferencing the first
    /// again should be valid (though it may require re-allocation):
    ///
    /// \code
    /// PropertiesDatabaseIterator A = primclex.properties<ValueType>().find_via_from(A_name);
    /// const MappedProperties& A_ref1 = *A;  // ok
    /// A_ref1.method();  // ok
    /// A->method(); // ok
    ///
    /// PropertiesDatabaseIterator<ValueType> B = primclex.properties<ValueType>().find_via_from(B_name);
    /// const MappedProperties& B_ref = *B;  // ok, but now A_ref1 may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    ///
    /// const MappedProperties& A_ref2 = *A;  // ok, but now B_ref may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A_ref2.method();  // ok
    /// B_ref.method(); // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    /// \endcode

    class PropertiesDatabaseIterator :

      public boost::iterator_facade <
      PropertiesDatabaseIterator,
      MappedProperties,
      std::forward_iterator_tag,
      const MappedProperties &,
      long > {

    public:

      /// Default constructor
      PropertiesDatabaseIterator() {}

      /// Construct iterator
      PropertiesDatabaseIterator(const PropertiesDatabaseIteratorBase &it) :
        m_ptr(notstd::clone(it)) {}

      PropertiesDatabaseIteratorBase *get() const {
        return m_ptr.unique().get();
      }

    private:

      friend boost::iterator_core_access;

      /// boost::iterator_facade implementation
      void increment() {
        m_ptr->increment();
      }

      /// boost::iterator_facade implementation
      const MappedProperties &dereference() const {
        return m_ptr->dereference();
      }

      /// boost::iterator_facade implementation
      bool equal(const PropertiesDatabaseIterator &B) const {
        return m_ptr->equal(*(B.m_ptr));
      }

      long distance_to(const PropertiesDatabaseIterator &B) const {
        return m_ptr->distance_to(*(B.m_ptr));
      }

      notstd::cloneable_ptr<PropertiesDatabaseIteratorBase> m_ptr;
    };


    class PropertiesDatabase : public DatabaseBase {

    public:

      /// \brief Compare type for std::set<std::string, Compare>, which sorts
      ///        to determine the 'from' configname for the best scoring MappedProperties
      ///        to a particular config
      class Compare {

      public:

        Compare(const PropertiesDatabase *_map, std::string _to_configname, const ScoreMappedProperties &_score) :
          m_map(_map),
          m_to(_to_configname),
          m_score(_score) {}

        /// \brief Compare mapped properties 'from_A' and 'from_B', preferring self-mapped results
        bool operator()(const std::string &from_A, const std::string &from_B) const;

        const ScoreMappedProperties &score_method() const {
          return m_score;
        }

      private:
        const PropertiesDatabase *m_map;
        std::string m_to;
        ScoreMappedProperties m_score;
      };

      typedef MappedProperties value_type;
      typedef PropertiesDatabaseIterator iterator;
      typedef Index size_type;

      PropertiesDatabase(const PrimClex &_primclex) :
        DatabaseBase(_primclex) {}


      /// \brief Begin iterator over pair<from_configname, pair<to_configname, data> >
      virtual iterator begin() const = 0;

      /// \brief End iterator over pair<from_configname, data>
      virtual iterator end() const = 0;

      virtual size_type size() const = 0;
      bool empty() const {
        return this->size() == 0;
      }

      /// \brief Return iterator to pair<from_configname, data> that is the best mapping to specified config
      ///
      /// - Prefers self-mapped, else best scoring
      virtual iterator find_via_to(std::string to_configname) const = 0;

      /// \brief Return iterator to pair<from_configname, data> that is from the specified config
      virtual iterator find_via_from(std::string from_configname) const = 0;


      /// \brief Name of the configuration that relaxed 'from'->'to'
      std::string relaxed_to(std::string from_configname) const {
        return find_via_from(from_configname)->to;
      }

      /// \brief Name of best scoring configuration that relaxed 'from'->'to'
      ///
      /// - Prefers self-mapped
      std::string relaxed_from(std::string to_configname) const {
        return find_via_to(to_configname)->from;
      }

      /// \brief Names of all configurations that relaxed 'from'->'to'
      ///
      /// - Empty set if none
      virtual std::set<std::string, Compare> relaxed_from_all(std::string to_configname) const = 0;

      /// \brief Change the score method for a single configuration
      virtual void set_score_method(std::string to_configname, const ScoreMappedProperties &score) = 0;

      /// \brief Name of ScoreMappedProperties method
      ScoreMappedProperties score_method(std::string to_configname) const {
        return relaxed_from_all(to_configname).value_comp().score_method();
      }

      /// \brief Best score of configurations that relaxed 'from'->'to'
      double best_score(std::string to_configname) const {
        return score(*find_via_to(to_configname));
      }

      /// \brief Score mapping 'from'->'to'
      double score(std::string from_configname) const {
        return score(*find_via_from(from_configname));
      }

      /// \brief Score mapping 'from'->'to'
      double score(const MappedProperties &value) const {
        return score_method(value.to)(value);
      }

      /// \brief Insert data
      std::pair<iterator, bool> insert(const MappedProperties &value);

      /// \brief Erase data
      iterator erase(iterator pos);

      /// \brief Erase data
      size_type erase_via_from(std::string from_configname) {
        auto it = find_via_from(from_configname);
        if(it == end()) {
          return 0;
        }
        erase(it);
        return 1;
      }


    private:

      /// \brief Private _insert MappedProperties, without modifying 'relaxed_from'
      virtual std::pair<iterator, bool> _insert(const MappedProperties &value) = 0;

      /// \brief Private _erase MappedProperties, without modifying 'relaxed_from'
      virtual iterator _erase(iterator pos) = 0;

      /// \brief Set sorted container of names of all configurations that relaxed 'from'->'to'
      virtual void _set_relaxed_from_all(
        std::string to_configname,
        const std::set<std::string, Compare> &_set) = 0;

    };

  }
}

#endif
