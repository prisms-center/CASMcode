#ifndef CASM_Database
#define CASM_Database

#include <iterator>
#include <memory>
#include <map>
#include <algorithm>
#include <boost/iterator/iterator_facade.hpp>
#include "casm/misc/cloneable_ptr.hh"
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  class PrimClex;

  namespace DB {

    /// Fully generic database interface for use by DatabaseHandler
    class DatabaseBase {

    public:

      DatabaseBase(const PrimClex &_primclex) :
        m_primclex(&_primclex) {}

      virtual ~DatabaseBase() {}

      virtual DatabaseBase &open() = 0;
      virtual void commit() = 0;
      virtual void close() = 0;

      const PrimClex &primclex() const {
        return *m_primclex;
      }

    private:

      const PrimClex *m_primclex;
    };

    template<typename ValueType>
    class DatabaseIterator;

    /// Dereferencing DatabaseIteratorBase only provides const references, whether
    /// the underlying resource is persistent (Supercell) or temporary
    /// (other types). Changing database entries must be done via copy then
    /// insert, or update.
    ///
    /// DatabaseIteratorBase should always be dereferenceable (except end or when
    /// default constructed), though the reference may be invalidated when a
    /// second DatabaseIteratorBase is dereferenced, dereferencing the first
    /// again should be valid (though it may require re-allocation):
    ///
    /// \code
    /// DerivedDatabaseIterator A = ... get from DerivedDatabase ...;
    /// const ValueType& A_ref1 = *A;  // ok
    /// A_ref1.method();  // ok
    /// A->method(); // ok
    ///
    /// DerivedDatabaseIterator B = ... get from DerivedDatabase ...;
    /// const ValueType& B_ref = *B;  // ok, but now A_ref1 may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    ///
    /// const ValueType& A_ref2 = *A;  // ok, but now B_ref may be invalidated
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
    /// - DatabaseIteratorBase *_clone() const
    ///
    template<typename ValueType>
    class DatabaseIteratorBase {

      friend DatabaseIterator<ValueType>;

    public:

      typedef ValueType value_type;
      typedef const value_type &reference;

      DatabaseIteratorBase() {}

      virtual std::string name() const = 0;

      std::unique_ptr<DatabaseIteratorBase> clone() const {
        return std::unique_ptr<DatabaseIteratorBase>(this->_clone());
      }

    private:

      virtual void increment() = 0;

      virtual reference dereference() const = 0;

      virtual bool equal(const DatabaseIteratorBase<ValueType> &other) const = 0;

      virtual DatabaseIteratorBase *_clone() const = 0;

      /* Do not implement this
      virtual long distance_to(const DatabaseIteratorBase<ValueType> &other) const = 0;
      */
    };

    /// \brief Wrapper class for specializations DatabaseIteratorBase
    ///
    /// - Gives all specialized Database<ValueType> the same iterator type
    ///
    /// Dereferencing DatabaseIterator only provides const references, whether
    /// the underlying resource is persistent (Supercell) or temporary
    /// (other types). Changing database entries must be done via copy then
    /// insert, or update.
    ///
    /// DatabaseIterator should always be dereferenceable (except end or when
    /// default constructed), though the reference may be invalidated when a
    /// second DatabaseIterator is dereferenced, dereferencing the first
    /// again should be valid (though it may require re-allocation):
    ///
    /// \code
    /// DatabaseIterator<ValueType> A = primclex.db<ValueType>().find(A_name);
    /// const ValueType& A_ref1 = *A;  // ok
    /// A_ref1.method();  // ok
    /// A->method(); // ok
    ///
    /// DatabaseIterator<ValueType> B = primclex.db<ValueType>().find(B_name);
    /// const ValueType& B_ref = *B;  // ok, but now A_ref1 may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    ///
    /// const ValueType& A_ref2 = *A;  // ok, but now B_ref may be invalidated
    /// A_ref1.method();  // maybe not ok
    /// A_ref2.method();  // ok
    /// B_ref.method(); // maybe not ok
    /// A->method(); // ok
    /// B->method(); // ok
    /// \endcode

    template<typename ValueType>
    class DatabaseIterator :

      public boost::iterator_facade <
      DatabaseIterator<ValueType>,
      ValueType,
      std::forward_iterator_tag,
      const ValueType &,
      long > {

    public:


      /// Default constructor
      DatabaseIterator() {}

      /// Construct iterator
      DatabaseIterator(const DatabaseIteratorBase<ValueType> &it) :
        m_ptr(notstd::clone(it)) {}

      std::string name() const {
        return m_ptr->name();
      }

      DatabaseIteratorBase<ValueType> *get() const {
        return m_ptr.unique().get();
      }

    private:

      friend boost::iterator_core_access;

      /// boost::iterator_facade implementation
      void increment() {
        m_ptr->increment();
      }

      /// boost::iterator_facade implementation
      const ValueType &dereference() const {
        return m_ptr->dereference();
      }

      /// boost::iterator_facade implementation
      bool equal(const DatabaseIterator &B) const {
        return m_ptr->equal(*(B.m_ptr));
      }

      /*
      long distance_to(const DatabaseIterator &B) const {
        Do not define this or boost::iterator_range<T>::size will compile but
        cause runtime errors. Use boost::distance instead.
      }
      */

      notstd::cloneable_ptr<DatabaseIteratorBase<ValueType> > m_ptr;
    };


    /// \brief Generic interface for database of a particular CASM type
    ///
    /// Currently all ValueType use the same interface, though it is possible to
    /// specialization this interface for a particular ValueType.
    ///
    /// ValueType must have:
    /// - std::string ValueType::name() const
    ///   - is unique, is constant, not necessarily deterministic
    /// - std::string ValueType::alias() const
    ///   - if ValueType::alias().empty(), no alias exists
    ///   - is unique, is not constant
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open() (remember to 'read_aliases')
    /// - void commit() (remember to 'write_aliases')
    /// - void close()
    /// - iterator begin()
    /// - iterator end()
    /// - size_type size() const
    /// - std::pair<iterator, bool> insert(const ValueType &obj)
    /// - std::pair<iterator, bool> insert(const ValueType &&obj)
    /// - iterator erase(iterator pos)
    /// - iterator find(const std::string &name_or_alias)
    ///
    /// Database insert methods by convention do not enforce canonical form or
    /// any other logic, they simply insert as is. By convention, a
    /// ValueType::insert method shold be implemented to enforce canonical form
    /// and any other constraints. But in cases where it is known that ValueType
    /// are being generated in the correct form for insertion in the database,
    /// the Database insert methods may be used directly. For example, see
    /// the method insert_unique_canon_configs.
    ///
    template<typename ValueType>
    class ValDatabase : public DatabaseBase {

    public:

      typedef DatabaseIterator<ValueType> iterator;
      typedef ValueType value_type;
      typedef Index size_type;

      ValDatabase(const PrimClex &_primclex) :
        DatabaseBase(_primclex) {}

      virtual ~ValDatabase() {}


      virtual iterator begin() const = 0;

      virtual iterator end() const = 0;

      virtual size_type size() const = 0;
      bool empty() const {
        return this->size() == 0;
      }

      /// For setting an alias
      std::pair<iterator, bool> set_alias(const std::string &name_or_alias, const std::string &alias);

      /// Get name from name_or_alias
      ///
      /// - Checks if name_or_alias is a known alias
      /// - If known alias, returns associated name
      /// - If not known alias, assumed to be a name and returns name_or_alias
      std::string name(const std::string &name_or_alias) const;

      /// Get alias from name_or_alias
      std::string alias(const std::string &name_or_alias) const;

      virtual std::pair<iterator, bool> insert(const ValueType &obj) = 0;
      virtual std::pair<iterator, bool> insert(const ValueType &&obj) = 0;

      template<typename InputIt>
      void insert(InputIt first, InputIt last) {
        std::for_each(first, last, [&](const ValueType & obj) {
          this->insert(obj);
        });
      }

      virtual iterator erase(iterator pos) = 0;
      size_type erase(const std::string &name_or_alias) {
        this->erase(this->find(name_or_alias));
        return 1;
      }
      size_type erase(const ValueType &obj) {
        return this->erase(obj.name());
      }

      virtual size_type count(const std::string &name_or_alias) const {
        if(this->find(name_or_alias) != this->end()) {
          return 1;
        }
        return 0;
      }

      virtual iterator find(const std::string &name_or_alias) const = 0;
      iterator find(const ValueType &obj) const {
        return find(obj.name());
      }

      virtual void commit() = 0;

    protected:

      void read_aliases();

      void write_aliases();

      /// Only ValDatabase<ValueType> is allowed to do a const id change
      template<typename _ValueType>
      void set_id(const _ValueType &obj, Index id) const {
        obj.set_id(id);
      }

    private:

      // enable lookup of name -> alias
      std::map<std::string, std::string> m_name_to_alias;

      // enable lookup of alias -> name
      std::map<std::string, std::string> m_alias_to_name;

    };

    template<typename ValueType>
    class Database;

  }
}

#endif
