#ifndef CASM_Database
#define CASM_Database

#include <iterator>
#include <memory>
#include "casm/external/boost.hh"
#include "casm/misc/cloneable_ptr.hh"


namespace CASM {

  class PrimClex;

  namespace DB {

    /// Fully generic database interface for use by DatabaseHandler
    class DatabaseBase {

    public:

      DatabaseBase(const PrimClex &_primclex) :
        m_primclex(&_primclex) {}

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
    /// Derived classes must implement private methods:
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

      std::unique_ptr<DatabaseIteratorBase> clone() const {
        return std::unique_ptr<DatabaseIteratorBase>(this->_clone());
      }

    private:

      virtual void increment() = 0;

      virtual reference dereference() const = 0;

      virtual bool equal(const DatabaseIteratorBase<ValueType> &other) const = 0;

      virtual DatabaseIteratorBase *_clone() const = 0;

      virtual long distance_to(const DatabaseIteratorBase<ValueType> &other) const = 0;

    };

    /// \brief Wrapper class for specializations DatabaseIteratorBase
    ///
    /// - Gives all specialized Database<ValueType> the same iterator type
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

      long distance_to(const DatabaseIterator &B) const {
        return m_ptr->distance_to(*(B.m_ptr));
      }

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
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    /// - iterator begin()
    /// - iterator end()
    /// - size_type size() const
    /// - std::pair<iterator, bool> set_alias(const std::string &name_or_alias, const std::string &alias)
    /// - std::pair<iterator, bool> insert(const ValueType &obj)
    /// - std::pair<iterator, bool> insert(const ValueType &&obj)
    /// - iterator erase(iterator pos)
    /// - iterator find(const std::string &name_or_alias)
    ///
    template<typename ValueType>
    class ValDatabase : public DatabaseBase {

    public:

      typedef DatabaseIterator<ValueType> iterator;
      typedef ValueType value_type;
      typedef Index size_type;


      virtual iterator begin() const = 0;

      virtual iterator end() const = 0;

      virtual size_type size() const = 0;
      bool empty() const {
        return this->size() == 0;
      }

      /// For setting an alias
      virtual std::pair<iterator, bool> set_alias(const std::string &name_or_alias, const std::string &alias) = 0;

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

      virtual size_type count(const std::string &name_or_alias) {
        if(this->find(name_or_alias) != this->end()) {
          return 1;
        }
        return 0;
      }

      virtual iterator find(const std::string &name_or_alias) const = 0;
      virtual iterator find(const ValueType &obj) const = 0;

      virtual void commit() = 0;

    };

    template<typename ValueType>
    class Database;

  }
}

#endif
