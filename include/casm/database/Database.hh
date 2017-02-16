#include "casm/misc/CASM_TMP.hh"

namespace CASM {

  namespace DB {

    /// Fully generic database interface
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

    /// Dereferencing DatabaseIterator only provides const references, whether
    /// the underlying resource is persistent (Supercell, PrimPeriodicDiffTransOrbit)
    /// or temporary (Configuration, DiffTransConfiguration). Changing database
    /// entries must be done via copy, then insert or update.
    ///
    /// Databases do not check canonical forms, etc.
    ///
    /// Derived classes must implement private methods:
    /// - name_type name() const
    /// - bool is_end() const
    /// - void increment()
    /// - reference dereference() const
    /// - DatabaseIteratorBase *_clone() const
    ///
    template<typename ValueType, typename NameType = std::string>
    class DatabaseIteratorBase {

      friend DatabaseIterator<ValueType, NameType>;

    public:

      typedef const value_type &reference;
      typedef NameType name_type;

      DatabaseIteratorBase() {}

      std::unique_ptr<DatabaseIteratorBase> clone() const {
        return std::unique_ptr<DatabaseIteratorBase>(this->_clone());
      }

    private:

      virtual name_type name() const = 0;

      virtual bool is_end() const = 0;

      virtual void increment() = 0;

      virtual reference dereference() const = 0;

      virtual bool equal(const DatabaseIteratorBase<ValueType, NameType> &B) const {
        bool this_is_end = this->is_end();
        bool other_is_end = other.is_end();

        if(this_is_end != other_is_end) {
          return false;
        }

        if(this_is_end) {
          return true;
        }

        return (this->name() == other.name());
      }

      virtual DatabaseIteratorBase *_clone() const = 0;

    };

    template<typename ValueType, typename NameType = std::string>
    class DatabaseIterator :

      public boost::iterator_facade <
      DatabaseIterator<ValueType, NameType>,
      ValueType,
      boost::forward_iterator_tag,
      const ValueType &,
      long > {

    public:

      typedef ValueType value_type;
      typedef const ValueType &reference;
      typedef NameType name_type;


      /// Default constructor
      DatabaseIterator() {}

      /// Construct iterator
      DatabaseIterator(const DatabaseIteratorBase<ValueType, NameType> &it) :
        m_ptr(notstd::clone(it)) {}


      name_type name() const {
        return m_ptr->name();
      }

      bool is_end() const {
        return m_ptr->is_end();
      }

      DatabaseIteratorBase<ValueType, NameType> *get() const {
        return m_ptr.unique().get();
      }

    private:

      friend boost::iterator_core_access;

      /// boost::iterator_facade implementation
      void increment() {
        m_ptr->_increment();
      }

      /// boost::iterator_facade implementation
      reference dereference() const {
        return m_ptr->_dereference();
      }

      /// boost::iterator_facade implementation
      bool equal(const DatabaseIterator &B) const {
        return m_ptr->_equal(*(B.m_ptr));
      }

      notstd::cloneable_ptr<DatabaseIteratorBase<ValueType, NameType> > m_ptr;
    };


    /// ValueType must have:
    /// - NameType ValueType::name() const
    ///
    /// Derived classes must implement public methods:
    /// - void DatabaseBase& open()
    /// - void commit()
    /// - void close()
    /// - iterator begin()
    /// - iterator end()
    /// - size_type size() const
    /// - std::pair<iterator, bool> insert(const ValueType &obj)
    /// - iterator erase(iterator pos)
    /// - iterator find(const name_type &name)
    ///
    template<typename ValueType, typename NameType = std::string>
    class Database {

    public:

      typedef DatabaseIterator<ValueType, NameType> iterator;
      typedef NameType name_type;
      typedef Index size_type;

      Database() {};

      virtual iterator begin() = 0;
      iterator begin() const {
        return static_cast<Database<ValueType>*>(this)->begin();
      }

      virtual iterator end() = 0;
      iterator end() const {
        return static_cast<Database<ValueType>*>(this)->end();
      }

      virtual size_type size() const = 0;
      bool empty() const {
        return this->size() == 0;
      }

      virtual std::pair<iterator, bool> insert(const ValueType &obj) = 0;

      virtual iterator erase(iterator pos) = 0;
      size_type erase(const name_type &name) {
        this->erase(this->find(name));
        return 1;
      }
      size_type erase(const ValueType &obj) {
        return this->erase(obj.name());
      }

      virtual size_type count(const name_type &name) {
        if(this->find(name) != this->end()) {
          return 1;
        }
        return 0;
      }

      virtual iterator find(const name_type &name) = 0;
      virtual const_iterator find(const name_type &name) const {
        return static_cast<>(this)->find(name);
      }

      virtual void commit() = 0;

    };


    /// Derived ScelDatabase must implement public methods:
    /// - iterator find(const Lattice &lat)
    /// - std::pair<iterator, bool> insert(const Lattice &lat)
    ///
    class ScelDatabase : public Database<Supercell> {

    public:

      virtual iterator find(const Lattice &lat);
      virtual std::pair<iterator, bool> insert(const Lattice &lat);

    private:

      std::map<std::string, Supercell> m_scel;

    };


    /// Derived ConfigDatabase must implement public methods:
    /// - std::pair<iterator, bool> rename(const name_type& old_name, const name_type& new_name)
    /// - std::pair<iterator, bool> update(const Configuration &config)
    /// - boost::iterator_range<iterator> scel_range(const name_type& scelname) const
    ///
    class ConfigDatabase : public Database<Configuration> {

    public:

      /// For setting an alias
      virtual std::pair<iterator, bool> set_alias(const name_type &name_or_alias, const name_type &alias);

      /// For updating properties
      virtual std::pair<iterator, bool> update(const Configuration &config);

      /// Range of Configuration in a particular supecell
      virtual boost::iterator_range<iterator> scel_range(const name_type &scelname) const;

    };


    class PrimPeriodicDiffTransOrbitDatabase : public Database<PrimPeriodicDiffTransOrbit> {

    public:

      /// For renaming
      virtual std::pair<iterator, bool> rename(const name_type &old_name, const name_type &new_name);

    };


    class DiffTransConfigDatabase : public Database<DiffTransConfiguration> {

    public:

      /// For renaming
      virtual std::pair<iterator, bool> rename(const name_type &old_name, const name_type &new_name);

      /// For updating properties
      virtual iterator update(const Configuration &config);

      /// Range of DiffTransConfiguration in a particular supecell
      virtual boost::iterator_range<iterator> scel_range(const name_type &scelname) const;

    };
  }
}
