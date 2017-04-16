#ifndef CASM_DatabaseSetIterator
#define CASM_DatabaseSetIterator

#include <memory>
#include <set>

namespace CASM {

  namespace DB {

    /// \brief DatabaseIterator for implementations using std::set<ValueType>
    ///
    /// DatabaseIterators must implement public methods:
    /// - Default constructor
    /// - std::string name() const
    /// - std::unique_ptr<DatabaseIteratorBase<ValueType> > clone() const
    ///
    /// DatabaseIterators must implement private methods:
    /// - bool is_end() const
    /// - void increment()
    /// - reference dereference() const
    /// - DatabaseIteratorBase *_clone() const
    ///
    template<typename ValueType, typename DatabaseType>
    class DatabaseSetIterator : public DatabaseIteratorBase<ValueType> {

    public:

      DatabaseSetIterator() {}

      std::string name() const override {
        return m_it->name();
      }

      std::unique_ptr<DatabaseIteratorBase<ValueType> > clone() const {
        return std::unique_ptr<DatabaseIteratorBase<ValueType> >(this->_clone());
      }

    private:

      friend DatabaseType;

      typedef typename std::set<ValueType>::iterator base_iterator;

      DatabaseSetIterator(base_iterator _it) :
        m_it(_it) {}

      base_iterator base() const {
        return m_it;
      }

      bool equal(const DatabaseIteratorBase<ValueType> &other) const override {
        return m_it == static_cast<const DatabaseSetIterator &>(other).m_it;
      }

      void increment() override {
        ++m_it;
      }

      const ValueType &dereference() const override {
        return *m_it;
      }

      long distance_to(const DatabaseIteratorBase<ValueType> &other) const override {
        return std::distance(m_it, static_cast<const DatabaseSetIterator &>(other).m_it);
      }

      DatabaseSetIterator *_clone() const override {
        return new DatabaseSetIterator(*this);
      }

      base_iterator m_it;
    };

  }
}

#endif
