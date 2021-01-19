#ifndef CASM_DatabaseNameIterator
#define CASM_DatabaseNameIterator

#include <memory>
#include <set>

namespace CASM {

namespace DB {

/// \brief DatabaseIterator for implementations using std::map<std::string,
/// ObjIterator>
///
/// For example, jsonDatabase<Configuration>, stores:
/// - a std::set<Configuration> to efficiently check for unique Configuration
/// - and a std::map<std::string, std::set<Configuration>::iterator> to
///   efficiently find Configuration by name.
/// By iterating over the map, configurations will be iterated in name order
/// rather than Configuration sort order.
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
template <typename ValueType, typename DatabaseType, typename ObjIterator>
class DatabaseNameIterator : public DatabaseIteratorBase<ValueType> {
 public:
  DatabaseNameIterator() {}

  std::string name() const override { return m_it->first; }

  std::unique_ptr<DatabaseIteratorBase<ValueType> > clone() const {
    return std::unique_ptr<DatabaseIteratorBase<ValueType> >(this->_clone());
  }

 private:
  friend DatabaseType;

  typedef
      typename std::map<std::string, ObjIterator>::const_iterator base_iterator;

  DatabaseNameIterator(base_iterator _it) : m_it(_it) {}

  base_iterator base() const { return m_it; }

  bool equal(const DatabaseIteratorBase<ValueType> &other) const override {
    return m_it == static_cast<const DatabaseNameIterator &>(other).m_it;
  }

  void increment() override { ++m_it; }

  const ValueType &dereference() const override { return *(m_it->second); }

  /*
  long distance_to(const DatabaseIteratorBase<ValueType> &other) const override
  { return std::distance(m_it, static_cast<const DatabaseNameIterator
  &>(other).m_it);
  }
  */

  DatabaseNameIterator *_clone() const override {
    return new DatabaseNameIterator(*this);
  }

  base_iterator m_it;
};

}  // namespace DB
}  // namespace CASM

#endif
