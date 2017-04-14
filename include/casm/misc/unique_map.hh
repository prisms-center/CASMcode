#ifndef CASM_unique_map_HH
#define CASM_unique_map_HH

#include <map>
#include <functional>
#include "casm/external/boost.hh"
#include "casm/misc/CASM_TMP.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace notstd {

  /// \brief An iterator over a UniqueMap
  ///
  /// UniqueType should be Unique or const Unique, and MapIteratorType should be std::map<int, std::unique_ptr<UniqueType> > with
  /// const-ness matching the UniqueType.
  ///
  template<typename TransformFunc, typename MapIteratorType>
  class UniqueMapIterator {

  public:

    typedef std::bidirectional_iterator_tag iterator_category;
    typedef typename MapIteratorType::difference_type difference_type;
    typedef typename std::result_of<TransformFunc(typename MapIteratorType::reference)>::type reference;
    typedef typename std::remove_reference<reference>::type value_type;
    typedef value_type *pointer;

    UniqueMapIterator() {}

    explicit UniqueMapIterator(MapIteratorType map_it) :
      m_it(map_it),
      m_transform(TransformFunc()) {}

    template<typename T, typename M>
    UniqueMapIterator(const UniqueMapIterator<T, M> &umap_it) :
      m_it(umap_it.base()),
      m_transform(TransformFunc()) {}

    template<typename T, typename M>
    bool operator==(const UniqueMapIterator<T, M> &B) const {
      return (!operator!=(B));
    }

    template<typename T, typename M>
    bool operator!=(const UniqueMapIterator<T, M> &B) const {
      return m_it != B.base();
    }

    int type() const {
      return m_it->first;
    }

    reference operator*() const {
      return m_transform(*m_it);
    }

    pointer operator->() const {
      return &(this->operator*());
    }

    // prefix
    UniqueMapIterator &operator++() {
      ++m_it;
      return *this;
    }

    // postfix
    UniqueMapIterator operator++(int) {
      UniqueMapIterator result(*this);
      ++m_it;
      return result;
    }

    // prefix
    UniqueMapIterator &operator--() {
      --m_it;
      return *this;
    }

    // postfix
    UniqueMapIterator operator--(int) {
      UniqueMapIterator result(*this);
      --m_it;
      return result;
    }

    MapIteratorType base() const {
      return m_it;
    }

  private:

    MapIteratorType m_it;
    TransformFunc m_transform;
  };


  template<typename MapType>
  struct GetSecond {

    typedef typename MapType::iterator::reference reference;
    typedef typename MapType::iterator::value_type::second_type &result_type;

    GetSecond() {}

    result_type operator()(reference pair) const {
      return pair.second;
    }

  };

  template<typename MapType>
  struct GetSecondConst {

    typedef typename MapType::const_iterator::reference reference;
    typedef const typename MapType::const_iterator::value_type::second_type &result_type;

    GetSecondConst() {}

    result_type operator()(reference pair) const {
      return pair.second;
    }

  };

  /// \brief std::map wrapper to enforce a 1-1 ValueType->KeyType relationship
  ///
  /// Template Parameters:
  /// - KeyType, type to use for key
  /// - ValueType, type to use as value
  /// - MapType, type of map to use internally.
  ///   - \code MapType::mapped_type(ValueType) \endcode must be a valid constructor.
  /// - TransformFunction and ConstTransformFunction, unary function types, to convert
  ///   *MapType::iterator to ValueType& and *MapType::const_iterator to const ValueType*
  ///
  template<typename KeyType,
           typename ValueType,
           typename MapType = std::map<KeyType, ValueType>,
           typename TransformFunction = GetSecond<MapType>,
           typename ConstTransformFunction = GetSecondConst<MapType> >
  class unique_map {

  public:

    typedef ValueType value_type;
    typedef KeyType key_type;
    typedef std::function<KeyType(const ValueType &)> KeyFuncType;
    typedef std::function<typename MapType::mapped_type(const ValueType &)> ConvertType;
    typedef ValueType &reference;
    typedef ValueType *pointer;

    typedef typename MapType::iterator MapIterator;
    typedef typename MapType::const_iterator ConstMapIterator;

    typedef UniqueMapIterator<TransformFunction, MapIterator> iterator;
    typedef UniqueMapIterator<ConstTransformFunction, ConstMapIterator> const_iterator;

    typedef typename MapType::size_type size_type;

    explicit unique_map(KeyFuncType keyfunc, ConvertType _converter) :
      m_keyfunc(keyfunc), m_converter(_converter) {}


    key_type key(const value_type &value) const {
      return m_keyfunc(value);
    }

    /// \brief Insert single value
    std::pair<iterator, bool> insert(const value_type &value) {
      return _insert(value);
    }

    /// \brief Insert single value
    iterator insert(const_iterator hint, const value_type &value) {
      return iterator(
               m_map.insert(
                 hint.base(),
                 std::make_pair(key(value), std::move(m_converter(value)))
               )
             );
    }

    /// \brief Insert single value
    iterator insert(iterator hint, const value_type &value) {
      return iterator(
               m_map.insert(
                 hint.base(),
                 std::make_pair(key(value), std::move(m_converter(value)))
               )
             );
    }

    /// \brief Variadic insert accepts as const UniqueMap& or const ValueType&
    template<typename... ValuesOrMaps>
    void insert(const ValuesOrMaps &... more) {
      CASM::CASM_TMP::ignore_returnvalues(_insert(more)...);
    }

    /// \brief Iterator range insert
    template<typename Iterator>
    void insert(Iterator begin, Iterator end, typename CASM::CASM_TMP::enable_if_iterator<Iterator>::type * = nullptr) {
      for(auto it = begin; it != end; ++it) {
        _insert(*it);
      }
    }

    value_type &operator[](const key_type &key) {
      // insert if not existing
      m_map[key];

      // return reference
      return *find(key);
    }

    iterator find(const key_type &key) {
      return iterator(m_map.find(key));
    }

    const_iterator find(const key_type &key) const {
      return const_iterator(m_map.find(key));
    }

    void clear() {
      m_map.clear();
    }

    iterator erase(const_iterator pos) {
      return iterator(m_map.erase(pos.base()));
    }

    iterator erase(const_iterator first, const_iterator last) {
      return iterator(m_map.erase(first.base(), last.base()));
    }

    size_type erase(const key_type &key) {
      return m_map.erase(key);
    }

    size_type size() const {
      return m_map.size();
    }

    bool empty() const {
      return m_map.empty();
    }

    iterator begin() {
      return iterator(m_map.begin());
    }

    const_iterator begin() const {
      return const_iterator(m_map.begin());
    }

    const_iterator cbegin() const {
      return const_iterator(m_map.cbegin());
    }

    iterator end() {
      return iterator(m_map.end());
    }

    const_iterator end() const {
      return const_iterator(m_map.end());
    }

    const_iterator cend() const {
      return const_iterator(m_map.cend());
    }


  private:

    /// \brief Copy insert
    std::pair<iterator, bool> _insert(const value_type &value) {
      auto result = m_map.insert(std::make_pair(key(value), std::move(m_converter(value))));
      return std::make_pair(iterator(result.first), result.second);
    }

    /// \brief Dictionary insert
    template<typename A, typename B, typename C, typename D, typename E>
    int _insert(const unique_map<A, B, C, D, E> &value) {
      insert(value.begin(), value.end());
      return 0;
    }

    MapType m_map;
    KeyFuncType m_keyfunc;
    ConvertType m_converter;
  };

}

#endif

