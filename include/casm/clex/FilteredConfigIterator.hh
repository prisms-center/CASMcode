#ifndef FilteredConfigIterator_HH
#define FilteredConfigIterator_HH

#include <iterator>

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  template <typename IteratorType>
  class FilteredConfigIterator;

  template<typename IteratorType>
  void swap(FilteredConfigIterator<IteratorType> &a, FilteredConfigIterator<IteratorType> &b);

  /// ConfigType bidirectional Iterator class
  ///   Can iterate over all ConfigType in all Supercells of the PrimClex,
  ///   where ConfigType = Configuration, const Configuration, Transition, or const Transition
  ///         PrimClexType = PrimClex or const PrimClex
  ///
  template <typename IteratorType>
  class FilteredConfigIterator : std::iterator<std::input_iterator_tag, typename IteratorType::value_type> {
  public:
    using pointer = typename IteratorType::pointer;
    using reference = typename IteratorType::reference;

    FilteredConfigIterator() {};

    FilteredConfigIterator(const IteratorType &_begin,
                           const IteratorType &_end,
                           const std::string &filter_expr,
                           const DataFormatterDictionary<Configuration> &_dict): m_it(_begin), m_end(_end) {
      m_filter = _dict.parse(filter_expr);
      ValueDataStream<bool> _stream;
      if(m_it != m_end) {
        _stream << m_filter(*m_it);
        if(!_stream.value())
          operator++();
      }
    }

    FilteredConfigIterator(const IteratorType &_begin,
                           const IteratorType &_end,
                           const std::vector<std::string> &filter_expr,
                           const DataFormatterDictionary<Configuration> &_dict): m_it(_begin), m_end(_end) {
      m_filter = _dict.parse(filter_expr);
      ValueDataStream<bool> _stream;
      if(m_it != m_end) {
        _stream << m_filter(*m_it);
        if(!_stream.value())
          operator++();
      }
    }

    FilteredConfigIterator(const IteratorType &_end):
      m_it(_end), m_end(_end), m_filter(ConstantValueFormatter<bool, Configuration>("invalid", false)) {

    }

    //FilteredConfigIterator &operator=(FilteredConfigIterator iter);

    reference operator*() const {
      return *m_it;
    }

    pointer operator->() const {
      return m_it.operator->();
    }

    bool operator==(const FilteredConfigIterator &iter) const {
      return m_it == iter.m_it;
    }

    bool operator!=(const FilteredConfigIterator &iter) const {
      return m_it != iter.m_it;
    }

    FilteredConfigIterator &operator++() {
      ++m_it;
      ValueDataStream<bool> _stream;
      _stream << m_filter(*m_it);
      while(m_it != m_end && !_stream.value()) {
        _stream << m_filter(*(++m_it));
      }
      return *this;
    }

    FilteredConfigIterator operator++(int) {
      FilteredConfigIterator<IteratorType> t_it(*this);
      operator++();
      return t_it;
    }

    friend void swap<>(FilteredConfigIterator &a, FilteredConfigIterator &b);

  private:
    IteratorType m_it, m_end;
    DataFormatter<Configuration> m_filter;
  };


  /// Definitions
  template<typename IteratorType>
  void swap(FilteredConfigIterator<IteratorType> &a, FilteredConfigIterator<IteratorType> &b) {
    using std::swap;

    swap(a.m_it, b.m_it);
    swap(a.m_end, b.m_end);
    swap(a.m_filter, b.m_filter);
  }

  template<typename IteratorType>
  FilteredConfigIterator<IteratorType> filter_begin(
    const IteratorType &it,
    const IteratorType &it_end,
    const std::vector<std::string> &filter_expr,
    const DataFormatterDictionary<Configuration> &_dict) {
    return FilteredConfigIterator<IteratorType>(it, it_end, filter_expr, _dict);
  }

  template<typename IteratorType>
  FilteredConfigIterator<IteratorType> filter_end(const IteratorType &it_end) {
    // Maybe this constructor should be private and make filter_end() a friend of FilteredConfigIterator
    return FilteredConfigIterator<IteratorType>(it_end);
  }


}
#endif
