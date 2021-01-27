#ifndef MULTICOUNTER_HH
#define MULTICOUNTER_HH

#include "casm/container/BaseCounter.hh"

namespace CASM {

/** \ingroup Counter

    @{
*/

template <typename SubCounterType>
class MultiCounter;

template <typename SubCounterType>
class CounterValueIterator<MultiCounter<SubCounterType> > {
 public:
  typedef MultiCounter<SubCounterType> CounterType;
  typedef typename CounterType::value_type value_type;
  typedef typename CounterType::size_type size_type;
  typedef typename SubCounterType::const_value_iterator sub_iterator;

  CounterValueIterator(CounterType const *_counter_ptr = nullptr,
                       size_type _ind = 0,
                       sub_iterator _sub_iter = sub_iterator())
      : m_counter_ptr(_counter_ptr), m_ind(_ind), m_sub_iter(_sub_iter) {}
  // CounterValueIterator(const CounterValueIterator &_it);

  bool operator==(const CounterValueIterator &_it) {
    return (m_counter_ptr == _it.m_counter_ptr && m_sub_iter == _it.m_sub_iter);
  }

  bool operator!=(const CounterValueIterator &_it) { return !((*this) == _it); }

  CounterValueIterator &operator++() {
    /* //assume user doesn't need safeguards
    if(m_counter_ptr == nullptr || m_ind < 0 || m_counter_ptr->size() <= m_ind)
    return *this;
    */

    ++m_sub_iter;
    while (m_sub_iter == (*m_counter_ptr)[m_ind].value_end() &&
           (m_ind + 1) < m_counter_ptr->size())
      m_sub_iter = (*m_counter_ptr)[++m_ind].value_begin();

    return *this;
  }

  CounterValueIterator operator++(int) {
    CounterValueIterator t_it(*this);
    ++(*this);
    return t_it;
  }

  const value_type &operator*() const { return *m_sub_iter; }
  value_type const *operator->() const { return &operator*(); }
  /*
  size_type index() {
    return m_ind;
  }
  */
 private:
  CounterType const *m_counter_ptr;
  size_type m_ind;
  sub_iterator m_sub_iter;
};

/// Close relative of the Counter and IsoCounter class.
/// MULTICOUNTER counts over multiple counter-type objects, enumerating all
/// allowed states of all allowed counters It has a template parameter
/// CounterType, which must be a counter class that has the same interface as
/// Counter and IsoCounter
template <typename CounterType>
class MultiCounter {
 protected:
  bool m_is_valid;
  Array<CounterType> m_counters;

 public:
  typedef typename CounterType::value_type value_type;
  typedef typename CounterType::size_type size_type;
  typedef CounterValueIterator<MultiCounter> const_value_iterator;

  MultiCounter();
  MultiCounter(const Array<CounterType> &_counter);
  ~MultiCounter() {}
  MultiCounter &operator++();
  void operator++(int) { ++(*this); }
  bool operator+=(Index steps);
  bool valid() const;

  void reset();

  operator const Array<CounterType> &() const { return m_counters; }

  Index size() const { return m_counters.size(); }

  const_value_iterator value_begin() const {
    if (size() != 0)
      return const_value_iterator(this, 0, m_counters[0].value_begin());
    else
      return const_value_iterator(this);
  }

  const_value_iterator value_end() const {
    if (size() != 0)
      return const_value_iterator(this, size() - 1,
                                  m_counters[size() - 1].value_end());
    else
      return const_value_iterator(this);
  }

  void push_back(const CounterType &new_counter) {
    m_counters.push_back(new_counter);
    reset();
  }

  CounterType &back() { return m_counters.back(); }

  const CounterType &back() const { return m_counters.back(); }

  CounterType &operator[](Index i) { return m_counters[i]; }

  const CounterType &operator[](Index i) const { return m_counters[i]; }

  const Array<CounterType> &current() const { return m_counters; }

  const Array<CounterType> &operator()() const { return m_counters; }
};

template <typename CounterType>
MultiCounter<CounterType>::MultiCounter() {
  m_is_valid = false;
}

template <typename CounterType>
MultiCounter<CounterType>::MultiCounter(const Array<CounterType> &_counters) {
  m_counters = _counters;
  reset();
}

template <typename CounterType>
void MultiCounter<CounterType>::reset() {
  m_is_valid = (size() != 0);
  for (Index i = 0; i < size(); i++) m_counters[i].reset();
}

template <typename CounterType>
bool MultiCounter<CounterType>::valid() const {
  return m_is_valid;
}

template <typename CounterType>
MultiCounter<CounterType> &MultiCounter<CounterType>::operator++() {
  Index index(0);

  while (valid() && index < size()) {
    if (!(++m_counters[index]).valid()) {
      index++;
      continue;
    }

    for (Index i = 0; i < index; i++) {
      m_counters[i].reset();
    }
    break;
  }

  if (index == size()) m_is_valid = false;

  return *this;
}

template <typename CounterType>
bool MultiCounter<CounterType>::operator+=(Index steps) {
  for (Index i = 0; i < steps && m_is_valid; i++) (*this)++;

  return m_is_valid;
}

template <typename CounterType>
std::ostream &operator<<(std::ostream &stream,
                         const MultiCounter<CounterType> &counter) {
  return stream << counter();
}

/** @} */
}  // namespace CASM

#endif
