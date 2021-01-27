#ifndef BASECOUNTER_HH
#define BASECOUNTER_HH

#include <iterator>

#include "casm/misc/CASM_TMP.hh"

namespace CASM {

/** \defgroup Counter

    @{
*/

template <typename CounterType>
class CounterValueIterator {
 public:
  typedef typename CounterType::value_type value_type;
  typedef typename CounterType::size_type difference_type;
  typedef typename CounterType::size_type size_type;
  typedef std::random_access_iterator_tag iterator_category;
  typedef const value_type &reference;
  typedef const value_type *pointer;

  CounterValueIterator(CounterType const *_counter_ptr = nullptr,
                       size_type _ind = 0)
      : m_counter_ptr(_counter_ptr), m_ind(_ind) {}

  // explicit CounterValueIterator(const CounterValueIterator &_it) :
  //  m_counter_ptr(_it.m_counter_ptr),
  //  m_ind(_it.m_ind) {}

  bool operator==(const CounterValueIterator &_it) {
    return (m_counter_ptr == _it.m_counter_ptr && m_ind == _it.m_ind);
  }

  bool operator!=(const CounterValueIterator &_it) { return !((*this) == _it); }

  CounterValueIterator &operator++() {
    ++m_ind;
    return *this;
  }

  CounterValueIterator operator++(int) {
    CounterValueIterator t_it(*this);
    ++(*this);
    return t_it;
  }

  reference operator*() const { return (*m_counter_ptr)[m_ind]; }
  pointer operator->() const { return &operator*(); }

 private:
  CounterType const *m_counter_ptr;
  size_type m_ind;
};

/// \brief A Counter allows looping over many incrementing variables in one
/// loop.
///
/// \tparam Container The type of container containing the variables begin
/// looped over.
///                   Container should be fully specified, so
///                   'std::vector<int>', not 'std::vector'.
/// \tparam value_type The type of variable contained by the Container. Must
/// have valid operators '+=' and '+'. \tparam size_type The type of the index
/// into the Container. \tparam Access A Functor or function with signature
/// 'value_type& Access::operator()(Container &, size_type);',
///                     which provides a reference to a element of the
///                     Container.
///
/// Default container access using the CASM_TMP::BracketAccess functor is
/// identical to 'container[index]'. Parentheses access, like
/// 'container(index)', can be obtained by setting the Access parameter using
/// the CASM_TMP::ParenthesesAccess functor.
///
/// Default value_type comparison is as 'value_type < value_type'. Custom
/// comparison can be provided using the Compare functor.
///
/// Several typedefs are provided:
/// - for Eigen::VectorXd: \ref EigenVectorXdCounter
/// - for Eigen::VectorXi: \ref EigenVectorXiCounter
/// - for Eigen::MatrixXd: \ref EigenMatrixXdCounter
/// - for Eigen::MatrixXi: \ref EigenMatrixXiCounter
///
/// The first element of the container is the inner loop, and the last element
/// of the container is the outer loop.
///
/// \b Example:
/// \code
///   int size = 3;
///   std::vector<int> initial(size, 0);
///   std::vector<int> final(size, 1);
///   std::vector<int> increment(size, 1);
///
///   casm::Counter<std::vector<int> > counter(initial, final, increment);
///
///   do {
///     for( int i=0; i<counter.size(); i++) {
///       std::cout << counter()[i] << " ";
///     }
///     std::cout << std::endl;
///   } while(counter++);
/// \endcode
///
/// \b Output:
/// \code
///   0 0 0
///   1 0 0
///   0 1 0
///   1 1 0
///   0 0 1
///   1 0 1
///   0 1 1
///   1 1 1
/// \endcode
///
/*template<typename Container,
         typename value_type = typename Container::value_type,
         typename size_type = typename Container::size_type,
         typename Access = CASM_TMP::BracketAccess<Container, value_type,
   size_type>, typename Compare = CASM_TMP::MuchLessThan<value_type> >
*/
template <typename DerivedCounter>
class BaseCounter {
 public:
  typedef typename traits<DerivedCounter>::Container Container;
  typedef typename traits<DerivedCounter>::value_type value_type;
  typedef typename traits<DerivedCounter>::size_type size_type;
  typedef typename traits<DerivedCounter>::Access Access;
  typedef typename traits<DerivedCounter>::Compare Compare;
  typedef CounterValueIterator<DerivedCounter> const_value_iterator;

  /// \brief Default construct a Counter
  BaseCounter() {}

  /// \brief Construct a Counter-type object
  ///
  /// \param _initial the initial Container values
  /// \param _final the final valid Container values
  /// \param _increment the amount to increment each Container value by
  BaseCounter(const Container &_initial, const Container &_final,
              const Container &_increment, Access _access = Access(),
              Compare _compare = Compare())
      : m_access(_access),
        m_compare(_compare),
        m_initial(_initial),
        m_final(_final),
        m_increment(_increment),
        m_lower(_initial),
        m_upper(_initial),
        m_current(_initial),
        m_valid(_initial.size() != 0) {
    //_init();
  }

  /// \returns false if the Counter::current value has been incremented past the
  /// Counter::final
  ///
  bool valid() const { return m_valid; };

  operator bool() const { return valid(); }

  /// const Access the current value of the Container (identical to current)
  ///
  /// \returns a const reference to the current value of the Container
  ///
  operator const Container &() const { return m_current; }

  const value_type &operator[](size_type index) const {
    return m_access(m_current, index);
  };

  /// Return the size of the Container
  ///
  /// \returns the size of the Container
  ///
  size_type size() const { return m_current.size(); }

  const_value_iterator value_begin() const {
    return const_value_iterator(static_cast<DerivedCounter const *>(this), 0);
  }

  const_value_iterator value_end() const {
    return const_value_iterator(static_cast<DerivedCounter const *>(this),
                                size());
  }

  /// Use internal Compare object to compare two values
  ///
  /// \returns bool that is result of Compare()(value A, value B)
  ///
  bool compare(const value_type &A, const value_type &B) {
    return m_compare(A, B);
  }

  /// const Access the current value of the Container
  ///
  /// \returns a const reference to the current value of the Container
  ///
  const Container &current() const { return m_current; }

  /// const Access the element 'index' of the current value of the Container
  ///
  /// \returns a const reference to an element of the Container
  ///
  const value_type &current(size_type index) const {
    return m_access(m_current, index);
  }

  /// const Access the current value of the Container (identical to current)
  ///
  /// \returns a const reference to the current value of the Container
  ///
  const Container &operator()() const { return m_current; }

  /// const Access the intial value of the Container
  ///
  /// \returns a const reference to the initial value of the Container
  ///
  const Container &initial() const { return m_initial; }

  /// const Access the final value of the Container
  ///
  /// \returns a const reference to the final value of the Container
  ///
  const Container & final() const { return m_final; }

  /// const Access the incrementing values of the Container
  ///
  /// \returns a const reference to the incrementing value of the Container
  ///
  const Container &increment() const { return m_increment; }

 protected:
  /// non-const reference bool validity flag
  bool &_valid() { return m_valid; }

  /// const Access to element in the current value of the Container
  ///
  /// \returns a const reference to the element 'index' in the current value of
  /// the Container
  ///
  const value_type &_current(size_type index) const {
    return m_access(m_current, index);
  }

  /// non-const Access to element in the current value of the Container
  ///
  /// \returns a non-const reference to the element 'index' in the current value
  /// of the Container
  ///
  value_type &_current(size_type index) { return m_access(m_current, index); }

  /// non-const Access the current value of the Container
  ///
  /// \returns a non-const reference to the current value of the Container
  ///
  Container &_current() { return m_current; }

  /// const Access to element in the initial value of the Container
  ///
  /// \returns a const reference to the element 'index' in the initial value of
  /// the Container
  ///
  const value_type &_initial(size_type index) const {
    return m_access(m_initial, index);
  }

  /// non-const Access to element in the initial value of the Container
  ///
  /// \returns a non-const reference to the element 'index' in the initial value
  /// of the Container
  ///
  value_type &_initial(size_type index) { return m_access(m_initial, index); }

  /// const Access to element in the final value of the Container
  ///
  /// \returns a const reference to the element 'index' in the final value of
  /// the Container
  ///
  const value_type &_final(size_type index) const {
    return m_access(m_final, index);
  }

  /// non-const Access to element in the final value of the Container
  ///
  /// \returns a non-const reference to the element 'index' in the final value
  /// of the Container
  ///
  value_type &_final(size_type index) { return m_access(m_final, index); }

  /// const Access to element in the increment value of the Container
  ///
  /// \returns a const reference to the element 'index' in the increment value
  /// of the Container
  ///
  const value_type &_increment(size_type index) const {
    return m_access(m_increment, index);
  }

  /// non-const Access to element in the increment value of the Container
  ///
  /// \returns a non-const reference to the element 'index' in the increment
  /// value of the Container
  ///
  value_type &_increment(size_type index) {
    return m_access(m_increment, index);
  }

  /// const Access to element in the upper value of the Container
  ///
  /// \returns a const reference to the element 'index' in the upper value of
  /// the Container
  ///
  const value_type &_upper(size_type index) const {
    return m_access(m_upper, index);
  }

  /// non-const Access to element in the upper value of the Container
  ///
  /// \returns a non-const reference to the element 'index' in the upper value
  /// of the Container
  ///
  value_type &_upper(size_type index) { return m_access(m_upper, index); }

  /// const Access to element in the lower value of the Container
  ///
  /// \returns a const reference to the element 'index' in the lower value of
  /// the Container
  ///
  const value_type &_lower(size_type index) const {
    return m_access(m_lower, index);
  }

  /// non-const Access to element in the lower value of the Container
  ///
  /// \returns a non-const reference to the element 'index' in the lower value
  /// of the Container
  ///
  value_type &_lower(size_type index) { return m_access(m_lower, index); }

 private:
  /// A functor that enables different ways of accessing the container, by
  /// default: container[index]
  Access m_access;

  /// A functor that enables custom comparison, by default: operator<
  Compare m_compare;

  /// Initial container values
  Container m_initial;

  /// Final container values
  Container m_final;

  /// Amount to increment each container value
  Container m_increment;

  /// The minimum of each value of m_initial and m_final
  Container m_lower;

  /// The maximum of each value of m_initial and m_final
  Container m_upper;

  /// The current state of the container
  Container m_current;

  /// True if m_current is within the allowed bounds,
  /// false once all values have been incremented to their limit
  bool m_valid;
};

/** @} */

}  // namespace CASM

#endif
