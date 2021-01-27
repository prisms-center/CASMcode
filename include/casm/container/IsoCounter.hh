#ifndef ISOCOUNTER_HH
#define ISOCOUNTER_HH

#include "casm/container/BaseCounter.hh"

namespace CASM {

/** \ingroup Counter

    @{
*/

/// \brief A IsoCounter allows looping over many incrementing variables in one
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
/// \tparam Compare A Functor or function with signature 'bool
/// Compare::operator()(const value_type& A, const value_type& B);',
///                     which provides a less than comparison, like 'operator<'.
///
/// Default container access using the CASM_TMP::BracketAccess functor is
/// identical to 'container[index]'. Parentheses access, like
/// 'container(index)', can be obtained by setting the Access parameter using
/// the CASM_TMP::ParenthesesAccess functor.
///
/// Default value_type comparison is as 'value_type < value_type'. Custom
/// comparison can be provided using the Compare functor.
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
template <typename _Container, typename _value_type, typename _size_type,
          typename _Access, typename _Compare>
class IsoCounter;

template <typename _Container, typename _value_type, typename _size_type,
          typename _Access, typename _Compare>
struct traits<
    IsoCounter<_Container, _value_type, _size_type, _Access, _Compare> > {
  typedef BaseCounter<
      IsoCounter<_Container, _value_type, _size_type, _Access, _Compare> >
      Base;
  typedef _Container Container;
  typedef _value_type value_type;
  typedef _size_type size_type;
  typedef _Access Access;
  typedef _Compare Compare;
};

template <typename _Container,
          typename _value_type = typename _Container::value_type,
          typename _size_type = typename _Container::size_type,
          typename _Access =
              CASM_TMP::BracketAccess<_Container, _value_type, _size_type>,
          typename _Compare = CASM_TMP::MuchLessThan<_value_type> >
class IsoCounter
    : public BaseCounter<
          IsoCounter<_Container, _value_type, _size_type, _Access, _Compare> > {
  typedef typename traits<IsoCounter>::Base Base;
  using Base::_current;
  using Base::_final;
  using Base::_increment;
  using Base::_initial;
  using Base::_lower;
  using Base::_upper;
  using Base::_valid;

 public:
  typedef _Container Container;
  typedef _value_type value_type;
  typedef _size_type size_type;
  typedef _Access Access;
  typedef _Compare Compare;

  using Base::compare;
  using Base::current;
  using Base::final;
  using Base::initial;
  using Base::size;
  using Base::valid;
  using Base::operator();
  using Base::operator[];
  using Base::operator bool;

  /// \brief Default construct a IsoCounter
  /// Container must be default constructible
  /// value_type must be constructible with zero
  IsoCounter() : IsoCounter(Container(), Container(), 0, 0) {}

  /// \brief Construct a IsoCounter-type object
  ///
  /// \param _initial the initial Container values
  /// \param _final the final valid Container values
  /// \param _increment the amount to increment each Container value by
  IsoCounter(const Container &_initial, const Container &_final,
             const value_type &_increment, const value_type &_sum,
             Access _access = Access(), Compare _compare = Compare())
      : Base(_initial, _final, Array<value_type>(1, _increment), _access,
             _compare),
        m_sum_constraint(_sum) {
    _init();
    reset();
  }

  const value_type &increment() const { return _increment(0); }

  value_type current_sum() const {
    if (size() == 0) return 0;
    value_type result(current(0));
    for (size_type i = 1; i < size(); i++) {
      result += current(i);
    }
    return result;
  }

  /// Increment the IsoCounter
  ///
  /// \returns true if the IsoCounter is in a valid state after incrementing
  ///
  IsoCounter &operator++() {
    if (!valid()) return (*this);
    size_type i(0), j(0);
    bool continue_flag = valid();
    value_type next_i = 0, next_j = 0;
    for (i = 0; i + 1 < size(); i++) {
      next_i = current(i) - increment();
      if (compare(_upper(i), next_i) || compare(next_i, _lower(i))) continue;
      j = i + 1;
      while (j < size() && !compare(_lower(j), _upper(j))) j++;
      if (j == size()) break;
      next_j = current(j) + increment();
      if (compare(_upper(j), next_j) || compare(next_j, _lower(j))) {
        i = j - 1;
        continue;
      }

      _current(i) = next_i;
      _current(j) = next_j;
      // next_i -= m_increment;

      continue_flag = false;
      break;
    }

    // if the loop finished without breaking, there were no more elements to
    // increment. that means we have reached the last permutation, so we
    // invalidate counter and return
    if (continue_flag || !valid()) {
      _valid() = false;
      return *this;
    }

    j = 0;
    next_j = current(j) + increment();
    // next_i already contains _current(i) from above
    next_i -= increment();
    while (j < i) {
      while (!(compare(_upper(j), next_j) || compare(next_j, _lower(j)) ||
               compare(_upper(i), next_i) || compare(next_i, _lower(i)))) {
        _current(i) = next_i;
        next_i -= increment();

        _current(j) = next_j;
        next_j += increment();
      }

      while (j < i &&
             (compare(_upper(i), next_i) || compare(next_i, _lower(i))))
        next_i = current(--i) - increment();

      while (j < i &&
             (compare(_upper(j), next_j) || compare(next_j, _lower(j))))
        next_j = current(++j) + increment();
    }
    return *this;
  }

  void operator++(int) { ++(*this); }

  void set_sum_constraint(const value_type &new_sum) {
    _sum_constraint() = new_sum;
    reset();
  }

  void set_current(const Container &new_current) {
    if (current().size() != new_current.size()) {
      std::cerr << "CRITICAL ERROR: In IsoCounter::set_current(), new state is "
                   "incompatible with this counter.\n"
                << "                Exiting...\n";
      exit(1);
    }
    _current() = new_current;
    _compute_validity();
  }

  /// Reset thse current value of the Counter to the initial value
  void reset() {
    _valid() = size() != 0;

    _current() = initial();

    // Find the residual sum: m_sum_constraint - m_current.sum()
    value_type sum_val(_sum_constraint() - current_sum());

    // if residual sum does not have same sign as increment, nothing can happen
    if (compare(sum_val + sum_val, sum_val) !=
        compare(increment() + increment(), increment())) {
      _valid() = false;
      return;
    }

    size_type i;

    for (i = 0; i < size(); i++) {
      // if(almost_equal(m_initial[i], m_final[i])) //avoid operating on fixed
      // elements continue;

      // if the i'th element isn't fixed, then increment it while decrementing
      // sum_val until one or the other cannot be incremented/decremented
      // further
      while (!(compare(_upper(i), _current(i) + increment()) ||
               compare(_current(i) + increment(), _lower(i)) ||
               almost_zero(sum_val))) {
        _current(i) += increment();
        sum_val -= increment();
      }
    }
    if (!almost_zero(sum_val)) {
      _valid() = false;
    }
    return;
  }

 private:
  value_type m_sum_constraint;

  value_type &_sum_constraint() { return m_sum_constraint; }

  const value_type &_sum_constraint() const { return m_sum_constraint; }

  bool _compute_validity() {
    value_type S(current_sum());
    _valid() =
        !(compare(_sum_constraint(), S) || compare(S, _sum_constraint()));
    for (size_type i = 0; i < size() && valid(); i++) {
      _valid() = valid() && compare(_current(i), _upper(i)) &&
                 compare(_lower(i), _current(i));
    }
    return valid();
  }

  /// \brief Called from the constructor to set m_lower and m_upper
  /// appropriately
  void _init() {
    for (size_type i = 0; i < initial().size(); i++) {
      if (compare(_initial(i), _final(i))) {
        _lower(i) = _initial(i);
        _upper(i) = _final(i);
      } else {
        _lower(i) = _final(i);
        _upper(i) = _initial(i);
      }
    }
  }

  value_type &_increment() { return Base::_increment(0); }

  const value_type &_increment() const { return Base::_increment(0); }
};

/** @} */

}  // namespace CASM

#endif
