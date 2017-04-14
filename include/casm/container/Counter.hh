#ifndef COUNTER_HH
#define COUNTER_HH

#include "casm/container/BaseCounter.hh"

namespace CASM {

  /** \defgroup Counter
   *
   *  \ingroup Container
   *  \brief Counters allow looping over many incrementing variables in one loop
   *
   *  @{
  */

  /// \brief A Counter allows looping over many incrementing variables in one loop.
  ///
  /// \tparam Container The type of container containing the variables begin looped over.
  ///                   Container should be fully specified, so 'std::vector<int>', not 'std::vector'.
  /// \tparam value_type The type of variable contained by the Container. Must have valid operators '+=' and '+'.
  /// \tparam size_type The type of the index into the Container.
  /// \tparam Access A Functor or function with signature 'value_type& Access::operator()(Container &, size_type);',
  ///                     which provides a reference to a element of the Container.
  /// \tparam Compare A Functor or function with signature 'bool Compare::operator()(const value_type& A, const value_type& B);',
  ///                     which provides a less than comparison, like 'operator<'.
  ///
  /// Default container access using the CASM_TMP::BracketAccess functor is identical to 'container[index]'. Parentheses access,
  /// like 'container(index)', can be obtained by setting the Access parameter using the CASM_TMP::ParenthesesAccess functor.
  ///
  /// Default value_type comparison is as 'value_type < value_type'. Custom comparison can be provided using
  /// the Compare functor.
  ///
  /// Several typedefs are provided:
  /// - for Eigen::VectorXd: \ref EigenVectorXdCounter
  /// - for Eigen::VectorXi: \ref EigenVectorXiCounter
  /// - for Eigen::MatrixXd: \ref EigenMatrixXdCounter
  /// - for Eigen::MatrixXi: \ref EigenMatrixXiCounter
  ///
  /// The first element of the container is the inner loop, and the last element of the container is the outer loop.
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
  template<typename _Container, typename _value_type, typename _size_type, typename _Access, typename _Compare>
  class Counter;

  namespace CASM_TMP {
    template<typename _Container, typename _value_type, typename _size_type, typename _Access, typename _Compare>
    struct traits<Counter<_Container, _value_type, _size_type, _Access, _Compare> > {
      typedef BaseCounter<Counter<_Container, _value_type, _size_type, _Access, _Compare> > Base;
      typedef _Container Container;
      typedef _value_type value_type;
      typedef _size_type size_type;
      typedef _Access Access;
      typedef _Compare Compare;
    };
  }

  template < typename _Container,
             typename _value_type = typename _Container::value_type,
             typename _size_type = typename _Container::size_type,
             typename _Access = CASM_TMP::BracketAccess<_Container, _value_type, _size_type>,
             typename _Compare = CASM_TMP::MuchLessThan<_value_type> >
  class Counter : public BaseCounter<Counter<_Container, _value_type, _size_type, _Access, _Compare> > {
    typedef typename CASM_TMP::traits<Counter>::Base Base;
    using Base::_valid;
    using Base::_current;
    using Base::_increment;
    using Base::_upper;
    using Base::_lower;
    using Base::_initial;
    using Base::_final;
  public:
    typedef _Container Container;
    typedef _value_type value_type;
    typedef _size_type size_type;
    typedef _Access Access;
    typedef _Compare Compare;

    using Base::valid;
    using Base::initial;
    using Base::current;
    using Base::size;
    using Base::compare;
    using Base::operator();
    using Base::operator[];
    using Base::operator bool;

    /// \brief Default construct a Counter
    Counter() {}


    /// \brief Construct a Counter-type object
    ///
    /// \param _initial the initial Container values
    /// \param _final the final valid Container values
    /// \param _increment the amount to increment each Container value by
    Counter(const Container &_initial,
            const Container &_final,
            const Container &_increment,
            Access _access = Access(),
            Compare _compare = Compare()) :
      Base(_initial, _final, _increment, _access, _compare) {
      _init();
    }

    /// Increment the Counter
    ///
    /// \returns true if the Counter is in a valid state after incrementing
    ///
    Counter &operator++() {
      if(!valid())
        return *this;
      size_type i = 0;
      while(compare(_upper(i), _current(i) + _increment(i)) ||
            compare(_current(i) + _increment(i), _lower(i))) {
        _current(i) = _initial(i);
        i++;
        if(i == this->size()) {
          _valid() = false;
          return *this;
        }
      }
      _current(i) += _increment(i);
      return *this;
    }

    void operator++(int) {
      ++(*this);
    }

    void set_current(const Container &new_current) {
      if(current().size() != new_current.size()) {
        std::cerr << "CRITICAL ERROR: In IsoCounter::set_current(), new state is incompatible with this counter.\n"
                  << "                Exiting...\n";
        exit(1);
      }
      _current() = new_current;
      _compute_validity();
    }


    /// Reset the current value of the Counter to the initial value
    void reset() {
      _current() = initial();
      _valid() = size() != 0;
    };

  private:
    /// \brief Called from the constructor to set m_lower and m_upper appropriately
    void _init() {
      for(size_type i = 0; i < initial().size(); i++) {
        if(compare(_initial(i), _final(i))) {
          _lower(i) = _initial(i);
          _upper(i) = _final(i);
        }
        else {
          _lower(i) = _final(i);
          _upper(i) = _initial(i);
        }
      }
    }

    bool _compute_validity() {
      _valid() = true;
      for(size_type i = 0; i < size() && valid(); i++) {
        _valid() = valid()
                   && compare(_current(i), _upper(i))
                   && compare(_lower(i), _current(i));
      }
      return valid();
    }

  };


  template <typename EigenType>
  using EigenCounter = Counter< EigenType, typename EigenType::Scalar, typename EigenType::Index,
        CASM_TMP::ParenthesesAccess< EigenType, typename EigenType::Scalar, typename EigenType::Index> >;

  /// \brief Counter for Eigen::VectorXd
  typedef EigenCounter<Eigen::VectorXd> EigenVectorXdCounter;
  /// \brief Counter for Eigen::MatrixXd
  typedef EigenCounter<Eigen::MatrixXd> EigenMatrixXdCounter;
  /// \brief Counter for Eigen::Matrix3d
  typedef EigenCounter<Eigen::Matrix3d> EigenMatrix3dCounter;
  /// \brief Counter for Eigen::Vector3d
  typedef EigenCounter<Eigen::Vector3d> EigenVector3dCounter;
  /// \brief Counter for Eigen::Vector3i
  typedef EigenCounter<Eigen::Vector3i> EigenVector3iCounter;
  /// \brief Counter for Eigen::VectorXi
  typedef EigenCounter<Eigen::VectorXi> EigenVectorXiCounter;
  /// \brief Counter for Eigen::MatrixXd
  typedef EigenCounter<Eigen::MatrixXi> EigenMatrixXiCounter;
  /// \brief Counter for Eigen::MatrixXd
  typedef EigenCounter<Eigen::Matrix3i> EigenMatrix3iCounter;
  /** @} */

}


#endif
