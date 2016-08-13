#ifndef CONFIGENUMITERATOR_HH
#define CONFIGENUMITERATOR_HH

#include <iterator>

namespace CASM {
  class jsonParser;

  template <typename ConfigType>
  class ConfigEnum;

  template <typename ConfigType>
  class ConfigEnumIterator {
  public:
    typedef typename ConfigEnum<ConfigType>::step_type difference_type;
    typedef typename ConfigEnum<ConfigType>::value_type value_type;
    typedef value_type const *pointer;
    typedef value_type const &reference;
    typedef std::input_iterator_tag iterator_category; //<-- only input_iterator functionality is guaranteed

  private:
    typedef ConfigEnum<ConfigType> enumerator_type;

    // pointer to the enumerator
    mutable enumerator_type *m_enum_ptr;
    difference_type m_step;

    difference_type _step() const {
      return m_step;
    }
    enumerator_type *_enum_ptr() const {
      return m_enum_ptr;
    }
  public:

    ConfigEnumIterator() : m_enum_ptr(nullptr), m_step(-1) {}

    ConfigEnumIterator(enumerator_type &_enumerator, difference_type _step) :
      m_enum_ptr(&_enumerator), m_step(_step) { }


    //ConfigEnumIterator(const ConfigEnumIterator<ConfigType> &B);

    //ConfigEnumIterator &operator=(const ConfigEnumIterator<ConfigType> &B);

    /// \brief Iterator comparison
    bool operator==(const ConfigEnumIterator<ConfigType> &B) const {
      return _step() == B._step() && _enum_ptr() == B._enum_ptr();
    }

    /// \brief Iterator comparison
    bool operator!=(const ConfigEnumIterator<ConfigType> &B) const {
      return !(operator==(B));
    }

    /// \brief Access the current config
    reference operator*() const {
      assert(_enum_ptr() && "ConfigEnumIterator does not point to valid ConfigEnumerator!");
      if(_step() != (_enum_ptr()->step()))
        _enum_ptr()->goto_step(_step());

      return _enum_ptr()->current();
    }

    /// \brief Access the current config
    pointer operator->() const {
      assert(_enum_ptr() && "ConfigEnumIterator does not point to valid ConfigEnumerator!");
      if(_step() != (_enum_ptr()->step()))
        _enum_ptr()->goto_step(_step());

      return &(_enum_ptr()->current());
    }


    /// \brief Prefix increment operator. Increment to next unique configuration.
    ConfigEnumIterator &operator++() {
      assert(_enum_ptr() && "ConfigEnumIterator does not point to valid ConfigEnumerator!");
      if(_step() == _enum_ptr()->step()) { // Case where we can guarantee dereferencing of Iterator
        _enum_ptr()->increment();
        m_step = _enum_ptr()->step();
        //std::cout << "Incrementing iterator inside if!\n";
      }
      else { // Can't guarantee deferencing, but preserve semantics if dereferencing is possible
        m_step++;
      }
      return *this;
    }

    /// \brief Postfix increment operator. Increment to next unique configurtion.
    ConfigEnumIterator operator++(int) {
      ConfigEnumIterator p_it(*this);
      operator++();
      return p_it; // <-- no guarantee that this can be dereferenced
    }

    /// \brief Shift the iterator by specified amount, without guarantee that it can be dereferenced afterwards.
    ConfigEnumIterator &operator+=(difference_type delta) {
      m_step += delta;
      return *this;
    }

    /// \brief Shift the iterator by specified amount, without guarantee that it can be dereferenced afterwards.
    ConfigEnumIterator &operator-=(difference_type delta) {
      m_step -= delta;
      return *this;
    }

    /// \brief Shift the iterator by specified amount, without guarantee that it can be dereferenced afterwards.
    ConfigEnumIterator operator+(difference_type delta) {
      ConfigEnumIterator<ConfigType> t_it(*this);
      return t_it += delta;
    }

    /// \brief Shift the iterator by specified amount, without guarantee that it can be dereferenced afterwards.
    ConfigEnumIterator operator-(difference_type delta) {
      ConfigEnumIterator<ConfigType> t_it(*this);
      return t_it -= delta;
    }

    const jsonParser &source() const {
      assert(_enum_ptr() && "ConfigEnumIterator does not point to valid ConfigEnumerator!");
      if(_step() != (_enum_ptr()->step()))
        _enum_ptr()->goto_step(_step());

      return _enum_ptr()->source();
    }

  };

}

#endif
