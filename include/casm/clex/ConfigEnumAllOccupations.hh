#ifndef CONFIGENUMALLOCCUPATIONS_HH
#define CONFIGENUMALLOCCUPATIONS_HH

#include "casm/clex/ConfigEnum.hh"
#include "casm/container/Counter.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  template <typename ConfigType>
  class ConfigEnumAllOccupations : public ConfigEnum<ConfigType> {
  public:
    typedef typename ConfigEnum<ConfigType>::step_type step_type;

    // ConfigType is either Configurations or ConfigDoF
    typedef typename ConfigEnum<ConfigType>::value_type value_type;

    typedef typename ConfigEnum<ConfigType>::iterator iterator;

    using ConfigEnum<ConfigType>::initial;
    using ConfigEnum<ConfigType>::final;
    using ConfigEnum<ConfigType>::current;
    using ConfigEnum<ConfigType>::num_steps;
    using ConfigEnum<ConfigType>::step;
  private:
    Counter<Array<int> > m_counter;
    PermuteIterator m_perm_begin, m_perm_end;
    using ConfigEnum<ConfigType>::_current;
    using ConfigEnum<ConfigType>::_step;
    using ConfigEnum<ConfigType>::_source;

    const PermuteIterator &_perm_begin() {
      return m_perm_begin;
    }
    const PermuteIterator &_perm_end() {
      return m_perm_end;
    }
  public:
    ConfigEnumAllOccupations(const value_type &_initial, const value_type &_final, PermuteIterator perm_begin, PermuteIterator perm_end);

    // **** Mutators ****
    // increment m_current and return a reference to it
    const value_type &increment();

    // set m_current to correct value at specified step and return a reference to it
    const value_type &goto_step(step_type _step);

  };

}

#include "casm/clex/ConfigEnumAllOccupations_impl.hh"

#endif
