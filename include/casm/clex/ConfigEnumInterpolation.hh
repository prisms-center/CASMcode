#ifndef CONFIGENUMINTERPOLATION_HH
#define CONFIGENUMINTERPOLATION_HH

namespace CASM {

  template <typename ConfigType>
  class ConfigEnumInterpolation : public ConfigEnum<ConfigType> {
  public:
    typedef typename ConfigEnum<ConfigType>::step_type step_type;

    // ConfigType is either Configuration or ConfigDoF
    typedef typename ConfigEnum<ConfigType>::value_type value_type;

    typedef typename ConfigEnum<ConfigType>::iterator iterator;

    using ConfigEnum<ConfigType>::initial;
    using ConfigEnum<ConfigType>::final;
    using ConfigEnum<ConfigType>::current;
    using ConfigEnum<ConfigType>::num_steps;
    using ConfigEnum<ConfigType>::step;
    using ConfigEnum<ConfigType>::source;
  private:
    typename value_type::displacement_matrix_t m_displacement_inc;
    Eigen::Matrix3d m_deformation_inc;

    using ConfigEnum<ConfigType>::_current;
    using ConfigEnum<ConfigType>::_step;
    using ConfigEnum<ConfigType>::_source;

  public:
    ConfigEnumInterpolation(const value_type &_initial, const value_type &_final, Index _N_steps);

    // **** Mutators ****
    // increment m_current and return a reference to it
    const value_type &increment();

    // set m_current to correct value at specified step and return a reference to it
    const value_type &goto_step(step_type _step);

  };

}

#include "casm/clex/ConfigEnumInterpolation_impl.hh"

#endif
