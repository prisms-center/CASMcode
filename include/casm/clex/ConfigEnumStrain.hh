#ifndef CONFIGENUMSTRAIN_HH
#define CONFIGENUMSTRAIN_HH

#include "casm/symmetry/PermuteIterator.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/clex/ConfigEnum.hh"

namespace CASM {

  class Supercell;

  template <typename ConfigType>
  class ConfigEnumStrain : public ConfigEnum<ConfigType> {
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
    using ConfigEnum<ConfigType>::source;

    ConfigEnumStrain(Supercell &scel, const value_type &_init, const std::vector<Index> &subspace_partitions, const std::vector<double> &magnitudes, std::string _mode);

    // **** Mutators ****
    // increment m_current and return a reference to it
    const value_type &increment();

    // set m_current to correct value at specified step and return a reference to it
    const value_type &goto_step(step_type _step);

  private:
    EigenCounter<Eigen::VectorXd> m_counter;
    StrainConverter m_strain_calc;
    Eigen::MatrixXd m_proj;
    PermuteIterator m_perm_begin, m_perm_end;
    Eigen::MatrixXd m_shape_factor;
    using ConfigEnum<ConfigType>::_current;
    using ConfigEnum<ConfigType>::_step;
    using ConfigEnum<ConfigType>::_source;

    const PermuteIterator &_perm_begin() {
      return m_perm_begin;
    }
    const PermuteIterator &_perm_end() {
      return m_perm_end;
    }
  };

}

#include "casm/clex/ConfigEnumStrain_impl.hh"

#endif
