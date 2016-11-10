#ifndef CASM_ConfigEnumInterpolation
#define CASM_ConfigEnumInterpolation

#include "casm/container/RandomAccessEnumerator.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  ENUMERATOR_INTERFACE_TRAITS(ConfigEnumInterpolation)

  /// Interpolate displacements and strains between two configurations with
  /// identical occupation
  ///
  /// \ingroup ConfigEnum
  ///
  class ConfigEnumInterpolation : public RandomAccessEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumInterpolation(const value_type &_initial, const value_type &_final, Index _size);

    ENUMERATOR_MEMBERS(ConfigEnumInterpolation)

  private:

    /// Implements goto_step
    Configuration *at_step(step_type n) override;


    // -- Unique members -------------------

    Configuration m_current;
    Configuration m_initial;
    Configuration m_final;

    typename value_type::displacement_matrix_t m_displacement_inc;
    Eigen::Matrix3d m_deformation_inc;

  };

}

#endif
