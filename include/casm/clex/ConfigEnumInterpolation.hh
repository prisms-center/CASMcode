#ifndef CASM_ConfigEnumInterpolation
#define CASM_ConfigEnumInterpolation

#include "casm/container/RandomAccessEnumerator.hh"
#include "casm/clex/Configuration.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumInterpolation_interface();
}

namespace CASM {

  /// Interpolate displacements and strains between two configurations with
  /// identical occupation
  ///
  /// \ingroup ConfigEnumGroup
  ///
  class ConfigEnumInterpolation : public RandomAccessEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    ConfigEnumInterpolation(const value_type &_initial, const value_type &_final, Index _size);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static const std::string interface_help;
    static int run(PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);

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
