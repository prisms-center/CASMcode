#include "casm/clex/ConfigEnumInterpolation.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumInterpolation_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumInterpolation>();
  }
}

namespace CASM {

  const std::string ConfigEnumInterpolation::enumerator_name = "ConfigEnumInterpolation";

  const std::string ConfigEnumInterpolation::interface_help =
    "ConfigEnumInterpolation: \n\n"

    "  ... include help documentation here ... \n\n";

  int ConfigEnumInterpolation::run(
    PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {
    throw std::runtime_error("ConfigEnumAllOccupations::run is not implemented");
  }

  /// \brief Constructor
  ///
  /// \param _initial,_final Initial and final configurations to interpolate between
  /// \param _size The total number of configurations to enumerate, including
  ///              the initial and final configurations
  ///
  /// - The `final` configuration is *not* pointed at by the end iterator,
  ///   which points past-the-final element, as is typical
  /// - `_size` will be equal to \code std::distance(this->begin(), this->end()) \endcode
  ConfigEnumInterpolation::ConfigEnumInterpolation(
    const value_type &_initial,
    const value_type &_final,
    Index _size):
    RandomAccessEnumeratorBase<Configuration>(_size),
    m_current(_initial),
    m_initial(_initial),
    m_final(_final) {

    reset_properties(m_current);
    reset_properties(m_initial);
    reset_properties(m_final);

    // interpolated configurations must be compatible in the following ways:
    if(m_initial.size() != m_final.size()
       || m_initial.occupation() != m_final.occupation()
       || m_initial.has_displacement() != m_final.has_displacement()
       || m_initial.has_deformation() != m_final.has_deformation()) {

      std::string s = "Error: " + this->name() + " attempting to"
                      "interpolate between incompatible Configurations";

      throw std::runtime_error(s);
    }

    if(m_initial.has_displacement()) {
      m_displacement_inc = (m_final.displacement() - m_initial.displacement()) / ((double) size() - 1);
    }

    if(m_initial.has_deformation()) {
      m_deformation_inc = (m_final.deformation() - m_initial.deformation()) / ((double) size() - 1);
    }

    this->_initialize(&m_current);
    _current().set_source(this->source(step()));
  }

  /// Set m_current to correct value at specified step and return a reference to it
  Configuration *ConfigEnumInterpolation::at_step(step_type n) {
    if(current().has_displacement())
      _current().set_displacement(m_initial.displacement() + ((double) n)*m_displacement_inc);

    if(current().has_deformation())
      _current().set_deformation(m_initial.deformation() + ((double) n)*m_deformation_inc);

    _current().set_source(this->source(step()));
    return &_current();
  }

}

