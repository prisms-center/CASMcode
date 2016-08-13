
namespace CASM {

  template <typename ConfigType>
  ConfigEnumInterpolation<ConfigType>::ConfigEnumInterpolation(const value_type &_initial, const value_type &_final, Index _N_steps) :
    ConfigEnum<ConfigType>(_initial, _final, _N_steps) {

    // interpolated configurations must be compatible in the following ways:
    if(initial().size() != final().size()
       || initial().occupation() != final().occupation()
       || initial().has_displacement() != final().has_displacement()
       || initial().is_strained() != final().is_strained()) {
      std::cerr << "CRITICAL ERROR: In ConfigEnumInterpolation constructor, attempting to interpolate between incompatible Configurations or ConfigDoFs.\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    if(initial().has_displacement())
      m_displacement_inc = (final().displacement() - initial().displacement()) / ((double)num_steps());

    if(initial().is_strained())
      m_deformation_inc = (final().deformation() - initial().deformation()) / ((double)num_steps());
  }

  //*******************************************************************************************

  // increment m_current and return a reference to it
  template <typename ConfigType>
  const typename ConfigEnumInterpolation<ConfigType>::value_type &ConfigEnumInterpolation<ConfigType>::increment() {
    if(current().has_displacement())
      _current().set_displacement(current().displacement() + m_displacement_inc);

    if(current().is_strained())
      _current().set_deformation(current().deformation() + m_deformation_inc);

    ++_step();
    _current().set_source(source());
    return current();
  }

  //*******************************************************************************************

  // set m_current to correct value at specified step and return a reference to it
  template <typename ConfigType>
  const typename ConfigEnumInterpolation<ConfigType>::value_type &ConfigEnumInterpolation<ConfigType>::goto_step(step_type new_step) {
    if(new_step == step())
      return current();

    if(current().has_displacement())
      _current().set_displacement(initial().displacement() + ((double)new_step)*m_displacement_inc);

    if(current().is_strained())
      _current().set_deformation(initial().deformation() + ((double)new_step)*m_deformation_inc);

    _step() = new_step;
    _current().set_source(source());
    return current();
  }
}

