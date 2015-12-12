#include "casm/clex/ConfigIterator.hh"

#include "casm/clex/PrimClex.hh"

namespace CASM {

  /// Specialize for Configuration, const Configuration, Transition, const Transition
  /*
  template <>
  Configuration &ConfigIterator<Configuration, PrimClex>::operator*() {
    return m_primclex->get_supercell(m_scel_index).get_config(m_config_index);
  }

  template <>
  const Configuration &ConfigIterator<const Configuration, const PrimClex>::operator*() {
    return m_primclex->get_supercell(m_scel_index).get_config(m_config_index);
  }
  */

  /*
  template <>
  Transition &ConfigIterator<Transition, PrimClex>::operator*() {
    return m_primclex->get_supercell(m_scel_index).get_transition(m_config_index);
  }

  template <>
  const Transition &ConfigIterator<const Transition, const PrimClex>::operator*() {
    return m_primclex->get_supercell(m_scel_index).get_transition(m_config_index);
  }
  */


  /// Specialize for Configuration, const Configuration, Transition, const Transition
  template<>
  void ConfigIterator<Configuration, PrimClex>::set_selected(bool _select) {
    m_primclex->get_supercell(m_scel_index).get_config(m_config_index).set_selected(_select);
  }

  template<>
  void ConfigIterator<const Configuration, const PrimClex>::set_selected(bool _select) {
    throw std::runtime_error("Attempting to alter const Configuration from ConfigIterator::set_selected. This is not possible!\n");
  }

  /// Specialize for Configuration, const Configuration, Transition, const Transition
  template<>
  int ConfigIterator<Configuration, PrimClex>::config_list_size() const {
    return m_primclex->get_supercell(m_scel_index).get_config_list().size();
  }

  template<>
  int ConfigIterator<const Configuration, const PrimClex>::config_list_size() const {
    return m_primclex->get_supercell(m_scel_index).get_config_list().size();
  }

  /*
  template<>
  int ConfigIterator<Transition, PrimClex>::config_list_size() const {
    return m_primclex()->get_supercell(m_scel_index).transition_list().size();
  }

  template<>
  int ConfigIterator<const Transition, const PrimClex>::config_list_size() const {
    return m_primclex()->get_supercell(m_scel_index).transition_list().size();
  }
  */

}

