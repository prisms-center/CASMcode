#include "casm/clex/ConfigEnumRandomLocal.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

namespace CASM {

  ConfigEnumRandomLocalParams::ConfigEnumRandomLocalParams(
    MTRand &_mtrand,
    DoFKey _dof_key,
    Index _n_config,
    double _mag,
    bool _normal_distribution):
    mtrand(_mtrand),
    dof_key(_dof_key),
    n_config(_n_config),
    mag(_mag),
    normal_distribution(_normal_distribution) {}

  ConfigEnumRandomLocal::ConfigEnumRandomLocal(
    ConfigEnumInput const &_in_config,
    ConfigEnumRandomLocalParams const &params):
    ConfigEnumRandomLocal(
      _in_config,
      params.dof_key,
      params.n_config,
      params.mag,
      params.normal_distribution,
      params.mtrand) {}

  ConfigEnumRandomLocal::ConfigEnumRandomLocal(ConfigEnumInput const &_in_config,
                                               DoFKey const &_dof_key,
                                               Index _n_config,
                                               double _mag,
                                               bool _normal,
                                               MTRand &_mtrand):
    m_n_config(_n_config),
    m_mtrand(_mtrand),
    m_mag(_mag),
    m_normal(_normal),
    m_unit_length(DoF::BasicTraits(_dof_key).unit_length()),
    m_site_selection(_in_config.sites().begin(), _in_config.sites().end()) {

    if(m_unit_length)
      m_normal = false;

    if(m_n_config < 0) {
      throw std::runtime_error("Error in ConfigEnumRandomLocal: n_config < 0");
    }
    if(m_n_config == 0) {
      this->_invalidate();
      return;
    }



    m_current = notstd::make_cloneable<Configuration>(_in_config.configuration());

    reset_properties(*m_current);
    this->_initialize(&(*m_current));

    m_dof_vals = &(m_current->configdof().local_dof(_dof_key));

    auto const &dof_info = m_dof_vals->info();
    for(Index l : m_site_selection)
      m_dof_dims.push_back(dof_info[m_current->sublat(l)].dim());

    // Make initial random config
    this->randomize();
    _set_step(0);
    m_current->set_source(this->source(step()));

    //std::cout << "Selection: " << m_site_selection <<"\n"
    //        << "dofkey: " << _dof_key << "\n"
    //        << "mag: " << m_mag << "\n"
    //        << "dof_dims: " << m_dof_dims << "\n"
    //        << "unit_length: " << m_unit_length << "\n"
    //        << "m_normal: " << m_normal << "\n";
  }

  std::string ConfigEnumRandomLocal::name() const {
    return enumerator_name;
  }

  const std::string ConfigEnumRandomLocal::enumerator_name = "ConfigEnumRandomLocal";

  /// Set m_current to correct value at specified step and return a reference to it
  void ConfigEnumRandomLocal::increment() {

    this->_increment_step();
    if(step() < m_n_config) {
      this->randomize();
      m_current->set_source(this->source(step()));
    }
    else {
      this->_invalidate();
    }

  }

  void ConfigEnumRandomLocal::randomize() {
    double tnorm;
    for(Index i = 0; i < m_site_selection.size(); ++i) {
      for(Index j = 0; j < m_dof_dims[i]; ++j) {
        m_dof_vals->site_value(m_site_selection[i])[j] = m_mtrand.randNorm(0., m_mag);
      }
      if(!m_normal) {
        tnorm = m_dof_vals->site_value(m_site_selection[i]).norm();
        if(!m_unit_length) {
          tnorm += -m_mag * std::log(m_mtrand.rand());
        }
        (m_dof_vals->site_value(m_site_selection[i])) /= tnorm;
      }
    }
    //std::cout << "Randomized: \n" << m_dof_vals->values() << "\n";
  }

}
