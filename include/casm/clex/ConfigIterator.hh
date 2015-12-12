#ifndef ConfigIterator_HH
#define ConfigIterator_HH

#include <iterator>

#include "casm/CASM_global_definitions.hh"

namespace CASM {

  template <typename ConfigType, typename PrimClexType>
  class ConfigIterator;

  template<typename ConfigType, typename PrimClexType>
  void swap(ConfigIterator<ConfigType, PrimClexType> &a, ConfigIterator<ConfigType, PrimClexType> &b);

  /// ConfigType bidirectional Iterator class
  ///   Can iterate over all ConfigType in all Supercells of the PrimClex,
  ///   where ConfigType = Configuration, const Configuration, Transition, or const Transition
  ///         PrimClexType = PrimClex or const PrimClex
  ///
  template <typename ConfigType, typename PrimClexType>
  class ConfigIterator : public std::iterator <std::bidirectional_iterator_tag, const ConfigType> {

    PrimClexType *m_primclex;
    Index m_scel_index;
    Index m_config_index;
    bool m_selected;
  public:

    ConfigIterator();

    //ConfigIterator(const ConfigIterator &iter);

    ConfigIterator(PrimClexType *primclex,
                   Index scel_index,
                   Index config_index,
                   bool _selected = false);

    //ConfigIterator &operator=(ConfigIterator iter);

    ConfigType &operator*() const;

    ConfigType *operator->() const;

    // could be const, but that seems weird
    void set_selected(bool _select);

    bool selected() const;

    bool operator==(const ConfigIterator &iter) const;

    bool operator!=(const ConfigIterator &iter) const;

    ConfigIterator &operator++();

    ConfigIterator operator++(int);

    ConfigIterator &operator--();

    ConfigIterator operator--(int);

    Index config_ind() const {
      return m_config_index;
    }
    Index scel_ind() const {
      return m_scel_index;
    }

    friend void swap<>(ConfigIterator &a, ConfigIterator &b);

  private:

    int config_list_size() const;
    void _next_config();
  };



  /// Definitions

  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType>::ConfigIterator() {}

  /*
  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType>::ConfigIterator(const ConfigIterator<ConfigType, PrimClexType> &iter)
    : m_primclex(iter.m_primclex), m_scel_index(iter.m_scel_index), m_config_index(iter.m_config_index) {
  }
  */

  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType>::ConfigIterator(PrimClexType *primclex,
                                                           Index scel_index,
                                                           Index config_index,
                                                           bool _selected)
    : m_primclex(primclex), m_scel_index(scel_index), m_config_index(config_index), m_selected(_selected) {
  }

  /*
  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType> &ConfigIterator<ConfigType, PrimClexType>::operator=(ConfigIterator<ConfigType, PrimClexType> iter) {
    swap(*this, iter);
    return *this;
    }
  */

  /// Specialize for Configuration, const Configuration, Transition, const Transition
  template <typename ConfigType, typename PrimClexType>
  ConfigType &ConfigIterator<ConfigType, PrimClexType>::operator*() const {
    return m_primclex->get_supercell(m_scel_index).get_config(m_config_index);
  }
  //template<> Configuration &ConfigIterator<Configuration, PrimClex>::operator*();
  //template<> const Configuration &ConfigIterator<const Configuration, const PrimClex>::operator*();
  //template<> Transition &ConfigIterator<Transition, const PrimClexType>::operator*();
  //template<> const Transition &ConfigIterator<const Transition, const PrimClexType>::operator*();

  template <typename ConfigType, typename PrimClexType>
  ConfigType *ConfigIterator<ConfigType, PrimClexType>::operator->() const {
    return &(operator*());
  }

  template <typename ConfigType, typename PrimClexType>
  bool ConfigIterator<ConfigType, PrimClexType>::selected() const {
    return (operator*()).selected();
  }

  template <typename ConfigType, typename PrimClexType>
  bool ConfigIterator<ConfigType, PrimClexType>::operator==(const ConfigIterator<ConfigType, PrimClexType> &iter) const {
    return (m_primclex == iter.m_primclex && m_scel_index == iter.m_scel_index && m_config_index == iter.m_config_index && m_selected == iter.m_selected);
  }

  template <typename ConfigType, typename PrimClexType>
  bool ConfigIterator<ConfigType, PrimClexType>::operator!=(const ConfigIterator<ConfigType, PrimClexType> &iter) const {
    return !(*this == iter);
  }


  template <typename ConfigType, typename PrimClexType>
  void ConfigIterator<ConfigType, PrimClexType>::_next_config() {
    m_config_index++;
    while(m_scel_index < m_primclex->get_supercell_list().size() && m_config_index >= config_list_size()) {
      m_scel_index++;
      m_config_index = 0;
    }
    //return *this;
  }

  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType> &ConfigIterator<ConfigType, PrimClexType>::operator++() {
    //std::cout << "m_scel_index = " << m_scel_index << "; m_config_index = " << m_config_index << ";";
    //if(m_scel_index <  m_primclex->get_supercell_list().size() && m_config_index < config_list_size())
    //  std::cout << " selected = " << m_primclex->get_supercell(m_scel_index).get_config(m_config_index).selected();
    //std::cout << "\n -increment-\n";
    _next_config();
    // don't increment past the end
    while(m_scel_index <  m_primclex->get_supercell_list().size() && (m_selected && !(m_primclex->get_supercell(m_scel_index).get_config(m_config_index).selected()))) {
      _next_config();
    }
    //std::cout << "m_scel_index = " << m_scel_index << "; m_config_index = " << m_config_index << ";";
    //if(m_scel_index <  m_primclex->get_supercell_list().size() && m_config_index < config_list_size())
    //  std::cout << " selected = " << m_primclex->get_supercell(m_scel_index).get_config(m_config_index).selected();
    //std::cout << "\n\n";
    return *this;
  }

  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType> ConfigIterator<ConfigType, PrimClexType>::operator++(int) {
    ConfigIterator<ConfigType, PrimClexType> cp(*this);
    ++(*this);
    return cp;
  }


  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType> &ConfigIterator<ConfigType, PrimClexType>::operator--() {
    // Can't always dereference result if supercell[0] has no configurations, but should be okay if you compare to
    // PrimClex::config_begin();
    while(m_config_index == 0 && (m_scel_index > 0 || config_list_size() == 0)) {
      m_scel_index--;
      m_config_index = config_list_size();
    }

    m_config_index--;
  }

  template <typename ConfigType, typename PrimClexType>
  ConfigIterator<ConfigType, PrimClexType> ConfigIterator<ConfigType, PrimClexType>::operator--(int) {
    ConfigIterator<ConfigType, PrimClexType> cp(*this);
    --(*this);
    return cp;
  }


  template<typename ConfigType, typename PrimClexType>
  void swap(ConfigIterator<ConfigType, PrimClexType> &a, ConfigIterator<ConfigType, PrimClexType> &b) {
    using std::swap;

    swap(a.m_primclex, b.m_primclex);
    swap(a.m_scel_index, b.m_scel_index);
    swap(a.m_config_index, b.m_config_index);
  }

  /// Specialize for Configuration, const Configuration, Transition, const Transition
  //template<> int ConfigIterator<Configuration, PrimClex>::config_list_size() const;
  //template<> int ConfigIterator<const Configuration, const PrimClex>::config_list_size() const;
  //template<> int ConfigIterator<Transition, PrimClex>::config_list_size() const;
  //template<> int ConfigIterator<const Transition, const PrimClex>::config_list_size() const;

}
#endif
