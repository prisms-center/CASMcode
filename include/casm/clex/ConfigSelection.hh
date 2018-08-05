#ifndef ConfigSelection_HH
#define ConfigSelection_HH

#include <limits>
#include "casm/casm_io/DataFormatter.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  template <bool IsConst, bool IsConstIterator>
  class ConfigSelectionIterator;

  template <bool IsConst>
  class ConfigSelection;

  typedef ConfigSelection<true> ConstConfigSelection;

  /** \defgroup Selection
   *  \brief Enables creating, reading, writing, and using selections of objects
   *  @{
   */

  template <bool IsConst, bool IsConstIterator>
  class ConfigSelectionIterator : public std::iterator <std::bidirectional_iterator_tag,
    typename std::conditional<IsConstIterator, const Configuration, CASM_TMP::ConstSwitch<IsConst, Configuration> >::type > {

    typedef std::iterator<std::bidirectional_iterator_tag,
            typename std::conditional<IsConstIterator, const Configuration,
            CASM_TMP::ConstSwitch<IsConst, Configuration> >::type > std_iterator_type;

  public:

    typedef typename std::conditional<IsConstIterator,
            std::map<std::string, bool>::const_iterator,
            std::map<std::string, bool>::iterator>::type MapIterator;
    typedef CASM_TMP::ConstSwitch<IsConst, PrimClex> PrimClexType;
    typedef typename std_iterator_type::reference reference;
    typedef typename std_iterator_type::pointer pointer;


    ConfigSelectionIterator();

    ConfigSelectionIterator(const MapIterator &it,
                            const MapIterator &begin,
                            const MapIterator &end,
                            PrimClexType *prim,
                            bool _selected_only = false);

    template<bool ArgIsConst, bool ArgIsConstIterator>
    ConfigSelectionIterator(const ConfigSelectionIterator<ArgIsConst, ArgIsConstIterator> &iter);

    bool selected() const;

    void set_selected(bool is_selected);

    const std::string &name() const;

    bool operator==(const ConfigSelectionIterator &_it) const;

    bool operator!=(const ConfigSelectionIterator &_it) const;

    reference operator*() const;

    pointer operator->() const;

    ConfigSelectionIterator &operator++();

    ConfigSelectionIterator operator++(int);

    ConfigSelectionIterator &operator--();

    ConfigSelectionIterator operator--(int);

  private:
    friend class ConfigSelectionIterator < IsConst, !IsConstIterator >;
    friend class ConfigSelection<IsConst>;

    PrimClexType *m_primclex;
    MapIterator m_it;
    MapIterator m_begin;
    MapIterator m_end;
    bool m_selected_only;
  };

  template <bool IsConst>
  class ConfigSelection {
  public:
    typedef CASM_TMP::ConstSwitch<IsConst, PrimClex> PrimClexType;
    typedef CASM_TMP::ConstSwitch<IsConst, Configuration> value_type;
    typedef ConfigSelectionIterator<IsConst, false> iterator;
    typedef ConfigSelectionIterator<IsConst, true> const_iterator;


    /// \brief Default constructor
    ConfigSelection() {};

    /// \brief Construct a configuration selection from all configurations
    //explicit ConfigSelection(PrimClexType &_primclex);

    /// \brief Reads a configuration selection from file at 'selection_path'
    /// - If selection_path=="MASTER", load master config_list as selection
    /// - If selection_path=="ALL", load all from config_list as selection
    /// - If selection_path=="NONE", load all from config_list as not selected
    /// - If selection_path=="CALCULATED", load configurations for which 'is_calculated' returns true
    /// - Else, assume path to file:
    ///   - Checks extension to determine file type:
    ///     - Ending in '.json' or '.JSON' for JSON formatted file
    ///     - Otherwise, CSV formatted file
    ///   - If column headers are detected, they are also stored, and can be accessed using 'ConfigSelection::col_headers()'
    ///     - Currently, the first two columns are always 'name' and 'selected' so these headers are assumed and not stored
    ///
    ConfigSelection(PrimClexType &_primclex, const fs::path &selection_path = "MASTER");

    //ConfigSelection(const ConfigSelection &) =default;

    ConfigSelection(const ConfigSelection<false> &RHS) :
      m_primclex(RHS.m_primclex), m_name(RHS.m_name), m_config(RHS.m_config), m_col_headers(RHS.m_col_headers) {
      //swap(m_col_headers,RHS.m_col_headers);
      //swap(m_config,RHS.m_config);
    }

    ConfigSelection &operator=(const ConfigSelection &) = default;

    Index size() const {
      return m_config.size();
    }

    iterator find(const std::string &configname) {
      return iterator(m_config.find(configname), m_config.begin(), m_config.end(), m_primclex);
    }

    const_iterator find(const std::string &configname) const {
      return const_iterator(m_config.find(configname), m_config.begin(), m_config.end(), m_primclex);
    }

    void read(std::istream &_input);

    const jsonParser &from_json(const jsonParser &_json);

    jsonParser &to_json(const DataFormatterDictionary<Configuration> &_dict,
                        jsonParser &_json,
                        bool only_selected = false) const;

    /// \brief check if configuration is selected (returns false if 'configname' cannot be found
    bool selected(const std::string &configname) const {
      auto find_it = m_config.find(configname);
      return find_it != m_config.end() && find_it->second;
    }

    bool selected(const Configuration &config) const {
      return selected(config.name());
    }

    void set_selected(const std::string &configname, bool is_selected) {
      m_config[configname] = is_selected;
    }

    void set_selected(const Configuration &config, bool is_selected) {
      set_selected(config.name(), is_selected);
    }

    iterator config_begin() {
      return iterator(m_config.begin(), m_config.begin(), m_config.end(), m_primclex);
    }

    const_iterator config_begin() const {
      return const_iterator(m_config.cbegin(), m_config.cbegin(), m_config.cend(), m_primclex);
    }

    iterator config_end() {
      return iterator(m_config.end(), m_config.begin(), m_config.end(), m_primclex);
    }

    const_iterator config_end() const {
      return const_iterator(m_config.cend(), m_config.cbegin(), m_config.cend(), m_primclex);
    }

    const_iterator config_cbegin() const {
      return const_iterator(m_config.cbegin(), m_config.cbegin(), m_config.cend(), m_primclex);
    }

    const_iterator config_cend() const {
      return const_iterator(m_config.cend(), m_config.cbegin(), m_config.cend(), m_primclex);
    }

    iterator selected_config_begin() {
      return iterator(m_config.begin(), m_config.begin(), m_config.end(), m_primclex, true);
    }

    const_iterator selected_config_begin() const {
      return const_iterator(m_config.cbegin(), m_config.cbegin(), m_config.cend(), m_primclex, true);
    }

    iterator selected_config_end() {
      return iterator(m_config.end(), m_config.begin(), m_config.end(), m_primclex, true);
    }

    const_iterator selected_config_end() const {
      return const_iterator(m_config.cend(), m_config.cbegin(), m_config.cend(), m_primclex, true);
    }

    const_iterator selected_config_cbegin() const {
      return const_iterator(m_config.cbegin(), m_config.cbegin(), m_config.cend(), m_primclex, true);
    }

    const_iterator selected_config_cend() const {
      return const_iterator(m_config.cend(), m_config.cbegin(), m_config.cend(), m_primclex, true);
    }

    std::pair<iterator, bool> insert(const std::pair<std::string, bool> &value) {
      auto res = m_config.insert(value);
      return std::make_pair(iterator(res.first,  m_config.begin(), m_config.end(), m_primclex), res.second);
    }

    const_iterator erase(const const_iterator &it) {
      bool selected_only = it.m_selected_only;
      return const_iterator(m_config.erase(it.m_it), m_config.cbegin(), m_config.cend(), m_primclex, selected_only);
    }

    int erase(const std::string &configname) {
      return m_config.erase(configname);
    }

    const std::vector<std::string> &col_headers() const {
      return m_col_headers;
    }

    const std::string &name() const {
      return m_name;
    }

    void print(const DataFormatterDictionary<Configuration> &_dict,
               std::ostream &_out,
               bool only_selected = false) const;

  private:
    friend class ConfigSelection < !IsConst >;

    PrimClexType *m_primclex;
    std::string m_name;
    std::map<std::string, bool> m_config;
    std::vector<std::string> m_col_headers;

  };

  template<bool IsConst, bool IsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator>::ConfigSelectionIterator() { }

  template<bool IsConst, bool IsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator>::ConfigSelectionIterator(
    const typename ConfigSelectionIterator<IsConst, IsConstIterator>::MapIterator &it,
    const typename ConfigSelectionIterator<IsConst, IsConstIterator>::MapIterator &begin,
    const typename ConfigSelectionIterator<IsConst, IsConstIterator>::MapIterator &end,
    typename ConfigSelectionIterator<IsConst, IsConstIterator>::PrimClexType *prim,
    bool _selected_only) :
    m_primclex(prim), m_it(it), m_begin(begin), m_end(end), m_selected_only(_selected_only) {

    while(m_selected_only && m_begin != m_end && !m_begin->second)
      ++m_begin;

    while(m_selected_only && m_it != m_end && !m_it->second)
      ++m_it;
  }

  /*
  template<>
  ConfigSelectionIterator<true, true>::ConfigSelectionIterator(const ConfigSelectionIterator<true, true> &iter) :
    m_it(iter.m_it),
    m_begin(iter.m_begin),
    m_end(iter.m_end),
    m_primclex(iter.m_primclex),
    m_selected_only(iter.m_selected_only) {}

  template<>
  ConfigSelectionIterator<true, true>::ConfigSelectionIterator(const ConfigSelectionIterator<true, false> &iter) :
    m_it(iter.m_it),
    m_begin(iter.m_begin),
    m_end(iter.m_end),
    m_primclex(iter.m_primclex),
    m_selected_only(iter.m_selected_only) {}

  template<>
  ConfigSelectionIterator<false, false>::ConfigSelectionIterator(const ConfigSelectionIterator<false, false> &iter) :
    m_it(iter.m_it),
    m_begin(iter.m_begin),
    m_end(iter.m_end),
    m_primclex(iter.m_primclex),
    m_selected_only(iter.m_selected_only) {}

  template<>
  ConfigSelectionIterator<false, true>::ConfigSelectionIterator(const ConfigSelectionIterator<false, false> &iter) :
    m_it(iter.m_it),
    m_begin(iter.m_begin),
    m_end(iter.m_end),
    m_primclex(iter.m_primclex),
    m_selected_only(iter.m_selected_only) {}
  */

  template<bool IsConst, bool IsConstIterator>
  template<bool ArgIsConst, bool ArgIsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator>::ConfigSelectionIterator(const ConfigSelectionIterator<ArgIsConst, ArgIsConstIterator> &iter) :
    m_it(iter.m_it),
    m_begin(iter.m_begin),
    m_end(iter.m_end),
    m_primclex(iter.m_primclex),
    m_selected_only(iter.m_selected_only) {}


  template<bool IsConst, bool IsConstIterator>
  bool ConfigSelectionIterator<IsConst, IsConstIterator>::selected() const {
    return m_it->second;
  }

  template<bool IsConst, bool IsConstIterator>
  void ConfigSelectionIterator<IsConst, IsConstIterator>::set_selected(bool is_selected) {
    m_it->second = is_selected;
  }

  template<bool IsConst, bool IsConstIterator>
  const std::string &ConfigSelectionIterator<IsConst, IsConstIterator>::name() const {
    return m_it->first;
  }

  template<bool IsConst, bool IsConstIterator>
  bool ConfigSelectionIterator<IsConst, IsConstIterator>::operator==(const ConfigSelectionIterator &_it) const {
    return (m_it == _it.m_it);
  }

  template<bool IsConst, bool IsConstIterator>
  bool ConfigSelectionIterator<IsConst, IsConstIterator>::operator!=(const ConfigSelectionIterator &_it) const {
    return !((*this) == _it);
  }

  template<bool IsConst, bool IsConstIterator>
  typename ConfigSelectionIterator<IsConst, IsConstIterator>::reference ConfigSelectionIterator<IsConst, IsConstIterator>::operator*() const {
    return m_primclex->configuration(m_it->first);
  }

  template<bool IsConst, bool IsConstIterator>
  typename ConfigSelectionIterator<IsConst, IsConstIterator>::pointer ConfigSelectionIterator<IsConst, IsConstIterator>::operator->() const {
    return &(m_primclex->configuration(m_it->first));
  }

  template<bool IsConst, bool IsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator> &ConfigSelectionIterator<IsConst, IsConstIterator>::operator++() {
    ++m_it;
    while(m_selected_only && m_it != m_end && !selected())
      ++m_it;
    return (*this);
  }

  template<bool IsConst, bool IsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator> ConfigSelectionIterator<IsConst, IsConstIterator>::operator++(int) {
    ConfigSelectionIterator t_it(*this);
    ++(*this);
    return t_it;
  }

  template<bool IsConst, bool IsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator> &ConfigSelectionIterator<IsConst, IsConstIterator>::operator--() {
    --m_it;
    while(m_selected_only && m_it != m_begin && !selected())
      --m_it;
    return (*this);
  }

  template<bool IsConst, bool IsConstIterator>
  ConfigSelectionIterator<IsConst, IsConstIterator> ConfigSelectionIterator<IsConst, IsConstIterator>::operator--(int) {
    ConfigSelectionIterator t_it(*this);
    --(*this);
    return t_it;
  }


  template<bool IsConst>
  std::ostream &operator<<(std::ostream &_stream, const ConfigSelection<IsConst> &selection) {
    selection.print(_stream);
    return _stream;
  }

  bool get_selection(const Array<std::string> &criteria, const Configuration &config, bool is_selected);

  namespace ConfigSelection_impl {

    bool is_operator(const std::string &q);

    std::string operate(const std::string &q, const std::string &A);

    std::string operate(const std::string &q, const std::string &A, const std::string &B);

    bool is_unary(const std::string &q);

    std::string convert_variable(const std::string &q, const Configuration &config);
  }

  /** @} */
}
#include "casm/clex/ConfigSelection_impl.hh"
#endif
