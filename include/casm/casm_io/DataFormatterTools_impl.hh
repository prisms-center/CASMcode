
namespace CASM {

  template<typename _Container,
           typename DataObject,
           typename _value_type,
           typename _size_type,
           typename Access>
  void Generic1DDatumFormatter<_Container,
       DataObject,
       _value_type,
       _size_type,
  Access>::init(const DataObject &_template_obj) const {
    if(_index_rules().size())
      return;

    Index size = m_size(m_evaluate(_template_obj));
    for(Index i = 0; i < size; i++) {
      _add_rule(std::vector<Index>({i}));
    }
  }

  //******************************************************************************

  template<typename _Container,
           typename DataObject,
           typename _value_type,
           typename _size_type,
           typename Access>
  std::string Generic1DDatumFormatter<_Container,
      DataObject,
      _value_type,
      _size_type,
  Access>::long_header(const DataObject &_template_obj) const {
    std::stringstream t_ss;
    auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
    Index s = max(8 - int(name().size()), 0);
    for(; it != end_it; ++it) {
      t_ss << " " << std::string(s, ' ') << name() << '(' << (*it)[0] << ')';
    }
    return t_ss.str();
  }

  //******************************************************************************

  template<typename _Container,
           typename DataObject,
           typename _value_type,
           typename _size_type,
           typename Access>
  void Generic2DDatumFormatter<_Container,
       DataObject,
       _value_type,
       _size_type,
  Access>::init(const DataObject &_template_obj) const {
    if(_index_rules().size())
      return;

    std::vector<Index> size = m_size(m_evaluate(_template_obj));
    Counter<std::vector<Index> > ind_count(std::vector<Index>(size.size(), 0), size, std::vector<Index>(size.size(), 1));
    for(; ind_count.valid(); ++ind_count) {
      _add_rule(ind_count());
    }

  }

  //******************************************************************************

  template<typename _Container,
           typename DataObject,
           typename _value_type,
           typename _size_type,
           typename Access>
  std::string Generic2DDatumFormatter<_Container,
      DataObject,
      _value_type,
      _size_type,
  Access>::long_header(const DataObject &_template_obj) const {
    std::stringstream t_ss;
    auto it(_index_rules().cbegin()), end_it(_index_rules().cend());
    Index s = max(10 - int(name().size()), 0);
    for(; it != end_it; ++it) {
      t_ss << " " << std::string(s, ' ') << name() << '(' << (*it)[0] << ',' << (*it)[1] << ')';
    }
    return t_ss.str();
  }

}
