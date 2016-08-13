#include <algorithm>
#include "casm/external/boost.hh"
namespace CASM {

  template<typename ValueType, typename ArgType, typename DataObject>
  bool DataFormatterOperator<ValueType,
       ArgType,
  DataObject>::parse_args(const std::string &_args) {

    if(_args.size() == 0)
      return true;
    if(!m_arg_formatter.empty())
      return false;

    std::vector<std::string> format_tags, subexprs;
    split_formatter_expression(_args, format_tags, subexprs);

    //std::cout << "_args: " << _args << "\n";
    //std::cout << "format_tabs: " << format_tags << "\n"
    //<< "subexprs: '" << subexprs << "'\n";
    for(Index i = 0; i < format_tags.size(); i++) {
      std::string ttag(format_tags[i].size(), ' ');
      std::transform(format_tags[i].cbegin(), format_tags[i].cend(), ttag.begin(), tolower);

      char ch = ttag[0];
      if(boost::is_any_of("-+.0123456789")(ch)) { // tag is a constant number
        if(std::any_of(ttag.cbegin(),
                       ttag.cend(),
                       boost::is_any_of(".e"))) {
          double val;
          try {
            val = std::stod(ttag);
          }
          catch(...) {
            throw std::runtime_error("Unable to parse '" + ttag + "' from subexpression " + _args + "\n");
          }
          m_arg_formatter.push_back(ConstantValueFormatter<double, DataObject>(ttag, val));
        }
        else {
          long val;
          try {
            val = std::stol(ttag);
          }
          catch(...) {
            throw std::runtime_error("Unable to parse '" + ttag + "' from subexpression " + _args + "\n");
          }
          m_arg_formatter.push_back(ConstantValueFormatter<long, DataObject>(ttag, val));
        }
      }
      else if(ttag == "false")
        m_arg_formatter.push_back(ConstantValueFormatter<bool, DataObject>(ttag, false));
      else if(ttag == "true")
        m_arg_formatter.push_back(ConstantValueFormatter<bool, DataObject>(ttag, true));
      else if((ch == '\'' || ch == '\"') && ch == ttag.back() && ttag.size() > 1) {
        m_arg_formatter.push_back(ConstantValueFormatter<std::string, DataObject>(format_tags[i], std::string((++format_tags[i].cbegin()), (--format_tags[i].cend()))));
      }
      else {
        //std::cout << "About to push_back parsed formatter!!\n";
        //can throw, should be handled by caller
        const BaseDatumFormatter<DataObject> &proto_format = *this->home().lookup(ttag);

        m_arg_formatter.push_back(proto_format, subexprs[i]);

      }
    }
    return true;
  }

  //******************************************************************************

  template<typename ValueType, typename ArgType, typename DataObject>
  std::string DataFormatterOperator<ValueType,
      ArgType,
  DataObject>::short_header(const DataObject &_obj)const {
    std::stringstream t_ss;
    m_arg_formatter.print_header(_obj, t_ss);

    std::vector<std::string> format_tags, subexprs;
    split_formatter_expression(t_ss.str(), format_tags, subexprs);
    t_ss.str("");
    t_ss.clear();
    t_ss << name();
    if(format_tags.size()) {
      t_ss << '(';
      for(Index i = 0; i < format_tags.size(); i++) {
        t_ss << format_tags[i];
        if(subexprs[i].size() > 0)
          t_ss << '(' << subexprs[i] << ')';
        if(i + 1 < format_tags.size())
          t_ss << ',';
      }
      t_ss << ')';
    }
    return t_ss.str();
  }

  // ******************************************************************************

  /*
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
  */

  // ******************************************************************************

  /*
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
  */

  /*
    // ******************************************************************************

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

    // ******************************************************************************

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
  */

}
