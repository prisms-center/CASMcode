#include "casm/casm_io/DataStream.hh"
#include "casm/container/Counter.hh"
#include "casm/casm_io/DataFormatterTools.hh"

namespace CASM {


  //****************************************************************************************
  template<typename DataObject>
  void BaseDatumFormatter<DataObject>::_parse_index_expression(const std::string &_expr) {
    //std::cout << "Parsing index expression: " << _expr << "\n";
    typedef boost::tokenizer<boost::char_separator<char> >
    tokenizer;
    boost::char_separator<char> sep1(","), sep2(" \t", ":");
    tokenizer tok1(_expr, sep1);
    std::vector<std::string> split_expr(tok1.begin(), tok1.end());
    std::vector<difference_type> ind_vec_begin(split_expr.size());
    std::vector<difference_type> ind_vec_end(split_expr.size());
    for(Index i = 0; i < split_expr.size(); i++) {
      Index rev_i = split_expr.size() - (i + 1);
      tokenizer tok2(split_expr[i], sep2);
      std::vector<std::string> ind_expr(tok2.begin(), tok2.end());
      if(ind_expr.size() == 1) {
        if(ind_expr[0][0] == ':')
          ind_vec_begin[rev_i] = -1;
        else
          ind_vec_begin[rev_i] = boost::lexical_cast<difference_type>(ind_expr[0]);
        ind_vec_end[rev_i] = ind_vec_begin[rev_i];
      }
      else if(ind_expr.size() == 3) {
        ind_vec_begin[rev_i] = boost::lexical_cast<difference_type>(ind_expr[0]);
        ind_vec_end[rev_i] = boost::lexical_cast<difference_type>(ind_expr[2]);
      }
      else
        throw std::runtime_error(std::string("In BaseDatumFormatter::_parse_index_expression(), invalid expression \"")
                                 + _expr + "\" passed as indices for format keyword '" + name() + "'\n");
    }
    Counter<std::vector<difference_type> > ind_count(ind_vec_begin, ind_vec_end, std::vector<difference_type>(split_expr.size(), 1));
    for(; ind_count.valid(); ++ind_count) {
      m_index_rules.push_back(std::vector<difference_type>(ind_count().rbegin(), ind_count().rend()));
    }

  }

  //******************************************************************************

  template<typename DataObject>
  bool DataFormatter<DataObject>::validate(const DataObject &_obj) const {
    if(!m_initialized)
      _initialize(_obj);
    for(Index i = 0; i < m_data_formatters.size(); i++)
      if(!m_data_formatters[i]->validate(_obj))
        return false;

    return true;
  }

  //******************************************************************************

  template<typename DataObject>
  void DataFormatter<DataObject>::inject(const DataObject &_obj, DataStream &_stream) const {
    if(!m_initialized)
      _initialize(_obj);

    Index num_pass(1), tnum;
    for(Index i = 0; i < m_data_formatters.size(); i++) {
      tnum = m_data_formatters[i]->num_passes(_obj);
      if(tnum == 1)
        continue;
      if(num_pass == 1 || tnum == num_pass)
        num_pass = tnum;
      else {
        std::cerr << "CRITICAL ERROR: Requesting to print formatted data elements that require different number of lines.\n"
                  << "                Exiting...\n";
        exit(1);
      }
    }
    for(Index np = 0; np < num_pass; np++) {
      for(Index i = 0; i < m_data_formatters.size(); i++) {
        m_data_formatters[i]->inject(_obj, _stream, np);
      }
      _stream.newline();
    }

    return;
  }

  //******************************************************************************

  template<typename DataObject>
  void DataFormatter<DataObject>::print(const DataObject &_obj, std::ostream &_stream) const {
    if(!m_initialized)
      _initialize(_obj);
    _stream << std::setprecision(m_prec) << std::fixed;
    Index num_pass(1), tnum;
    for(Index i = 0; i < m_data_formatters.size(); i++) {
      tnum = m_data_formatters[i]->num_passes(_obj);
      if(tnum == 1)
        continue;
      if(num_pass == 1 || tnum == num_pass)
        num_pass = tnum;
      else {
        std::cerr << "CRITICAL ERROR: Requesting to print formatted data elements that require different number of lines.\n"
                  << "                Exiting...\n";
        exit(1);
      }
    }
    std::stringstream t_ss;
    int depad_request;
    for(Index np = 0; np < num_pass; np++) {
      depad_request = 0;
      for(Index i = 0; i < m_data_formatters.size(); i++) {
        t_ss.clear();
        t_ss.str(std::string());
        m_data_formatters[i]->print(_obj, t_ss, np);

        // two space fixed separator
        _stream << "  ";
        //variable separator
        if(depad_request + 2 < m_col_sep[i])
          _stream << std::string(m_col_sep[i] - depad_request - 2, ' ');
        _stream << t_ss.str();
        depad_request = m_col_sep[i] + int(t_ss.str().size()) - m_col_width[i];

      }
      _stream << std::endl;
    }

    return;
  }

  //******************************************************************************

  template<typename DataObject>
  jsonParser &DataFormatter<DataObject>::to_json(const DataObject &_obj, jsonParser &json) const {
    if(!m_initialized)
      _initialize(_obj);
    for(Index i = 0; i < m_data_formatters.size(); i++) {
      m_data_formatters[i]->to_json(_obj, json[m_data_formatters[i]->short_header(_obj)]);
    }

    return json;
  }

  //******************************************************************************

  template<typename DataObject>
  jsonParser &DataFormatter<DataObject>::to_json_arrays(const DataObject &_obj, jsonParser &json) const {
    if(!m_initialized) {
      _initialize(_obj);
    }

    jsonParser::iterator it;
    jsonParser::iterator end = json.end();

    for(Index i = 0; i < m_data_formatters.size(); i++) {
      jsonParser tmp;
      m_data_formatters[i]->to_json(_obj, tmp);
      it = json.find(m_data_formatters[i]->short_header(_obj));
      if(it == end) {
        json[m_data_formatters[i]->short_header(_obj)].put_array().push_back(tmp);
      }
      else {
        it->push_back(tmp);
      }
    }

    return json;
  }

  //******************************************************************************

  template<typename DataObject>
  void DataFormatter<DataObject>::print_header(const DataObject &_template_obj, std::ostream &_stream) const {
    _stream << m_comment;
    if(!m_initialized)
      _initialize(_template_obj);
    int header_size, twidth;
    for(Index i = 0; i < m_data_formatters.size(); i++) {
      std::stringstream t_ss;
      m_data_formatters[i]->print(_template_obj, t_ss);
      header_size = (m_data_formatters[i]->long_header(_template_obj)).size();
      twidth = m_sep;
      if(i == 0)
        twidth += max(t_ss.str().size(), header_size + m_comment.size());
      else
        twidth += max(int(t_ss.str().size()), header_size);
      m_col_width[i] = twidth;
      m_col_sep[i] = twidth - int(t_ss.str().size());
      if(i == 0)
        twidth -= m_comment.size();
      _stream << std::string(twidth - header_size, ' ') << m_data_formatters[i]->long_header(_template_obj);
    }
    _stream <<  std::endl;
    return;
  }

  //******************************************************************************

  template<typename DataObject>
  void DataFormatter<DataObject>::_initialize(const DataObject &_template_obj) const {
    for(Index i = 0; i < m_data_formatters.size(); i++)
      m_data_formatters[i]->init(_template_obj);
    m_initialized = true;
    return;
  }

  //******************************************************************************


  /// \brief Equivalent to find, but throw error with suggestion if _name not found
  template<typename DataObject, typename DatumFormatterType>
  typename DataFormatterDictionary<DataObject, DatumFormatterType>::const_iterator
  DataFormatterDictionary<DataObject, DatumFormatterType>::lookup(
    const key_type &_name) const {

    auto res = this->find(_name);
    if(res != this->end()) {
      return res;
    }
    else {

      // If no match, try to use demerescau-levenshtein distance to make a helpful suggestion
      auto it = this->begin();
      int min_dist(-1);
      for(; it != this->end(); ++it) {
        int dist = dl_string_dist(_name, it->name());
        if(min_dist < 0 || dist < min_dist) {
          min_dist = dist;
          //std::cout << "New best: \"" << it->first << "\" aka \"" << it->second->name() << "\"\n";
          res = it;
        }
      }

      throw std::runtime_error("CRITICAL ERROR: Invalid format flag \"" + _name + "\" specified.\n"
                               + "                Did you mean \"" + res->name() + "\"?\n");

    }

  }

  /// \brief Generates formatted help using the 'name' and 'description' of all
  ///        contained BaseDatumFormatter
  ///
  template<typename DataObject, typename DatumFormatterType>
  void DataFormatterDictionary<DataObject, DatumFormatterType>::print_help(
    std::ostream &_stream,
    typename BaseDatumFormatter<DataObject>::FormatterType ftype,
    int width, int separation) const {

    const_iterator it_begin(this->cbegin()), it_end(this->cend());
    std::string::size_type len(0);
    for(auto it = it_begin; it != it_end; ++it) {
      if(ftype == it->type())
        len = max(len, it->name().size());
    }
    for(auto it = it_begin; it != it_end; ++it) {
      if(ftype != it->type())
        continue;
      _stream << std::string(5, ' ') << it->name() << std::string(len - it->name().size() + separation, ' ');
      std::string::size_type wcount(0);
      std::string::const_iterator str_end(it->description().cend());
      for(std::string::const_iterator str_it = it->description().cbegin(); str_it != str_end; ++str_it) {
        if(wcount >= width && isspace(*str_it)) {
          _stream << std::endl << std::string(5 + len + separation, ' ');
          wcount = 0;
        }
        else {
          _stream << *str_it;
          wcount++;
        }
      }
      _stream << std::endl << std::endl;
    }
  }


  /// \brief Use the vector of strings to build a DataFormatter<DataObject>
  ///
  /// Expects vector of "formattername(argument1,argument2,...)"
  ///
  /// Uses DataFormatterDictionary<DataObject>::lookup to suggest alternatives if
  /// exact request not found.
  ///
  template<typename DataObject, typename DatumFormatterType>
  DataFormatter<DataObject> DataFormatterDictionary<DataObject, DatumFormatterType>::parse(
    const std::vector<std::string> &input) const {

    DataFormatter<DataObject> formatter;
    std::vector<std::string> format_tags, format_args;
    for(Index i = 0; i < input.size(); i++) {
      split_formatter_expression(input[i], format_tags, format_args);
    }
    for(Index i = 0; i < format_tags.size(); i++) {
      formatter.push_back(*lookup(format_tags[i]), format_args[i]);
    }
    return formatter;
  }

  /// \brief Use a single string to build a DataFormatter<DataObject>
  ///
  /// Expects string of "formattername(argument1,argument2,...) formattername(argument1,argument2,...) ..."
  ///
  /// Uses DataFormatterDictionary<DataObject>::lookup to suggest alternatives if
  /// exact request not found.
  ///
  template<typename DataObject, typename DatumFormatterType>
  DataFormatter<DataObject> DataFormatterDictionary<DataObject, DatumFormatterType>::parse(
    const std::string &input) const {

    DataFormatter<DataObject> formatter;
    std::vector<std::string> format_tags, format_args;
    split_formatter_expression(input, format_tags, format_args);

    for(Index i = 0; i < format_tags.size(); i++) {
      formatter.push_back(*lookup(format_tags[i]), format_args[i]);
    }
    return formatter;
  }

  //****************************************************************************************

  template<typename DataObject>
  void DataFormatterParser<DataObject>::load_aliases(const fs::path &alias_path) {
    if(!fs::exists(alias_path)) {
      return;
    }
    jsonParser mjson(alias_path);

    auto it(mjson.cbegin()), it_end(mjson.cend());
    for(; it != it_end; ++it) {
      add_custom_formatter(datum_formatter_alias<DataObject>(it.name(), it->get<std::string>()));
    }
  }

}
