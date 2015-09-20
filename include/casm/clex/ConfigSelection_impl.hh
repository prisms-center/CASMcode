#include <boost/algorithm/string.hpp>
#include "casm/clex/ConfigIO.hh"
namespace CASM {



  template <bool IsConst>
  ConfigSelection<IsConst>::ConfigSelection(typename ConfigSelection<IsConst>::PrimClexType &_primclex, const fs::path &selection_path)
    : m_primclex(&_primclex) {
    if(selection_path.extension() == ".json" || selection_path.extension() == ".JSON")
      from_json(jsonParser(selection_path));
    else {
      std::ifstream select_file(selection_path.string().c_str());
      read(select_file);
      select_file.close();
    }
  }

  //******************************************************************************

  template <bool IsConst>
  void ConfigSelection<IsConst>::read(std::istream &_input) {
    std::string tname;
    bool tselect;
    _input >> std::ws;
    if(_input.peek() == '#') {
      _input.get();
      // discard first two columns
      _input >> tname;
      _input >> tname;
      std::getline(_input, tname, '\n');
      m_col_headers.clear();
      boost::split(m_col_headers, tname, boost::is_any_of(" \t"), boost::token_compress_on);
      //_input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    while(_input >> tname >> tselect) {
      m_config[tname] = tselect;
      _input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
  }

  //******************************************************************************

  template <bool IsConst>
  const jsonParser &ConfigSelection<IsConst>::from_json(const jsonParser &_json) {
    if(!_json.is_array()) {
      std::cerr << "CRITICAL ERROR: Unable to initialize a ConfigSelection from passed JSON input." << std::endl
                << "                JSON input must be an array of Configuration records." << std::endl
                << "                Exiting..." << std::endl;
      exit(1);
    }
    bool tselect;
    std::map<std::string, bool> prop_map;
    std::string tname;
    bool tselected;
    bool contains_name;
    m_config.clear();
    for(Index i = 0; i < _json.size(); i++) {
      auto it(_json[i].cbegin()), it_end(_json[i].cend());


      contains_name = false;
      tselected = true;  // default to selected

      for(; it != it_end; ++it) {
        if(it.name() == "name") {
          tname = it->get<std::string>();
          contains_name = true;
        }
        else if(it.name() == "selected") {
          tselected = it->get<bool>();
        }
        else
          prop_map[it.name()] = true;
      }

      if(!contains_name) {
        throw std::runtime_error(
          std::string("CRITICAL ERROR: Field 'name' is missing from Configuration ") +
          std::to_string(i) + " of json Array passed to ConfigSelection::from_json()." +
          " This field is required.");
      }

      m_config[tname] = tselected;

    }

    std::map<std::string, bool>::iterator it(prop_map.begin()), it_end(prop_map.end());
    for(; it != it_end; ++it)
      m_col_headers.push_back(it->first);
    return _json;
  }

  //******************************************************************************

  template <bool IsConst>
  jsonParser &ConfigSelection<IsConst>::to_json(jsonParser &_json, bool only_selected) const {
    _json.put_array();

    DataFormatter<Configuration> tformat(ConfigIOParser::parse(m_col_headers));

    const_iterator it, it_end;

    if(only_selected) {
      it = selected_config_cbegin();
      it_end = selected_config_cend();
    }
    else {
      it = config_cbegin();
      it_end = config_cend();
    }

    for(; it != it_end; ++it) {
      _json.push_back(tformat(*it));
      _json[_json.size() - 1]["name"] = it.name();
      _json[_json.size() - 1]["selected"] = it.selected();
    }
    return _json;
  }
  //******************************************************************************
  template <bool IsConst>
  void ConfigSelection<IsConst>::print(std::ostream &_out, bool only_selected) const {
    DataFormatter<Configuration> tformat(ConfigIOParser::parse(m_col_headers));
    tformat.set_header_prefix("                  name     selected");
    _out << FormatFlag(_out).print_header(false);

    const_iterator it, it_end;

    if(only_selected) {
      it = selected_config_cbegin();
      it_end = selected_config_cend();
    }
    else {
      it = config_cbegin();
      it_end = config_cend();
    }

    if(it == it_end)
      return;

    tformat.print_header(*it, _out);
    std::ostream::streampos pos;
    tformat.set_indent(-30);
    for(; it != it_end; ++it) {


      pos = _out.tellp();
      _out << "    " << it.name() << "        " << it.selected() << "  ";

      _out << tformat(*it);
    }
  }
  //******************************************************************************

  template< bool IsConst>
  typename ConfigSelection<IsConst>::value_type &ConfigSelection<IsConst>::operator[](const std::string &configname) const {
    return m_primclex->configuration(configname);
  }
}

