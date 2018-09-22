#include "casm/external/boost.hh"
#include "casm/clex/ConfigIO.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigIOSelected.hh"
namespace CASM {



  template <bool IsConst>
  ConfigSelection<IsConst>::ConfigSelection(typename ConfigSelection<IsConst>::PrimClexType &_primclex, const fs::path &selection_path)
    : m_primclex(&_primclex), m_name(selection_path.string()) {
    if(selection_path == "MASTER") {
      for(auto it = _primclex.config_begin(); it != _primclex.config_end(); ++it) {
        set_selected(it->name(), it->selected());
      }
    }
    else if(selection_path == "NONE") {
      for(auto it = _primclex.config_cbegin(); it != _primclex.config_cend(); ++it) {
        set_selected(it->name(), false);
      }
    }
    else if(selection_path == "ALL") {
      for(auto it = _primclex.config_cbegin(); it != _primclex.config_cend(); ++it) {
        set_selected(it->name(), true);
      }
    }
    else if(selection_path == "CALCULATED") {
      for(auto it = _primclex.config_cbegin(); it != _primclex.config_cend(); ++it) {
        set_selected(it->name(), is_calculated(*it));
      }
    }
    else {
      if(!fs::exists(selection_path)) {
        std::stringstream ss;
        ss << "ERROR in parsing configuation selection name. \n"
           << "  Expected <filename>, 'ALL', 'NONE', 'CALCULATED', or 'MASTER' <--default \n"
           << "  Received: '" << selection_path << "'\n"
           << "  No file named '" << selection_path << "'.";
        throw std::runtime_error(ss.str());
      }
      m_name = fs::absolute(selection_path).string();
      if(selection_path.extension() == ".json" || selection_path.extension() == ".JSON")
        from_json(jsonParser(selection_path));
      else {
        std::ifstream select_file(selection_path.string().c_str());
        read(select_file);
        select_file.close();
      }
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
        if(it.name() == "configname") {
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
          std::string("CRITICAL ERROR: Field 'configname' is missing from Configuration ") +
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
  jsonParser &ConfigSelection<IsConst>::to_json(
    const DataFormatterDictionary<Configuration> &_dict,
    jsonParser &_json,
    bool only_selected) const {
    _json.put_array();

    DataFormatter<Configuration> tformat(ConfigIO::configname(), datum_formatter_alias("selected", ConfigIO::selected_in(*this)));

    tformat.append(_dict.parse(m_col_headers));

    if(only_selected) {
      _json = tformat(selected_config_cbegin(),
                      selected_config_cend());
    }
    else {
      _json = tformat(config_cbegin(),
                      config_cend());
    }

    return _json;
  }
  //******************************************************************************
  template <bool IsConst>
  void ConfigSelection<IsConst>::print(
    const DataFormatterDictionary<Configuration> &_dict,
    std::ostream &_out,
    bool only_selected) const {
    DataFormatter<Configuration> tformat(ConfigIO::configname(), datum_formatter_alias("selected", ConfigIO::selected_in(*this)));

    tformat.append(_dict.parse(m_col_headers));

    if(only_selected) {
      _out << tformat(selected_config_cbegin(),
                      selected_config_cend());
    }
    else {
      _out << tformat(config_cbegin(),
                      config_cend());
    }

  }
  //******************************************************************************
  /*
    template< bool IsConst>
    typename ConfigSelection<IsConst>::value_type &ConfigSelection<IsConst>::operator[](const std::string &configname) const {
      return m_primclex->configuration(configname);
    }
  */
}
