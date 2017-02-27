#include "casm/database/Selection.hh"

namespace CASM {

  namespace DB {

    template<typename ObjType>
    Selection<ObjType>::Selection(Database<ObjType> &_db, fs::path selection_path) :
      m_db(&_db),
      m_primclex(&m_db->primclex()),
      m_name(selection_path.string()) {

      if(selection_path == "MASTER" || selection_path.empty()) {
        fs::path master_selection_path = primclex().dir().template master_selection<ObjType>();
        if(fs::exists(master_selection_path)) {
          fs::ifstream select_file(master_selection_path);
          read(select_file);
          select_file.close();
        }
        else {
          for(const auto &obj : db()) {
            m_data.insert(std::make_pair(obj.name(), false));
          }
        }
      }
      else if(selection_path == "NONE") {
        for(const auto &obj : db()) {
          m_data.insert(std::make_pair(obj.name(), false));
        }
      }
      else if(selection_path == "ALL") {
        for(const auto &obj : db()) {
          m_data.insert(std::make_pair(obj.name(), true));
        }
      }
      else if(selection_path == "CALCULATED") {
        for(const auto &obj : db()) {
          m_data.insert(std::make_pair(obj.name(), is_calculated(obj)));
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
        if(selection_path.extension() == ".json" || selection_path.extension() == ".JSON") {
          from_json(jsonParser(selection_path));
        }
        else {
          fs::ifstream select_file(selection_path);
          read(select_file);
          select_file.close();
        }
      }
    }

    /// \brief True if obj is in Selection and is selected; false otherwise
    template<typename ObjType>
    bool Selection<ObjType>::selected(const ObjType &obj) const {
      if(!obj.alias().empty()) {
        auto it = m_data.find(obj.alias());
        if(it != m_data.end() && it->second) {
          return true;
        }
      }
      auto it = m_data.find(obj.name());
      if(it != m_data.end() && it->second) {
        return true;
      }
      return false;
    }

    /// \brief Ensure obj is in Selection and set selected to specified value
    template<typename ObjType>
    void Selection<ObjType>::set_selected(const ObjType &obj, bool selected) const {
      if(!obj.alias().empty()) {
        auto it = m_data.find(obj.alias());
        if(it != m_data.end()) {
          it->second = selected;
          return;
        }
      }
      auto it = m_data.find(obj.name());
      if(it != m_data.end()) {
        it->second = selected;
      }
      return;
    }

    //******************************************************************************

    template<typename ObjType>
    void Selection<ObjType>::read(std::istream &_input) {
      std::string tname;
      bool tselect;
      _input >> std::ws;
      if(_input.peek() == '#') {
        _input.get();
        // discard first two columns (name, selected)
        _input >> tname;
        _input >> tname;
        std::getline(_input, tname, '\n');
        m_col_headers.clear();
        boost::split(m_col_headers, tname, boost::is_any_of(" \t"), boost::token_compress_on);
        //_input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      while(_input >> tname >> tselect) {
        m_data[tname] = tselect;
        _input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    }

    //******************************************************************************

    template<typename ObjType>
    const jsonParser &Selection<ObjType>::from_json(const jsonParser &_json) {

      if(!_json.is_array()) {
        std::cerr << "CRITICAL ERROR: Unable to initialize a Selection from passed JSON input." << std::endl
                  << "                JSON input must be an array of records." << std::endl
                  << "                Exiting..." << std::endl;
        exit(1);
      }

      std::string tname;
      bool tselected;
      bool contains_name;

      std::set<std::string> prop_set;
      m_data.clear();

      for(Index i = 0; i < _json.size(); i++) {
        auto it(_json[i].cbegin()), it_end(_json[i].cend());

        contains_name = false;
        tselected = true;  // default to selected

        for(; it != it_end; ++it) {

          // for compatibility, include "configname" as "name"
          if(it.name() == "name" || it.name() == "configname") {
            tname = it->get<std::string>();
            contains_name = true;
          }
          else if(it.name() == "selected") {
            tselected = it->get<bool>();
          }
          else {
            prop_set.insert(it.name());
          }
        }

        if(!contains_name) {
          throw std::runtime_error(
            std::string("CRITICAL ERROR: Field 'name' is missing from ") +
            std::to_string(i) + " of json Array passed to Selection::from_json()." +
            " This field is required.");
        }

        m_data[tname] = tselected;

      }

      m_col_headers = std::vector<std::string>(prop_set.begin(), prop_set.end());

      return _json;
    }

    //******************************************************************************

    template<typename ObjType>
    jsonParser &Selection<ObjType>::to_json(
      const DataFormatterDictionary<ObjType> &_dict,
      jsonParser &_json,
      bool only_selected) const {
      _json.put_array();

      DataFormatter<ObjType> tformat(
        name_or_alias<ObjType>(),
        datum_formatter_alias("selected", selected_in(*this)));

      tformat.append(_dict.parse(m_col_headers));

      if(only_selected) {
        _json = tformat(selected().begin(), selected().end());
      }
      else {
        _json = tformat(all().begin(), all().end());
      }

      return _json;
    }

    //******************************************************************************

    template<typename ObjType>
    void Selection<ObjType>::print(
      const DataFormatterDictionary<ObjType> &_dict,
      std::ostream &_out,
      bool only_selected) const {

      DataFormatter<ObjType> tformat(
        name_or_alias<ObjType>(),
        datum_formatter_alias("selected", selected_in(*this)));

      tformat.append(_dict.parse(m_col_headers));

      if(only_selected) {
        _out << tformat(selected().begin(), selected().end());
      }
      else {
        _out << tformat(all().begin(), all().end());
      }

    }

  }
}

