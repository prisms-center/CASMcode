#include "casm/database/Selection_impl.hh"

#include "casm/clex/PrimClex_impl.hh"
#include "casm/casm_io/DataFormatter_impl.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/casm_io/stream_io/container.hh"
#include "casm/database/Selected_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"

// explicit template instantiations
#define INST_Selection(r, data, type) \
template class SelectionIterator<type,Selection<type>::base_iterator>; \
template class SelectionIterator<type,Selection<type>::base_const_iterator>; \
template class Selection<type>; \

namespace CASM {
  namespace DB {

    namespace {

      // Initialization for "CALCULATED" depends on type:

      /// Use SFINAE to implement if ObjType is a calculable config type:
      ///
      /// - Set selected to value of 'is_calculated(const ObjType&)'
      template<typename ObjType, typename IfConfigType<ObjType>::type * = nullptr>
      void init_calculated(
        typename Selection<ObjType>::map_type &m_data,
        Database<ObjType> &db) {
        for(const auto &obj : db) {
          m_data.insert(std::make_pair(obj.name(), is_calculated(obj)));
        }
      }

      /// Use SFINAE to implement if ObjType is not a calculable config type
      ///
      /// - throw runtime_error
      template<typename ObjType, typename IfNotConfigType<ObjType>::type * = nullptr>
      void init_calculated(
        typename Selection<ObjType>::map_type &m_data,
        Database<ObjType> &db) {
        std::stringstream msg;
        msg << "Selection \"CALCULATED\" is not allowed for type: " << traits<ObjType>::short_name;
        throw std::runtime_error(msg.str());
      }
    }


    /// boost::iterator_facade implementation
    template<typename ObjType, typename BaseIterator>
    const ObjType &SelectionIterator<ObjType, BaseIterator>::dereference() const {
      return *(m_list->db().find(m_it->first));
    }

    /// \brief Use default ObjType database
    template<typename ObjType>
    Selection<ObjType>::Selection(const PrimClex &_primclex, fs::path selection_path) :
      Selection(_primclex.db<ObjType>(), selection_path) {}

    template<typename ObjType>
    Selection<ObjType>::Selection(Database<ObjType> &_db, fs::path selection_path) :
      Logging(_db.primclex()),
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
      else if(selection_path == "EMPTY") {

      }
      else if(selection_path == "ALL") {
        for(const auto &obj : db()) {
          m_data.insert(std::make_pair(obj.name(), true));
        }
      }
      else if(selection_path == "CALCULATED") {
        init_calculated(m_data, db());
      }
      else {
        if(!fs::exists(selection_path)) {
          std::stringstream ss;
          ss << "ERROR in parsing configuration selection name. \n"
             << "  Expected <filename>, 'ALL', 'NONE', 'EMPTY', 'CALCULATED', or 'MASTER' <--default \n"
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

    template<typename ObjType>
    Index Selection<ObjType>::selected_size() const {
      return boost::distance(selected());
    }

    /// \brief True if obj is in Selection and is selected; false otherwise
    template<typename ObjType>
    bool Selection<ObjType>::is_selected(const std::string &name_or_alias) const {
      auto it = m_data.find(db().name(name_or_alias));
      if(it == m_data.end()) {
        return false;
      }
      return it->second;
    }

    /// \brief Set selected objects to value of criteria
    template<typename ObjType>
    void Selection<ObjType>::set(const DataFormatterDictionary<ObjType> &dict,
                                 const std::string &criteria) {
      try {
        if(criteria.size()) {
          DataFormatter<ObjType> tformat(dict.parse(criteria));
          auto it = all().begin();
          auto end = all().end();
          for(; it != end; ++it) {
            ValueDataStream<bool> select_stream;
            if(select_stream.fail()) {
              err_log() << "Warning: Unable to apply criteria \"" << criteria
                        << "\" to " << traits<ObjType>::name << " " << it.name()
                        << "\n";
              continue;
            }
            select_stream << tformat(*it);
            it.is_selected() = select_stream.value();
          }
        }
      }
      catch(std::exception &e) {
        throw std::runtime_error(
          std::string("Failure to select using criteria \"") + criteria +
          "\" for " + traits<ObjType>::name + "\n"
          "    Reason:  " + e.what());
      }
    }

    /// \brief Set selected objects to value, if criteria true
    template<typename ObjType>
    void Selection<ObjType>::set(const DataFormatterDictionary<ObjType> &dict,
                                 const std::string &criteria,
                                 bool value) {
      try {
        if(criteria.size()) {
          DataFormatter<ObjType> tformat(dict.parse(criteria));
          auto it = all().begin();
          auto end = all().end();
          for(; it != end; ++it) {
            if(it.is_selected() == value) {
              continue;
            }
            ValueDataStream<bool> select_stream;
            if(select_stream.fail()) {
              err_log() << "Warning: Unable to apply criteria \"" << criteria
                        << "\" to " << traits<ObjType>::name << " " <<  it.name()
                        << "\n";
              continue;
            }

            select_stream << tformat(*it);
            if(select_stream.value()) {
              it.is_selected() = value;
            }
          }
        }
        else {
          auto it = all().begin();
          auto end = all().end();
          for(; it != end; ++it) {
            it.is_selected() = value;
          }
        }
      }
      catch(std::exception &e) {
        throw std::runtime_error(
          std::string("Failure to select using criteria \"") + criteria +
          "\" for " + traits<ObjType>::name + "\n"
          "    Reason:  " + e.what());
      }
    }

    /// \brief Read csv selection from stream
    template<typename ObjType>
    void Selection<ObjType>::read(std::istream &_input) {
      std::string tname_or_alias;
      bool tselect;
      _input >> std::ws;
      if(_input.peek() == '#') {
        _input.get();
        // discard first two columns (name_or_alias, selected)
        _input >> tname_or_alias;
        _input >> tname_or_alias;
        std::getline(_input, tname_or_alias, '\n');
        m_col_headers.clear();
        boost::split(m_col_headers, tname_or_alias, boost::is_any_of(" \t"), boost::token_compress_on);
        //_input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
      while(_input >> tname_or_alias >> tselect) {
        // skip unknown quietly (not sure what is best)
        // this typically arises in cases an object was deleted
        if(db().find(db().name(tname_or_alias)) == db().end()) {
          continue;
        }
        m_data[db().name(tname_or_alias)] = tselect;
        _input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }
    }

    /// \brief Read selection from JSON
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

          // for compatibility, include "configname" also
          if(it.name() == "name" || it.name() == "alias" || it.name() == "alias_or_name" || it.name() == "configname") {
            tname = db().name(it->get<std::string>());
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

    /// \brief Write selection to JSON
    template<typename ObjType>
    jsonParser &Selection<ObjType>::to_json(
      const DataFormatterDictionary<ObjType> &_dict,
      jsonParser &_json,
      bool only_selected) const {
      _json.put_array();

      DataFormatter<ObjType> tformat(
        alias_or_name<ObjType>(),
        datum_formatter_alias("selected", Selected<ObjType>(*this)));

      tformat.append(_dict.parse(m_col_headers));

      if(only_selected) {
        _json = tformat(selected().begin(), selected().end());
      }
      else {
        _json = tformat(all().begin(), all().end());
      }

      return _json;
    }

    /// \brief Print csv selection to stream
    template<typename ObjType>
    void Selection<ObjType>::print(
      const DataFormatterDictionary<ObjType> &_dict,
      std::ostream &_out,
      bool only_selected) const {

      DataFormatter<ObjType> tformat(
        alias_or_name<ObjType>(),
        datum_formatter_alias("selected", Selected<ObjType>(*this)));

      tformat.append(_dict.parse(m_col_headers));

      if(only_selected) {
        _out << tformat(selected().begin(), selected().end());
      }
      else {
        _out << tformat(all().begin(), all().end());
      }

    }

    /// \brief Write selection to file
    ///
    /// \returns 0 if successful, ERR_EXISTING_FILE if file already exists and !force
    template<typename ObjType>
    bool Selection<ObjType>::write(const DataFormatterDictionary<ObjType> &dict,
                                   bool force,
                                   const fs::path &_out_path,
                                   bool write_json,
                                   bool only_selected) const {

      fs::path out_path(_out_path);
      if(out_path.string() == "MASTER") {
        out_path = primclex().dir().template master_selection<ObjType>();
        force = true;
      }

      if(fs::exists(out_path) && !force) {
        err_log() << "File " << out_path
                  << " already exists. Use --force to force overwrite." << std::endl;
        return ERR_EXISTING_FILE;
      }

      if(write_json || out_path.extension() == ".json" || out_path.extension() == ".JSON") {
        jsonParser json;
        this->to_json(dict, json, only_selected);
        SafeOfstream sout;
        sout.open(out_path);
        json.print(sout.ofstream());
        sout.close();
      }
      else {
        SafeOfstream sout;
        sout.open(out_path);
        this->print(dict, sout.ofstream(), only_selected);
        sout.close();
      }
      return 0;
    }

    BOOST_PP_SEQ_FOR_EACH(INST_Selection, _, CASM_DB_TYPES)
  }
}
