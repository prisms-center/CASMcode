#include "casm/app/QueryHandler.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/DataFormatterTools.hh"


namespace CASM {

  template<typename DataObject>
  QueryHandler<DataObject>::QueryHandler(const ProjectSettings &set) :
    m_set(&set),
    m_dict(make_dictionary<DataObject>()) {

    // add plugins
    load_query_plugins(
      set,
      std::inserter(m_dict, m_dict.end()),
      std::inserter(m_lib, m_lib.end()));
  };

  template<typename DataObject>
  QueryHandler<DataObject>::~QueryHandler() {
    // order of deletion matters
    m_dict.clear();
    m_lib.clear();
  }

  template<typename DataObject>
  DataFormatterDictionary<DataObject> &QueryHandler<DataObject>::dict() {
    return m_dict;
  }

  template<typename DataObject>
  const DataFormatterDictionary<DataObject> &QueryHandler<DataObject>::dict() const {
    return m_dict;
  }

  /// \brief Set the selection to be used for the 'selected' column
  ///
  /// - ToDo: generalize ConfigIO::Selected
  template<typename DataObject>
  void QueryHandler<DataObject>::set_selected(const typename QueryTraits<DataObject>::Selected &selection) {
    if(m_dict.find("selected") != m_dict.end()) {
      m_dict.erase("selected");
    }
    m_dict.insert(
      datum_formatter_alias(
        "selected",
        selection,
        "Returns true if configuration is specified in the input selection"
      )
    );
  }

  /// \brief Set the selection to be used for the 'selected' column
  ///
  /// - ToDo: generalize ConstConfigSelection
  template<typename DataObject>
  void QueryHandler<DataObject>::set_selected(const typename QueryTraits<DataObject>::Selection &selection) {
    typedef typename QueryTraits<DataObject>::Selected selected_type;
    set_selected(selected_type(selection));
  }

  /// \brief Add user-defined query alias
  ///
  /// - Aliases are added to the dictionary and ProjectSettings, but not saved
  ///   to file
  template<typename DataObject>
  void QueryHandler<DataObject>::add_alias(const std::string &alias_name, const std::string &alias_command) {

    auto new_formatter = datum_formatter_alias<DataObject>(alias_name, alias_command, m_dict);
    auto key = m_dict.key(new_formatter);

    // if not in dictionary (includes operator dictionary), add
    if(m_dict.find(key) == m_dict.end()) {
      m_dict.insert(new_formatter);
    }
    // if a user-created alias, over-write with message
    else if(aliases().find(alias_name) != aliases().end()) {
      m_set->err_log() << "WARNING: I already know '" << alias_name << "' as:\n"
                       << "             " << _aliases()[alias_name] << "\n"
                       << "         I will forget it and learn '" << alias_name << "' as:\n"
                       << "             " << alias_command << std::endl;
      m_dict.insert(new_formatter);
    }
    // else do not add, throw error
    else {
      std::stringstream ss;
      ss << "Error: Attempted to over-write standard CASM query name with user alias.\n";
      throw std::runtime_error(ss.str());
    }

    // save alias
    _aliases()[alias_name] = alias_command;

  }

  /// \brief Load enumerator plugins from a CASM project
  template<typename DataFormatterDictInserter, typename RuntimeLibInserter>
  std::pair<DataFormatterDictInserter, RuntimeLibInserter>
  load_query_plugins(
    const ProjectSettings &set,
    DataFormatterDictInserter dict_it,
    RuntimeLibInserter lib_it) {

    typedef typename DataFormatterDictInserter::container_type dict_type;
    typedef typename dict_type::DataObject DataObject;
    typedef typename dict_type::DatumFormatterType formatter_type;
    typedef formatter_type *formatter_type_ptr;
    typedef formatter_type_ptr(signature)();

    const DirectoryStructure &dir = set.dir();

    if(dir.root_dir().empty()) {
      return std::make_pair(dict_it, lib_it);
    }

    if(fs::is_directory(dir.query_plugins<DataObject>())) {

      // loop over custom query files *.cc
      for(auto &entry : boost::make_iterator_range(fs::directory_iterator(dir.query_plugins<DataObject>()), {})) {

        fs::path p = entry.path();
        std::string p_s = p.string();
        auto p_size = p_s.size();

        if(fs::is_regular_file(p) && p_s.compare(p_size - 3, p_size, ".cc") == 0) {

          fs::path f = p.filename();
          std::string f_s = f.string();
          auto f_size = f_s.size();

          std::string msg = "compiling new custom query: " + f_s.substr(0, f_size - 3);

          // '-L$CASM_PREFIX/.libs' is a hack so 'make check' works
          auto lib_ptr = std::make_shared<RuntimeLibrary>(
                           p_s.substr(0, p_size - 3),
                           set.compile_options() + " " + include_path(dir.query_plugins<DataObject>()),
                           set.so_options() + " -lcasm ",
                           msg,
                           set);

          auto make_formatter = lib_ptr->template get_function<signature>(
            "make_" + f_s.substr(0, f_size - 3) + "_formatter");

          std::unique_ptr<formatter_type> ptr(make_formatter());

          // will clone on insert
          *dict_it++ = *ptr;
          *lib_it++ = std::make_pair(ptr->name(), lib_ptr);
        }
      }
    }

    return std::make_pair(dict_it, lib_it);
  }

}
