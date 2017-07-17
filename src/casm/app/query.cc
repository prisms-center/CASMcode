#include "casm/app/query.hh"
#include "casm/app/DBInterface.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/DatabaseTypeDefs.hh"
#include "casm/database/Selection.hh"

namespace CASM {

  namespace Completer {

    QueryOption::QueryOption(): OptionHandlerBase("query") {};

    void QueryOption::initialize() {
      add_general_help_suboption();
      add_selection_suboption();
      add_db_type_suboption(traits<Configuration>::short_name, DB::types_short());
      //add_db_type_suboption(traits<Supercell>::short_name, DB::types_short());
      //add_db_type_suboption(traits<PrimPeriodicDiffTransOrbit>::short_name, DB::types_short());
      add_output_suboption();
      add_gzip_suboption();

      m_desc.add_options()
      ("columns,k", po::value<std::vector<std::string> >(&m_columns_vec)->multitoken()->zero_tokens()->value_name(ArgHandler::query()), "List of values you want printed as columns")
      ("json,j", "Print in JSON format (CSV otherwise, unless output extension is .json/.JSON)")
      ("verbatim,v", "Print exact properties specified, without prepending 'name' and 'selected' entries")
      ("all,a", "Print results all objects in input selection, whether or not they are selected.")
      ("no-header,n", "Print without header (CSV only)")
      ("alias", po::value<std::vector<std::string> >(&m_new_alias_vec)->multitoken(),
       "Create an alias for a query that will persist within this project. "
       "Ex: 'casm query --alias is_Ni_dilute = lt(atom_frac(Ni),0.10001)'")
      ("write-pos", "Write POS file for each configuration");

      return;
    }

    bool QueryOption::verbatim_flag() const {
      return vm().count("verbatim");
    }

    const std::vector<std::string> &QueryOption::columns_vec() const {
      return m_columns_vec;
    }

    const std::vector<std::string> &QueryOption::new_alias_vec() const {
      return m_new_alias_vec;
    }

  }


  // -- QueryCommandImplBase --------------------------------------------

  /// Defaults used if DataObject type doesn't matter or not given
  class QueryCommandImplBase : public Logging {
  public:

    QueryCommandImplBase(const QueryCommand &cmd);

    virtual ~QueryCommandImplBase() {}

    virtual int help() const;

    virtual int desc() const;

    virtual int run() const;

  protected:

    const QueryCommand &m_cmd;
  };


  QueryCommandImplBase::QueryCommandImplBase(const QueryCommand &cmd) :
    Logging(cmd),
    m_cmd(cmd) {}

  int QueryCommandImplBase::help() const {
    log() << std::endl << m_cmd.opt().desc() << std::endl;

    log() <<
          "Prints the properties for the objects currently selected in the MASTER\n"
          "selection or the specifed selection file via --selection (-c).\n\n"

          "Property values are output in column-separated (default) or JSON format.  By default, \n"
          "entries for 'name' and 'selected' values are included in the output. \n\n"

          "The type of objects acted on is specified via --type (-t).\n\n";

    m_cmd.print_names(log());
    return 0;
  }

  int QueryCommandImplBase::desc() const {
    return help();
  }

  int QueryCommandImplBase::run() const {
    throw CASM::runtime_error("Unknown error in 'casm query'.", ERR_UNKNOWN);
  }


  // -- QueryCommandImpl -----------------

  /// 'casm query' implementation, templated by type
  ///
  /// This:
  /// - holds a DB::InterfaceData object which stores dictionaries and selections
  /// - provides the implementation for 'help' (i.e. print allowed query commands)
  /// - provides the implementation for 'run' (i.e. perform query)
  ///
  template<typename DataObject>
  class QueryCommandImpl : public QueryCommandImplBase {
  public:
    QueryCommandImpl(const QueryCommand &cmd);

    int help() const override;

    int desc() const override;

    int run() const override;

  private:

    int _alias() const;

    int _write_pos() const;

    int _query() const;

    int _count(std::string s) const {
      return m_cmd.vm().count(s);
    }

    const fs::path &_selection_path() const {
      return m_cmd.opt().selection_path();
    }

    const std::vector<std::string> &_columns_vec() const {
      return m_cmd.opt().columns_vec();
    }

    std::vector<std::string> _all_columns() const;

    fs::path _output_path() const {
      return m_cmd.opt().output_path();
    }

    bool _write_json() const;

    bool _write_gz() const;

    DataFormatterDictionary<DataObject> &_dict() {
      return m_data.dict();
    }

    const DataFormatterDictionary<DataObject> &_dict() const {
      return m_data.dict();
    }

    std::string _sel_str() const {
      return m_data.sel_str();
    }

    double _sel_size() const {
      return m_data.sel_size();
    }

    DB::Selection<DataObject> &_sel(Index i = 0) {
      return m_data.sel(i);
    }

    const DB::Selection<DataObject> &_sel(Index i = 0) const {
      return m_data.sel(i);
    }

  private:

    // access dictionary and selections
    DB::InterfaceData<DataObject> m_data;
  };


  template<typename DataObject>
  QueryCommandImpl<DataObject>::QueryCommandImpl(const QueryCommand &cmd) :
    QueryCommandImplBase(cmd),
    m_data(cmd) {}

  template<typename DataObject>
  int QueryCommandImpl<DataObject>::help() const {

    if(!m_cmd.opt().help_opt_vec().size()) {
      return QueryCommandImplBase::help();
    }

    m_cmd.log() <<
                "Prints the properties for the objects currently selected in the MASTER\n"
                "selection or the specifed selection file.\n\n"

                "Property values are output in column-separated (default) or JSON format.  By default, \n"
                "entries for 'name' and 'selected' values are included in the output. \n\n";

    for(const std::string &str : m_cmd.opt().help_opt_vec()) {
      if(str.empty()) {
        continue;
      }

      if(str[0] == 'o') {
        m_cmd.log() << "Available operators for use within queries:" << std::endl;
        _dict().print_help(log(), BaseDatumFormatter<DataObject>::Operator);
      }
      else if(str[0] == 'p') {
        m_cmd.log() << "Available property tags are currently:" << std::endl;
        _dict().print_help(log(), BaseDatumFormatter<DataObject>::Property);
      }
      log() << std::endl;
    }
    log() << std::endl;
    return 0;
  }

  template<typename DataObject>
  int QueryCommandImpl<DataObject>::desc() const {
    return help();
  }

  template<typename DataObject>
  int QueryCommandImpl<DataObject>::run() const {

    if(_count("alias")) {
      return _alias();
    }
    else if(_count("write-pos")) {
      return _write_pos();
    }
    else if(_count("columns")) {
      return _query();
    }

    throw runtime_error("unknown error in 'casm query'", ERR_UNKNOWN);
  }

  template<typename DataObject>
  int QueryCommandImpl<DataObject>::_alias() const {
    ProjectSettings &set = m_cmd.primclex().settings();

    // get user input
    std::string new_alias_str;
    for(auto const &substr : m_cmd.opt().new_alias_vec()) {
      new_alias_str += substr;
    }

    // parse new_alias_str to create formatter
    auto it = std::find(new_alias_str.cbegin(), new_alias_str.cend(), '=');
    std::string alias_name = boost::trim_copy(std::string(new_alias_str.cbegin(), it));
    std::string alias_command = boost::trim_copy(std::string(++it, new_alias_str.cend()));

    try {
      set.query_handler<DataObject>().add_alias(alias_name, alias_command);
      set.commit();
      return 0;
    }
    catch(std::runtime_error &e) {
      err_log() << "Unable to learn alias\n"
                << "   \"" << alias_name << " = " << alias_command << "\"\n"
                << e.what() << std::endl;
      return ERR_UNKNOWN;
    }
  }

  template<typename DataObject>
  int QueryCommandImpl<DataObject>::_write_pos() const {
    for(const auto &obj : _sel().selected()) {
      obj.write_pos();
    }
    return 0;
  }

  namespace {
    bool _check_gz(fs::path p) {
      if(p.extension() == ".gz" || p.extension() == ".GZ") {
        return true;
      }
      return false;
    }

    bool _check_json(fs::path p) {
      if(p.extension() == ".json" || p.extension() == ".JSON") {
        return true;
      }
      return false;
    }
  }

  template<typename DataObject>
  int QueryCommandImpl<DataObject>::_query() const {

    // set output_stream: where the query results are written
    std::unique_ptr<std::ostream> uniq_fout;
    std::ostream &output_stream = make_ostream_if(_count("output"), log(), uniq_fout, _output_path(), _write_gz());
    output_stream << FormatFlag(output_stream).print_header(!_count("no-header"));

    // set status_stream: where query settings and PrimClex initialization messages are sent
    Log &status_log = (_output_path().string() == "STDOUT") ? err_log() : log();
    if(_output_path().string() == "STDOUT") {
      log().set_verbosity(0);
    }

    // Print info
    status_log << "Print:" << std::endl;
    for(int p = 0; p < _columns_vec().size(); p++) {
      status_log << "   - " << _columns_vec()[p] << std::endl;
    }
    if(_count("output")) {
      if(_output_path().string() == "STDOUT") {
        status_log << "to " << _output_path() << std::endl;
      }
      else {
        status_log << "to " << fs::absolute(_output_path()) << std::endl;
      }
    }
    status_log << std::endl;

    // construct formatter
    DataFormatter<DataObject> formatter = _dict().parse(_all_columns());

    auto begin = _count("all") ? _sel().all().begin() : _sel().selected().begin();
    auto end = _count("all") ? _sel().all().end() : _sel().selected().end();

    if(_write_json()) {
      jsonParser json;
      json = formatter(begin, end);
      output_stream << json;
    }
    else { // write CSV
      output_stream << formatter(begin, end);
    }

    if(!uniq_fout) {
      status_log << "\n   -Output printed to terminal, since no output file specified-\n";
    }

    status_log << "  DONE." << std::endl << std::endl;

    return 0;
  }

  template<typename DataObject>
  std::vector<std::string> QueryCommandImpl<DataObject>::_all_columns() const {
    // Construct DataFormatter
    std::vector<std::string> all_columns;
    if(!m_cmd.opt().verbatim_flag()) {
      all_columns.push_back("name");
      all_columns.push_back("selected");
    }
    all_columns.insert(all_columns.end(), _columns_vec().begin(), _columns_vec().end());
    return all_columns;
  }

  template<typename DataObject>
  bool QueryCommandImpl<DataObject>::_write_json() const {

    if(_count("json")) {
      return true;
    }
    // Checks for: X.json.gz / X.json / X.gz  (also accepts .JSON or .GZ)
    else if(_check_gz(_output_path())) {
      return _check_json(_output_path().stem());
    }
    else {
      return _check_json(_output_path());
    }
  }

  template<typename DataObject>
  bool QueryCommandImpl<DataObject>::_write_gz() const {

    if(m_cmd.opt().gzip_flag()) {
      return true;
    }
    // Checks for: X.json.gz / X.json / X.gz  (also accepts .JSON or .GZ)
    else {
      return _check_gz(_output_path());
    }
  }


  // -- class QueryCommand ----------------------------------------------------

  const std::string QueryCommand::name = "query";

  QueryCommand::QueryCommand(const CommandArgs &_args, Completer::QueryOption &_opt) :
    APICommand<Completer::QueryOption>(_args, _opt) {}

  QueryCommand::~QueryCommand() {}

  int QueryCommand::vm_count_check() const {
    if(!in_project()) {
      err_log().error("No casm project found");
      err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    std::string cmd;
    std::vector<std::string> allowed_cmd = {"alias", "columns", "write-pos"};

    Index num_cmd(0);
    for(const std::string &cmd_str : allowed_cmd) {
      if(vm().count(cmd_str)) {
        num_cmd++;
        cmd = cmd_str;
      }
    }

    if(num_cmd != 1) {
      err_log() << "Error in 'casm query'. Exactly one of the following must be used: "
                << allowed_cmd << std::endl;
      return ERR_INVALID_ARG;
    }
    return 0;
  }

  int QueryCommand::help() const {
    return impl().help();
  }

  int QueryCommand::desc() const {
    return impl().desc();
  }

  int QueryCommand::run() const {
    return impl().run();
  }

  QueryCommandImplBase &QueryCommand::impl() const {
    if(!m_impl) {
      if(in_project()) {
        if(!opt().db_type_opts().count(opt().db_type())) {
          std::stringstream msg;
          msg << "--type " << opt().db_type() << " is not allowed for 'casm " << name << "'.";
          print_names(err_log());
          throw CASM::runtime_error(msg.str(), ERR_INVALID_ARG);
        }

        DB::for_type_short(opt().db_type(), DB::ConstructImpl<QueryCommand>(m_impl, *this));
      }
      else {
        m_impl = notstd::make_unique<QueryCommandImplBase>(*this);
      }
    }
    return *m_impl;
  }

  void QueryCommand::print_names(std::ostream &sout) const {
    sout << "The allowed types are:\n";

    for(const auto &db_type : opt().db_type_opts()) {
      sout << "  " << db_type << std::endl;
    }
  }

}
