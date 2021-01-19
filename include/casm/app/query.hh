#ifndef CASM_query
#define CASM_query

#include "casm/app/APICommand.hh"
#include "casm/app/DBInterface.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {
/*namespace Completer {

  class QueryOption : public OptionHandlerBase {

  public:

    using OptionHandlerBase::selection_path;
    using OptionHandlerBase::output_path;
    using OptionHandlerBase::gzip_flag;
    using OptionHandlerBase::help_opt_vec;
    using OptionHandlerBase::db_type;
    using OptionHandlerBase::db_type_opts;

    QueryOption();

    bool verbatim_flag() const;

    const std::vector<std::string> &columns_vec() const;

    const std::vector<std::string> &new_alias_vec() const;

  private:

    void initialize() override;

    std::vector<std::string> m_columns_vec;

    std::vector<std::string> m_new_alias_vec;

  };

}*/
}

namespace CASM {

class QueryCommandImplBase;
template <typename T>
class QueryCommandImpl;

/// 'casm query' command
///
/// This class wraps a pointer to an implementation class
/// QueryCommandImpl<DataObject> which does most of the work. The implementation
/// classes are defined in the source file.
///
class QueryCommand : public APICommand<Completer::QueryOption> {
 public:
  using ImplBase = QueryCommandImplBase;
  template <typename ConfigType>
  using Impl = QueryCommandImpl<ConfigType>;

  static const std::string name;

  QueryCommand(const CommandArgs &_args, Completer::QueryOption &_opt);

  ~QueryCommand();

  int vm_count_check() const override;

  int help() const override;

  int desc() const override;

  int run() const override;

  // -- custom --

  QueryCommandImplBase &impl() const;

  void print_names(std::ostream &sout) const;

 private:
  mutable std::unique_ptr<QueryCommandImplBase> m_impl;
};

}  // namespace CASM

#endif
