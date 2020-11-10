#ifndef CASM_app_enum_stream_io
#define CASM_app_enum_stream_io

#include <iostream>
#include <string>

namespace CASM {

  struct DoFSpace;
  class Log;

  /// Print DoFSpace information
  template<typename PermuteIteratorIt>
  void print_dof_space(Log &log,
                       std::string const &name,
                       DoFSpace const &dof_space,
                       PermuteIteratorIt permute_begin,
                       PermuteIteratorIt permute_end,
                       bool sym_axes,
                       bool calc_wedges);

  template<typename NamedInitialStatesType>
  void print_initial_states(Log &log, NamedInitialStatesType const &named_initial_states);

  void print_options(Log &log, ConfigEnumOptions const &options);
}

#endif
