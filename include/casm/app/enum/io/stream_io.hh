#ifndef CASM_app_enum_stream_io
#define CASM_app_enum_stream_io

#include <iostream>
#include <string>

namespace CASM {

  class ConfigEnumInput;
  class DoFSpace;
  class Log;

  /// Print DoFSpace information
  template<typename PermuteIteratorIt>
  void print_dof_space(Log &log,
                       DoFSpace const &dof_space,
                       std::string const &identifier,
                       ConfigEnumInput const &input_state,
                       PermuteIteratorIt permute_begin,
                       PermuteIteratorIt permute_end,
                       bool sym_axes,
                       bool calc_wedges);

  template<typename NamedInitialStatesType>
  void print_initial_states(Log &log, NamedInitialStatesType const &named_initial_states);

  void print_options(Log &log, ConfigEnumOptions const &options);
}

#endif
