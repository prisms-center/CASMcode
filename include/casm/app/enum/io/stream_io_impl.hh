#ifndef CASM_app_enum_stream_io_impl
#define CASM_app_enum_stream_io_impl

#include "casm/app/enum/io/stream_io.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

namespace CASM {

/// Print DoFSpace information
///
/// \param log Log to print to
/// \param name Enumeration state name
/// \param dof_space DoFSpace to print information for.
/// \param sym_axes If true, calculate and print symmetry-adapted axes. If
/// false, print only the
///        axis glossary
/// \param calc_wedges If this and sym_axes are true, calculate and print
/// SubWedge
///
/// Notes:
/// - If !sym_axes, only prints axis glossary
/// - If sym_axes==true, construct and print VectorSpaceSymReport
///   - TODO: Currently this prints the VectorSpaceSymReport using the JSON
///   method. Could be
///     specialized for nicer formatting.
template <typename PermuteIteratorIt>
void print_dof_space(Log &log, DoFSpace const &dof_space,
                     std::string const &identifier,
                     ConfigEnumInput const &input_state,
                     PermuteIteratorIt permute_begin,
                     PermuteIteratorIt permute_end, bool sym_axes,
                     bool calc_wedges) {
  log.subsection();
  log.begin_section();
  log << std::endl;
  if (!sym_axes) {
    log.begin(std::string("DoF space: ") + identifier);
    jsonParser json;
    to_json(dof_space, json);
    log << json << std::endl << std::endl;
  } else {
    log.begin(std::string("DoF space and symmetry report: ") + identifier);
    log << "For large spaces this may be slow... basis shape = ["
        << dof_space.basis().rows() << ", " << dof_space.basis().cols() << "]"
        << std::endl;

    VectorSpaceSymReport sym_report = vector_space_sym_report(
        dof_space, input_state, permute_begin, permute_end, calc_wedges);

    log << "Symmetry report complete" << std::endl;

    // TODO: specialized print to stream instead of JSON to stream?
    jsonParser json;
    to_json(dof_space, json, identifier, input_state, sym_report);
    log << json << std::endl << std::endl;
  }
  log.end_section();
}

template <typename NamedInitialStatesType>
void print_initial_states(Log &log,
                          NamedInitialStatesType const &named_initial_states) {
  log.indent() << "# of initial enumeration states: "
               << named_initial_states.size() << std::endl;
  log.subsection().begin_section<Log::verbose>();
  log.indent() << "initial enumeration states:" << std::endl;
  log.increase_indent();
  for (auto const &named_initial_state : named_initial_states) {
    log.indent() << named_initial_state.first << std::endl;
  }
  log.decrease_indent();
  log.end_section();
}
}  // namespace CASM

#endif
