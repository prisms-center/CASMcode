// #include "casm/app/enum/io/stream_io.hh"
// #include "casm/casm_io/Log.hh"
// #include "casm/casm_io/json/jsonParser.hh"
// #include "casm/enumerator/DoFSpace.hh"
// #include "casm/symmetry/SymRepTools.hh"
// #include "casm/symmetry/io/json/SymRepTools.hh"
//
// namespace CASM {
//
//   /// Print DoFSpace information
//   ///
//   /// \param log Log to print to
//   /// \param name Enumeration state name
//   /// \param dof_space DoFSpace to print information for.
//   /// \param sym_axes If true, calculate and print symmetry-adapted axes. If false, print only the
//   ///        axis glossary
//   /// \param calc_wedges If this and sym_axes are true, calculate and print SubWedge
//   ///
//   /// Notes:
//   /// - If !sym_axes, only prints axis glossary
//   /// - If sym_axes==true, construct and print VectorSpaceSymReport
//   ///   - TODO: Currently this prints the VectorSpaceSymReport using the JSON method. Could be
//   ///     specialized for nicer formatting.
//   template<typename PermuteIteratorIt>
//   void print_dof_space(Log &log,
//                        std::string const &name,
//                        DoFSpace const &dof_space,
//                        PermuteIteratorIt permute_begin,
//                        PermuteIteratorIt permute_end,
//                        bool sym_axes,
//                        bool calc_wedges) {
//     if(!sym_axes) {
//       log.begin(std::string("DoF Space Axes: ") + name);
//       log << "Note: column and site indexing begin with 1" << std::endl;
//       auto axis_glossary = make_axis_glossary(dof_space.dof_key,
//                                               dof_space.config_region.configuration(),
//                                               dof_space.config_region.sites());
//       for(Index column_index = 0; column_index != axis_glossary.size(); ++column_index) {
//         log << "column: " << column_index + 1 << " DoF: " << axis_glossary[column_index] << std::endl;
//       }
//     }
//     else {
//       log.begin(std::string("DoF Vector Space Symmetry Report: ") + name);
//       log << "For large spaces this may be slow... total dimension = "
//           << dof_space.dof_subspace.cols() << std::endl;
//
//       VectorSpaceSymReport sym_report = vector_space_sym_report(dof_space,
//                                                                 permute_begin,
//                                                                 permute_end,
//                                                                 calc_wedges);
//       // TODO: specialized print to stream instead of JSON to stream?
//       log << jsonParser {sym_report} << std::endl << std::endl;
//       log.end_section();
//     }
//   }
// }
