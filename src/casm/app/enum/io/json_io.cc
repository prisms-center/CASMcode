// #include "casm/app/enum/io/json_io.hh"
// #include "casm/casm_io/container/json_io.hh"
// #include "casm/casm_io/json/InputParser_impl.hh"
// #include "casm/misc/CASM_Eigen_math.hh"
//
// namespace CASM {
//
//   /// Parse a number or array value and return a vector
//   ///
//   /// \param parser Generic KwargsParser
//   /// \param attribute_name std::string containing the name of the attribute to be parsed
//   /// \param dimension If a number is given, return a constant vector of this dimension and the
//   ///        given value.
//   /// \param default_value Optional, pointer to default value to return if attribute does not
//   ///        exist or is null. If `default_value==nullptr`, the attribute is a required input and an
//   ///        error will be inserted.
//   ///
//   /// Options:
//   /// - If attribute is not present or null and `default_value!=nullptr`: return *default_value;
//   /// - If attribute is not present or null and `default_value==nullptr`: insert error message and return Eigen::VectorXd{}
//   /// - If attribute is number: return a Eigen::VectorXd::Constant(dimension, number)
//   /// - If attribute is an array: return array as a Eigen::VectorXd
//   /// - If attribute is any other type: insert error message and return Eigen::VectorXd{}
//   Eigen::VectorXd parse_vector_from_number_or_array(
//     KwargsParser &parser,
//     std::string attribute_name,
//     Index dimension,
//     Eigen::VectorXd const *default_value) {
//
//     Eigen::VectorXd result;
//     if(!parser.self.contains(attribute_name) || parser.self[attribute_name].is_null()) {
//       if(default_value != nullptr) {
//         return *default_value;
//       }
//       else {
//         parser.require(result, attribute_name);
//         return result;
//       }
//     }
//     if(parser.self[attribute_name].is_number()) {
//       return Eigen::VectorXd::Constant(dimension, parser.self[attribute_name].get<double>());
//     }
//     else if(parser.self[attribute_name].is_array()) {
//       parser.optional(result, attribute_name);
//       if(result.size() != dimension) {
//         std::stringstream msg;
//         std::string option_type = (default_value != nullptr) ? "(optional)" : "(required)";
//         msg << "Error: \"" << attribute_name << "\" " << option_type
//             << " must be a number or array of number. Expected array of size " << dimension << ".";
//         parser.error.insert(msg.str());
//       }
//       return result;
//     }
//     else {
//       std::stringstream msg;
//       std::string option_type = (default_value != nullptr) ? "(optional)" : "(required)";
//       msg << "Error: \"" << attribute_name << "\" " << option_type
//           << " must be a number or array of number.";
//       parser.error.insert(msg.str());
//       return result;
//     }
//   }
//
//   /// Parse DoF space axes from a JSON object
//   ///
//   /// Expected format: (this parses "axes" only)
//   ///
//   /// Full space case:
//   /// {
//   ///   "axes": {
//   ///     "q1": [q10, q11, q12, ...],            /// note: starts with "q1", not "q0"
//   ///     "q2": [q20, q21, q22, ...],
//   ///     "q3": [q30, q31, q32, ...],
//   ///     ...
//   ///   },
//   ///   "min": [q1min, q2min, q3min, ...],        /// optional, default is zeros vector
//   ///   "max": [q1max, q2max, q3max, ...],        /// required
//   ///   "increment": [q1inc, q2inc, q3inc, ...],  /// required
//   /// }
//   ///
//   /// Subspace case:
//   /// - if some qi (i in range [1, dof_space_dimension]) are missing, then use the subspace
//   ///   specified by the axes that are provided
//   /// - min.size() / max.size() / increment.size() must equal axes.size()
//   /// - example:
//   /// {
//   ///   "axes": {
//   ///     "q1": [q10, q11, q12, ...],
//   ///     "q2": [q20, q21, q22, ...],
//   ///     "q3": [q30, q31, q32, ...],
//   ///     "q5": [q50, q51, q52, ...],
//   ///   },
//   ///   "min": [q1min, q2min, q3min, q5min],        /// optional, default is zeros vector
//   ///   "max": [q1max, q2max, q3max, q5max],        /// required
//   ///   "increment": [q1inc, q2inc, q3inc, q5inc],  /// required
//   /// }
//   ///
//   /// notes:
//   /// - if total dof_space_dimension is double digit, then use "q01", "q02", ... etc.
//   /// - if total dof_space_dimension is triple digit, then use "q001", "q002", ... etc.
//   ///
//   void parse_axes_from_object(KwargsParser &parser,
//                               Eigen::MatrixXd &axes,
//                               Index dof_space_dimension) {
//
//     Eigen::MatrixXd inaxes = Eigen::MatrixXd::Zero(dof_space_dimension, dof_space_dimension);
//
//     std::set<Index> found;
//     for(Index i = 0; i < dof_space_dimension; ++i) {
//       std::string axis_name = "q" + to_sequential_string(i + 1, dof_space_dimension);
//       fs::path location = fs::path {"axes"} / axis_name;
//       auto value_ptr = parser.optional_at<Eigen::VectorXd>(location);
//
//       if(value_ptr != nullptr) {
//         found.insert(i);
//         Eigen::VectorXd const &value = *value_ptr;
//         if(value.size() != dof_space_dimension) {
//           std::stringstream msg;
//           msg << "Error reading axis vector \"" << axis_name
//               << "\": expected size=" << dof_space_dimension
//               << " found size=" << value.size();
//           parser.error.insert(msg.str());
//         }
//         else {
//           inaxes.col(found.size()) = value;
//         }
//       }
//     }
//
//     Index subspace_dimension = found.size();
//     axes = inaxes.leftCols(subspace_dimension);
//
//   }
//
//   /// Read "axes" from a row-vector JSON matrix, store in `axes` argument as a column vector matrix
//   void parse_axes_from_array(KwargsParser &parser,
//                              Eigen::MatrixXd &axes,
//                              Index dof_space_dimension) {
//     Eigen::MatrixXd row_vector_axes;
//     parser.require(row_vector_axes, "axes");
//     axes = row_vector_axes.transpose();
//
//     // check axes dimensions:
//     if(axes.rows() != dof_space_dimension) {
//       // Note: message "columns" refers to JSON input axes, transpose of params.axes
//       std::stringstream msg;
//       msg << "Number of columns of \"axes\" must be equal to site DoF space dimension ("
//           << dof_space_dimension << "). Size as parsed: " << axes.rows();
//       parser.error.insert(msg.str());
//     }
//     if(axes.cols() > dof_space_dimension) {
//       // Note: message "rows" refers to JSON input axes, transpose of params.axes
//       std::stringstream msg;
//       msg << "Number of coordinate axes (number of rows of \"axes\") must be less than or equal to "
//           "site DoF space dimension (" << dof_space_dimension << "). Number of axes parsed: "
//           << axes.cols();
//       parser.error.insert(msg.str());
//     }
//   }
//
//
//
//   /// Parse "axes", "min", "max", and "increment" for continuous DoF enumeration
//   ///
//   /// Accepts a JSON object of axes or a row vector matrix of axes. May be rank deficient. For JSON
//   /// object format see `parse_axes_from_object`
//   ///
//   /// \param parser JSON to be parsed
//   /// \param axes Set to be column vector matrix of axes
//   /// \param min_val Set to be minimum counter value. Must be same size as axes.cols().
//   /// \param max_val Set to be maximum counter value. Must be same size as axes.cols().
//   /// \param inc_val Set to be counter increment value. Must be same size as axes.cols().
//   /// \param dof_space_dimension Full DoF space dimension. If no `"axes"` in the JSON, the default
//   ///        axes are identity matrix of size dof_space_dimension. Errors are inserted into parser
//   ///        if axes vectors are not of length equal to dof_space_dimension.
//   ///
//   void parse_dof_space_axes(KwargsParser &parser,
//                             Eigen::MatrixXd &axes,
//                             Eigen::VectorXd &min_val,
//                             Eigen::VectorXd &max_val,
//                             Eigen::VectorXd &inc_val,
//                             Index dof_space_dimension) {
//
//     // "axes" -> axes
//     if(!parser.self.contains("axes")) {
//       axes = Eigen::MatrixXd::Identity(dof_space_dimension, dof_space_dimension);
//     }
//     else if(parser.self.contains("axes") && parser.self["axes"].is_obj()) {
//       parse_axes_from_object(parser, axes, dof_space_dimension);
//     }
//     else if(parser.self.contains("axes") && parser.self["axes"].is_array()) {
//       parse_axes_from_array(parser, axes, dof_space_dimension);
//     }
//     else {
//       std::stringstream msg;
//       msg << "The \"axes\" must be a JSON object of axis vectors (named \"q1\", \"q2\", etc.), or a row-vector matrix";
//       parser.error.insert(msg.str());
//     }
//
//     // "min" -> min_val  (number or array of number, else zeros vector)
//     Eigen::VectorXd default_value = Eigen::VectorXd::Zero(axes.cols());
//     min_val = parse_vector_from_number_or_array(parser, "min", axes.cols(), &default_value);
//
//     // "max" -> max_val (number or array of number, required)
//     max_val = parse_vector_from_number_or_array(parser, "max", axes.cols());
//
//     // "increment" -> inc_val (number or array of number, required)
//     inc_val = parse_vector_from_number_or_array(parser, "increment", axes.cols());
//   }
//
//
// }
