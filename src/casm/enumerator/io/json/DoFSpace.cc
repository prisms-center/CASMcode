#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"


namespace CASM {

  jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json, std::string name) {
    json["dof"] = dofspace.dof_key;
    {
      jsonParser &cjson = json["initial_configuration"];

      cjson["identifier"] = name;

      SimpleStructure sstruc = make_simple_structure(dofspace.config_region.configuration());

      cjson["lattice_vectors"] = sstruc.lat_column_mat.transpose();

      jsonParser &sjson = cjson["sites"];

      for(Index i = 0; i < sstruc.mol_info.size(); ++i) {
        to_json_array(sstruc.mol_info.cart_coord(i), sjson[to_sequential_string(i + 1, sstruc.mol_info.size())][sstruc.mol_info.names[i]]);
      }

      cjson["selected_sites"].put_array();
      for(Index s : dofspace.config_region.sites()) {
        cjson["selected_sites"].push_back(s + 1);
      }

    }

    json["glossary"] = make_axis_glossary(dofspace.dof_key,
                                          dofspace.config_region.configuration(),
                                          dofspace.config_region.sites());

    json["dof_subspace"] = dofspace.dof_subspace.transpose();

    return json;
  }

  /// Parse a number or array value and return a vector
  ///
  /// \param parser AxesCounterParams parser
  /// \param attribute_name std::string containing the name of the attribute to be parsed
  /// \param dimension If a number is given, return a constant vector of this dimension and the
  ///        given value.
  /// \param default_value Optional, pointer to default value to return if attribute does not
  ///        exist or is null. If `default_value==nullptr`, the attribute is a required input and an
  ///        error will be inserted.
  ///
  /// Options:
  /// - If attribute is not present or null and `default_value!=nullptr`: return *default_value;
  /// - If attribute is not present or null and `default_value==nullptr`: insert error message and return Eigen::VectorXd{}
  /// - If attribute is number: return a Eigen::VectorXd::Constant(dimension, number)
  /// - If attribute is an array: return array as a Eigen::VectorXd
  /// - If attribute is any other type: insert error message and return Eigen::VectorXd{}
  Eigen::VectorXd parse_vector_from_number_or_array(
    InputParser<AxesCounterParams> &parser,
    std::string attribute_name,
    Index dimension,
    Eigen::VectorXd const *default_value = nullptr) {

    std::string option_type = (default_value != nullptr) ? "(optional)" : "(required)";
    Eigen::VectorXd no_result;
    if(!parser.self.contains(attribute_name) || parser.self[attribute_name].is_null()) {
      if(default_value != nullptr) {
        return *default_value;
      }
      else {
        std::stringstream msg;
        msg << "Error: missing required option '" << attribute_name << "'. "
            << "Must be a number or array of numbers of size " << dimension << ".";
        parser.insert_error(attribute_name, msg.str());
        return no_result;
      }
    }
    if(parser.self[attribute_name].is_number()) {
      return Eigen::VectorXd::Constant(dimension, parser.self[attribute_name].get<double>());
    }
    else if(parser.self[attribute_name].is_array()) {
      auto subparser = parser.subparse<Eigen::VectorXd>(attribute_name);
      if(subparser->value == nullptr) {
        return no_result;
      }
      if(subparser->value->size() != dimension) {
        std::stringstream msg;
        msg << "Error: '" << attribute_name << "' " << option_type
            << " must be a number or array of numbers of size " << dimension << ".";
        subparser->error.insert(msg.str());
      }
      return *subparser->value;
    }
    else {
      std::stringstream msg;
      msg << "Error: '" << attribute_name << "' " << option_type
          << " must be a number or array of numbers of size " << dimension << ".";
      parser.insert_error(attribute_name, msg.str());
      return no_result;
    }
  }

  /// Parse DoF space axes from a JSON object
  ///
  /// Expected format: (this parses "axes" only)
  ///
  /// Full space case:
  /// {
  ///   "axes": {
  ///     "q1": [q10, q11, q12, ...],            /// note: starts with "q1", not "q0"
  ///     "q2": [q20, q21, q22, ...],
  ///     "q3": [q30, q31, q32, ...],
  ///     ...
  ///   },
  ///   "min": [q1min, q2min, q3min, ...],        /// optional, default is zeros vector
  ///   "max": [q1max, q2max, q3max, ...],        /// required
  ///   "increment": [q1inc, q2inc, q3inc, ...],  /// required
  /// }
  ///
  /// Subspace case:
  /// - if some qi (i in range [1, dof_space_dimension]) are missing, then use the subspace
  ///   specified by the axes that are provided
  /// - min.size() / max.size() / increment.size() must equal axes.size()
  /// - example:
  /// {
  ///   "axes": {
  ///     "q1": [q10, q11, q12, ...],
  ///     "q2": [q20, q21, q22, ...],
  ///     "q3": [q30, q31, q32, ...],
  ///     "q5": [q50, q51, q52, ...],
  ///   },
  ///   "min": [q1min, q2min, q3min, q5min],        /// optional, default is zeros vector
  ///   "max": [q1max, q2max, q3max, q5max],        /// required
  ///   "increment": [q1inc, q2inc, q3inc, q5inc],  /// required
  /// }
  ///
  /// notes:
  /// - if total dof_space_dimension is double digit, then use "q01", "q02", ... etc.
  /// - if total dof_space_dimension is triple digit, then use "q001", "q002", ... etc.
  ///
  void parse_axes_from_object(InputParser<AxesCounterParams> &parser,
                              Eigen::MatrixXd &axes,
                              Index dof_space_dimension) {

    Eigen::MatrixXd inaxes = Eigen::MatrixXd::Zero(dof_space_dimension, dof_space_dimension);

    std::set<Index> found;
    for(Index i = 0; i < dof_space_dimension; ++i) {
      std::string axis_name = "q" + to_sequential_string(i + 1, dof_space_dimension);
      auto subparser = parser.subparse_if<Eigen::VectorXd>(fs::path {"axes"} / axis_name);
      if(subparser->value != nullptr) {
        found.insert(i);
        if(subparser->value->size() != dof_space_dimension) {
          std::stringstream msg;
          msg << "Error reading axis vector '" << axis_name
              << "': expected size=" << dof_space_dimension
              << " found size=" << subparser->value->size();
          subparser->error.insert(msg.str());
        }
        else {
          inaxes.col(found.size()) = *subparser->value;
        }
      }
    }

    Index subspace_dimension = found.size();
    axes = inaxes.leftCols(subspace_dimension);

  }

  /// Read "axes" from a row-vector JSON matrix, store in `axes` argument as a column vector matrix
  void parse_axes_from_array(InputParser<AxesCounterParams> &parser,
                             Eigen::MatrixXd &axes,
                             Index dof_space_dimension) {
    // Eigen::MatrixXd row_vector_axes;
    auto subparser = parser.subparse<Eigen::MatrixXd>("axes");
    if(subparser->value == nullptr) {
      return;
    }
    axes = subparser->value->transpose();

    // check axes dimensions:
    if(axes.rows() != dof_space_dimension) {
      // Note: message "columns" refers to JSON input axes, transpose of params.axes
      std::stringstream msg;
      msg << "Number of columns of 'axes' must be equal to site DoF space dimension ("
          << dof_space_dimension << "). Size as parsed: " << axes.rows();
      subparser->error.insert(msg.str());
    }
    if(axes.cols() > dof_space_dimension) {
      // Note: message "rows" refers to JSON input axes, transpose of params.axes
      std::stringstream msg;
      msg << "Number of coordinate axes (number of rows of 'axes') must be less than or equal to "
          "site DoF space dimension (" << dof_space_dimension << "). Number of axes parsed: "
          << axes.cols();
      subparser->error.insert(msg.str());
    }
  }


  /// Parse DoF space "axes"
  ///
  /// Accepts a JSON object of axes or a row vector matrix of axes. May be rank deficient.
  ///
  /// \param parser JSON to be parsed
  /// \param axes Set to be column vector matrix of axes
  /// \param dof_space_dimension Full DoF space dimension. If no `"axes"` in the JSON, the default
  ///        axes are identity matrix of size dof_space_dimension. Errors are inserted into parser
  ///        if axes vectors are not of length equal to dof_space_dimension.
  ///
  /// TODO: document "axes" here
  ///
  void parse_dof_space_axes(InputParser<AxesCounterParams> &parser,
                            Eigen::MatrixXd &axes,
                            Index dof_space_dimension) {
    // "axes" -> axes
    if(!parser.self.contains("axes")) {
      axes = Eigen::MatrixXd::Identity(dof_space_dimension, dof_space_dimension);
    }
    else if(parser.self.contains("axes") && parser.self["axes"].is_obj()) {
      parse_axes_from_object(parser, axes, dof_space_dimension);
    }
    else if(parser.self.contains("axes") && parser.self["axes"].is_array()) {
      parse_axes_from_array(parser, axes, dof_space_dimension);
    }
    else {
      std::stringstream msg;
      msg << "The 'axes' must be a JSON object of axis vectors (named 'q1', 'q2', etc.), or a row-vector matrix";
      parser.error.insert(msg.str());
    }
  }


  /// Parse DoF space counter values from "min", "max", and "increment"
  ///
  /// \param parser JSON to be parsed
  /// \param dof_subspace DoFSpace subspace (column vector matrix)
  /// \param min_val Set to be minimum counter value. Must be same size as axes.cols().
  /// \param max_val Set to be maximum counter value. Must be same size as axes.cols().
  /// \param inc_val Set to be counter increment value. Must be same size as axes.cols().
  ///
  /// Errors are inserted into parser if counter vector sizes are not equal to dof_subspace.cols().
  ///
  /// TODO: document "min", "max", "increment" here
  ///
  void parse_dof_space_counter(InputParser<AxesCounterParams> &parser,
                               Eigen::MatrixXd const &dof_subspace,
                               Eigen::VectorXd &min_val,
                               Eigen::VectorXd &max_val,
                               Eigen::VectorXd &inc_val) {

    // "min" -> min_val  (number or array of number, else zeros vector)
    Eigen::VectorXd default_value = Eigen::VectorXd::Zero(dof_subspace.cols());
    min_val = parse_vector_from_number_or_array(parser, "min", dof_subspace.cols(), &default_value);

    // "max" -> max_val (number or array of number, required)
    max_val = parse_vector_from_number_or_array(parser, "max", dof_subspace.cols());

    // "increment" -> inc_val (number or array of number, required)
    inc_val = parse_vector_from_number_or_array(parser, "increment", dof_subspace.cols());
  }

  /// Parse DoFSpace subspace from  "axes" and normal coordinate grid counter from "min", "max", and "increment"
  void parse(InputParser<AxesCounterParams> &parser,
             Index dof_space_dimension) {
    AxesCounterParams params;
    parse_dof_space_axes(parser, params.axes, dof_space_dimension);
    parse_dof_space_counter(parser, params.axes, params.min_val, params.max_val, params.inc_val);
    if(parser.valid()) {
      parser.value = notstd::make_unique<AxesCounterParams>(params);
    }
  }

}
