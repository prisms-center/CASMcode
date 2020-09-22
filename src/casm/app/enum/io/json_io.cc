#include "casm/app/enum/io/json_io.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  /// Parse a number or array value and return a vector
  ///
  /// \param parser Generic KwargsParser
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
    KwargsParser &parser,
    std::string attribute_name,
    Index dimension,
    Eigen::VectorXd const *default_value) {

    Eigen::VectorXd result;
    if(!parser.self.contains(attribute_name) || parser.self[attribute_name].is_null()) {
      if(default_value != nullptr) {
        return *default_value;
      }
      else {
        parser.require(result, attribute_name);
        return result;
      }
    }
    if(parser.self[attribute_name].is_number()) {
      return Eigen::VectorXd::Constant(dimension, parser.self[attribute_name].get<double>());
    }
    else if(parser.self[attribute_name].is_array()) {
      parser.optional(result, attribute_name);
      if(result.size() != dimension) {
        std::stringstream msg;
        std::string option_type = (default_value != nullptr) ? "(optional)" : "(required)";
        msg << "Error: \"" << attribute_name << "\" " << option_type
            << " must be a number or array of number. Expected array of size " << dimension << ".";
        parser.error.insert(msg.str());
      }
      return result;
    }
    else {
      std::stringstream msg;
      std::string option_type = (default_value != nullptr) ? "(optional)" : "(required)";
      msg << "Error: \"" << attribute_name << "\" " << option_type
          << " must be a number or array of number.";
      parser.error.insert(msg.str());
      return result;
    }
  }

  // if(result.size() != subspace_dimension) {
  //   std::stringstream msg;
  //   msg << "Error: \"" << attribute_name << "\" array size must match dimension of \"axes\" rows");
  //   parser.error.insert(msg.str());
  // }


}
