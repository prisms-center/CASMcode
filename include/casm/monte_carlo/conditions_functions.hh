#ifndef CASM_Monte_conditions_functions
#define CASM_Monte_conditions_functions

#include <optional>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace Monte {
namespace conditions {

/// \brief If lhs and rhs have values, do `*lhs += *rhs`
template <typename MatrixType>
void incr(std::optional<MatrixType> &lhs,
          std::optional<MatrixType> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error incrementing conditions: existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  if (lhs->cols() != rhs->cols() || lhs->rows() != rhs->rows()) {
    throw std::runtime_error("Error incrementing conditions: shape mismatch");
  }
  *lhs += *rhs;
}

/// \brief If lhs and rhs have values, do `*lhs -= *rhs`
template <typename MatrixType>
void decr(std::optional<MatrixType> &lhs,
          std::optional<MatrixType> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error decrementing conditions: existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  if (lhs->cols() != rhs->cols() || lhs->rows() != rhs->rows()) {
    throw std::runtime_error("Error decrementing conditions: shape mismatch");
  }
  *lhs -= *rhs;
}

/// \brief If lhs and rhs have values, do `almost_equal(*lhs, *rhs, tol)`
template <typename MatrixType>
bool almost_equal(std::optional<MatrixType> const &lhs,
                  std::optional<MatrixType> const &rhs, double tol) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error(
        "Error checking conditions equality: existence mismatch");
  }
  if (!lhs.has_value()) {
    return true;
  }
  if (lhs->cols() != rhs->cols() || lhs->rows() != rhs->rows()) {
    throw std::runtime_error(
        "Error checking conditions equality: shape mismatch");
  }
  return almost_equal(*lhs, *rhs, tol);
}

/// \brief If lhs and rhs have values, update max_division if round((*lhs)(i) /
/// (*rhs)(i)) > max_division for any i, and ignoring divide by 0
template <typename MatrixType>
void find_max_division(int &max_division, std::optional<MatrixType> const &lhs,
                       std::optional<MatrixType> const &rhs) {
  if (lhs.has_value() != rhs.has_value()) {
    throw std::runtime_error("Error in find_max_division: existence mismatch");
  }
  if (!lhs.has_value()) {
    return;
  }
  if (lhs->cols() != rhs->cols() || lhs->rows() != rhs->rows()) {
    throw std::runtime_error("Error in find_max_division: shape mismatch");
  }
  for (Index i = 0; i < lhs->size(); i++) {
    if (almost_zero((*rhs)(i))) {
      continue;
    }
    int temp_division = round((*lhs)(i) / (*rhs)(i));
    if (temp_division > max_division) {
      max_division = temp_division;
    }
  }
}

/// \brief If lhs and rhs have values, update max_division if round((*lhs)(i) /
/// (*rhs)(i)) > max_division for any i, and ignoring divide by 0
template <typename MatrixType>
void print_param(std::ostream &sout, std::string name,
                 std::optional<MatrixType> const &param) {
  if (!param.has_value()) {
    return;
  }
  if (param->cols() == 1) {
    jsonParser json;
    sout << name << ": " << to_json_array(*param, json) << "\n";
  } else {
    jsonParser json;
    sout << name << ": " << to_json(*param, json) << "\n";
  }
}

}  // namespace conditions
}  // namespace Monte
}  // namespace CASM

#endif
