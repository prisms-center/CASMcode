#ifndef CASM_monte2_Definitions
#define CASM_monte2_Definitions

#include <map>
#include <string>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace Monte2 {

typedef long CountType;
typedef double TimeType;

/// Map of value name to vector value
typedef std::map<std::string, Eigen::VectorXd> VectorValueMap;

}  // namespace Monte2
}  // namespace CASM

#endif
