#ifndef CASM_ClusterSpecsParser_impl
#define CASM_ClusterSpecsParser_impl

#include "casm/clusterography/ClusterSpecsParser.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/CASM_global_definitions.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  template<typename RequiredType>
  bool OrbitBranchSpecsParser::require_previous(
    jsonParser::const_iterator it,
    std::string option) {

    auto prev_it = previous(it);
    if(prev_it == self_it->end()) {
      std::stringstream msg;
      msg << "Error: "
          << "branch '" << branch_to_string(branch_to_int(it) - 1) << "' is required because branch '"
          << it.name() << "' is included.";
      error.insert(msg.str());
      return false;
    }
    return bool(require<RequiredType>(it, option));
  }

  template<typename RequiredType>
  bool OrbitBranchSpecsParser::require_nonincreasing(
    jsonParser::const_iterator it,
    std::string option) {

    if(!require<RequiredType>(it, option) || !require_previous<RequiredType>(it, option)) {
      return false;
    }

    if(it->find(option)->get<double>() > previous(it)->find(option)->get<double>()) {
      std::stringstream msg;
      msg << "Error: "
          << "'" << option << "' increases from branch '" << previous(it).name()
          << "' to branch '" << it.name() << "'";
      error.insert(msg.str());
      return false;
    }
    return true;
  }

}

#endif
