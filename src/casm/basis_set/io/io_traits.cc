#include "casm/basis_set/io/io_traits.hh"
#include "casm/basis_set/BasisFunctionSpecs.hh"

namespace CASM {

  const std::string traits<PARAM_PACK_TYPE>::name = "param_pack_type";

  const std::multimap<PARAM_PACK_TYPE, std::vector<std::string> > traits<PARAM_PACK_TYPE>::strval = {
    {PARAM_PACK_TYPE::DEFAULT, {"DEFAULT", "Default", "default"} },
    {PARAM_PACK_TYPE::DIFF, {"DIFF", "Diff", "diff"} }
  };

}
