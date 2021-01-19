#include "casm/clex/ConfigurationTraits.hh"

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <vector>

#include "casm/casm_io/Log.hh"
#include "casm/global/definitions.hh"

namespace CASM {

const std::string traits<Configuration>::name = "Configuration";

const std::string traits<Configuration>::short_name = "config";

/// Tokenizes 'SCELV_A_B_C_D_E_F/I' to integral values {V, A, B, C, D, E, F, I}
/// and does lexicographical comparison
bool traits<Configuration>::name_compare(std::string A, std::string B) {
  try {
    std::vector<std::string> splt_vec_A;
    boost::split(splt_vec_A, A, boost::is_any_of("L_/"),
                 boost::token_compress_on);
    std::vector<std::string> splt_vec_B;
    boost::split(splt_vec_B, B, boost::is_any_of("L_/"),
                 boost::token_compress_on);
    for (int i = 1; i < splt_vec_A.size(); ++i) {
      Index i_A = boost::lexical_cast<Index>(splt_vec_A[i]);
      Index i_B = boost::lexical_cast<Index>(splt_vec_B[i]);
      if (i_A != i_B) {
        return i_A < i_B;
      }
    }
    return false;
  } catch (std::exception &e) {
    err_log().error("In traits<Configuration>::name_compare");
    err_log() << "A: " << A << std::endl;
    err_log() << "B: " << B << std::endl;
    err_log() << e.what() << std::endl;
    throw e;
  }
};
}  // namespace CASM
