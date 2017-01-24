#ifndef UNIT_COMMON_HH
#define UNIT_COMMON_HH

#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"

namespace test {
  /// \brief Check expected JSON vs calculated JSON using BOOST_CHECK_EQUAL
  bool check(std::string test,
             const jsonParser &expected,
             const jsonParser &calculated,
             fs::path test_cases_path,
             bool quiet,
             double tol = 0.0);

}

#endif
