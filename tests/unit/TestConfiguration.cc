#include "TestConfiguration.hh"
#include "casm/casm_io/jsonParser.hh"

namespace test {

  TestConfiguration::TestConfiguration(
    const PrimClex &primclex,
    const Eigen::Matrix3i &T,
    const std::vector<int> &_occupation) :
    TestConfiguration(primclex, make_supercell(primclex.prim().lattice(), T), _occupation) {}

  TestConfiguration::TestConfiguration(
    const PrimClex &primclex,
    const Lattice &lat,
    const std::vector<int> &_occupation) :
    TestSupercell(primclex, lat),
    config(this->scel, jsonParser(), ConfigDoF(_occupation, primclex.crystallography_tol())),
    config_permute_fg(config.factor_group()),
    config_sym_fg(make_sym_group(config_permute_fg.begin(), config_permute_fg.end())) {}


}

