#include "TestCompN.hh"
#include "casm/clex/Configuration.hh"

extern "C" {
  CASM::BaseDatumFormatter<CASM::Configuration> *make_TestCompN_formatter() {
    return new CASM::ConfigIO::TestCompN();
  }
}

namespace CASM {

  namespace ConfigIO {
    // --- CompN implementations -----------

    const std::string TestCompN::Name = "test_comp_n";

    const std::string TestCompN::Desc =
      "Number of each species per unit cell, including vacancies. "
      "No argument prints all available values. Ex: test_comp_n, test_comp_n(Au), test_comp_n(Pt), etc.";

    /// \brief Returns the number of each species per unit cell
    Eigen::VectorXd TestCompN::evaluate(const Configuration &config) const {
      return comp_n(config);
    }
  }
}
