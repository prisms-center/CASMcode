#include "casm/clex/Reference.hh"

#include "casm/clex/Supercell.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  // --- ConstantReference implementations -----------

  const std::string ConstantReference::Name = "ConstantReference";

  const std::string ConstantReference::Desc = "Returns a constant scalar.";


  // --- HyperPlaneReferenceBase implementations -----------

  /// \brief Return the reference hyperplane used for a particular configuration
  ///
  /// Returns the Configuration specific hyperplane if it exists, else the
  /// Supercell specific hyperplane if it exists, else the global hyperplane.
  Eigen::VectorXd HyperPlaneReferenceBase::hyperplane(const Configuration &config) const {

    auto c_it = m_config_ref.find(config.name());
    if(c_it != m_config_ref.end()) {
      return c_it->second;
    }

    auto s_it = m_supercell_ref.find(config.get_supercell().get_name());
    if(s_it != m_supercell_ref.end()) {
      return s_it->second;
    }

    return m_global_ref;
  }


  // --- HyperPlaneReference implementations -----------

  const std::string HyperPlaneReference::Name = "HyperPlaneReference";

  const std::string HyperPlaneReference::Desc = "Returns a reference value based on the value of a hyperplane.";




}
