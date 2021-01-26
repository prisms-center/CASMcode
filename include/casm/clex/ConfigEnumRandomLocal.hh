#ifndef CASM_ConfigEnumRandomLocal
#define CASM_ConfigEnumRandomLocal

#include <functional>

#include "casm/clex/Configuration.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/misc/cloneable_ptr.hh"

class MTRand;

namespace CASM {

class ConfigEnumInput;

/** \defgroup ConfigEnumGroup Configuration Enumerators
 *  \ingroup Configuration
 *  \ingroup Enumerator
 *  \brief Enumerates Configuration
 *  @{
 */

/// Parameters controlling ConfigEnumRandomLocal
struct ConfigEnumRandomLocalParams {
  ConfigEnumRandomLocalParams(MTRand &_mtrand, DoFKey _dof_key, Index _n_config,
                              double _mag, bool _normal_distribution);

  /// Random number generator
  MTRand &mtrand;

  /// Name of site degree of freedom for which normal coordinates are to be
  /// generated.
  ///
  /// DoFKey is a typedef for std::string
  DoFKey dof_key;

  /// Number of random configurations to generate
  Index n_config;

  /// Magnitude used to scale random vector at each site
  double mag;

  /// True if using "normal" distribution, else using "uniform" distribution
  bool normal_distribution;
};

/// Enumerate random values for continuous degrees of freedom
class ConfigEnumRandomLocal : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  ConfigEnumRandomLocal(ConfigEnumInput const &_in_config,
                        ConfigEnumRandomLocalParams const &params);

  ConfigEnumRandomLocal(ConfigEnumInput const &_in_config,
                        DoFKey const &_dof_key, Index _n_config, double _mag,
                        bool _normal, MTRand &_mtrand);

  std::string name() const override;

  static const std::string enumerator_name;

 private:
  /// Implements increment
  void increment() override;

  // -- Unique -------------------

  void randomize();

  // Key of DoF being perturbed
  DoFKey m_dof_key;

  // Number of configurations to be enumerated
  Index m_n_config;

  LocalContinuousConfigDoFValues *m_dof_vals;

  // Pseudo-random number generator
  MTRand &m_mtrand;

  // std deviation of normal distribution
  // max magnitude if uniform distribution
  double m_mag;

  // true if normally distributed
  // false if uniformly distributed
  bool m_normal;

  // true if dof is unit vector
  // false otherwise
  bool m_unit_length;

  // Sites on which enumeration is being performed
  std::vector<Index> m_site_selection;

  // Dimension of DoF at each selected site
  std::vector<Index> m_dof_dims;

  notstd::cloneable_ptr<Configuration> m_current;
};

/** @}*/
}  // namespace CASM

#endif
