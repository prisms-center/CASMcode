#ifndef CASM_ConfigEnumRandomOccupations
#define CASM_ConfigEnumRandomOccupations

#include <functional>

#include "casm/app/enum/EnumInterface.hh"
#include "casm/clex/Configuration.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/misc/cloneable_ptr.hh"

extern "C" {
CASM::EnumInterfaceBase *make_ConfigEnumRandomOccupations_interface();
}

class MTRand;

namespace CASM {

class ConfigEnumInput;

/** \defgroup ConfigEnumGroup Configuration Enumerators
 *  \ingroup Configuration
 *  \ingroup Enumerator
 *  \brief Enumerates Configuration
 *  @{
 */

/// Parameters controlling ConfigEnumRandomOccupations
struct ConfigEnumRandomOccupationsParams {
  ConfigEnumRandomOccupationsParams(MTRand &_mtrand, Index _n_config);

  /// Random number generator
  MTRand &mtrand;

  /// Number of random configurations to generate
  Index n_config;
};

/// \brief Enumerate n random occupations in a particular Supercell
///
class ConfigEnumRandomOccupations : public InputEnumeratorBase<Configuration> {
  // -- Required members -------------------

 public:
  ConfigEnumRandomOccupations(ConfigEnumInput const &_in_config,
                              ConfigEnumRandomOccupationsParams const &params);

  ConfigEnumRandomOccupations(ConfigEnumInput const &_in_config,
                              Index _n_config, MTRand &_mtrand);

  std::string name() const override;

  static const std::string enumerator_name;

 private:
  /// Implements increment
  void increment() override;

  // -- Unique -------------------

  void randomize();

  Index m_n_config;
  MTRand &m_mtrand;
  std::vector<int> m_max_allowed;
  std::vector<Index> m_site_selection;
  notstd::cloneable_ptr<Configuration> m_current;
};

/** @}*/
}  // namespace CASM

#endif
