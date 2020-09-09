#ifndef CASM_ConfigEnumRandomLocal
#define CASM_ConfigEnumRandomLocal

#include "casm/app/enum/EnumInterface.hh"
#include "casm/clex/Configuration.hh"
#include "casm/container/Counter.hh"
#include "casm/enumerator/InputEnumerator.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/misc/cloneable_ptr.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumRandomLocal_interface();
}

class MTRand;

namespace CASM {

  /** \defgroup ConfigEnumGroup Configuration Enumerators
   *  \ingroup Configuration
   *  \ingroup Enumerator
   *  \brief Enumerates Configuration
   *  @{
  */

  /// \brief Enumerate random values for continuous degrees of freedom
  ///
  class ConfigEnumRandomLocal : public InputEnumeratorBase<Configuration> {

    // -- Required members -------------------

  public:

    /// \brief Construct with a Supercell, using all permutations
    ConfigEnumRandomLocal(
      ConfigEnumInput const &_in_config,
      DoFKey const &_dof_key,
      Index _n_config,
      double _mag,
      bool _normal,
      MTRand &_mtrand);

    std::string name() const override {
      return enumerator_name;
    }

    static const std::string enumerator_name;
    static std::string interface_help();
    static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt, EnumeratorMap const *interface_map = nullptr);

  private:
    /// Implements increment
    void increment() override;

    // -- Unique -------------------

    void randomize();

    // Key of DoF being perturbed
    DoFKey m_dof_key;

    // Number of configurations to be enumerated
    Index m_n_config;

    ConfigDoF::LocalDoFContainerType *m_dof_vals;

    // Pointer to pseudo-random number generator
    MTRand *m_mtrand;

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
}

#endif
