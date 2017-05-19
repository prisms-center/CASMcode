#ifndef CASM_DiffTransConfigInterpolation
#define CASM_DiffTransConfigInterpolation

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/clex/ConfigEnumInterpolation.hh"
#include "casm/kinetics/DiffTransconfiguration.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
// extern "C" {
//   CASM::EnumInterfaceBase *make_DiffTransConfigInterpolation_interface();
// }

namespace CASM {

  namespace Kinetics {
    /** \defgroup ConfigEnumGroup Configuration Enumerators
     *  \ingroup Configuration
     *  \ingroup Enumerator
     *  \brief Enumerates Configuration
     *  @{
     */

    /// \brief Enumerate over all possible occupations in a particular Supercell
    ///
    class DiffTransConfigInterpolation : public InputEnumeratorBase<Configuration> {

      // -- Required members -------------------

    public:

      /// \brief Construct with a Supercell, using all permutations
      DiffTransConfigInterpolation(const DiffTransConfiguration &_diff_trans_config);

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;
      static const std::string interface_help;
      
      static int run(const std::string config_path, const int n_images, DiffTransConfiguration diff_trans_config);
      
    private:
      // Returns to and from configs that are edited to have same occupants on the diff_trans_orbit
      // output used to intepolate rest of the configuration other that diff trans orbit 

      DiffTransConfiguration m_diff_trans_config;

      // /// Implements increment
      // void increment() override;


      // // -- Unique -------------------
      
      // /// Returns true if current() is primitive and canonical
      // bool _check_current() const;
      
      // Counter<std::vector<int> > m_counter;
      // notstd::cloneable_ptr<Configuration> m_current;
    };
    Configuration make_config_diff_trans_free(const Configuration &config, const DiffusionTransformation &diff_trans);
  }
}

#endif
