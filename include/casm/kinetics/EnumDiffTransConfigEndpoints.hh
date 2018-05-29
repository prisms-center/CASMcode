#ifndef CASM_EnumDiffTransConfigEndpoints
#define CASM_EnumDiffTransConfigEndpoints

#include "casm/container/Counter.hh"
#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
// extern "C" {
//   CASM::EnumInterfaceBase *make_EnumDiffTransConfigEndpoints_interface();
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
    class EnumDiffTransConfigEndpoints : public InputEnumeratorBase<Configuration> {

      // -- Required members -------------------

    public:

      /// \brief Construct with a DiffTransConfiguration
      /// Adds the end points of an in project DiffTransConfiguration to the CASM project
      EnumDiffTransConfigEndpoints(const DiffTransConfiguration &_diff_trans_config);

      std::string name() const override {
        return enumerator_name;
      }

      static const std::string enumerator_name;
      static const std::string interface_help;

      /// Command line interface parsing to give instructions on how to enumerate the endpoints
      static int run(const PrimClex &primclex, const jsonParser &kwargs, const Completer::EnumOption &enum_opt);

    private:
      Configuration m_current;
      DiffTransConfiguration m_source;
    };
  }
}
#endif
