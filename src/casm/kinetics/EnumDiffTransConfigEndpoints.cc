#include "casm/kinetics/EnumDiffTransConfigEndpoints.hh"

#include "casm/container/Enumerator_impl.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/database/DiffTransConfigDatabase.hh"


// extern "C" {
//   CASM::EnumInterfaceBase *make_EnumDiffTransConfigEndpoints_interface() {
//     return new CASM::EnumInterface<CASM::EnumDiffTransConfigEndpoints>();
//   }
// }

namespace CASM {

  namespace Kinetics {

    /// \brief Construct with a Supercell, using all permutations
    ///
    /// \param _diff_trans_config abcdcd
    /// \param n_images number of images excluding end points
    EnumDiffTransConfigEndpoints::EnumDiffTransConfigEndpoints(
      const DiffTransConfiguration &_diff_trans_config):
      InputEnumeratorBase<Configuration>(),
      m_current(_diff_trans_config.sorted().from_config()),
      m_source(_diff_trans_config) {
      this->_initialize(&m_current);
      m_current.set_source(this->source(step()));
    }


    const std::string EnumDiffTransConfigEndpoints::enumerator_name = "EnumDiffTransConfigEndpoints";

    const std::string EnumDiffTransConfigEndpoints::interface_help =
      "EnumDiffTransConfigEndpoints: \n\n"
      "  selection: string (optional, default=MASTER) \n"
      "    The name of a selection of diff_trans_configs for which to enumerate endpoints \n"
      "  names: JSON array of strings (optional, default=[]) \n"
      "    The names of diff_trans_configs for which to enumerate endpoints \n"
      "  NOTE: although both of these are considered optional one of them must be specified.\n\n"
      "  Example:\n"
      "  {\n"
      "    \"names\": [\"diff_trans/0/SCEL8_2_2_2_0_0_0/2\",\"diff_trans/0/SCEL8_2_2_2_0_0_0/4\"],\n"
      "    \"selection\": \"dtc_to_calculate\"\n"
      "  }\n\n";

    int EnumDiffTransConfigEndpoints::run(const PrimClex &primclex,
                                          const jsonParser &kwargs,
                                          const Completer::EnumOption &enum_optconst) {

      // get selection filename from json/enumoption // do json for now
      // Constrct a DB selection of DiffTransConfiguration from json and enumoption inputs
      DB::Selection<DiffTransConfiguration> sel = make_selection<DiffTransConfiguration>(primclex, kwargs, "names", "selection", enumerator_name, OnError::THROW);

      Log &log = primclex.log();
      auto &db_config = primclex.db<Configuration>();

      Index Ninit = db_config.size();
      log << "# configurations in this project: " << Ninit << "\n" << std::endl;
      for(const auto &config : sel.selected()) {

        db_config.insert(config.from_config().canonical_form());
        db_config.insert(config.to_config().canonical_form());
        config.cache_clear();

      }

      log << "  DONE." << std::endl << std::endl;

      Index Nfinal = db_config.size();

      log << "# new configurations: " << Nfinal - Ninit << "\n";
      log << "# configurations in this project: " << Nfinal << "\n" << std::endl;

      log << "Writing config database..." << std::endl;
      db_config.commit();
      log << "  DONE" << std::endl;
      // setup error methods
      return 0;
    }

  }
}
