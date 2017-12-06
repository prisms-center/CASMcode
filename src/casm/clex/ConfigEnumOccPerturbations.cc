/*
#include "casm/clex/ConfigEnumOccPerturbations_impl.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/container/Enumerator_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_ConfigEnumOccPerturbations_interface() {
    return new CASM::EnumInterface<CASM::ConfigEnumOccPerturbations>();
  }
}

namespace CASM {

  const std::string ConfigEnumOccPerturbations::enumerator_name = "ConfigEnumOccPerturbations";

  const std::string ConfigEnumOccPerturbations::interface_help =
    "ConfigEnumOccPerturbations: \n\n"

    "  background_configs: JSON array of strings \n "
    "    Indicate which configurations will be the background structures to be"
    "    perturbed. The JSON array of strings \"background_configs\" should be"
    "    names of Configurations that exist in your CASM project.\n\n"

    "  local_cspecs: JSON object (optional,default= {}) \n "
    "    Specify the clusters that should be perturbed. The string \"local_cspecs\"
    "    should be a local cspecs style initialization used in casm bset "
    "    enumeration. This option takes precedence over the following option.\n\n"

    "  local_cspecs_filepath: string (optional,default=\"\") \n "
    "    Indicate the local cspecs file that specifies the clusters that should "
    "    be perturbed. The string \"local_cspecs_filepath\" should be the file "
    "    path to a JSON file containing the local cspecs.\n\n"

    "  dry_run: bool (optional, default=false)\n"
    "    Perform dry run.\n\n"

    "  Example:\n"
    "  {\n"
    "   \"background_configs\":[\"SCEL8_2_2_2_0_0_0/2\"],\n"
    "    \"local_cspecs\":{\n"
    "       \"basis_functions\" : {\n"
    "        \"site_basis_functions\" : \"occupation\"\n"
    "      },\n"
    "      \"orbit_branch_specs\" : { \n"
    "       \"1\" : {\"cutoff_radius\" : 6.0},\n"
    "       \"2\" : {\"max_length\" : 6.01,\"cutoff_radius\" : 6.0},\n"
    "       \"3\" : {\"max_length\" : 4.01,\"cutoff_radius\" : 5.0}\n"
    "      }\n"
    "    }\n"
    "  }\n\n";

  int ConfigEnumOccPerturbations::run(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {}

  /// \brief Construct with a Supercell, using all permutations
  ConfigEnumOccPerturbations::ConfigEnumOccPerturbations(const Configuration &_background_config){}

  /// Implements _increment over all occupations
  void ConfigEnumOccPerturbations::increment() {}

}
*/
