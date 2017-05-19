#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/container/Enumerator_impl.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransConfigInterpolation.hh"
// extern "C" {
//   CASM::EnumInterfaceBase *make_DiffTransConfigInterpolation_interface() {
//     return new CASM::EnumInterface<CASM::DiffTransConfigInterpolation>();
//   }
// }

namespace CASM {

  namespace Kinetics {

    /// \brief Construct with a Supercell, using all permutations
    DiffTransConfigInterpolation::DiffTransConfigInterpolation(const DiffTransConfiguration &_diff_trans_config):m_diff_trans_config(_diff_trans_config) {}

    const std::string DiffTransConfigInterpolation::enumerator_name = "DiffTransConfigInterpolation";

    const std::string DiffTransConfigInterpolation::interface_help =
      "DiffTransConfigInterpolation: \n\n"
      
      "      }' \n\n";

    int DiffTransConfigInterpolation::run(const std::string config_path,
                                          const int n_images,
                                          const DiffTransConfiguration diff_trans_config) {
      Configuration from_config = diff_trans_config.from_config();
      Configuration to_config = diff_trans_config.to_config();
      DiffusionTransformation diff_trans  = diff_trans_config.diff_trans();
      Configuration from_config_edited = make_config_diff_trans_free(from_config, diff_trans);
      Configuration to_config_edited = make_config_diff_trans_free(to_config, diff_trans);
      ConfigEnumInterpolation e(from_config_edited, to_config_edited, n_images);

      Index i = 0;
      for(const auto &tconfig : e) {
        Configuration con = tconfig;
        ++i;
      }
      // auto lambda = [&](const DiffTransConfiguration & diff_trans_config) {
      //   return notstd::make_unique<DiffTransConfigInterpolation>(diff_trans_config);
      // };
      
      // int returncode = make_config_interpolation(
      //                    config_path,
      //                    n_images);
      return 0;
    }
    Configuration make_config_diff_trans_free(const Configuration &config,
                                              const DiffusionTransformation &diff_trans) {
      Configuration result = config;
      for(auto traj : diff_trans.specie_traj()) {
        Index l = config.supercell().linear_index(traj.from.uccoord);
        result.set_occ(l, traj.from.occ);
      }
      return result;
    }
  }
}
