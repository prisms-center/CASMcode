#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/container/Enumerator_impl.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransconfiguration.hh"
#include "casm/kinetics/DiffTransConfigInterpolation.hh"
#include "casm/database/Selection.hh"
// extern "C" {
//   CASM::EnumInterfaceBase *make_DiffTransConfigInterpolation_interface() {
//     return new CASM::EnumInterface<CASM::DiffTransConfigInterpolation>();
//   }
// }

namespace CASM {

  namespace Kinetics {

    /// \brief Construct with a Supercell, using all permutations
    DiffTransConfigInterpolation::DiffTransConfigInterpolation(const DiffTransConfiguration &_diff_trans_config,
                                                               int n_images):m_diff_trans_config(_diff_trans_config) {
      Configuration from_config = m_diff_trans_config.from_config();
      Configuration to_config = m_diff_trans_config.to_config();

      // for testing purposes
      // include displacements & deformation
      Eigen::Vector3d zero(0., 0., 0.);
      Eigen::Vector3d dx(0.001, 0., 0.);
      Eigen::Vector3d dy(0., 0.001, 0.);
      Eigen::Vector3d dz(0., 0., 0.001);
      to_config.set_disp(5, dx*4);
      to_config.set_disp(6, dx*3);
      to_config.set_disp(26, dy*3);
      to_config.set_disp(30, dx*1 + dy*4);
      to_config.set_disp(31, dz*4 + dy*2);


      DiffusionTransformation diff_trans  = m_diff_trans_config.diff_trans();
      Configuration from_config_edited = make_config_diff_trans_free(from_config, diff_trans);
      Configuration to_config_edited = make_config_diff_trans_free(to_config, diff_trans);
      ConfigEnumInterpolation e(from_config_edited, to_config_edited, n_images);
      // from_config.write_pos(std::cout);
      // to_config.write_pos(std::cout);
      // from_config_edited.write_pos(std::cout);
      // to_config_edited.write_pos(std::cout);
      
      Index i = 0;
      for(const auto &tconfig : e) {
        Configuration con = tconfig;
        // con.write_pos(std::cout);
        Configuration::displacement_matrix_t config_disp = con.displacement();
        std::cout << i << "\n";
        std::cout << config_disp << "\n";
        ++i;
      }
    }

    const std::string DiffTransConfigInterpolation::enumerator_name = "DiffTransConfigInterpolation";

    const std::string DiffTransConfigInterpolation::interface_help =
      "DiffTransConfigInterpolation: \n\n"
      
      "      }' \n\n";

    int DiffTransConfigInterpolation::run(const PrimClex &primclex,
                                          const jsonParser &kwargs,
                                          const Completer::EnumOption &enum_optconst) {

      // Constrct a DB selection of DiffTransConfiguration from json and enumoption inputs
      DB::Selection<DiffTransConfiguration> sel(primclex);
      int n_images = 4; // For Testing purposes
      for (const auto &config : sel.selected()){
        // Create a interpolation object  
        DiffTransConfigInterpolation enumerator(config, n_images);
        for (const auto &thing: enumerator){
          //do a print of things
        }
      }
      // setup error methods
      return 0;
    }

    // In namespace Kinetics 
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
