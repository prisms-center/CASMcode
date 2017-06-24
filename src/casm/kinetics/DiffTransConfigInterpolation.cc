#include "casm/casm_io/Log.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
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
    ///
    /// \param _diff_trans_config abcdcd
    /// \param n_images number of images excluding end points
    DiffTransConfigInterpolation::DiffTransConfigInterpolation(
      const DiffTransConfiguration &_diff_trans_config,
      const int n_images):
      RandomAccessEnumeratorBase<Configuration>(n_images+2),
      m_current(_diff_trans_config.from_config()),
      m_diff_trans_config(_diff_trans_config) {
      
      Configuration from_config = m_diff_trans_config.from_config();
      // set displacements from the properties.
      jsonParser json = from_config.calc_properties();
      // std::vector<std::vector<double>> deform = json["relaxed_deformation"];
      // Eigen::Matrix3d deform = json["relaxed_deformation"];
      // from_config.set_deformation(deform);
      // from_config.set_displacement(json["relaxed_displacement"]);

      Configuration to_config = get_relaxed_to_config(m_diff_trans_config);
      DiffusionTransformation diff_trans  = m_diff_trans_config.diff_trans();
      Configuration to_config_mutated = prepare_to_config(to_config, diff_trans);
      m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(from_config,
                                                                            to_config_mutated, n_images+2); // +2 for end states
      
      // Index i = 0;
      // for(const auto &tconfig : *m_config_enum_interpol) {
      //   Configuration con = tconfig;
      //   // con.write_pos(std::cout);
      //   Configuration::displacement_matrix_t config_disp = con.displacement();
      //   std::cout << "Image #:" << i << "\n";
      //   std::cout << config_disp << "\n";
      //   ++i;
      // }

      this->_initialize(&m_current);
      m_current.set_source(this->source(step()));
    }

    const std::string DiffTransConfigInterpolation::enumerator_name = "DiffTransConfigInterpolation";

    const std::string DiffTransConfigInterpolation::interface_help =
      "DiffTransConfigInterpolation: \n\n"
      
      "      }' \n\n";

    int DiffTransConfigInterpolation::run(const PrimClex &primclex,
                                          const jsonParser &kwargs,
                                          const Completer::EnumOption &enum_optconst) {

      // Constrct a DB selection of DiffTransConfiguration from json and enumoption inputs
      // DB::Selection<DiffTransConfiguration> sel(primclex);
      // int n_images = 4; // For Testing purposes
      // for (const auto &config : sel.selected()){
      //   // Create a interpolation object  
      //   DiffTransConfigInterpolation enumerator(config, n_images);
      //   for (const auto &thing: enumerator){
      //     //do a print of things
      //   }
      // }
      // setup error methods
      return 0;
    }

    Configuration DiffTransConfigInterpolation::prepare_to_config(const Configuration &config,
                                                                  const DiffusionTransformation &diff_trans){
      Configuration result = config;
      for(auto traj : diff_trans.specie_traj()) {
        Index k = config.supercell().linear_index(traj.from.uccoord);
        Index l = config.supercell().linear_index(traj.to.uccoord);
        result.set_occ(k, traj.from.occ);
        Eigen::Vector3d displacement = config.disp(l);
        const Eigen::Vector3d from_pos = config.supercell().coord(k).frac();
        const Eigen::Vector3d to_pos = config.supercell().coord(l).frac();
        Eigen::Vector3d ideal_pos_inc = to_pos - from_pos;
        Eigen::Vector3d final_disp = displacement + ideal_pos_inc;
        result.set_disp(k,final_disp);
      }
      
      return result;
    }

    const Configuration *DiffTransConfigInterpolation::at_step(step_type n) {
      Configuration result = (*m_config_enum_interpol)[n];
      m_current.set_displacement(result.displacement());
      return &m_current;
    }

    // In namespace Kinetics 
    Configuration get_relaxed_to_config(const DiffTransConfiguration &dfc) {
      Configuration to_config = dfc.to_config();
      
      // find the configname
      // extract the deformations and displacemts of the cannonical config
      // find transformation operation from cannonical to to_config
      // apply the transformation to deformations and displacements
      // apply the transformed deformations and displacements to the to_config

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
      
      return to_config;
    }
  }
}
