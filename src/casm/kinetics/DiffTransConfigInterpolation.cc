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
      
      Configuration from_config = get_relaxed_config(m_diff_trans_config, m_diff_trans_config.from_config());
      Configuration to_config = get_relaxed_config(m_diff_trans_config, m_diff_trans_config.to_config());
      DiffusionTransformation diff_trans  = m_diff_trans_config.diff_trans();
      Configuration to_config_mutated = prepare_to_config(to_config, diff_trans);
      m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(from_config,
                                                                            to_config_mutated,
                                                                            n_images+2); // +2 for end states
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
        const Eigen::Vector3d from_pos = config.supercell().coord(k).const_cart();
        const Eigen::Vector3d to_pos = config.supercell().coord(l).const_cart();
        Eigen::Vector3d ideal_pos_inc = to_pos - from_pos;
        Eigen::Vector3d final_disp = displacement + ideal_pos_inc;
        result.set_disp(k,final_disp);
      }
      return result;
    }

    const Configuration *DiffTransConfigInterpolation::at_step(step_type n) {
      Configuration result = (*m_config_enum_interpol)[n];
      m_current.set_displacement(result.displacement());
      m_current.set_deformation(result.deformation());
      return &m_current;
    }

    // In namespace Kinetics
    // move it to DiffTransConfig
    Configuration get_relaxed_config(const DiffTransConfiguration &dfc, Configuration config) {
      config.init_deformation();
      config.init_displacement();
      return config;
      
      // auto it = to_config.from_canonical();
      // Configuration to_config_canonical = to_config.canonical_form();
      // to_config_canonical.init_deformation();
      // to_config_canonical.init_displacement();
      // bool is_data, is_data_complete;
      // jsonParser json;
      // std::tie(json, is_data, is_data_complete) = to_config_canonical.read_calc_properties();
      // auto it_1 = json.find("relaxation_deformation");
      // if(it_1 != json.end()) {
      //   to_config_canonical.set_deformation(it_1->get<Eigen::Matrix3d>());
      // }
      // auto it_2 = json.find("relaxation_displacement");
      // if(it_2 != json.end()) {
      //   to_config_canonical.set_displacement(it_2->get<Eigen::MatrixXd>());
      // }
      // Configuration to_config_final = copy_apply(it, to_config_canonical);
      // return to_config_final;
    }
  }
}
