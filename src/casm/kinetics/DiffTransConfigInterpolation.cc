#include "casm/kinetics/DiffTransConfigInterpolation.hh"
#include "casm/kinetics/DiffTransConfiguration_impl.hh"

#include "casm/container/Enumerator_impl.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Selection_impl.hh"


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
      const int n_images, std::string calctype):
      RandomAccessEnumeratorBase<Configuration>(n_images + 2),
      m_current(_diff_trans_config.sorted().from_config()),
      m_diff_trans_config(_diff_trans_config.sorted()) {
      auto configs = get_relaxed_endpoints(m_diff_trans_config, calctype);
      Configuration from_config = configs.first;
      Configuration to_config = configs.second;
      DiffusionTransformation diff_trans  = m_diff_trans_config.diff_trans();
      Configuration to_config_mutated = prepare_to_config(to_config, diff_trans);
      m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(from_config,
                                                                            to_config_mutated,
                                                                            n_images + 2); // +2 for end states
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
                                                                  const DiffusionTransformation &diff_trans) {
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
        result.set_disp(k, final_disp);
      }
      return result;
    }

    const Configuration *DiffTransConfigInterpolation::at_step(step_type n) {
      Configuration result = (*m_config_enum_interpol)[n];
      m_current.set_displacement(result.displacement());
      m_current.set_deformation(result.deformation());
      return &m_current;
    }

    /// Returns copies of from and to config updated such they reflect the relaxed structures from the properties database
    std::pair<Configuration, Configuration> get_relaxed_endpoints(const DiffTransConfiguration &dfc, std::string calctype) {
      Configuration rlx_frm = copy_apply_properties(make_configuration(dfc.primclex(), dfc.from_configname()), calctype);

      Configuration rlx_to = copy_apply_properties(make_configuration(dfc.primclex(), dfc.to_configname()), calctype);

      return std::make_pair(copy_apply(dfc.from_config_from_canonical(), rlx_frm),
                            copy_apply(dfc.to_config_from_canonical(), rlx_to));

    }
  }
}
