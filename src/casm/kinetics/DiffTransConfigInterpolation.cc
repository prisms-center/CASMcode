#include "casm/kinetics/DiffTransConfigInterpolation.hh"
#include "casm/kinetics/DiffTransConfiguration_impl.hh"

#include "casm/enumerator/Enumerator_impl.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/database/DiffTransConfigDatabase.hh"


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
      m_current(_diff_trans_config.sorted().from_config()) {
      auto diff_trans_config = _diff_trans_config.sorted();
      auto configs = get_relaxed_endpoints(diff_trans_config, calctype);
      Configuration from_config = configs.first;
      Configuration to_config = configs.second;
      DiffusionTransformation diff_trans  = diff_trans_config.diff_trans();
      //if(!from_config.has_displacement()) {
      //from_config.init_displacement();
      //}
      //if(!from_config.has_deformation()) {
      //from_config.init_deformation();
      //}
      //if(!to_config.has_displacement()) {
      //to_config.init_displacement();
      //}
      //if(!to_config.has_deformation()) {
      //to_config.init_deformation();
      //}
      Configuration to_config_mutated = prepare_to_config(to_config, diff_trans);
      m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(from_config,
                                                                            to_config_mutated,
                                                                            n_images + 2); // +2 for end states
      this->_initialize(&m_current);
      m_current.set_source(this->source(step()));
    }

    /// \brief Construct with a Supercell, using all permutations
    ///
    /// \param _diff_trans_config abcdcd
    /// \param n_images number of images excluding end points
    DiffTransConfigInterpolation::DiffTransConfigInterpolation(
      const DiffusionTransformation &diff_trans,
      const Configuration from_config,
      const Configuration to_config,
      const int n_images):
      RandomAccessEnumeratorBase<Configuration>(n_images + 2),
      m_current(make_attachable(diff_trans, from_config)) {
      Configuration to_config_mutated = prepare_to_config(to_config, diff_trans);
      m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(make_attachable(diff_trans, from_config),
                                                                            to_config_mutated,
                                                                            n_images + 2); // +2 for end states
      this->_initialize(&m_current);
      m_current.set_source(this->source(step()));
    }

    const std::string DiffTransConfigInterpolation::enumerator_name = "DiffTransConfigInterpolation";

    std::string DiffTransConfigInterpolation::interface_help() {
      return
        "DiffTransConfigInterpolation: \n\n"
        "  n_images: integer \n"
        "    The number of images to interpolate for each diff_trans_config\n\n"

        "  selection: string (optional, default="") \n"
        "    The names of a selection of diff_trans_configs to interpolate images \n"
        "    within.\n\n"
        "  names: JSON array of strings (optional, default=[]) \n"
        "    The names of a selection of diff_trans_configs to interpolate images \n"
        "    within.\n\n"

        "  calctype: string (optional, default=$current_calctype)\n"
        "    The name of the calctype to obtain the fixed lattice calculations for  \n"
        "    the endpoints of the diff_trans_configs. \n\n"

        "  Example:\n"
        "  {\n"
        "    \"n_images\": 4,\n"
        "    \"selection\": \"low_barrier_diff_trans\",\n"
        "    \"calctype\": \"fixed_lattice\""
        "  }\n\n";
    }

    int DiffTransConfigInterpolation::run(const PrimClex &primclex,
                                          const jsonParser &kwargs,
                                          const Completer::EnumOption &enum_optconst) {

      // get selection filename from json/enumoption // do json for now
      // selcection by names
      // Constrct a DB selection of DiffTransConfiguration from json and enumoption inputs
      DB::Selection<DiffTransConfiguration> dtc_sel = make_selection<DiffTransConfiguration>(
                                                        primclex, kwargs, "names", "selection", enumerator_name, OnError::THROW);
      int n_images = kwargs["n_images"].get<int>(); // set defaults with get_else
      std::string calctype = kwargs["calctype"].get<std::string>();
      Index i;
      for(const auto &config : dtc_sel.selected()) {
        // Create a interpolation object
        DiffTransConfigInterpolation enumerator(config, n_images, calctype);
        i = 0;
        for(const auto &img_config : enumerator) {
          // file_path = $project_dir/training_data/diff_trans/$diff_trans_name/$scelname/$configid/$image_number/POSCAR
          int n_images = kwargs[config.name()]["n_images"].get<int>(); // set defaults with get_else
          std::string calctype = kwargs[config.name()]["calctype"].get<std::string>();
          auto file_path = primclex.dir().configuration_calc_dir(config.name(), calctype);
          file_path += "N_images_" + std::to_string(n_images) + "/0" + std::to_string(i) + "/POSCAR";
          fs::ofstream file(file_path);
          img_config.write_pos(file);
          i++;
        }
      }
      // setup error methods
      return 0;
    }

    Configuration DiffTransConfigInterpolation::prepare_to_config(const Configuration &config,
                                                                  const DiffusionTransformation &diff_trans) {
      Configuration result = config;
      for(auto traj : diff_trans.species_traj()) {
        Index k = config.supercell().linear_index(traj.from.uccoord);
        Index l = config.supercell().linear_index(traj.to.uccoord);

        result.set_occ(k, traj.from.occ);

        //Eigen::Vector3d displacement = config.disp(l);

        const Eigen::Vector3d from_pos = traj.from.uccoord.coordinate().const_cart();
        const Eigen::Vector3d to_pos = traj.to.uccoord.coordinate().const_cart();
        Eigen::Vector3d ideal_pos_inc = to_pos - from_pos;
        //Eigen::Vector3d final_disp = displacement + ideal_pos_inc;

        //result.set_disp(k, final_disp);

      }
      return result;
    }

    const Configuration *DiffTransConfigInterpolation::at_step(step_type n) {
      Configuration result = (*m_config_enum_interpol)[n];
      //m_current.set_displacement(result.displacement());
      //m_current.set_deformation(result.deformation());
      return &m_current;
    }

    /// Returns copies of from and to config updated such they reflect the relaxed structures from the properties database
    std::pair<Configuration, Configuration> get_relaxed_endpoints(const DiffTransConfiguration &dfc, std::string calctype) {
      Configuration ret_frm = dfc.from_config();
      Configuration ret_to = dfc.to_config();
      if(dfc.from_configname().find("none") == std::string::npos) {
        Configuration rlx_frm = copy_apply_properties(make_configuration(dfc.primclex(), dfc.from_configname()), calctype);
        ret_frm = copy_apply(dfc.from_config_from_canonical(), rlx_frm);
      }
      if(dfc.to_configname().find("none") == std::string::npos) {
        Configuration rlx_to = copy_apply_properties(make_configuration(dfc.primclex(), dfc.to_configname()), calctype);
        ret_to = copy_apply(dfc.to_config_from_canonical(), rlx_to);
      }
      return std::make_pair(ret_frm, ret_to);
    }

    /// applies deformation from input config onto outpur_config and prints it to a file.
    void apply_deformation(const PrimClex primclex, std::string output_configname, std::string output_path,
                           std::string input_configname, std::string calctype) {
      Configuration input_config = copy_apply_properties(make_configuration(primclex, input_configname), calctype);
      //Eigen::Matrix3d deformation = input_config.deformation();
      Configuration final_config = make_configuration(primclex, output_configname);
      //final_config.set_deformation(deformation);
      fs::ofstream file(output_path);
      final_config.write_pos(file);
      return;
    }
  }
}
