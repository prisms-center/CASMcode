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
#include "casm/casm_io/jsonParser.hh"


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
      const int n_images, std::string calctype,
      bool override_mirrors):
      RandomAccessEnumeratorBase<Configuration>(n_images + 2),
      m_current(_diff_trans_config.sorted().from_config()) {
      auto diff_trans_config = _diff_trans_config.sorted();
      if(diff_trans_config.from_configname() == diff_trans_config.to_configname() && n_images > 1) {
        std::cout << diff_trans_config.name() << " has symmetrically equivalent endpoints. \n"
                  << " A full NEB calculation is not required to determine barrier height." << std::endl
                  << " We suggest interpolating a single image and calculating each of the two structures individually." << std::endl;

      }
      auto configs = get_relaxed_endpoints(diff_trans_config, calctype);
      Configuration from_config = configs.first;

      Configuration to_config = configs.second;
      if(!override_mirrors && diff_trans_config.from_configname() == diff_trans_config.to_configname()) {
        std::cout << "I am going to interpolate 1 image with the same deformed \n" <<
                  "lattice as endpoints but ideally interpolated coordinates" << std::endl;
        from_config.clear_displacement();
        to_config.clear_displacement();

      }
      m_current = from_config;
      if(from_config.supercell().lattice().lat_column_mat() != to_config.supercell().lattice().lat_column_mat()) {
        throw std::runtime_error("Attempting to interpolate between configs with different lattices.\n You will have a bad time.");
      }
      DiffusionTransformation diff_trans  = diff_trans_config.diff_trans();
      if(!from_config.has_displacement()) {
        from_config.init_displacement();
      }
      if(!from_config.has_deformation()) {
        from_config.init_deformation();
      }
      if(!to_config.has_displacement()) {
        to_config.init_displacement();
      }
      if(!to_config.has_deformation()) {
        to_config.init_deformation();
      }
      Configuration to_config_mutated = prepare_to_config(to_config, diff_trans);
      if(!override_mirrors && diff_trans_config.from_configname() == diff_trans_config.to_configname()) {
        m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(from_config,
                                                                              to_config_mutated,
                                                                              3); // +2 for end states
      }
      else {
        m_config_enum_interpol = notstd::make_unique<ConfigEnumInterpolation>(from_config,
                                                                              to_config_mutated,
                                                                              n_images + 2); // +2 for end states
      }
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
      std::string endpt_calctype = kwargs["endstate_calctype"].get<std::string>();
      std::string calctype = kwargs["calctype"].get<std::string>();
      bool override_mirrors;
      kwargs.get_else(override_mirrors, "override_mirrors", false);
      Index i;
      for(const auto &config : dtc_sel.selected()) {
        // Create a interpolation object
        DiffTransConfigInterpolation enumerator(config, n_images, endpt_calctype, override_mirrors);

        int length = std::distance(enumerator.m_config_enum_interpol->begin(), enumerator.m_config_enum_interpol->end()) - 2;
        for(i = 0; i < length + 2; ++i) {
          // file_path = $project_dir/training_data/diff_trans/$diff_trans_name/$scelname/$configid/$image_number/POSCAR
          //int n_images = kwargs[config.name()]["n_images"].get<int>(); // set defaults with get_else
          auto file_path = primclex.dir().configuration_calc_dir(config.name(), calctype);
          file_path += "/N_images_" + std::to_string(length) + "/poscars/0" + std::to_string(i) + "/POSCAR";
          fs::create_directories(file_path.parent_path());
          fs::ofstream file(file_path);
          enumerator.at_step(i)->write_pos(file);
        }
        jsonParser endpts_json;
        fs::ofstream endpts(primclex.dir().configuration_calc_dir(config.name(), calctype) / ("/N_images_" + std::to_string(length)) / "endpoint_props.json");
        std::vector<std::string> tokens;
        std::string name = config.from_config().name();
        boost::split(tokens, name, boost::is_any_of("."), boost::token_compress_on);
        std::string canon_config_name = tokens[0];
        //lol
        canon_config_name = config.from_config().canonical_form().name();

        if(canon_config_name.find("none") == std::string::npos) {
          if(tokens.size() == 4) {
            Configuration canon_config = *primclex.db<Configuration>().find(config.from_config().canonical_form().name());
            canon_config.calc_properties(endpt_calctype);
            endpts_json["0"] = apply(config.from_config_from_canonical(), canon_config).print_properties(endpt_calctype);
          }
          else {
            endpts_json["0"] = make_configuration(primclex, config.from_config().name()).print_properties(endpt_calctype);
          }
        }
        else {
          std::cout << "Found 'none' in your from configname" << std::endl;
        }

        std::vector<std::string> tokens2;
        std::string name2 = config.to_config().name();
        boost::split(tokens2, name2, boost::is_any_of("."), boost::token_compress_on);
        std::string canon_config_name2 = tokens2[0];
        //lol
        canon_config_name2 = config.to_config().canonical_form().name();

        if(canon_config_name2.find("none") == std::string::npos) {
          if(tokens2.size() == 4) {
            Configuration canon_config = *primclex.db<Configuration>().find(config.to_config().canonical_form().name());
            canon_config.calc_properties(endpt_calctype);
            endpts_json[std::to_string(i - 1)] = apply(config.to_config_from_canonical(), canon_config).print_properties(endpt_calctype);
          }
          else {
            endpts_json[std::to_string(i - 1)] = make_configuration(primclex, config.to_config().name()).print_properties(endpt_calctype);
          }
        }
        else {
          std::cout << "Found 'none' in your to configname" << std::endl;
        }
        endpts_json.print(endpts);
        endpts.close();
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

        Eigen::Vector3d displacement = config.disp(l);

        const Eigen::Vector3d from_pos = traj.from.uccoord.coordinate().const_cart();
        const Eigen::Vector3d to_pos = traj.to.uccoord.coordinate().const_cart();
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
      Configuration ret_frm = dfc.from_config();
      Configuration ret_to = dfc.to_config();
      if(dfc.from_configname().find("none") == std::string::npos) {
        Configuration rlx_frm = copy_apply_properties(make_configuration(dfc.primclex(), dfc.from_configname()), calctype);
        ret_frm = copy_apply(dfc.from_config_from_canonical(), rlx_frm);
      }
      else {
        dfc.primclex().log() << ">>>>>>!!!!!WARNING: FROM (INITIAL) CONFIGURATION NOT PRESENT IN CONFIG LIST!!!!!" << std::endl
                             << "I suggest enumerating and calculating it before interpolating!!!!!<<<<<<" << std::endl << std::endl;
      }
      if(dfc.to_configname().find("none") == std::string::npos) {
        Configuration rlx_to = copy_apply_properties(make_configuration(dfc.primclex(), dfc.to_configname()), calctype);
        ret_to = copy_apply(dfc.to_config_from_canonical(), rlx_to);
      }
      else {
        dfc.primclex().log() << ">>>>>>!!!!!WARNING: TO (FINAL) CONFIGURATION NOT PRESENT IN CONFIG LIST!!!!!" << std::endl
                             << "I suggest enumerating and calculating it before interpolating!!!!!<<<<<<" << std::endl << std::endl;
      }
      return std::make_pair(ret_frm, ret_to);
    }

    /// applies deformation from input config onto outpur_config and prints it to a file.
    void apply_deformation(const PrimClex primclex, std::string output_configname, std::string output_path,
                           std::string input_configname, std::string calctype) {
      Configuration input_config = copy_apply_properties(make_configuration(primclex, input_configname), calctype);
      Eigen::Matrix3d deformation = input_config.deformation();
      Configuration final_config = make_configuration(primclex, output_configname);
      final_config.set_deformation(deformation);
      fs::ofstream file(output_path);
      final_config.write_pos(file);
      return;
    }
  }
}
