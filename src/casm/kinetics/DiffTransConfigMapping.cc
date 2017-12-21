#include "casm/kinetics/DiffTransConfigMapping.hh"
#include "casm/clex/ConfigMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/completer/Handlers.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/kinetics/DiffTransConfigInterpolation.hh"
#include "casm/crystallography/jsonStruc.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/OrbitDecl.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/kinetics/DiffusionTransformationTraits.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"
namespace CASM {
  namespace Kinetics {
    //*******************************************************************************************

    DiffTransConfigMapper::DiffTransConfigMapper(const PrimClex &_pclex,
                                                 double _lattice_weight,
                                                 double _max_volume_change/*=0.5*/,
                                                 int options/*=robust*/,
                                                 double _tol/*=TOL*/) :
      m_pclex(&_pclex),
      m_lattice_weight(_lattice_weight),
      m_max_volume_change(_max_volume_change),
      m_min_va_frac(0.),
      m_max_va_frac(1.),
      m_robust_flag(options & robust),
      m_strict_flag(options & strict),
      m_rotate_flag(options & rotate),
      m_tol(max(1e-9, _tol)) {
      //squeeze lattice_weight into (0,1] if necessary
      m_lattice_weight = max(min(_lattice_weight, 1.0), 1e-9);
      ParamComposition param_comp(_pclex.prim());
      m_fixed_components = param_comp.fixed_components();
      m_max_volume_change = max(m_tol, _max_volume_change);
    }

    //*******************************************************************************************

    DiffTransConfigMapperResult DiffTransConfigMapper::import_structure_occupation(const fs::path &pos_path) const {
      return import_structure_occupation(pos_path, nullptr);
    }

    //*******************************************************************************************

    DiffTransConfigMapperResult DiffTransConfigMapper::import_structure_occupation(
      const fs::path &pos_path,
      const Kinetics::DiffTransConfiguration *hint_ptr) const {

      DiffTransConfigMapperResult result;
      result.structures = _get_structures(pos_path);
      ConfigMapper mapper(primclex(),
                          lattice_weight(),
                          m_max_volume_change,
                          m_robust_flag || m_rotate_flag || m_strict_flag,
                          primclex().crystallography_tol());
      //Find out which species are moving from which basis site to the other

      ConfigDoF from_dof;
      Lattice from_lat;
      mapper.struc_to_configdof(result.structures[0], from_dof, from_lat);
      auto scel_ptr = std::make_shared<Supercell>(&(mapper.primclex()), from_lat);
      Configuration from_config(scel_ptr, jsonParser(), from_dof);
      std::vector<UnitCellCoord> from_uccoords;
      std::vector<UnitCellCoord> to_uccoords;
      ConfigMapperResult from_res;
      if(hint_ptr != nullptr) {
        Configuration tmp = hint_ptr->from_config().canonical_form();
        from_res = mapper.import_structure_occupation(result.structures[0], &(tmp));
      }
      else {
        from_res = mapper.import_structure_occupation(result.structures[0]);
      }
      //This rigid rotation and rigid shift seems unnecessary surprisingly
      /*
         SymOp op(from_res.cart_op);
          Coordinate rigid_shift = copy_apply(op,result.structures[0].basis[0]) - Coordinate(from_config.uccoord(from_res.best_assignment[0]));
      Coordinate t_shift(primclex().prim().lattice());
      t_shift.cart() = rigid_shift.cart();
      */
      //Maybe check coordinate similarity after applying deformations
      //Check the unitcell coordinate within a tolerance of the maxium displacement of any atom in the from config
      //This max_displacement is not considering rigid translational shifts of the structures's basis to the primclex's basis
      double max_displacement = from_config.displacement().colwise().norm().maxCoeff() * 2 + primclex().crystallography_tol();
      // For image 00 set reference of POSCAR index to  basis site linear index
      for(auto &site : result.structures[0].basis) {
        from_uccoords.emplace_back(primclex().prim(), site, max_displacement);
      }

      // For last image  find POSCAR index to basis site linear index
      for(auto &site : result.structures[result.structures.size() - 1].basis) {
        to_uccoords.emplace_back(primclex().prim(), site, max_displacement);
      }
      //ConfigMapperResult to_config_result = mapper.import_structure_occupation(result.structures[result.structures.size()-1]);
      std::vector<Index> moving_atoms;
      for(int i = 0 ; i < from_uccoords.size(); i++) {
        if(from_uccoords[i] != to_uccoords[i]) {
          moving_atoms.push_back(i);
        }
      }
      ////if this isn't a closed loop one of the species is a vacancy
      std::set<UnitCellCoord> vacancy_from;
      std::set<UnitCellCoord> vacancy_to;
      for(int i = 0 ; i < moving_atoms.size() ; i++) {
        if(vacancy_from.find(to_uccoords[moving_atoms[i]]) == vacancy_from.end()) {
          vacancy_from.insert(to_uccoords[moving_atoms[i]]);
        }
        else {
          vacancy_from.erase(vacancy_from.find(to_uccoords[moving_atoms[i]]));
        }
        if(vacancy_to.find(from_uccoords[moving_atoms[i]]) == vacancy_to.end()) {
          vacancy_to.insert(from_uccoords[moving_atoms[i]]);
        }
        else {
          vacancy_to.erase(vacancy_to.find(from_uccoords[moving_atoms[i]]));
        }
      }
      Kinetics::DiffusionTransformation diff_trans = _make_hop(result.structures[0],
                                                               from_uccoords,
                                                               to_uccoords,
                                                               vacancy_from,
                                                               vacancy_to,
                                                               moving_atoms);
      //diff_trans constructed from end points might not be shortest path
      //in order to set up shortest path create eigen counter that represents all adjacent supercells
      //to which the unitcell coords will be reachable from the from uccoords
      // for each item in counter take item wise product with maximum uccoord in from configs scel
      // + 1 to account for starting at 0
      // make temp storage "final_diff_trans" replace if max length is smaller for current diff_trans
      //
      Kinetics::DiffusionTransformation final_diff_trans = _shortest_hop(diff_trans, from_config.supercell());
      //THIS IS THE FIRST CASE IN WHICH WE DON'T WANT TO SORT DIFFTRANS ON CONSTRUCTION -speak with brian about removing sorting from prepare
      //or stick with prepareless workaround. Alternatively check if sorted, if not then sort diff trans and flip from/to then create.
      //Need to somehow only use occupation from the config to construct diff_trans_config
      diff_trans = final_diff_trans;
      //Attach hop to ideal from config in same orientation
      from_config.clear_deformation();
      from_config.clear_displacement();
      from_config.init_deformation();
      from_config.init_displacement();
      result.config = notstd::make_unique<Kinetics::DiffTransConfiguration>(from_config, diff_trans);
      //NEED TO SET ORBIT NAME OF DTC SOMEHOW FOR NEW DIFF TRANS
      PrimPeriodicDiffTransSymCompare sym_c {primclex().crystallography_tol()};
      OrbitGenerators<PrimPeriodicDiffTransOrbit> generators(primclex().prim().factor_group(), sym_c);
      generators.insert(sym_c.prepare(diff_trans));
      std::vector<PrimPeriodicDiffTransOrbit> dt_orbits;
      generators.make_orbits(std::back_inserter(dt_orbits), primclex());
      auto insert_res = primclex().db<PrimPeriodicDiffTransOrbit>().insert(dt_orbits.front());
      primclex().db<PrimPeriodicDiffTransOrbit>().commit();
      result.config->set_orbit_name(insert_res.first.name());
      //use this to interpolate same amount of images
      Kinetics::DiffTransConfigInterpolation interpolater(result.config->diff_trans(), result.config->from_config(), result.config->to_config(), result.structures.size() - 2); //<- using current calctype here
      int image_no = 0;
      for(auto it = interpolater.begin(); it != interpolater.end(); ++it) {
        result.relaxation_properties.push_back(jsonParser());
        result.relaxation_properties[image_no].put_obj();
        Structure pseudoprim = make_deformed_struc(*it);
        Structure img = Structure(result.structures[image_no]);
        ConfigMapperResult tmp_result = ConfigMapping::structure_mapping(pseudoprim, img, lattice_weight());
        result.relaxation_properties[image_no]["lattice_deformation"] = tmp_result.relaxation_properties["best_mapping"]["lattice_deformation"];
        result.relaxation_properties[image_no]["basis_deformation"] = tmp_result.relaxation_properties["best_mapping"]["basis_deformation"];
        image_no++;
      }

      //Structure config.supercell().superstructure(config) //<---how to get structure from ideal config
      //calculate strain scores and basis scores for every image and sum/average/sumsq
      // set relaxation properties and indicate successful mapping or not
      if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
        jsonParser all_strucs;
        to_json(pos_path, all_strucs);
        int count = 0;
        std::vector<double> energies;
        for(auto &img : all_strucs) {
          energies.push_back(img["relaxed_energy"].get<double>());
          result.relaxation_properties[count]["relaxed_energy"] = img["relaxed_energy"];
          count++;
        }
        if(!all_strucs.contains("kra")) {
          all_strucs["kra"] = *(std::max_element(energies.begin(), energies.end())) - (energies.front() + energies.back()) / 2.0;
        }
        result.kra = all_strucs["kra"].get<double>();
      }
      result.success = true;
      return result;
    }


    Kinetics::DiffusionTransformation DiffTransConfigMapper::_shortest_hop(Kinetics::DiffusionTransformation &diff_trans, const Supercell &scel) const {
      Kinetics::DiffusionTransformation final_diff_trans = diff_trans;
      EigenCounter<Eigen::Vector3l> counter(Eigen::Vector3l::Constant(-1), Eigen::Vector3l::Constant(1), Eigen::Vector3l::Ones());
      for(int i = 0 ; i < diff_trans.occ_transform().size(); i++) {
        while(counter.valid()) {
          Eigen::Vector3l scelvec = (scel.uccoord(scel.num_sites() - 1).unitcell() + Eigen::Vector3l::Constant(1));
          UnitCell shift = counter().cwiseProduct(scelvec);
          Kinetics::DiffusionTransformation tmp = diff_trans;
          UnitCellCoord replace_this = tmp.occ_transform()[i].uccoord;
          tmp.occ_transform()[i].uccoord = replace_this + shift;
          for(auto &traj : tmp.specie_traj()) {
            if(traj.from.uccoord == replace_this) {
              traj.from.uccoord = replace_this + shift;
            }
            if(traj.to.uccoord == replace_this) {
              traj.to.uccoord = replace_this + shift;
            }
          }
          if(tmp.max_length() <= final_diff_trans.max_length()) {
            final_diff_trans = tmp;
          }
          counter++;
        }
      }
      return final_diff_trans;
    }

    Kinetics::DiffusionTransformation DiffTransConfigMapper::_make_hop(BasicStructure<Site> &from_struc,
                                                                       std::vector<UnitCellCoord> &from_coords,
                                                                       std::vector<UnitCellCoord> &to_coords,
                                                                       std::set<UnitCellCoord> &vacancy_from,
                                                                       std::set<UnitCellCoord> &vacancy_to,
                                                                       std::vector<Index> &moving_atoms) const {
      //From the moving species and basis sites, should be able to create hop
      Kinetics::DiffusionTransformation diff_trans(primclex().prim());
      for(int i = 0; i < moving_atoms.size(); i++) {
        diff_trans.occ_transform().emplace_back(from_coords[moving_atoms[i]], 0, 0);
      }
      if(vacancy_from.size() && vacancy_to.size()) {
        diff_trans.occ_transform().emplace_back(*vacancy_from.begin(), 0, 0);
      }
      for(int i = 0; i < moving_atoms.size(); i++) {
        std::vector<std::string> allowed_from_occs = primclex().prim().basis[from_coords[moving_atoms[i]].sublat()].allowed_occupants();
        Index from_occ_index = std::distance(allowed_from_occs.begin(), std::find(allowed_from_occs.begin(), allowed_from_occs.end(), from_struc.basis[moving_atoms[i]].occ_name()));
        //for now pos is 0 because Molecules are hard
        Kinetics::SpecieLocation from_loc(from_coords[moving_atoms[i]], from_occ_index, 0);
        std::vector<std::string> allowed_to_occs = primclex().prim().basis[to_coords[moving_atoms[i]].sublat()].allowed_occupants();
        Index to_occ_index = std::distance(allowed_to_occs.begin(), std::find(allowed_to_occs.begin(), allowed_to_occs.end(), from_struc.basis[moving_atoms[i]].occ_name()));
        //for now pos is 0 because Molecules are hard
        Kinetics::SpecieLocation to_loc(to_coords[moving_atoms[i]], to_occ_index, 0);
        diff_trans.specie_traj().emplace_back(from_loc, to_loc);
        for(auto &occ_trans : diff_trans.occ_transform()) {
          if(occ_trans.uccoord == from_coords[moving_atoms[i]]) {
            occ_trans.from_value = from_occ_index;
          }
          if(occ_trans.uccoord == to_coords[moving_atoms[i]]) {
            occ_trans.to_value = to_occ_index;
          }
        }
      }
      if(vacancy_from.size() && vacancy_to.size()) {
        std::vector<std::string> allowed_from_occs = primclex().prim().basis[vacancy_from.begin()->sublat()].allowed_occupants();
        Index from_occ_index = std::distance(allowed_from_occs.begin(), std::find(allowed_from_occs.begin(), allowed_from_occs.end(), "Va"));
        Kinetics::SpecieLocation from_loc(*vacancy_from.begin(), from_occ_index, 0);
        std::vector<std::string> allowed_to_occs = primclex().prim().basis[vacancy_to.begin()->sublat()].allowed_occupants();
        Index to_occ_index = std::distance(allowed_to_occs.begin(), std::find(allowed_to_occs.begin(), allowed_to_occs.end(), "Va"));
        Kinetics::SpecieLocation to_loc(*vacancy_to.begin(), to_occ_index, 0);
        diff_trans.specie_traj().emplace_back(from_loc, to_loc);
        for(auto &occ_trans : diff_trans.occ_transform()) {
          if(occ_trans.uccoord == *vacancy_from.begin()) {
            occ_trans.from_value = from_occ_index;
          }
          if(occ_trans.uccoord == *vacancy_to.begin()) {
            occ_trans.to_value = to_occ_index;
          }
        }
      }
      return diff_trans;
    }


    std::vector<BasicStructure<Site>> DiffTransConfigMapper::_get_structures(const fs::path &pos_path) const {
      std::map<Index, BasicStructure<Site>> bins;
      std::vector<BasicStructure<Site>> images;
      if(pos_path.extension() == ".json" || pos_path.extension() == ".JSON") {
        jsonParser all_strucs;
        to_json(pos_path, all_strucs);
        for(auto it =  all_strucs.begin(); it != all_strucs.end(); ++it) {
          BasicStructure<Site> struc;
          try {
            int img_no = std::stoi(it.name());
            from_json(simple_json(struc, "relaxed_"), *it);
            bins.insert(std::make_pair(img_no, struc));
          }
          catch(std::invalid_argument) {
          }
        }
      }
      /* else if(fs::exists(pos_path / "properties.calc.json")) {
         jsonParser all_strucs;
         to_json(pos_path / "properties.calc.json", all_strucs);
         int count = 0;
         for(auto &img : all_strucs) {
           BasicStructure<Site> struc;
           from_json(simple_json(struc, "relaxed_"), img);
           bins.insert(std::make_pair(count, struc));
           count++;
         }
       }*/
      else {
        for(auto &dir_path : fs::directory_iterator(pos_path)) {
          try {
            int img_no = std::stoi(dir_path.path().filename().string());
            if(fs::is_directory(dir_path)) {
              if(fs::is_regular(dir_path / "CONTCAR")) {
                bins.insert(std::make_pair(img_no, BasicStructure<Site>(dir_path / "CONTCAR")));
              }
              else if(fs::is_regular(dir_path / "POSCAR")) {
                bins.insert(std::make_pair(img_no, BasicStructure<Site>(dir_path / "POSCAR")));
              }
              else {
                std::cerr << "NO POSCAR OR CONTCAR FOUND IN " << dir_path << std::endl;
              }
            }
          }
          catch(...) {
          }
        }
      }
      for(int i = 0 ; i < bins.size(); i++) {
        try {
          images.push_back(bins[i]);
        }
        catch(...) {
          std::cerr << "IMAGE NUMBERS NOT CONSECUTIVE IN " << pos_path << std::endl;
        }
      }
      return images;
    }

  }
}
//*******************************************************************************************
