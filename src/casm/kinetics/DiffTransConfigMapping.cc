#include "casm/kinetics/DiffTransConfigMapping.hh"
#include "casm/clex/ConfigMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/completer/Handlers.hh"
#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"
#include "casm/kinetics/DiffTransConfigInterpolation.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/OrbitDecl.hh"
#include "casm/symmetry/OrbitGeneration.hh"
#include "casm/symmetry/OrbitGeneration_impl.hh"
#include "casm/kinetics/DiffusionTransformationTraits.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"
namespace CASM {
  namespace Kinetics {
    //*******************************************************************************************

    DiffTransConfigMapper::DiffTransConfigMapper(const PrimClex &_pclex,
                                                 ConfigMapping::Settings const &_settings,
                                                 double _tol/*=TOL*/) :
      m_pclex(&_pclex),
      m_settings(_settings),
      m_tol(max(1e-9, _tol)) {
      //squeeze lattice_weight into (0,1] if necessary
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
                          this->settings(),
                          primclex().crystallography_tol());
      //Find out which species are moving from which basis site to the other
      if(hint_ptr != nullptr) {
        Configuration tmp = hint_ptr->from_config().canonical_form();
        std::vector<std::string> scelname;
        scelname.push_back(tmp.supercell().name());
        mapper.clear_allowed_lattices();
        mapper.add_allowed_lattices({scelname});
        //	*result.config =*hint_ptr;
      }
      //      else {
      SimpleStructure from_child = make_simple_structure(result.structures[0]);
      std::set<MappingNode> from_set = mapper.struc_mapper().map_deformed_struc(from_child);

      SimpleStructure to_child = make_simple_structure(result.structures.back());
      std::set<MappingNode> to_set = mapper.struc_mapper().map_deformed_struc(to_child);

      if(from_set.empty() || to_set.empty()) {
        return result;
      }
      MappingNode from_node = *from_set.begin();
      MappingNode to_node = *to_set.begin();

      auto scel_ptr = std::make_shared<Supercell>(&(mapper.primclex()), from_node.lat_node.parent.scel_lattice());
      std::vector<UnitCellCoord> from_uccoords;
      std::vector<UnitCellCoord> to_uccoords;

      std::set<UnitCellCoord> vacancy_from;
      std::set<UnitCellCoord> vacancy_to;
      //This rigid rotation and rigid shift seems unnecessary surprisingly
      Coordinate first_site = Coordinate(result.structures[0].basis(0).const_frac(), scel_ptr->lattice(), FRAC);

      Coordinate t_shift =
        copy_apply(SymOp::point_op(from_node.lat_node.isometry), first_site)
        - Coordinate(scel_ptr->uccoord(find_index(from_node.atom_permutation, 0)).coordinate(this->primclex().prim()));

      t_shift.set_lattice(primclex().prim().lattice(), CART);

      //std::cout << "t shift " << t_shift.const_frac() << std::endl;
      //Maybe check coordinate similarity after applying deformations
      //Check the unitcell coordinate within a tolerance of the maxium displacement of any atom in the from config
      //This max_displacement is not considering rigid translational shifts of the structures's basis to the primclex's basis
      BasicStructure<Site> from = result.structures[0];
      BasicStructure<Site> to = result.structures[result.structures.size() - 1];
      //std::cout << "unconditioned from struc" << std::endl;
      //from.print_xyz(std::cout,true);
      //std::cout << "unconditioned to struc" << std::endl;
      //to.print_xyz(std::cout,true);
      _precondition_from_and_to(from_node.lat_node.isometry,
                                from_node.lat_node.stretch,
                                t_shift.const_cart(),
                                from,
                                to);

      //std::cout << "trans" << t_shift.const_cart() << std::endl;
      //std::cout << "from struc" << std::endl;
      //from.print_xyz(std::cout,true);
      //std::cout << "to struc" << std::endl;
      //to.print_xyz(std::cout,true);
      double uccoord_tol = 1.1 * max(from_node.atom_displacement.colwise().norm().maxCoeff(), to_node.atom_displacement.colwise().norm().maxCoeff()) + cuberoot(abs(scel_ptr->lattice().vol())) * primclex().crystallography_tol();
      //std::cout << "uccoord_tol " << uccoord_tol <<std::endl;
      std::vector<Index> moving_atoms = _analyze_atoms_ideal(from,
                                                             to,
                                                             *scel_ptr,
                                                             uccoord_tol,
                                                             from_uccoords,
                                                             to_uccoords,
                                                             vacancy_from,
                                                             vacancy_to);
      //std::cout << " MOVING ATOM IS GOING FROM " << from_uccoords[moving_atoms[0]]<< " TO " << to_uccoords[moving_atoms[0]]    << std::endl;
      Kinetics::DiffusionTransformation diff_trans = _make_hop(result.structures[0],
                                                               from_uccoords,
                                                               to_uccoords,
                                                               vacancy_from,
                                                               vacancy_to,
                                                               moving_atoms);
      //std::cout << "possibly not shortest hop" << diff_trans << std::endl;
      //diff_trans constructed from end points might not be shortest path
      //in order to set up shortest path create eigen counter that represents all adjacent supercells
      //to which the unitcell coords will be reachable from the from uccoords
      // for each item in counter take item wise product with maximum uccoord in from configs scel
      // + 1 to account for starting at 0
      // make temp storage "final_diff_trans" replace if max length is smaller for current diff_trans
      //

      {
        SimpleStructure oriented_struc = mapper.struc_mapper().calculator().resolve_setting(from_node, from_child);

        auto tdof = to_configdof(oriented_struc, *scel_ptr);
        Configuration from_config(scel_ptr, jsonParser(), tdof.first);

        //Configuration from_config(Configuration::zeros(scel_ptr));
        //from_config.set_occupation(mapper.occupation(make_simple_structure(result.structures[0]),
        //from_node));
        Kinetics::DiffusionTransformation final_diff_trans = _shortest_hop(diff_trans, from_config.supercell());
        //THIS IS THE FIRST CASE IN WHICH WE DON'T WANT TO SORT DIFFTRANS ON CONSTRUCTION -speak with brian about removing sorting from prepare
        //or stick with prepareless workaround. Alternatively check if sorted, if not then sort diff trans and flip from/to then create.
        //Need to somehow only use occupation from the config to construct diff_trans_config
        diff_trans = final_diff_trans;
        //std::cout << "diff_trans is " << std::endl;
        //std::cout << diff_trans << std::endl;
        //Attach hop to ideal from config in same orientation

        //std::cout << "from config is " << std::endl << from_config << std::endl;
        result.config = notstd::make_unique<Kinetics::DiffTransConfiguration>(from_config, diff_trans);
        if(!result.config->has_valid_from_occ()) {
          throw std::runtime_error("Moving forward with invalid diff_trans_config is a bad idea.");
        }
      }
      //std::cout << "to config is " << std::endl << non_canon_to_config << std::endl;
      //NEED TO SET ORBIT NAME OF DTC SOMEHOW FOR NEW DIFF TRANS

      PrimPeriodicDiffTransSymCompare sym_c {this->primclex().shared_prim(), primclex().crystallography_tol()};
      OrbitGenerators<PrimPeriodicDiffTransOrbit> generators(primclex().prim().factor_group(), sym_c);
      generators.insert(sym_c.prepare(diff_trans));
      std::vector<PrimPeriodicDiffTransOrbit> dt_orbits;
      generators.make_orbits(std::back_inserter(dt_orbits), primclex());
      auto insert_res = primclex().db<PrimPeriodicDiffTransOrbit>().insert(dt_orbits.front());
      primclex().db<PrimPeriodicDiffTransOrbit>().commit();
      result.config->set_orbit_name(insert_res.first.name());
      //    }

      //use this to interpolate same amount of images

      Kinetics::DiffTransConfigInterpolation interpolater(result.config->diff_trans(), result.config->from_config(), result.config->to_config(), result.structures.size() - 2); //<- using current calctype here
      int image_no = 0;
      for(auto it = interpolater.begin(); it != interpolater.end(); ++it) {
        //result.relaxation_properties.push_back(jsonParser());
        //result.relaxation_properties[image_no].put_obj();
        //Structure pseudoprim = make_deformed_struc(*it);
        //Structure img = Structure(result.structures[image_no]);
        // ConfigMapperResult tmp_result = ConfigMapping::structure_mapping(pseudoprim, img, lattice_weight());
        //result.relaxation_properties[image_no]["lattice_deformation"] = tmp_result.relaxation_properties["best_mapping"]["lattice_deformation"];
        //result.relaxation_properties[image_no]["lattice_deformation"] = -1.0;
        //result.relaxation_properties[image_no]["basis_deformation"] = tmp_result.relaxation_properties["best_mapping"]["basis_deformation"];
        //result.relaxation_properties[image_no]["basis_deformation"] = -1.0;
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
          //result.relaxation_properties[count]["relaxed_energy"] = img["relaxed_energy"];
          count++;
        }
        if(!all_strucs.contains("kra")) {
          result.kra = *(std::max_element(energies.begin(), energies.end())) - (energies.front() + energies.back()) / 2.0;
        }
        else {
          result.kra = all_strucs["kra"].get<double>();
        }
      }
      result.success = true;
      return result;
    }

    std::vector<Index> DiffTransConfigMapper::_analyze_atoms_ideal(const BasicStructure<Site> &from,
                                                                   const BasicStructure<Site> &to,
                                                                   const Supercell &scel,
                                                                   double uccoord_tol,
                                                                   std::vector<UnitCellCoord> &from_uccoords,
                                                                   std::vector<UnitCellCoord> &to_uccoords,
                                                                   std::set<UnitCellCoord> &vacancy_from,
                                                                   std::set<UnitCellCoord> &vacancy_to) const {
      // For image 00 set reference of POSCAR index to  basis site linear index
      // tolerance for UnitCellCoord mapping has 20% wiggle room from max displacement
      // instead of introducing wiggle room maybe take max disp between from map and to map
      int from_count = 0;
      for(auto &site : from.basis()) {
        from_uccoords.push_back(scel.prim_grid().within(_site_to_uccoord(site, primclex(), uccoord_tol)));
        from_count++;
      }

      // For last image  find POSCAR index to basis site linear index
      int to_count = 0;
      for(auto &site : to.basis()) {
        to_uccoords.push_back(scel.prim_grid().within(_site_to_uccoord(site, primclex(), uccoord_tol)));
        to_count++;
      }
      std::vector<Index> moving_atoms;
      for(int i = 0 ; i < from_uccoords.size(); i++) {
        if(from_uccoords[i] != to_uccoords[i]) {
          moving_atoms.push_back(i);
        }
      }
      ////if this isn't a closed loop one of the species is a vacancy
      for(int i = 0 ; i < moving_atoms.size() ; i++) {
        if(vacancy_from.find(to_uccoords[moving_atoms[i]]) == vacancy_from.end()) {
          vacancy_from.insert(to_uccoords[moving_atoms[i]]);
        }
        else {
          std::cout << "There should be only 1 vacancy in hop!" << std::endl;
          vacancy_from.erase(vacancy_from.find(to_uccoords[moving_atoms[i]]));
        }
        if(vacancy_to.find(from_uccoords[moving_atoms[i]]) == vacancy_to.end()) {
          vacancy_to.insert(from_uccoords[moving_atoms[i]]);
        }
        else {
          std::cout << "There should be only 1 vacancy in hop!" << std::endl;
          vacancy_to.erase(vacancy_to.find(from_uccoords[moving_atoms[i]]));
        }
      }
      return moving_atoms;
    }

    UnitCellCoord DiffTransConfigMapper::_site_to_uccoord(const Site &site, const PrimClex &pclex, double tol) const {
      return UnitCellCoord::from_coordinate(pclex.prim(), site, tol);
    }

    void DiffTransConfigMapper::_precondition_from_and_to(const Eigen::Matrix3d &cart_op, const Eigen::Matrix3d &strain, const Eigen::Vector3d &trans, BasicStructure<Site> &from, BasicStructure<Site> &to) const {
      from.set_lattice(Lattice(strain.inverse() * (cart_op.transpose()*from.lattice().lat_column_mat())), FRAC);
      from += Coordinate(trans, from.lattice(), CART);
      for(auto &site : from.set_basis()) {
        site.within();
      }
      to.set_lattice(Lattice(strain.inverse() * (cart_op.transpose()*to.lattice().lat_column_mat())), FRAC);
      to += Coordinate(trans, to.lattice(), CART);
      for(auto &site : to.set_basis()) {
        site.within();
      }
      return;
    }

    Kinetics::DiffusionTransformation DiffTransConfigMapper::_shortest_hop(const Kinetics::DiffusionTransformation &diff_trans, const Supercell &scel) const {
      Kinetics::DiffusionTransformation final_diff_trans = diff_trans;
      EigenCounter<Eigen::Vector3l> counter(Eigen::Vector3l::Constant(-1), Eigen::Vector3l::Constant(1), Eigen::Vector3l::Ones());
      for(int i = 0 ; i < diff_trans.occ_transform().size(); i++) {
        Eigen::Matrix3d lat_mat = scel.lattice().lat_column_mat();
        while(counter.valid()) {
          UnitCell shift = lround(scel.prim().lattice().inv_lat_column_mat() * lat_mat * counter().cast<double>());
          Kinetics::DiffusionTransformation tmp = diff_trans;
          UnitCellCoord replace_this = tmp.occ_transform()[i].uccoord;
          tmp.occ_transform()[i].uccoord = replace_this + shift;
          for(auto &traj : tmp.species_traj()) {
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

    Kinetics::DiffusionTransformation DiffTransConfigMapper::_make_hop(const BasicStructure<Site> &from_struc,
                                                                       const std::vector<UnitCellCoord> &from_coords,
                                                                       const std::vector<UnitCellCoord> &to_coords,
                                                                       const std::set<UnitCellCoord> &vacancy_from,
                                                                       const std::set<UnitCellCoord> &vacancy_to,
                                                                       const std::vector<Index> &moving_atoms) const {
      //From the moving species and basis sites, should be able to create hop
      Kinetics::DiffusionTransformation diff_trans(primclex().prim());
      for(int i = 0; i < moving_atoms.size(); i++) {
        diff_trans.occ_transform().emplace_back(from_coords[moving_atoms[i]], 0, 0);
      }
      if(vacancy_from.size() && vacancy_to.size()) {
        diff_trans.occ_transform().emplace_back(*vacancy_from.begin(), 0, 0);
      }
      for(int i = 0; i < moving_atoms.size(); i++) {
        std::vector<std::string> allowed_from_occs = primclex().prim().basis()[from_coords[moving_atoms[i]].sublattice()].allowed_occupants();
        Index from_occ_index = std::distance(allowed_from_occs.begin(),
                                             std::find(allowed_from_occs.begin(),
                                                       allowed_from_occs.end(),
                                                       from_struc.basis()[moving_atoms[i]].occ_name()));
        //for now pos is 0 because Molecules are hard
        Kinetics::SpeciesLocation from_loc(from_coords[moving_atoms[i]], from_occ_index, 0);
        std::vector<std::string> allowed_to_occs = primclex().prim().basis()[to_coords[moving_atoms[i]].sublattice()].allowed_occupants();
        Index to_occ_index = std::distance(allowed_to_occs.begin(),
                                           std::find(allowed_to_occs.begin(),
                                                     allowed_to_occs.end(),
                                                     from_struc.basis()[moving_atoms[i]].occ_name()));

        //for now pos is 0 because Molecules are hard
        Kinetics::SpeciesLocation to_loc(to_coords[moving_atoms[i]], to_occ_index, 0);
        diff_trans.species_traj().emplace_back(from_loc, to_loc);
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
        std::vector<std::string> allowed_from_occs = primclex().prim().basis()[vacancy_from.begin()->sublattice()].allowed_occupants();
        Index from_occ_index = std::distance(allowed_from_occs.begin(), std::find(allowed_from_occs.begin(), allowed_from_occs.end(), "Va"));
        Kinetics::SpeciesLocation from_loc(*vacancy_from.begin(), from_occ_index, 0);
        std::vector<std::string> allowed_to_occs = primclex().prim().basis()[vacancy_to.begin()->sublattice()].allowed_occupants();
        Index to_occ_index = std::distance(allowed_to_occs.begin(), std::find(allowed_to_occs.begin(), allowed_to_occs.end(), "Va"));
        Kinetics::SpeciesLocation to_loc(*vacancy_to.begin(), to_occ_index, 0);
        diff_trans.species_traj().emplace_back(from_loc, to_loc);
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
            throw std::runtime_error("DiffTransConfigMapper must be re-implemented to use SimpleStructure");
            //int img_no = std::stoi(it.name());
            //from_json(simple_json(struc, "relaxed_"), *it);
            //bins.insert(std::make_pair(img_no, struc));
          }
          catch(std::invalid_argument) {
          }
        }
      }
      else {
        fs::directory_iterator end;
        for(fs::directory_iterator begin(pos_path); begin != end; ++begin) {
          try {
            int img_no = std::stoi(begin->path().filename().string());
            if(fs::is_directory(*begin)) {
              if(fs::is_regular(*begin / "CONTCAR")) {
                bins.insert(std::make_pair(img_no, BasicStructure<Site>(*begin / "CONTCAR")));
              }
              else if(fs::is_regular(*begin / "POSCAR")) {
                bins.insert(std::make_pair(img_no, BasicStructure<Site>(*begin / "POSCAR")));
              }
              else {
                std::cerr << "NO POSCAR OR CONTCAR FOUND IN " << *begin << std::endl;
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
