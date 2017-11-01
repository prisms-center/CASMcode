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

namespace CASM {

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
    result.relaxation_properties.put_obj();
    ConfigMapper mapper(primclex(),
                        lattice_weight(),
                        m_max_volume_change,
                        m_robust_flag || m_rotate_flag || m_strict_flag,
                        primclex().crystallography_tol());
    //Find out which species are moving from which basis site to the other

    ConfigDoF from_dof;
    Lattice from_lat;
    mapper.struc_to_configdof(result.structures[0], from_dof, from_lat);
    Supercell scel(&(mapper.primclex()), from_lat);
    Configuration from(scel, jsonParser(), from_dof);
    std::vector<UnitCellCoord> from_uccoords;
    std::vector<UnitCellCoord> to_uccoords;
    //Maybe check coordinate similarity after applying deformations
    //Check the unitcell coordinate within a tolerance of the maxium displacement of any atom in the from config
    double max_displacement = primclex().crystallography_tol();
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
      std::cout << from_uccoords[i] << to_uccoords[i] << std::endl;
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
    if(vacancy_from.size() && vacancy_to.size()) {
      std::cout << "Vacancy path" << std::endl;
      std::cout << *vacancy_from.begin() << *vacancy_to.begin() << std::endl;
    }
    //From the moving species and basis sites, should be able to create hop

    //Map first and last structure using ConfigMapping
    ////first is the from config last is the to config <--(This may get flipped flopped upon sorting of DiffTransConfiguration)
    //Attach hop to ideal from config in same orientation
    //use this to interpolate same amount of images
    //calculate strain scores and basis scores for every image and sum/average/sumsq
    // set relaxation properties and indicate successful mapping or not



    return result;
  }

  std::vector<BasicStructure<Site>> DiffTransConfigMapper::_get_structures(const fs::path &pos_path) const {
    std::map<Index, BasicStructure<Site>> bins;
    std::vector<BasicStructure<Site>> images;
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
//*******************************************************************************************
