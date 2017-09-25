#include "casm/kinetics/DiffTransConfigMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/crystallography/Lattice.hh"
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
    //result.structures = _get_structures(pos_path);
    result.relaxation_properties.put_obj();

    //Find out which species are moving from which basis site to the other
    ////if this isn't a closed loop one of the species is a vacancy
    //From the moving species and basis sites, should be able to create hop
    //Map first and last structure using ConfigMapping
    ////first is the from config last is the to config <--(This may get flipped flopped upon sorting of DiffTransConfiguration)
    //Attach hop to ideal from config in same orientation
    //use this to interpolate same amount of images
    //calculate strain scores and basis scores for every image and sum/average/sumsq
    // set relaxation properties and indicate successful mapping or not



    return result;
  }

}
//*******************************************************************************************
