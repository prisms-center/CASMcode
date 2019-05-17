#include "casm/clex/ConfigMapping.hh"

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/clex/PrimClex.hh"
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


  namespace ConfigMapping {
    ConfigMapperResult structure_mapping(Structure &host, Structure &other, double lattice_weight) {
      const PrimClex &pclex = PrimClex(host, null_log());
      ConfigMapper tmp_mapper(pclex, lattice_weight, 0.0);
      //tmp_mapper.set_max_va_frac(0.0);
      return tmp_mapper.import_structure(SimpleStructure(other));
    }
  }

  PrimStrucMapCalculator::PrimStrucMapCalculator(BasicStructure<Site> const &_prim,
                                                 StrucMapping::SpeciesMode _species_mode/*=StrucMapping::ATOM*/) :
    StrucMapCalculatorInterface(SimpleStructure(_prim)),
    m_prim(_prim),
    m_species_mode(_species_mode) {

    SymGroup fg;
    _prim.generate_factor_group(fg);
    set_point_group(fg.copy_no_trans());
    ParamComposition param_comp(m_prim);
    m_fixed_components = param_comp.fixed_components();
  }
  //*******************************************************************************************

  ConfigMapper::ConfigMapper(const PrimClex &_pclex,
                             double _lattice_weight,
                             double _max_volume_change/*=0.5*/,
                             int options/*=robust*/,
                             double _tol/*=-1.*/) :
    m_pclex(&_pclex),
    m_struc_mapper(PrimStrucMapCalculator(_pclex.prim()),
                   _lattice_weight,
                   1,//_Nbest,
                   //std::vector<SymOp>({SymOp()}),
                   _max_volume_change,
                   options,
                   _tol > 0. ? _tol : _pclex.crystallography_tol()) {

  }

}
