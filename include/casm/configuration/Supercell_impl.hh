#ifndef CASM_config_Supercell_impl
#define CASM_config_Supercell_impl

#include "casm/configuration/Supercell.hh"
#include "casm/configuration/SupercellSymInfo_impl.hh"
#include "casm/cry.hh"

namespace CASM {
namespace config {

Supercell::Supercell(std::shared_ptr<Prim const> const &_prim,
                     Superlattice const &_superlattice)
    : prim(_prim),
      superlattice(_superlattice),
      unitcell_index_converter(_superlattice.transformation_matrix_to_super()),
      unitcellcoord_index_converter(
          _superlattice.transformation_matrix_to_super(), _prim->n_sublat),
      sym_info(_prim, _superlattice) {}

}  // namespace config
}  // namespace CASM

#endif
