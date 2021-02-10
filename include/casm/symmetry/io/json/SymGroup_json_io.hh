#ifndef CASM_symmetry_io_json_SymGroup
#define CASM_symmetry_io_json_SymGroup

#include "casm/global/definitions.hh"

namespace CASM {

class SymGroup;
class SymGroupRepHandle;
class SymGroupRepID;
class jsonParser;

void write_symop(SymGroup const &grp, Index i, jsonParser &j);

void write_symgroup(SymGroup const &grp, jsonParser &json);

void write_basis_permutation_rep(SymGroup const &grp, jsonParser &json,
                                 SymGroupRepID symgrouprep_id);

void write_occ_permutation_rep(SymGroup const &grp, jsonParser &json,
                               std::vector<SymGroupRepID> occupant_symrep_IDs);

void write_matrix_rep(SymGroupRepHandle const &grp, jsonParser &json);

}  // namespace CASM

#endif
