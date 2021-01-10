#ifndef CASM_DoFSpace_impl
#define CASM_DoFSpace_impl

#include "casm/crystallography/Structure.hh"
#include "casm/clex/Supercell.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {

  /// Make VectorSpaceSymReport
  ///
  /// \param dof_space DoFSpace to make VectorSpaceSymReport for
  /// \param input_state ConfigEnumInput, whose invariant group w.r.t. [group_begin, group_end) is
  ///        is used to specify the symmetry of the vector space to be analyzed
  /// \param group_begin, group_end Range of PermuteIterator
  /// \param calc_wedges If true, calculate the irreducible wedges for the vector space. This may
  ///        take a long time.
  template<typename PermuteIteratorIt>
  VectorSpaceSymReport vector_space_sym_report(DoFSpace const &dof_space,
                                               ConfigEnumInput const &input_state,
                                               PermuteIteratorIt group_begin,
                                               PermuteIteratorIt group_end,
                                               bool calc_wedges) {

    if(!is_valid_dof_space(input_state.configuration(), dof_space)) {
      std::stringstream msg;
      msg << "Error in vector_space_sym_report: DoFSpace is not valid for given input state.";
      throw std::runtime_error(msg.str());
    }
    if(!AnisoValTraits(dof_space.dof_key()).global()) {
      if(input_state.sites() != dof_space.sites()) {
        std::stringstream msg;
        msg << "Error in vector_space_sym_report: Mismatch between input state sites and dof_space sites.";
        throw std::runtime_error(msg.str());
      }
    }

    // We need a temporary mastersymgroup to manage the symmetry representation for the DoF
    MasterSymGroup g;
    SymGroupRepID id;
    DoFKey dof_key = dof_space.dof_key();
    SupercellSymInfo const &sym_info = input_state.configuration().supercell().sym_info();
    xtal::BasicStructure const &prim_struc = input_state.configuration().prim().structure();
    std::vector<PermuteIterator> invariant_subgroup = make_invariant_subgroup(input_state,
                                                                              group_begin,
                                                                              group_end);

    // Populate temporary objects for two cases
    // CASE 1: DoF is global (use prim_struc's list of global DoFs, rather than relying on val_traits.is_global()
    if(prim_struc.global_dofs().count(dof_key)) {

      // Global DoF, use point group only
      SymGroup pointgroup = make_point_group(invariant_subgroup, sym_info.supercell_lattice());
      g = make_master_sym_group(pointgroup,
                                sym_info.supercell_lattice());

      id = g.allocate_representation();
      SymGroupRep const &rep = *(sym_info.global_dof_symrep(dof_key).rep_ptr());
      for(Index i = 0; i < pointgroup.size(); ++i) {
        Index fg_ix = pointgroup[i].index();
        g[i].set_rep(id, *rep[fg_ix]);
      }

    }
    // CASE 2: DoF is local
    else {
      auto group_and_ID = collective_dof_symrep(input_state.sites().begin(),
                                                input_state.sites().end(),
                                                sym_info,
                                                dof_key,
                                                invariant_subgroup);
      g = group_and_ID.first;
      g.is_temporary_of(group_and_ID.first);
      id = group_and_ID.second;
    }

    SymGroupRep const &rep = g.representation(id);

    // Generate report, based on constructed inputs
    VectorSpaceSymReport result = vector_space_sym_report(rep,
                                                          g,
                                                          dof_space.basis(),
                                                          calc_wedges);
    result.axis_glossary = dof_space.axis_glossary();

    return result;
  }

}

#endif
