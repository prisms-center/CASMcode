#ifndef CASM_DoFSpace_impl
#define CASM_DoFSpace_impl

#include "casm/enumerator/DoFSpace.hh"

namespace CASM {

  template<typename PermuteIteratorIt>
  VectorSpaceSymReport vector_space_sym_report(DoFSpace const &dof_space,
                                               PermuteIteratorIt group_begin,
                                               PermuteIteratorIt group_end,
                                               bool calc_wedges) {

    // We need a temporary mastersymgroup to manage the symmetry representation for the DoF
    MasterSymGroup g;
    SymGroupRepID id;
    DoFKey dof_key = dof_space.dof_key;

    ConfigEnumInput const &config_region = dof_space.config_region;
    SupercellSymInfo const &sym_info = config_region.configuration().supercell().sym_info();
    xtal::BasicStructure const &prim_struc = config_region.configuration().prim().structure();

    std::vector<PermuteIterator> invariant_subgroup = make_invariant_subgroup(config_region,
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
      auto group_and_ID = collective_dof_symrep(config_region.sites().begin(),
                                                config_region.sites().end(),
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
                                                          dof_space.dof_subspace,
                                                          calc_wedges);
    result.axis_glossary = make_axis_glossary(
                             dof_key,
                             config_region.configuration(),
                             config_region.sites());

    return result;
  }

}

#endif
