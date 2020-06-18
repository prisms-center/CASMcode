#include "casm/enumerator/DoFSpace.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SupercellSymInfo_impl.hh"
#include "casm/symmetry/SymRepTools.hh"
namespace CASM {

  //Initializes DoF subspace to identity, if empty matrix is provided
  DoFSpace::DoFSpace(ConfigEnumInput _config_region,
                     DoFKey _dof_key,
                     Eigen::MatrixXd _dof_subspace) :
    config_region(std::move(_config_region)),
    dof_key(_dof_key) {
    Index dofdim = 0;
    if(AnisoValTraits(dof_key).global()) {
      dofdim = config_region.config().configdof().global_dof(dof_key).dim();
    }
    else if(dof_key == "occ") {
      std::vector<int> max_occ = config_region.supercell().max_allowed_occupation();
      for(Index i : config_region.sites()) {
        dofdim += max_occ[i] + 1;
      }
    }
    else {
      dofdim = config_region.config().configdof().local_dof(dof_key).dim() * config_region.sites().size();
    }

    if(dof_subspace.size() == 0)
      dof_subspace.setIdentity(dofdim, dofdim);
    if(dof_subspace.rows() != dofdim) {
      throw std::runtime_error("Attempting to construct DoFSpace with " + std::to_string(dof_subspace.rows())
                               + " but " + std::to_string(dofdim) + " rows are required.");
    }

  }


  VectorSpaceSymReport vector_space_sym_report(DoFSpace const &_space,
                                               bool calc_wedges) {

    // We need a temporary mastersymgroup to manage the symmetry representation for the DoF
    // The particular group used is specified by DoFSpace, and defaults to factor group of configuration/supercell
    MasterSymGroup g;
    SymGroupRepID id;
    DoFKey dof_key = _space.dof_key;
    AnisoValTraits val_traits(dof_key);

    ConfigEnumInput const &config_region = _space.config_region;
    SupercellSymInfo const &sym_info = config_region.supercell().sym_info();
    xtal::BasicStructure const &prim_struc = config_region.config().prim().structure();

    // The axis_glossary gives names un-symmetrized coordinate system
    // It will be populated based on whether the DoF is global or local
    std::vector<std::string> axis_glossary;

    // Populate temporary objects for two cases
    // CASE 1: DoF is global (use prim_struc's list of global DoFs, rather than relying on val_traits.is_global()
    if(prim_struc.global_dofs().count(dof_key)) {

      // Global DoF, use point group only
      SymGroup pointgroup = make_point_group(config_region.group(), sym_info.supercell_lattice());
      g = make_master_sym_group(pointgroup,
                                sym_info.supercell_lattice());

      id = g.allocate_representation();
      SymGroupRep const &rep = *(sym_info.global_dof_symrep(dof_key).rep_ptr());
      for(Index i = 0; i < pointgroup.size(); ++i) {
        Index fg_ix = pointgroup[i].index();
        g[i].set_rep(id, *rep[fg_ix]);
      }

      // Global DoF, axis_glossary comes straight from the DoF
      axis_glossary = component_descriptions(prim_struc.global_dof(dof_key));
    }
    // CASE 2: DoF is local
    else {
      auto group_and_ID = collective_dof_symrep(config_region.sites().begin(),
                                                config_region.sites().end(),
                                                sym_info,
                                                dof_key,
                                                config_region.group());
      g = group_and_ID.first;
      g.is_temporary_of(group_and_ID.first);
      id = group_and_ID.second;

      // Generate full axis_glossary for all active sites of the config_region
      for(Index l : config_region.sites()) {
        Index b = config_region.config().sublat(l);
        if(!prim_struc.basis()[b].dofs().count(dof_key))
          continue;
        std::vector<std::string> tdescs = component_descriptions(prim_struc.basis()[b].dof(dof_key));
        for(std::string const &desc : tdescs) {
          axis_glossary.push_back(desc + "[" + std::to_string(l + 1) + "]");
        }
      }

    }

    SymGroupRep const &rep = g.representation(id);

    // Generate report, based on constructed inputs
    VectorSpaceSymReport result = vector_space_sym_report(rep,
                                                          g,
                                                          _space.dof_subspace,
                                                          calc_wedges);
    result.axis_glossary = std::move(axis_glossary);
    return result;

  }


}
