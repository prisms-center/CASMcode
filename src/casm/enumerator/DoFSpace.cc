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

    Index dof_space_dimension = get_dof_space_dimension(
                                  _dof_key,
                                  config_region.configuration(),
                                  config_region.sites());

    if(dof_subspace.size() == 0)
      dof_subspace.setIdentity(dof_space_dimension, dof_space_dimension);
    if(dof_subspace.rows() != dof_space_dimension) {
      throw std::runtime_error("Attempting to construct DoFSpace with " + std::to_string(dof_subspace.rows())
                               + " but " + std::to_string(dof_space_dimension) + " rows are required.");
    }

  }

  Index get_dof_space_dimension(DoFKey dof_key,
                                Configuration const &configuration,
                                std::set<Index> const &sites) {

    Index dof_space_dimension = 0;
    auto const &supercell = configuration.supercell();
    auto const &configdof = configuration.configdof();

    if(AnisoValTraits(dof_key).global()) {
      dof_space_dimension = configdof.global_dof(dof_key).dim();
    }
    else if(dof_key == "occ") {
      std::vector<int> max_occ = supercell.max_allowed_occupation();
      for(Index i : sites) {
        dof_space_dimension += max_occ[i] + 1;
      }
    }
    else {
      dof_space_dimension = configdof.local_dof(dof_key).dim() * sites.size();
    }
    return dof_space_dimension;
  }

  /// Return a vector<std::string> with DoF component descriptions
  ///
  /// Note:
  /// - Site DoF components are written with site index using a "beginning from 1" convention
  ///   - For example: if single DoF component per site with name "dx", and "sites" is {0, 2}, the
  ///     output is: {"dx[1]", "dx[3]"}
  std::vector<std::string> make_axis_glossary(DoFKey dof_key,
                                              Configuration const &configuration,
                                              std::set<Index> const &sites) {

    xtal::BasicStructure const &prim_struc = configuration.prim().structure();

    // The axis_glossary gives names un-symmetrized coordinate system
    // It will be populated based on whether the DoF is global or local
    std::vector<std::string> axis_glossary;

    if(prim_struc.global_dofs().count(dof_key)) {
      // Global DoF, axis_glossary comes straight from the DoF
      axis_glossary = component_descriptions(prim_struc.global_dof(dof_key));
    }
    else {
      // Generate full axis_glossary for all active sites of the config_region
      for(Index site_index : sites) {
        Index sublattice_index = configuration.sublat(site_index);
        xtal::Site const &site = prim_struc.basis()[sublattice_index];
        if(!site.dofs().count(dof_key))
          continue;
        std::vector<std::string> tdescs = component_descriptions(site.dof(dof_key));
        for(std::string const &desc : tdescs) {
          axis_glossary.push_back(desc + "[" + std::to_string(site_index + 1) + "]");
        }
      }
    }
    return axis_glossary;
  }

}
