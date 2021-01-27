#include "casm/app/ProjectBuilder.hh"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "casm/app/AppIO.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {

void build_project(ProjectSettings const &project_settings,
                   Structure const &prim) {
  ProjectSettings const &set = project_settings;

  if (!set.has_dir()) {
    throw std::runtime_error(std::string("Error in build_project:") +
                             " Project settings do not have a root dir.");
  }
  auto checkroot = find_casmroot(set.root_dir());
  if (checkroot == set.root_dir()) {
    throw std::runtime_error(std::string("Error in build_project:") +
                             " Project already exists at '" +
                             checkroot.string() + "'");
  }
  DirectoryStructure const &dir = set.dir();

  // Create project directories -------------
  log().indent() << "Creating CASM project directory tree at: "
                 << dir.root_dir() << std::endl;
  create_all_directories(set);

  // Write project_settings.json file -------
  commit(set);

  // Write prim.json file (into .casm directory) -------------------
  log().indent() << "Writing prim file: " << dir.prim() << std::endl;
  write_prim(prim, dir.prim(), FRAC);

  // Calculate symmetry  --------------------
  SymGroup lattice_point_grp{SymGroup::lattice_point_group(prim.lattice())};

  // Write symmetry info files

  // Write lattice point group
  {
    fs::ofstream outfile;
    jsonParser json;
    outfile.open(dir.lattice_point_group());
    write_symgroup(lattice_point_grp, json);
    json.print(outfile);
    outfile.close();
  }

  // Write factor group
  {
    fs::ofstream outfile;
    jsonParser json;
    outfile.open(dir.factor_group());
    write_symgroup(prim.factor_group(), json);
    json.print(outfile);
    outfile.close();
  }

  // Write crystal point group
  {
    fs::ofstream outfile;
    jsonParser json;
    outfile.open(dir.crystal_point_group());
    write_symgroup(prim.point_group(), json);
    json.print(outfile);
    outfile.close();
  }

  // Generate empty composition_axes.json --------------------
  CompositionAxes().write(dir.composition_axes());
}

ProjectSettings make_default_project_settings(xtal::BasicStructure const &prim,
                                              std::string project_name) {
  ProjectSettings settings{project_name};
  auto W = default_nlist_weight_matrix(prim, settings.crystallography_tol());
  settings.set_nlist_weight_matrix(W);
  settings.set_nlist_sublat_indices(default_nlist_sublat_indices(prim));
  settings.set_default_clex(default_configuration_clex());
  settings.set_required_properties("Configuration", "default",
                                   {"relaxed_energy"});
  return settings;
}

ProjectSettings make_default_project_settings(xtal::BasicStructure const &prim,
                                              std::string project_name,
                                              fs::path root_dir) {
  ProjectSettings settings = make_default_project_settings(prim, project_name);
  settings.set_root_dir(root_dir);
  return settings;
}

Eigen::Matrix3l default_nlist_weight_matrix(xtal::BasicStructure const &prim,
                                            double tol) {
  return PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(),
                                              10, tol);
}

std::set<int> default_nlist_sublat_indices(xtal::BasicStructure const &prim) {
  std::set<int> sublat_indices;
  for (int b = 0; b < prim.basis().size(); ++b) {
    if (prim.basis()[b].occupant_dof().size() >= 2 ||
        prim.basis()[b].dof_size() > 0) {
      sublat_indices.insert(b);
    }
  }
  return sublat_indices;
}
}  // namespace CASM
