#ifndef CASM_app_clex_file_io
#define CASM_app_clex_file_io

namespace CASM {

  class Configuration;
  class DirectoryStructure;
  class Supercell;

  /// Write supercell "LAT" file (supercell lattice vectors as column vectors) to standard location
  void write_lat(Supercell const &supercell, DirectoryStructure const &dir);

  /// Write configuration "POS" file (VASP POSCAR) to standard location
  void write_pos(Configuration const &configuration, DirectoryStructure const &dir);

  /// Write configuration "structure.json" file (structure that results from appling DoF) to standard location
  void write_structure_json(Configuration const &configuration, DirectoryStructure const &dir);

  /// Write configuration "config.json" file (DoF values in standard basis) to standard location
  void write_config_json(Configuration const &configuration, DirectoryStructure const &dir);

}

#endif
