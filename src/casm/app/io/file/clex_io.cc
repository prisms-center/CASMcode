#ifndef CASM_app_clex_file_io
#define CASM_app_clex_file_io

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/app/DirectoryStructure.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/clex/io/stream/Configuration_stream_io.hh"

namespace CASM {

  /// Write supercell "LAT" file (supercell lattice vectors as column vectors) to standard location
  void write_lat(Supercell const &supercell, DirectoryStructure const &dir) {
    fs::create_directories(dir.configuration_dir(supercell.name()));
    fs::ofstream outfile {dir.LAT(supercell.name())};
    outfile << supercell.lattice().lat_column_mat() << std::endl;
  }

  /// Write configuration "POS" file (VASP POSCAR) to standard location
  void write_pos(Configuration const &configuration, DirectoryStructure const &dir) {
    fs::create_directories(dir.configuration_dir(configuration.name()));
    fs::ofstream outfile {dir.POS(configuration.name())};
    print_poscar(configuration, outfile);
  }

  /// Write configuration "structure.json" file (structure that results from appling DoF) to standard location
  void write_structure_json(Configuration const &configuration, DirectoryStructure const &dir) {
    fs::create_directories(dir.configuration_dir(configuration.name()));
    fs::ofstream outfile {dir.structure_json(configuration.name())};
    jsonParser json;
    to_json(make_simple_structure(configuration), json);
    json.print(outfile);
  }

  /// Write configuration "config.json" file (DoF values in standard basis) to standard location
  void write_config_json(Configuration const &configuration, DirectoryStructure const &dir) {
    fs::create_directories(dir.configuration_dir(configuration.name()));
    fs::ofstream outfile {dir.config_json(configuration.name())};
    jsonParser json;
    to_json(configuration, json);
    json.print(outfile);
  }

}

#endif
