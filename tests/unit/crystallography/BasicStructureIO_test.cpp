#include "casm/crystallography/io/BasicStructureIO.hh"

#include "autotools.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "gtest/gtest.h"

using namespace CASM;

// There should be more extensive BasicStructureIO unit tests...

TEST(SelectiveDynamicsTest, ConstructorTest1) {
  using namespace xtal;

  AnisoValTraits selectivedynamics("selectivedynamics");

  xtal::SpeciesProperty selective_dynamics_on(selectivedynamics,
                                              Eigen::Vector3d(1., 1., 1.));
  xtal::SpeciesProperty selective_dynamics_off(selectivedynamics,
                                               Eigen::Vector3d(0., 0., 0.));
  // std::map<std::string, xtal::SpeciesProperty> prop_map;

  Molecule Zr = Molecule::make_atom("Zr");
  Zr.set_properties({{selective_dynamics_off.name(), selective_dynamics_off}});
  SiteDoFSet disp_xyz{AnisoValTraits::disp()};

  Molecule H = Molecule::make_atom("H");
  H.set_properties({{selective_dynamics_on.name(), selective_dynamics_on}});

  // SiteDoFSet disp_xyz {
  //     AnisoValTraits::disp(),  // AnisoVal type
  //     {"dx", "dy", "dz"},      // axes names
  //     Eigen::Matrix3d::Identity(),
  //     {}  // excluded_occs
  // };

  Lattice lat{Eigen::Vector3d{0., 2.4106965, 2.4106965},
              Eigen::Vector3d{2.4106965, 0., 2.4106965},
              Eigen::Vector3d{2.4106965, 2.4106965, 0.}};

  BasicStructure struc{lat};
  struc.set_basis({Site{Coordinate{0.0, 0.0, 0.0, lat, FRAC}, {Zr}, {disp_xyz}},
                   Site{Coordinate{0.25, 0.25, 0.25, lat, FRAC}, {H}},
                   Site{Coordinate{0.75, 0.75, 0.75, lat, FRAC}, {H}}});
  struc.set_global_dofs({AnisoValTraits::strain("GL")});

  EXPECT_EQ(struc.basis().size(), 3);

  jsonParser json;
  to_json(struc, json, FRAC);
  std::cout << json << std::endl;
}

TEST(SelectiveDynamicsTest, JsonIOTest1) {
  using namespace xtal;

  jsonParser prim_json = CASM::jsonParser::parse(std::string(R"({
    "basis" : [
      {
        "coordinate" : [ 0.0000000, 0.0000000, 0.0000000 ],
        "occupant_dof" : [ "Zr" ],
        "dofs": {
          "disp": {}
        }
      },
      {
        "coordinate" : [ 0.2500000, 0.2500000, 0.2500000 ],
        "occupant_dof" : [ "H" ]
      },
      {
        "coordinate" : [ 0.7500000, 0.7500000, 0.7500000 ],
        "occupant_dof" : [ "H" ]
      }
    ],
    "species" : {
      "H": {
        "properties": {
          "selectivedynamics": {
            "value": [1, 1, 1]
          }
        }
      },
      "Zr": {
        "properties": {
          "selectivedynamics": {
            "value": [0, 0, 0]
          }
        }
      }
    },
    "dofs" : {
        "GLstrain" : {}
    },
    "coordinate_mode" : "Fractional",
    "description" : "Cubic ZrH_{2}",
    "lattice_vectors" : [
      [0.       , 2.4106965, 2.4106965],
      [2.4106965, 0.       , 2.4106965],
      [2.4106965, 2.4106965, 0.       ]
    ],
    "title" : "ZrH2"
  })"));

  BasicStructure prim = read_prim(prim_json, TOL);
  EXPECT_EQ(prim.basis().size(), 3);

  COORD_TYPE coordtype = FRAC;
  Eigen::Matrix3d inv_lat_column_mat = Eigen::Matrix3d::Identity();
  // change this to use FormatFlag
  if (coordtype == FRAC) {
    inv_lat_column_mat = prim.lattice().inv_lat_column_mat();
  }

  // jsonParser tjson;
  // for (auto const &site : prim.basis()) {
  //   for (auto const &mol : site.occupant_dof()) {
  //     for (auto const &atom_position : mol.atoms()) {
  //       for (auto const &prop_name_value_pair : atom_position.properties()) {
  //         std::cout << "atom_position prop_name_value_pair" << std::endl;
  //         to_json(prop_name_value_pair.second, tjson);
  //         std::cout << tjson << std::endl;
  //       }
  //       std::cout << "atom_position" << std::endl;
  //       to_json(atom_position, tjson, inv_lat_column_mat);
  //       std::cout << tjson << std::endl;
  //     }
  //     for (auto const &prop_name_value_pair : mol.properties()) {
  //       std::cout << "mol prop_name_value_pair" << std::endl;
  //       to_json(prop_name_value_pair.second, tjson);
  //       std::cout << tjson << std::endl;
  //     }
  //     std::cout << "molecule" << std::endl;
  //     to_json(mol, tjson, inv_lat_column_mat);
  //     std::cout << tjson << std::endl;
  //   }
  //   std::cout << "site" << std::endl;
  //   to_json(site, tjson, coordtype);
  //   std::cout << tjson << std::endl;
  // }
  // std::cout << "basic_structure" << std::endl;
  jsonParser json_again;
  to_json(prim, json_again, coordtype);
  // std::cout << json_again << std::endl;
  EXPECT_EQ(json_again["basis"].size(), 3);
}
