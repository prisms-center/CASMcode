#ifndef CASMtest_TestStructures
#define CASMtest_TestStructures

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Molecule.hh"

// using namespace CASM;

namespace test {

  inline CASM::xtal::BasicStructure empty_prim() {

    using namespace CASM;
    using namespace CASM::xtal;

    BasicStructure struc { Lattice {Eigen::Matrix3d::Identity()} };
    struc.set_title("empty");
    Molecule A = Molecule::make_atom("A");
    struc.push_back(Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), FRAC), std::vector<Molecule> {A}));

    return struc;
  }

  inline CASM::xtal::BasicStructure ZrO_prim() {

    using namespace CASM;
    using namespace CASM::xtal;

    // lattice vectors as rows
    Eigen::Matrix3d lat;
    lat << 3.233986860000, 0.000000000000, 0.000000000000,
        -1.616993430000, 2.800714770000, 0.000000000000,
        0.000000000000, 0.000000000000, 5.168678340000;

    BasicStructure struc(Lattice(lat.transpose()));
    struc.set_title("ZrO");

    Molecule O = Molecule::make_atom("O");
    Molecule Zr = Molecule::make_atom("Zr");
    Molecule Va = Molecule::make_vacancy();

    struc.push_back(Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), FRAC), std::vector<Molecule> {Zr}));
    struc.push_back(Site(Coordinate(Eigen::Vector3d(2. / 3., 1. / 3., 1. / 2.), struc.lattice(), FRAC), std::vector<Molecule> {Zr}));
    struc.push_back(Site(Coordinate(Eigen::Vector3d(1. / 3., 2. / 3., 1. / 4.), struc.lattice(), FRAC), std::vector<Molecule> {Va, O}));
    struc.push_back(Site(Coordinate(Eigen::Vector3d(1. / 3., 2. / 3., 3. / 4.), struc.lattice(), FRAC), std::vector<Molecule> {Va, O}));

    return struc;
  }

  inline CASM::xtal::BasicStructure FCC_ternary_prim() {

    using namespace CASM;
    using namespace CASM::xtal;

    // lattice vectors as cols
    Eigen::Matrix3d lat;
    lat << 0.0, 2.0, 2.0,
        2.0, 0.0, 2.0,
        2.0, 2.0, 0.0;

    BasicStructure struc {Lattice{lat}};
    struc.set_title("FCC_ternary");

    Molecule A = Molecule::make_atom("A");
    Molecule B = Molecule::make_atom("B");
    Molecule C = Molecule::make_atom("C");

    struc.push_back(Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), CART), std::vector<Molecule> {A, B, C}));

    return struc;

  }

  inline CASM::xtal::BasicStructure SimpleCubicGLstrain_prim() {

    using namespace CASM;
    using namespace CASM::xtal;

    // lattice vectors as cols
    Eigen::Matrix3d lat;
    lat << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;

    BasicStructure struc {Lattice{lat}};
    struc.set_title("SimpleCubicGLstrain");

    Molecule A = Molecule::make_atom("A");

    struc.push_back(Site(Coordinate(Eigen::Vector3d::Zero(), struc.lattice(), CART), std::vector<Molecule> {A}));

    // Add global DoF
    // GLstrain: Green-Lagrange strain
    struc.set_global_dofs({ AnisoValTraits::strain("GL") });

    return struc;

  }

  inline CASM::xtal::BasicStructure FCC_ternary_strain_disp_prim() {

    using namespace CASM;
    using namespace CASM::xtal;

    BasicStructure struc = FCC_ternary_prim();

    // Add global DoF
    // GLstrain: Green-Lagrange strain
    struc.set_global_dofs({ AnisoValTraits::strain("GL") });

    // Update basis to add displacement
    SiteDoFSet disp_dofset {AnisoValTraits::disp()};
    std::vector<Site> new_basis;
    for(auto &site : struc.basis()) {
      new_basis.emplace_back(Coordinate {site}, site.occupant_dof(), std::vector<SiteDoFSet> {disp_dofset});
    }
    struc.set_basis(new_basis);

    return struc;

  }
}

#endif
