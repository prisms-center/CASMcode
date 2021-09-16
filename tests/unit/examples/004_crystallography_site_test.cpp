#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/crystallography/Site.hh"
#include "gtest/gtest.h"

TEST(ExampleCrystallographySite, SiteConstructor) {
  using CASM::AnisoValTraits;
  using CASM::xtal::Coordinate;
  using CASM::xtal::Lattice;
  using CASM::xtal::Molecule;
  using CASM::xtal::Site;
  using CASM::xtal::SiteDoFSet;
  using Eigen::Vector3d;

  // First, construct a Lattice
  Vector3d a{1.0, 0.0, 0.0};
  Vector3d b{0.0, 2.0, 0.0};
  Vector3d c{0.0, 0.0, 3.0};
  Lattice lattice{a, b, c};

  // The Site is a Coordinate (via inheritance) and it also contains
  //   information on the degrees of freedom (DoF) allowed on a particular
  //   sublattice of the infinte crystal.
  //
  // The allowed DoF may be:
  // - discrete "occupants", represented using Molecule
  // - continuous valued properties (displacement, magspin, etc.), represented
  // using SiteDoFSet.
  //
  // Molecules are composed of atoms, which may be isotropic, or may have
  // properties which are
  //    anisotropic.
  //
  //      Molecules have:
  //      - a name string
  //      - a vector of atoms, represented as AtomPosition
  //      - a map of <attribute name string>:<SpeciesProperty>
  //      - can be made 'divisible' or 'indivisible', for kinetics purposes
  //
  //      AtomPosition have:
  //      - a name string
  //      - a coordinate
  //      - a map of <attribute name string>:<SpeciesProperty>
  //
  //      SpeciesProperty have:
  //      - AnisoValTraits, which provides the attribute type name, a standard
  //      coordinate system
  //        (the "standard basis"), and specifies how values transform under
  //        application of symmetry.
  //      - a Eigen::VectorXd value, with meaning determined by the type of
  //      SpeciesProperty
  //
  //   Molecules can be compared using the `Molecule::is_identical` function,
  //   which requires that
  //     all components and subcomponents of the two molecules are identical up
  //     to a specified tolerance.
  //
  //   The function `Molecule::make_atom` constructs a Molecule with a single
  //   isotropic atom. The molecule's
  //     name and the atom's name are equivalent.
  //
  //   The function `Molecule::make_vacancy` constructs an atom-like Molecule
  //   with the special name "Va" to
  //     indicate that it represents a vacancy.
  //
  // SiteDoFSet represent continuous degrees of freedom (DoF) which are "local"
  // to a particular
  //   sublattice (i.e. displacement).
  //
  //      SiteDoFSet have:
  //      - AnisoValTraits, which provides the DoF type name, a standard
  //      coordinate system (the
  //        "standard basis"), and specifies how values transform under
  //        application of symmetry.
  //      - a "DoF basis", a set of named basis vectors which are denoted
  //      relative to the standard
  //        basis, allowing the user to specify the DoFSet components, name
  //        them, and restrict DoF values to a particular subspace
  //      - a list of site occupants for which the DoF does not apply
  //      ("excluded_occupants",
  //        std::unordered_set<std::string>)
  //
  // Examples of standard basis specified by AnisoValTraits:
  // - "disp" -> (dx, dy, dz) -> displacement components relative to fixed
  // laboratory frame
  // - "strain" -> (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy)
  // -> tensor elements

  auto make_atom = &Molecule::make_atom;  // function pointer to static member

  // Construct a Site with coordinate and occupants, no continuous DoF:
  // - binary occupation: "A" and "B" atoms
  Site site_1{Coordinate{Vector3d{0., 0., 0.}, lattice, CASM::FRAC},
              std::vector<Molecule>{make_atom("A"), make_atom("B")}};
  EXPECT_EQ(site_1.occupant_dof().size(), 2);
  EXPECT_EQ(site_1.dofs().size(), 0);

  // Construct displacement SiteDoFSet
  SiteDoFSet disp_dofset{AnisoValTraits::disp()};

  // Construct a Site with coordinate, fixed occupation, and displacement
  // (explicit types)
  Site site_2{Coordinate{Eigen::Vector3d{0., 0., 0.}, lattice, CASM::FRAC},
              std::vector<Molecule>{make_atom("A")},
              std::map<std::string, SiteDoFSet>{
                  {disp_dofset.type_name(), disp_dofset}}};
  EXPECT_EQ(site_2.occupant_dof().size(), 1);
  EXPECT_EQ(site_2.dofs().size(), 1);
}
