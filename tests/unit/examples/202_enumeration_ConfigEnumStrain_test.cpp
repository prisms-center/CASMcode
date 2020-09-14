#include "gtest/gtest.h"
#include "autotools.hh"
#include "Common.hh"

#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/SymRepTools.hh"
#include "casm/symmetry/io/json/SymRepTools.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// This introduces how to use CASM to analyze and sample the space of global continuous degrees of
// freedom (DoF), using strain as an example.

// A brief review of how CASM represents an infinite crystal and the space of allowed perturbations
// (the DoF space):
//
//   Each type of continuous DoF has:
//   - AnisoValTraits, which provides the DoF type name, a standard coordinate system (the "standard
//     basis"), and specifies how values transform under application of symmetry. Examples of the
//     "standard basis" specified by AnisoValTraits include:
//     - "disp" -> (dx, dy, dz) -> displacement components relative to fixed laboratory frame
//     - "strain" -> (e_xx, e_yy, e_zz, sqrt(2)*e_yz, sqrt(2)*e_xz, sqrt(2)*e_xy) -> tensor elements
//
//   A BasicStructure represents an infinite crystal and its allowed perturbations by specifying one
//   unit cell. BasicStructure has:
//   - a Lattice
//   - global continuous DoF (std::map<DoFKey, xtal::DoFSet>)
//   - a basis (std::vector<Site>). Each Site may have:
//     - local discrete DoF ("occupant_dof", std::vector<Molecule>)
//     - local continuous DoF ("dof", std::map<DoFKey, xtal::SiteDoFSet>)
//
//     Note: DoFKey is a typedef for std::string, matching AnisoValTraits::name() for a particular
//           type of DoF.
//
//   Each xtal::DoFSet or xtal::SiteDoFSet has:
//   - a "DoF basis" specifying a basis which is a linear combination of the standard basis defined
//     by AnisoValTraits.
//
//     Note: The DoF basis may be a subspace of the standard basis. This allows specifying that for
//           a particular crystal only a subset of strain perturbations are allowed, or for
//           particular sites only particular displacements are allowed. These subspaces are taken
//           into account when determining the symmetry group of the BasicStructure.
//
//   A Structure is an object that holds a (shared, const) BasicStructure and its symmetry
//   information. A Structure has:
//   - a BasicStructure (std::shared_ptr<const xtal::BasicStructure>)
//   - the basic structure's factor group (MasterSymGroup)
//   - keys (SymGroupRepID) used to access symmetry group representations (SymGroupRep)
//
//     Note: In CASM, a Structure is often referred to as the "primitive crystal structure" or "prim",
//           especially in the context where it is generally expected to the primitive unit cell.
//
//
// The DoF Space
// -------------
//
// The total crystal degree of freedom space corresponding to a particular prim structure is the
// combination of the global continuous DoFs, local discrete DoF, and local continuous DoFs.
//
// - Essential symmetry information (irreducible representations) can be obtained for various
//   subspaces of the total crystal degree of freedom space and subgroups of the prim structure's
//   factor group.
// - The degree of freedom subspaces are primarily specified by a combination of DoF type (e.g.,
//   "disp", "GLstrain", "occ") and, for local DoF, which sites to include.
// - The symmetry of a degree of freedom subspace can not be higher than that of the prim factor
//   group, but it may have lower symmetry. It may have lower translational symmetry corresponding
//   to a particular supercell lattice periodicity. It may have lower symmetry corresponding to
//   order/disorder of local degrees of freedom within a supercell. Thus, this symmetry constraint
//   may be specified by using the factor group of a particular Configuration.
//
// The DoFSpace class is used to specify a particular degree of freedom space and the symmetry
// constraints on it. DoFSpace have:
//   - dof_key (DoFKey): A string indicating which DoF type (e.g., "disp", "Hstrain", "occ")
//   - config_region (ConfigEnumInput): Specifies both the symmetry group for the DoFSpace (from
//     "config_region.group()") and, for local DoF, which DoF to include (DoF on sites listed
//     "config_region.sites()"). The symmetry group "config_region.group()" is the factor group
//     of the ConfigEnumInput's configuration, excluding any operation that causes permutation
//     between those sites included in "config_region.sites()" and those sites excluded from
//     "config_region.sites()".
//   - dof_subspace (Eigen::MatrixXd): Allows specifying a subspace of the space determined from
//     config_region and dof_key. Examples:
//
//       For dof_key=="disp", and config_region.sites().size()==4:
//         The default DoF space has dimension 12 corresponding to (x,y,z) for each of the 4 sites.
//         Then, dof_subspace is a matrix with 12 rows and 12 or fewer columns.
//
//       For dof_key=="occ", and config_region.sites().size()==4, with the particular sites selected
//       having allowed occupations ["A", "B"], ["A", "B", "C"], ["A, B"], and ["A", "B", "C"] :
//         The default DoF space has dimension 10 = 2 + 3 + 2 + 3.
//         Then, dof_subspace is a matrix with 10 rows and 10 or fewer columns.
//
//       For dof_key=="GLstrain":
//         The default DoF space has dimension 6 for the six independent GLstrain components.
//         Then, dof_subspace is a matrix with 6 rows and 6 or fewer columns.
//
//     By default, dof_subspace is the full DoF space (identity matrix with dimensions matching
//     the full space for the particular combination of config_region and dof_key).
//
//
// Sampling the DoF Space
// ----------------------
//
// For sampling strains within the full GLstrain space, one way to proceed is to find symmetrically
// unique "wedges", symmetrically unique portions of space that fall between high symmetry
// directions. Note that wedges represent an infinite volume of space. For example, in strain space
// a wedge begins at the origin, representing zero strain, and comprises the region of space between
// high symmetry strain directions at any magnitude.
//
// There are multiple types of wedges:
// - The SymRepTools::IrrepWedge class represents a wedge in an irreducible vector space (a subspace
//   of the full DoF space that is unchanged under transformation by the relevant symmetry group).
// - The SymRepTools::SubWedge class represents a vector of SymRepTools::IrrepWedge, one from each
//   irreducible subspace of the full DoF space. Together, the IrrepWedge comprising a Subwedge span
//   the full DoF space. While its axes do span the full space, the orbit of a SubWedge may not fill
//   the full DoF space.
// - The VectorSpaceSymReport::irreducible_wedge is a vector of SubWedge whose orbits combined
//   together fill the full DoF space.
//
// Sampling of degree of freedom values in the full DoF space can be performed by sampling, in turn,
// values in each SubWedge comprising the irreducible wedge for the full space. The standard
// approach in CASM to enumerate degree of freedom values in a SubWedge is to construct a grid of
// points aligned with the subwedge axes. Especially for subwedges whose axes form acute angles,
// this grid may include many unwanted points with large magnitudes. So CASM also provides a
// "trim_corners" option which excludes any points lying outside the ellipsoid that passes through
// the maximum enumerated value lying on each wedge axis.
//
// The following examples show how to specify a DoF space, create a VectorSpaceSymReport to analyze
// it, and construct and use wedges to sample the DoF space.

namespace enumeration_test_impl {

  // Helper functions to make the examples easier to read

  SymGroupRep const &global_dof_symrep(Structure const &prim, DoFKey global_dof_key);

  SymGroupRep const &global_dof_symrep(ConfigEnumInput const &config_input, DoFKey global_dof_key);

  SymGroup make_point_group(ConfigEnumInput const &config_input);

  // check if there is a solution X, for:  B = A * X
  bool is_equivalent_column_space(Eigen::MatrixXcd const &A, Eigen::MatrixXcd const &B, double tol);
}

// This test fixture class constructs a CASM PrimClex and Supercells for enumeration examples
class ExampleEnumerationSimpleCubicConfigEnumStrain : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  // values to check
  std::vector<Eigen::MatrixXcd> expected_irrep_subspace;
  double tol;

  ExampleEnumerationSimpleCubicConfigEnumStrain();

  // Check that the IrrepInfo trans_mat value matches expectations:
  void check_irrep(SymRepTools::IrrepInfo const &irrep);

};

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, VectorSpaceSymReport) {

  // This example demonstrates specifying a DoFSpace and constructing its VectorSpaceSymReport

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // The VectorSpaceSymReport provides symmetry information about a particular crystal "degree of
  // freedom space" under particular symmetry constraints. VectorSpaceSymReport has:
  // - TODO: more description ...

  // Construct the SimpleCubic GLstrain DoF space.
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space {config_input, dof_key};

  // Construct the VectorSpaceSymReport for the SimpleCubic GLstrain space.
  bool calc_wedges = true;  // explanation TODO
  VectorSpaceSymReport sym_report = vector_space_sym_report(dof_space, calc_wedges);

  // Uncomment to print dof_space:
  // jsonParser dof_space_json;
  // to_json(dof_space, dof_space_json);
  // std::cout << "DoFSpace:\n" << dof_space_json << std::endl;

  // Uncomment to print sym_report:
  // jsonParser sym_report_json;
  // to_json(sym_report, sym_report_json);
  // std::cout << "VectorSpaceSymReport:\n" << sym_report_json << std::endl;

  // Expect three irreducible representations
  EXPECT_EQ(sym_report.irreps.size(), 3);

  // Check the calculated irreducible representations
  for(SymRepTools::IrrepInfo const &irrep : sym_report.irreps) {
    this->check_irrep(irrep);
  }

}

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, FullSpaceIrreps) {

  // Example constructing the full GLstrain space's irreducible representations directly

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // Use full prim structure symmetry, specified by using the volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";

  std::vector<CASM::Configuration> configurations;
  SymGroup point_group = make_point_group(config_input);
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);

  // Construct irreducible representations
  bool allow_complex = true;
  std::vector<SymRepTools::IrrepInfo> irreps = irrep_decomposition(
                                                 strain_symrep, point_group, allow_complex);

  // Expect three irreducible representations
  EXPECT_EQ(irreps.size(), 3);

  // Check the calculated irreducible representations
  for(SymRepTools::IrrepInfo const &irrep : irreps) {
    this->check_irrep(irrep);
  }
}


TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, MakeFullSpaceIrrepWedges) {

  // This example demonstrates constructing the unique SymRepTools::IrrepWedge for the full GLstrain
  //   space

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // Use full prim structure symmetry, specified by using the volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";

  SymGroup point_group = make_point_group(config_input);
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);

  // Make the unique IrrepWedge in full GLstrain space
  int dim = strain_symrep.dim();
  Eigen::MatrixXd subspace = Eigen::MatrixXd::Identity(dim, dim);
  std::vector<SymRepTools::IrrepWedge> irrep_wedges = SymRepTools::irrep_wedges(
                                                        strain_symrep, point_group, subspace);

  // TODO: checks that demonstrate what the IrrepWedge are and that they are what they should be
  EXPECT_EQ(irrep_wedges.size(), 3);

}


TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, MakeFullSpaceSubWedges) {

  // This example demonstrates constructing the unique SymRepTools::SubWedge for the full GLstrain
  //   space

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // Use full prim structure symmetry, specified by using the volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";

  SymGroup point_group = make_point_group(config_input);
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);

  // Make the unique SubWedge in full GLstrain space
  std::vector<SymRepTools::SubWedge> subwedges = SymRepTools::symrep_subwedges(strain_symrep, point_group);

  // TODO: checks that demonstrate what the SubWedge are and that they are what they should be
  EXPECT_EQ(subwedges.size(), 6);

}

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, FullSpaceWedgesEnum) {

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // This example demonstrates constructing SymRepTools::SubWedge for the full GLstrain space and
  //   using them to enumerate strain configurations

  // Use full prim structure symmetry, specified by using the volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";

  SymGroup point_group = make_point_group(config_input);
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);

  // Make wedges in GLstrain space
  std::vector<SymRepTools::SubWedge> subwedges = SymRepTools::symrep_subwedges(strain_symrep, point_group);

  // ** Enumerate degree of freedom values by sampling in each SubWedge **

  // SubWedge dimension == full GLstrain space dimension (6)
  int dim = strain_symrep.dim();

  // Specify a grid with three points along each SubWedge axis, with minimum value of 0.0
  //   and increment value of 0.5.
  double _min = 0.0;
  double _incr = 0.5;
  int n_increment = 3;
  Eigen::VectorXd min_value = Eigen::VectorXd::Constant(dim, _min);
  Eigen::VectorXd max_value = Eigen::VectorXd::Constant(dim, _min + (n_increment - 1) * _incr + _incr / 10.);
  Eigen::VectorXd increment_value = Eigen::VectorXd::Constant(dim, _incr);

  // Construct the CASM::ConfigEnumStrain enumerator
  bool auto_range = false;  // explanation TODO
  bool trim_corners = false; // explanation TODO
  CASM::ConfigEnumStrain enumerator {config_input, subwedges, min_value, max_value, increment_value,
                                     dof_key, auto_range, trim_corners};

  std::vector<CASM::Configuration> configurations;
  std::copy(enumerator.begin(), enumerator.end(), std::back_inserter(configurations));

  // Check the number of configurations:
  // - Each subwedge has a grid with n_increment points along each dimension, so the total number
  //   of enumerated configurations is pow(n_increment, dim) * subwedges.size()
  EXPECT_EQ(configurations.size(), pow(n_increment, dim)*subwedges.size());

}

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, SubSpaceIrreps) {

  // Example demonstrating generating irreducible representations for the irreducible subspaces
  //   containing particular directions

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // Use full prim structure symmetry, specified by using the volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";

  SymGroup point_group = make_point_group(config_input);
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);

  // Example 1:
  // A subspace that includes a single direction will always have a single irrep
  for(int i = 0; i < expected_irrep_subspace.size(); ++i) {

    // Construct subspace consisting of a single direction in an irreducible subspace
    Eigen::MatrixXd subspace {6, 1};
    subspace.col(0) = expected_irrep_subspace[i].col(0).real();

    // Construct irreducible representations spanning the subspace.
    bool allow_complex = true;
    std::vector<SymRepTools::IrrepInfo> irreps = irrep_decomposition(
                                                   strain_symrep, point_group, subspace, allow_complex);

    EXPECT_EQ(irreps.size(), 1);
    this->check_irrep(irreps[0]);
  }

  // Example 2:
  // A subspace that spans two irreducible subspaces will produce two irreps
  for(int i = 0; i < expected_irrep_subspace.size(); ++i) {
    for(int j = 0; j < i; ++j) {

      // Construct subspace consisting of one direction from each of two irreducible subspaces
      Eigen::MatrixXd subspace {6, 2};
      subspace.col(0) = expected_irrep_subspace[i].col(0).real();
      subspace.col(1) = expected_irrep_subspace[j].col(0).real();

      // Construct irreducible representations spanning the subspace.
      bool allow_complex = true;
      std::vector<SymRepTools::IrrepInfo> irreps = irrep_decomposition(
                                                     strain_symrep, point_group, subspace, allow_complex);

      EXPECT_EQ(irreps.size(), 2);

      // TODO: check it contains the expected two irreps
      for(auto const &irrep : irreps) {
        this->check_irrep(irrep);
      }
    }
  }

}

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, SubSpaceWedgesEnum) {

  // Example enumerating GLstrain in each of the known irreducible subspaces

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  // Use full prim structure symmetry, specified by using the volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";
  SymGroup point_group = make_point_group(config_input);
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);
  int dim = strain_symrep.dim();

  // Loop over the known irreducible subspaces
  for(int i = 0; i < expected_irrep_subspace.size(); ++i) {
    int irrep_dim = expected_irrep_subspace[i].cols();

    // Construct subspace consisting of a single direction from the irreducible subspace
    Eigen::MatrixXd subspace {dim, 1};
    subspace.col(0) = expected_irrep_subspace[i].col(0).real();

    // Make SubWedges in the irreducible subspace
    std::vector<SymRepTools::SubWedge> subwedges = SymRepTools::symrep_subwedges(
                                                     strain_symrep, point_group, subspace);

    // Check single SubWedge generated from single IrrepWedge
    EXPECT_EQ(subwedges.size(), 1);
    EXPECT_EQ(subwedges[0].irrep_wedges().size(), 1);
    EXPECT_EQ(subwedges[0].trans_mat().rows(), dim);
    EXPECT_EQ(subwedges[0].trans_mat().cols(), irrep_dim);

    // Enumerate by sampling a portion of the subwedges:
    // - Using min=0.2, max=1.21, increment=0.5 results in a sampling grid with three points along
    //   each irreducible subspace dimension
    Eigen::VectorXd min_value = Eigen::VectorXd::Constant(irrep_dim, 0.2);
    Eigen::VectorXd max_value = Eigen::VectorXd::Constant(irrep_dim, 1.21);
    Eigen::VectorXd increment_value = Eigen::VectorXd::Constant(irrep_dim, 0.5);

    // Construct the CASM::ConfigEnumStrain enumerator
    bool auto_range = false;  // explanation TODO
    bool trim_corners = false; // explanation TODO
    CASM::ConfigEnumStrain enumerator {config_input, subwedges, min_value, max_value, increment_value,
                                       dof_key, auto_range, trim_corners};


    // Enumerate strained configurations
    std::vector<CASM::Configuration> configurations;
    std::copy(enumerator.begin(), enumerator.end(), std::back_inserter(configurations));

    // Check number of configurations
    EXPECT_EQ(configurations.size(), pow(3, irrep_dim)*subwedges.size());

    // Check first value matches initial grid point
    Configuration initial_config = configurations[0];
    auto const &GLstrain = initial_config.configdof().global_dof(dof_key).values();
    auto expected_GLstrain = subwedges[0].trans_mat() * min_value;
    EXPECT_EQ(almost_equal(GLstrain, expected_GLstrain, tol), true);

  }
}

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, FullSpaceDirectGridStrainEnum) {

  // It is also possible to use ConfigEnumStrain to directly sample GLstrain on any user-specified
  //   grid. This is done by creating a "dummy subwedge" which specifies the axes of the sampling
  //   grid irrespective of any symmetry considerations. This example demonstrates how to do this.

  using namespace enumeration_test_impl; // for global_dof_symrep

  // Enumerate GLstrain applied to the default volume 1 supercell
  Supercell const &supercell = *primclex.db<Supercell>().begin();
  CASM::ConfigEnumInput config_input {supercell};
  DoFKey dof_key = "GLstrain";

  // To sample strains in the entire GLstrain space, we want the "dummy subwedge" to span the entire
  //   6-dimensional GLstrain space. For any global DoF we can get the dimensionality of its space
  //   from the SymGroupRep.
  SymGroupRep const &strain_symrep = global_dof_symrep(config_input, dof_key);
  int dim = strain_symrep.dim();

  // To sample GLstrain on a grid spanning the entire GLstrain space, the sampled_space_rank is the
  //   full dimensionality of the GLstrain space.
  int sampled_space_rank = dim;

  // In this example, use the reference GLstrain space axes for the sampling grid axes. To rotate
  //   the sampling grid, change "axes" so that sampling grid axes are column vectors:
  //
  //     GLstrain_value = axes * sampling_grid_value
  //
  Eigen::MatrixXd axes = Eigen::MatrixXd::Identity(dim, sampled_space_rank);

  // We'll "trick" the enumerator by creating a "dummy subwedge" which defines the axes of the space
  //   we want to sample.
  std::vector<SymRepTools::SubWedge> subwedges;
  subwedges.push_back(SymRepTools::SubWedge({SymRepTools::IrrepWedge::make_dummy_irrep_wedge(axes)}));

  EXPECT_EQ(almost_equal(subwedges[0].trans_mat(), axes, tol), true);
  EXPECT_EQ(subwedges.size(), 1);

  // Specify a grid with three points along each SubWedge axis, with minimum value of -0.5
  //   and increment value of 0.5.
  double _min = -0.5;
  double _incr = 0.5;
  int n_increment = 3;
  Eigen::VectorXd min_value = Eigen::VectorXd::Constant(sampled_space_rank, _min);
  Eigen::VectorXd max_value = Eigen::VectorXd::Constant(sampled_space_rank, _min + (n_increment - 1) * _incr + _incr / 10.);
  Eigen::VectorXd increment_value = Eigen::VectorXd::Constant(sampled_space_rank, _incr);

  // Construct the CASM::ConfigEnumStrain enumerator
  bool auto_range = false;  // explanation TODO
  bool trim_corners = false; // explanation TODO
  CASM::ConfigEnumStrain enumerator {config_input, subwedges, min_value, max_value, increment_value,
                                     dof_key, auto_range, trim_corners};

  std::vector<CASM::Configuration> configurations;
  std::copy(enumerator.begin(), enumerator.end(), std::back_inserter(configurations));

  // Check the number of configurations:
  // - Each subwedge has a grid with n_increment points along each dimension, so the total number
  //   of enumerated configurations is pow(n_increment, dim) * subwedges.size()
  EXPECT_EQ(configurations.size(), pow(n_increment, sampled_space_rank)*subwedges.size());
};


namespace enumeration_test_impl {

  // Helper functions to make the examples easier to read
  SymGroupRep const &global_dof_symrep(Structure const &prim, DoFKey global_dof_key) {
    return prim.factor_group().representation(prim.global_dof_symrep_ID(global_dof_key));
  }

  SymGroupRep const &global_dof_symrep(ConfigEnumInput const &config_input, DoFKey global_dof_key) {
    return global_dof_symrep(config_input.configuration().supercell().prim(), global_dof_key);
  }

  SymGroup make_point_group(ConfigEnumInput const &config_input) {
    return make_point_group(
             make_invariant_group(config_input),
             config_input.configuration().supercell().sym_info().supercell_lattice());
  }

  // check if there is a solution X, for:  B = A * X
  bool is_equivalent_column_space(Eigen::MatrixXcd const &A, Eigen::MatrixXcd const &B, double tol) {
    Eigen::MatrixXcd X = A.fullPivHouseholderQr().solve(B);
    return (A * X).isApprox(B, tol);
  }
}

ExampleEnumerationSimpleCubicConfigEnumStrain::ExampleEnumerationSimpleCubicConfigEnumStrain():
  title("ExampleEnumerationSimpleCubicConfigEnumStrain"),
  shared_prim(std::make_shared<CASM::Structure const>(test::SimpleCubicGLstrain())),
  project_settings(make_default_project_settings(*shared_prim, title)),
  primclex(project_settings, shared_prim) {

  int begin_volume {1};
  int end_volume {2};
  std::string dirs {"abc"};
  Eigen::Matrix3i generating_matrix {Eigen::Matrix3i::Identity()};
  CASM::xtal::ScelEnumProps enumeration_params {begin_volume, end_volume, dirs, generating_matrix};
  bool existing_only = false;

  // The ScelEnumByProps variant that accepts a PrimClex in the constructor inserts Supercells into
  // the Supercell database available at `primclex.db<Supercell>()` as it constructs them.
  CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

  // Increments the enumerator iterators to construct all Supercell
  int count = std::distance(enumerator.begin(), enumerator.end());
  EXPECT_EQ(count, 1);
  EXPECT_EQ(primclex.db<Supercell>().size(), 1);

  // Expected symmetry adapted GLstrain axes for a simple cubic structure
  Eigen::MatrixXcd Q {6, 6};

  // "q1"
  Q.col(0) << 0.577350269190, 0.577350269190, 0.577350269190, 0.000000000000, 0.000000000000, 0.000000000000;
  expected_irrep_subspace.push_back(Q.block(0, 0, 6, 1));

  // "q2, q3"
  Q.col(1) << -0.408248290464, -0.408248290464, 0.816496580928, 0.000000000000, 0.000000000000, 0.000000000000;
  Q.col(2) << 0.707106781187, -0.707106781187, -0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000;
  expected_irrep_subspace.push_back(Q.block(0, 1, 6, 2));

  // "q4, q5, q6"
  Q.col(3) << 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000;
  Q.col(4) <<  0.000000000000, 0.000000000000, 0.000000000000, -1.000000000000, 0.000000000000, 0.000000000000;
  Q.col(5) << 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000;
  expected_irrep_subspace.push_back(Q.block(0, 3, 6, 3));

  tol = 1e-10;

}

// Check that the IrrepInfo trans_mat value matches expectations:
void ExampleEnumerationSimpleCubicConfigEnumStrain::check_irrep(SymRepTools::IrrepInfo const &irrep) {

  using namespace enumeration_test_impl; // for is_equivalent_column_space

  // The SymRepTools::IrrepInfo::trans_mat is a matrix that multiplies a vector in the original
  //   reference space and transforms it into irreducible representation vector space:
  //
  //      v_irrep_space = irrep.trans_mat * v_reference_space
  //
  // The reference space has dimensions (vector_dim  x vector_dim)
  // An irreducible representation vector space has dimensions (vector_dim  x irrep_dim)
  // The IrrepInfo::trans_mat has dimensions (irrep_dim x vector_dim)
  //
  // For this example, GLstrain has vector_dim=6, and the irrep_dim depends on the particular
  //   irreducible representation. The irreducible representations are named with two indices:
  //
  //      irrep name: irrep_<index>_<index>  // TODO: explain mult/indices
  //
  //   irrep_1_1: irrep_dim=1, irrep_axes=["q1"]
  //     q1 = 0.577350269190, 0.577350269190, 0.577350269190, 0.000000000000, 0.000000000000, 0.000000000000
  //   irrep_2_1: irrep_dim=1, irrep_axes=["q2", "q3"]
  //     q2 = -0.408248290464, -0.408248290464, 0.816496580928, 0.000000000000, 0.000000000000, 0.000000000000
  //     q3 = 0.707106781187, -0.707106781187, -0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000
  //   irrep_2_1: irrep_dim=1, irrep_axes=["q4", "q5", "q6"]
  //     q4 = 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000
  //     q5 = 0.000000000000, 0.000000000000, 0.000000000000, -1.000000000000, 0.000000000000, 0.000000000000
  //     q6 = 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 1.000000000000, 0.000000000000

  EXPECT_EQ(irrep.trans_mat.rows(), irrep.irrep_dim());
  EXPECT_EQ(irrep.trans_mat.cols(), irrep.vector_dim());

  // The simple cubic GLstrain irreps are all distinct by irrep_dim,
  //   so use those to identify and then check the subspaces
  if(irrep.irrep_dim() == 1) {  // irrep_1_1: ["q1"]
    EXPECT_EQ(irrep.complex, false);
    EXPECT_EQ(irrep.pseudo_irrep, false);
    EXPECT_EQ(irrep.index, 0);
    EXPECT_EQ(is_equivalent_column_space(irrep.trans_mat.transpose(), expected_irrep_subspace[0], tol), true);
    // TODO: IrrepInfo::directions, not currently checked
    // TODO: IrrepInfo::characters, not currently checked
  }
  else if(irrep.irrep_dim() == 2) { // irrep_2_1: ["q2", "q3"]
    EXPECT_EQ(irrep.complex, false);
    EXPECT_EQ(irrep.pseudo_irrep, false);
    EXPECT_EQ(irrep.index, 0);
    EXPECT_EQ(is_equivalent_column_space(irrep.trans_mat.transpose(), expected_irrep_subspace[1], tol), true);
    // TODO: IrrepInfo::directions, not currently checked
    // TODO: IrrepInfo::characters, not currently checked
  }
  else if(irrep.irrep_dim() == 3) { // irrep_3_1: ["q4", "q5", "q6"]
    EXPECT_EQ(irrep.complex, false);
    EXPECT_EQ(irrep.pseudo_irrep, false);
    EXPECT_EQ(irrep.index, 0);
    EXPECT_EQ(is_equivalent_column_space(irrep.trans_mat.transpose(), expected_irrep_subspace[2], tol), true);
    // TODO: IrrepInfo::directions, not currently checked
    // TODO: IrrepInfo::characters, not currently checked
  }
  else {
    EXPECT_EQ(true, false) << " for SimpleCubic GLstrain space, no irrep should have irrep_dim > 3";
  }
}
