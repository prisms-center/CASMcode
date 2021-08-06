#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Strain.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/io/stream/SymInfo_stream_io.hh"
#include "gtest/gtest.h"

/// Crystal properties
/// ------------------
///
/// CASM needs to know how to use crystal properties in several places:
/// - Site and global properties in xtal::SimpleStructure
/// - Degrees of freedoms in xtal::BasicStructure through xtal::DoFSet and
/// xtal::SiteDoFSet
/// - Molecule properties in xtal::SpeciesProperty
/// - Calculated properties in MappedProperties
///
/// Examples of various property types, and their name strings, include:
/// - occupation ("occ"): which atom or molecule occupies a particular crystal
/// site
/// - displacement ("disp"): displacement of crystal sites from an ideal
/// location
/// - strain ("Ustrain", "GLstrain", etc.): how lattice vectors are transformed
/// from the ideal
///   lattice vectors, under various metrics ("Ustrain": stretch tensor,
///   "GLstrain": Green-Lagrange strain metric, etc.)
/// - energy ("energy"): the crystal energy, relative to user chosen reference
/// states
/// - force ("force"): forces on atoms present without/with allowing atomic and
/// lattice relaxation
///
/// The AnisoValTraits class contains all the necessary information for CASM to
/// deal with different types of properties. For each property type there must
/// exist one AnisoValTraits object which controls how CASM treats the property,
/// including:
/// - the property name
/// - whether the property is a local (site) property or global property
/// - whether the property is continuous or discrete
/// - the dimension of the property value vector
/// - the standard basis for the property value (what each component means)
/// - how symmetry operations transform the property value
/// - etc.
///
/// The AnisoValTraits class contains two constructors. One allows defining a
/// new property type which gets placed in a single static container. The other,
/// with just a "name" argument, returns a copy of the already existing
/// AnisoValTraits object for the property type with the given name.
///
/// A particular instance of a property, such as in BasicStructure or
/// SimpleStructure, is labeled with a property name string. Property name
/// strings must end with the property type name (i.e. "occ", "disp", "Ustrain",
/// "energy", etc.) and may also include a modifier describing the particular
/// property (i.e. "formation" in "formation_energy"). If a modifier exists, it
/// must be separated from the property type by an underscore character ('_').
/// The name of a custom property type may not include an underscore.
///
/// Thus, whenever a particular property name is encountered, the property type
/// name can be determined and used to lookup the corresponding AnisoValTraits
/// object.
///
/// Examples:
/// - property name: "energy" -> property type (AnisoValTraits) name: "energy"
/// - property name: "relaxed_energy" -> property type (AnisoValTraits) name:
/// "energy"
/// - property name: "formation_energy" -> property type (AnisoValTraits) name:
/// "energy"
/// - property name: "relaxedenergy" -> Error (because there is no underscore)
/// - property name: "Ustrain" -> property type (AnisoValTraits) name: "Ustrain"
/// - property name: "strain" ->  Error (because there no strain type prefix)
///

namespace {

/// Adapts CASM::xtal::SymOp to more complicated CASM::SymOp type before
/// printing description
std::string brief_description(CASM::xtal::SymOp const &symop,
                              CASM::xtal::Lattice const &lattice) {
  return brief_description(
      CASM::adapter::Adapter<CASM::SymOp, CASM::xtal::SymOp>()(symop), lattice);
}
}  // namespace

TEST(ExampleCrystallographyProperties, ApplySymmetry) {
  // First, construct a Lattice
  Eigen::Vector3d a{1.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 1.0, 0.0};
  Eigen::Vector3d c{0.0, 0.0, 1.0};
  CASM::xtal::Lattice lattice{a, b, c};

  // Construct the lattice point group
  std::vector<CASM::xtal::SymOp> point_group =
      CASM::xtal::make_point_group(lattice);

  // Construct a map of properties: pairs of <property name>: <property value
  // vector>
  Eigen::VectorXd dx{3};
  dx << 0.01, 0.0, 0.0;

  // Construct a strain metric vector
  //
  // Deformation tensor, F:
  //     after_strain_lat_column_mat = F * before_strain_lat_column_mat
  Eigen::Matrix3d F;
  F << 1.01, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;

  CASM::StrainConverter strain_converter{"U"};
  Eigen::VectorXd Ustrain_Exx_unrolled =
      strain_converter.unrolled_strain_metric(F);

  // Property vectors
  std::map<std::string, Eigen::VectorXd> properties{
      {CASM::AnisoValTraits::disp().name(), dx},  // a displacement along x
      {CASM::AnisoValTraits::strain("U").name(),
       Ustrain_Exx_unrolled}  // Right-stretch strain tensor Exx
  };

  // // Uncomment to print properties
  // std::cout << "Original properties:\n";
  // for(auto const &property_pair : properties) {
  //   std::string property_name = property_pair.first;
  //   Eigen::VectorXd property_value = property_pair.second;
  //   std::cout << "- " << property_name << ": " << property_value.transpose()
  //   << std::endl;
  // }
  // std::cout << std::endl;

  // For each symmetry operation in the point group, transform the properties
  int i = 0;
  for (CASM::xtal::SymOp symop : point_group) {
    // // Uncomment to print point group operations
    // std::cout << "Lattice point group operation " << i++ << ": " <<
    // brief_description(symop, lattice) << std::endl;
    std::map<std::string, Eigen::VectorXd> transformed_properties;

    for (auto const &property_pair : properties) {
      std::string property_name = property_pair.first;
      Eigen::VectorXd property_value = property_pair.second;
      CASM::AnisoValTraits traits{property_name};
      auto matrix = traits.symop_to_matrix(symop.matrix, symop.translation,
                                           symop.is_time_reversal_active);
      transformed_properties[property_name] = matrix * property_value;
    }

    // // Uncomment to print properties
    // std::cout << "Transformed properties:\n";
    // for(auto const &property_pair : transformed_properties) {
    //   std::string property_name = property_pair.first;
    //   Eigen::VectorXd property_value = property_pair.second;
    //   std::cout << "- " << property_name << ": " <<
    //   property_value.transpose() << std::endl;
    // }
    // std::cout << std::endl;
  }
}
