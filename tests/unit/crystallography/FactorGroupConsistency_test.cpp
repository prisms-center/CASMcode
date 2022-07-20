#include <algorithm>
#include <fstream>
#include <ostream>
#include <vector>

#include "Common.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/SymTypeComparator.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/global/eigen.hh"
#include "gtest/gtest.h"
using namespace CASM;

namespace {
void print_poscar(const xtal::BasicStructure &printable,
                  std::ostream &outstream) {
  CASM::VaspIO::PrintPOSCAR p(xtal::make_simple_structure(printable));
  p.sort();
  p.print(outstream);
  return;
}

// TODO: Should this go somewhere?
xtal::SymOpVector make_factor_group_via_mapping(
    const xtal::BasicStructure &struc, double tol) {
  xtal::SimpleStructure sstruc = xtal::make_simple_structure(struc);
  xtal::StrucMapper mapper((xtal::SimpleStrucMapCalculator(sstruc)));
  auto sym_set = mapper.map_deformed_struc_impose_lattice(
      sstruc, xtal::Lattice(sstruc.lat_column_mat), 0,
      xtal::StrucMapping::big_inf(), tol);

  return adapter::Adapter<xtal::SymOpVector, decltype(sym_set)>()(sym_set);
}
}  // namespace

class CrystalGroupTest : public testing::Test {
 protected:
  void add_primitive_structure(std::string prim_file_in_data_dir,
                               int factor_group_size) {
    auto test_files_dir = test::data_dir("crystallography");
    fs::ifstream s;

    s.open(test_files_dir / prim_file_in_data_dir);
    primitive_structures.emplace_back(
        xtal::BasicStructure::from_poscar_stream(s));
    expected_primitive_factor_group_size.emplace_back(factor_group_size);
    s.close();

    return;
  }

  void SetUp() override {
    add_primitive_structure("POS1.txt", 16);
    add_primitive_structure("PRIM4.txt", 48);
    add_primitive_structure("hcp_mg.vasp", 24);

    Eigen::Matrix3l super_mat_0;
    super_mat_0 << 1, 0, 0, 0, 1, 0, 0, 0, 4;
    prim_to_super_trasnformation_matrices.emplace_back(std::move(super_mat_0));

    Eigen::Matrix3l super_mat_1;
    super_mat_1 << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    prim_to_super_trasnformation_matrices.emplace_back(std::move(super_mat_1));

    for (const xtal::BasicStructure &prim : primitive_structures) {
      superstructures.push_back({});
      for (const Eigen::Matrix3l &mat : prim_to_super_trasnformation_matrices) {
        superstructures.back().emplace_back(
            xtal::make_superstructure(prim, mat));
      }
    }

    return;
  }

  std::vector<xtal::BasicStructure> primitive_structures;
  std::vector<int> expected_primitive_factor_group_size;

  // Each inner vector is all the primitive structures with the same
  // transformation applied to each primitive structure
  std::vector<std::vector<xtal::BasicStructure>> superstructures;
  std::vector<Eigen::Matrix3l> prim_to_super_trasnformation_matrices;

  double tol = 1e-5;

  void compare_factor_groups(const xtal::BasicStructure &struc) {
    xtal::SymOpVector mapping_factor_group =
        ::make_factor_group_via_mapping(struc, tol);
    xtal::SymOpVector factor_group = xtal::make_factor_group(struc, tol);

    EXPECT_EQ(mapping_factor_group.size(), factor_group.size());

    for (const xtal::SymOp &mapping_op : mapping_factor_group) {
      UnaryCompare_f<xtal::SymOpPeriodicCompare_f> matches_mapping_op(
          mapping_op, struc.lattice(), tol);
      EXPECT_TRUE(std::find_if(factor_group.begin(), factor_group.end(),
                               matches_mapping_op) != factor_group.end());
    }
  }
};

TEST_F(CrystalGroupTest, PrimitiveFactorGroupSizes) {
  for (int i = 0; i < primitive_structures.size(); ++i) {
    auto factor_group = xtal::make_factor_group(primitive_structures[i], tol);
    int expected_size = expected_primitive_factor_group_size[i];
    EXPECT_EQ(factor_group.size(), expected_size);
  }
}

TEST_F(CrystalGroupTest, PrimitiveFactorGroupViaMapping) {
  for (const xtal::BasicStructure &prim : primitive_structures) {
    compare_factor_groups(prim);
  }
}

TEST_F(CrystalGroupTest, SuperstructureFactorGroupViaMapping) {
  for (int i = 0; i < superstructures.size(); ++i) {
    for (const xtal::BasicStructure &superstruc : superstructures[i]) {
      compare_factor_groups(superstruc);
    }
  }
}
