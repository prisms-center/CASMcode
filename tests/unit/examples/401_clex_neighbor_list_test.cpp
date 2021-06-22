#include "casm/app/ProjectBuilder.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/crystallography/LinearIndexConverter.hh"
#include "crystallography/TestStructures.hh"  // for test::ZrO_prim
#include "gtest/gtest.h"

// Neighbor lists
// --------------
//
// Neighbor lists are used primarily by Clexulator (the basis function
// evaluating functor class) for efficient lookup of the DoF values that are the
// basis function arguments. The Clexulator is evaluated for a single unit cell
// and written in terms of the DoF values in that unit cell and its neighbors.
// The neighbor list enables quickly finding the neighboring site DoF values, in
// a consistent order, for any unit cell. Neighbor lists may also be useful for
// evaluating other functions of local environment.
//
// Neighbor list index naming conventions:
// - "unitcell_index", a linear index given to all unit cell "within" a
// supercell. Convert
//   between unitcell_index and xtal::UnitCell with
//   xtal::UnitCellIndexConverter.
// - "site_index", a linear index given to all sites "within" a supercell.
// Convert between
//   site_index and xtal::UnitCellCoord with xtal::UnitCellCoordIndexConverter.
//
// The PrimNeighborList provides a list of unit cell that are neighbors of the
// origin unit cell. The PrimNeighborList is constructed with a weight matrix,
// W, that defines the shape of the neighborhood. The PrimNeighborList can be
// expanded as needed to increase the range of the neighborhood without
// affecting the order of the neighbors. The PrimNeighborList has:
// - a "weight_matrix" (Eigen::Matrix3l), W, an integer matrix that defines an
// integer distance
//    metric from the origin unit cell to a neighboring UnitCell.
//   - The canonical order of neighboring UnitCell is obtained by
//   lexicographically sorting
//     [r, i, j, k], where r = (i,j,k).transpose() * W * (i,j,k).
//   - The canonical order of UnitCellCoord is obtained by lexicographically
//   sorting [r, i, j, k, b].
// - "sublat_indices" (std::set<int>), containing the indices of sublattices
// that will be included
//   in the neighbor list of "site_index". For purposes of evaluating basis
//   functions, sublattices that have no DoF can be excluded.
//
// Note:
// - If the lattice vectors are very skewed a careful choice of the weight
// matrix may reduce the
//   size of the neighbor list. In practice, using the weight matrix generated
//   by the function `default_nlist_weight_matrix` is nearly always sufficient.
// - The `default_nlist_sublat_indices` function can be used to generate the set
// of
//  "sublat_indices" that have any DoF (2 or more allowed occupant, 1 or more
//  continuous DoF).
// - The size of ConfigDoF value vectors and matrices is always determined by
// the total number of
//   sublattices, irrespective of the particular "sublat_indices" specified for
//   the PrimNeighborList.
//
//
// The SuperNeighborList takes the ordering of unit cells neighboring the origin
// unit cell from the PrimNeighborList and uses it to construct lists of linear
// indices of neighboring unit cells and neighboring sites for each unit cell in
// a particular supercell. SuperNeighborList has:
// - "unitcells", (std::vector<Index>, one for each unit cell in the supercell),
// provides a vector
//   of "unitcell_index" specifying the unit cells that are neighbors to a
//   particular unit cell
// - "sites", (std::vector<Index>, one for each unit cell in the supercell),
// provides a vector of
//   "site_index" specifying sites that are neighbors to a particular unit celln
//

TEST(ExampleCrystallography, NeighborList) {
  // example generating PrimNeighborList and SuperNeighborList

  // Construct a BasicStructure, representing the infinite crystal and allowed
  // perturbations.
  CASM::xtal::BasicStructure basic_structure = test::ZrO_prim();

  // Constuct and expand a PrimNeighborList:

  // PrimNeighborList constructor arguments
  double tol = basic_structure.lattice().tol();
  Eigen::Matrix3l weight_matrix =
      CASM::default_nlist_weight_matrix(basic_structure, tol);
  std::set<int> nlist_sublat_indices =
      CASM::default_nlist_sublat_indices(basic_structure);

  // Construct a PrimNeighborList
  CASM::PrimNeighborList prim_nlist{weight_matrix, nlist_sublat_indices.begin(),
                                    nlist_sublat_indices.end(),
                                    basic_structure.basis().size()};

  // Initial size == 1, origin is always included
  EXPECT_EQ(prim_nlist.size(), 1);

  // Expand, adding a nearest neighbor unit cells (adds the 6 in-plane nearest
  // neighbors)
  prim_nlist.expand({1, 0, 0});
  EXPECT_EQ(prim_nlist.size(), 7);

  // // Uncomment to print neighboring unit cells of the origin unit cell
  // std::cout << "PrimNeighborList (expanded with {0, 0, 0}):\n";
  // for(CASM::xtal::UnitCell unitcell : prim_nlist) {
  //   std::cout << unitcell.transpose() << std::endl;
  // }
  // std::cout << std::endl;

  // Expand, adding an existing neighbor does nothing
  prim_nlist.expand({1, 0, 0});
  EXPECT_EQ(prim_nlist.size(), 7);

  // Expand to include farther neighbors. Will include 27 unit cells (check as
  // an exercise):
  // - 13 in the plane z==0
  // - 7 each in the planes z==+c, z==-c
  prim_nlist.expand({1, 1, 1});
  EXPECT_EQ(prim_nlist.size(), 27);

  // // Uncomment to print neighboring unit cells of the origin unit cell
  // std::cout << "PrimNeighborList (expanded with {1, 1, 1}):\n";
  // for(CASM::xtal::UnitCell unitcell : prim_nlist) {
  //   std::cout << unitcell.transpose() << std::endl;
  // }
  // std::cout << std::endl;

  // Construct a SuperNeighborList based on the PrimNeighborList

  // The prim to super lattice transformation matrix, T:
  //   super_lattice_column_matrix = T * prim_lattice_column_matrix
  Eigen::Matrix3l T;
  T << 2, 0, 0, 0, 2, 0, 0, 0, 2;

  // Construct index converters for the supercell
  CASM::xtal::UnitCellIndexConverter unitcell_index_converter{T};
  CASM::xtal::UnitCellCoordIndexConverter unitcellcoord_index_converter{
      T, (int)basic_structure.basis().size()};

  // Construct the SuperNeighborList
  CASM::SuperNeighborList supercell_nlist{T, prim_nlist};

  // Find neighbors for each unit cell in the supercell
  for (int i = 0; i < unitcell_index_converter.total_sites(); ++i) {
    // SuperNeighborList::unitcells provides a vector of "unitcell_index"
    // specifying the unit cells that are neighbors to a particular unit cell
    std::vector<CASM::Index> unitcell_index_neighbors_to_i =
        supercell_nlist.unitcells(i);
    EXPECT_EQ(unitcell_index_neighbors_to_i.size(), prim_nlist.size());

    // // Uncomment to print neighboring unit cells of unit cell i
    // std::cout << "\nSuperNeighborList of " <<
    // unitcell_index_converter(i).transpose() << ":\n"; for(CASM::Index
    // unitcell_index : unitcell_index_neighbors_to_i) {
    //   std::cout << unitcell_index_converter(unitcell_index).transpose() <<
    //   std::endl;
    // }
    // std::cout << std::endl;

    // SuperNeighborList::sites provides a vector of "site_index" specifying
    // sites that are neighbors to a particular unit cell
    std::vector<CASM::Index> site_index_neighbors_to_i =
        supercell_nlist.sites(i);
    EXPECT_EQ(site_index_neighbors_to_i.size(),
              prim_nlist.size() * prim_nlist.sublat_indices().size());

    // // Uncomment to print neighboring sites of unit cell i
    // std::cout << "\nSuperNeighborList of " <<
    // unitcellcoord_index_converter(i) << ":\n"; for(CASM::Index site_index :
    // site_index_neighbors_to_i) {
    //   std::cout << unitcellcoord_index_converter(site_index) << std::endl;
    // }
    // std::cout << std::endl;
  }
}
