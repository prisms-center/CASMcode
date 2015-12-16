#ifndef CASM_CLASSES_HH
#define CASM_CLASSES_HH

//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE std::size_t


// /////////////////////////////////////////
// Class declarations


// /////////////////////////////////////////
// Header files

// json IO
#include "casm/casm_io/jsonParser.hh"

// system
#include "casm/system/RuntimeLibrary.hh"

// Global variables & functions
#include "casm/CASM_global_definitions.hh"

// Containers
#include "casm/container/Array.hh"					// template //
#include "casm/container/LinearAlgebra.hh"			// template //
#include "casm/container/Counter.hh"					// template //
#include "casm/container/IsoCounter.hh"					// template //
#include "casm/container/MultiCounter.hh"					// template //
#include "casm/container/Template_Algorithms.hh"		// template functions //
#include "casm/container/Permutation.hh"			  // Depends on Array //
//#include "casm/container/PolyTrie.hh"			  // Depends on Array //

// I/O
#include "casm/casm_io/FormatFlag.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataStream.hh"
#include "casm/casm_io/EigenDataStream.hh"
#include "casm/casm_io/Args.hh"


// Misc
#include "casm/misc/HierarchyID.hh"					// template //

// Crystallography - Coordinate
#include "casm/crystallography/CoordinateSystems.hh"
#include "casm/crystallography/Coordinate.hh"		// contains Lattice*
#include "casm/crystallography/UnitCellCoord.hh"

// Symmetry
#include "casm/symmetry/SymOpRepresentation.hh"			// contains MasterSymGroup*
#include "casm/symmetry/SymPermutation.hh"			        // inherits SymOpRepresentation
#include "casm/symmetry/SymBasisPermute.hh"			        // inherits SymOpRepresentation
#include "casm/symmetry/SymMatrixXd.hh"			        // inherits SymOpRepresentation
#include "casm/symmetry/SymOp.hh"					// contains Coordinate, Lattice*; inherits SymOpRepresentation
#include "casm/symmetry/SymGroup.hh"					// contains SymOp
#include "casm/symmetry/SymGroupRep.hh"				// contains SymOpRepresentation

// Containers - Tensors
#include "casm/container/Tensor.hh"					// template // args SymOp, SymGroup

// Crystallography - Molecule
#include "casm/crystallography/Molecule.hh"			// contains Lattice*, SymOp, Coordinate

// Basis Sets
//#include "casm/basis_set/MathExpressions.hh"
#include "casm/basis_set/DoF.hh"						// template //
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/BasisFunction.hh"			// template // contains DoF, Tensor, SparseTensor; args SymOp
#include "casm/basis_set/OccupantFunction.hh"
#include "casm/basis_set/BasisSet.hh"			// template // contains DoF; args SymOp

// Crystallography - Site, Lattice
#include "casm/crystallography/Site.hh"				// contains Molecule, Coordinate, (Lattice*), BasisSet, DoF*
#include "casm/crystallography/Lattice.hh"			// contains Lattice*, MasterSymGroup, GridPoint<double>; source uses GridPoint<double>

// Clusterography
#include "casm/clusterography/Cluster.hh"			// template // contains Lattice*, SymGroup, (SymOp), Tensor, TensorBasis; args Structure, Lattice, GridPoint<T>
#include "casm/clusterography/SiteCluster.hh"		// contains GenericCluster<Site>, BasisSet
#include "casm/clusterography/Orbit.hh"				// template // contains SymOp, Cluster<T>; args Coordinate, Lattice, Structure
#include "casm/clusterography/OrbitBranch.hh"		// template // contains GenericOrbit< Cluster<T> >, SiteCluster
#include "casm/clusterography/Orbitree.hh"			// template // contains GenericOrbitBranch< Cluster<T> >, Lattice; args GridPoint<T>
#include "casm/clusterography/jsonClust.hh"

#include "casm/clusterography/ClusterFunctions.hh"

// Crystallography - Structure
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/BasicStructure.hh"			// contains Lattice, Specie, Molecule, Site, SiteOrbitBranch, MasterSymGroup, Coordinate
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/crystallography/Structure.hh"			// contains Lattice, Specie, Molecule, Site, SiteOrbitBranch, MasterSymGroup, Coordinate

//Clex-dependent symmetry
#include "casm/symmetry/PermuteIterator.hh"

// Clex                                         // contains things for making configurations and correlations.
#include "casm/clex/Properties.hh"
#include "casm/clex/Correlation.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/ConfigEnum.hh"
#include "casm/clex/ConfigEnumIterator.hh"
#include "casm/clex/ConfigEnumInterpolation.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/clex/CompositionConverter.hh"
#include "casm/clex/DoFManager.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "clex/ConfigSelection.hh"
#include "clex/ConfigIO.hh"
#include "casm/clex/ConfigIOHull.hh"
#include "casm/clex/ConfigIONovelty.hh"
#include "casm/clex/ConfigIOStrucScore.hh"
#include "casm/clex/ConfigIOStrain.hh"
#include "casm/clex/ConfigIOSelected.hh"


// Hull
#include "casm/hull/Hull.hh"

#endif
