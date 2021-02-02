#ifndef CASM_StrucMapping
#define CASM_StrucMapping

#include <unordered_set>
#include <vector>

#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {
namespace xtal {
class Lattice;
class LatticeMap;
class SimpleStructure;
class StrucMapCalculatorInterface;

// In this file:
struct LatticeNode;
struct AssignmentNode;
struct MappingNode;
class StrucMapper;

namespace StrucMapping {

typedef std::vector<std::vector<Index>> PermuteOpVector;
/// \brief Very large value used to denote invalid or impossible mapping
inline double big_inf() { return 10E20; }

/// \brief use as default value to initialize mapping costs. Does not indicate
/// ivalidity
inline double small_inf() { return 10E10; }

inline bool is_inf(double _val) { return _val > small_inf() / 2.; }

/// \brief Calculate the basis cost function of a MappingNode as the normalized
/// mean-square displacement of its atoms The displacement vectors are deformed
/// to the CHILD structure's coordinate system before calculating \param
/// mapped_result a proposed mapping; both lattice_node and atomic_node of
/// mapped_result must be initialized \param Nsites number of atoms (excluding
/// vacancies) in the relaxed structure, for proper normalization result is
/// dimensionless, having been normalized by the squared radius of a sphere
/// having the same atomic volume of CHILD structure
double atomic_cost_child(const MappingNode &mapped_result, Index Nsites);

/// \brief Calculate the basis cost function of a MappingNode as the normalized
/// mean-square displacement of its atoms The displacement vectors are deformed
/// to the PARENT structure's coordinate system before calculating \param
/// mapped_result a proposed mapping; both lattice_node and atomic_node of
/// mapped_result must be initialized \param Nsites number of atoms (excluding
/// vacancies) in the relaxed structure, for proper normalization result is
/// dimensionless, having been normalized by the squared radius of a sphere
/// having the same atomic volume of PARENT structure
double atomic_cost_parent(const MappingNode &mapped_result, Index Nsites);

/// \brief Calculate the basis cost function of a MappingNode as the average of
/// atomic_cost_child and atomic_cost_parent
double atomic_cost(const MappingNode &mapped_config, Index Nsites);

/// \brief Calculate the symmetry breaking atomic cost of a MappingNode
double atomic_cost(
    const MappingNode &basic_mapping_node, SymOpVector &factor_group,
    const std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,
                                               Index>> &permutation_group,
    Index Nsites);

}  // namespace StrucMapping

/// \brief Class describing the lattice-mapping portion of a particular mapping
/// A general map for child_struc onto parent_struc may require forming a
/// supercell of parent_struc (most commonly) and/or of child_struc. As such,
/// the LatticeNode is specified in terms of superlattices of both parent_struc
/// and child_struc, as well as deformation and rotation information sufficient
/// to fully define the lattice map
struct LatticeNode {
  /// \brief stretch tensor that takes child superlattice from its de-rotated,
  /// deformed state to its ideal, parent-mapped state We utilize a convention
  /// in which 'stretch' is applied *after* 'isometry', so 'stretch' is the
  /// *left* stretch tensor
  Eigen::Matrix3d stretch;

  /// \brief cartesian rotation/reflection/rotoreflection that rotates the child
  /// superlattice to its de-rotated, deformed state
  Eigen::Matrix3d isometry;

  /// \brief PrimGrid for supercell of parent structure
  /// The parent structure defines the ideal strain state, and the child
  /// structure must undergo an 'idealization' deformation to perfectly tile the
  /// parent structure supercell.
  ///   Define:
  ///      - 'Sp(i)' as the ideal parent supercell
  ///      - 'Sc(i)' as the ideal child supercell
  ///      - 'Sc(d)' as the deformed child supercell
  ///      - 'R' as cartesian rotation/reflection/rotoreflection ('isometry')
  ///      - 'U' as stretch tensor ('stretch')
  ///   Then, the ideal and deformed superlattices satisfy
  ///      Sp(i) == Sc(i)
  ///   and
  ///      Sc(i) == U * R * Sc(d)
  ///   and
  ///      Sc(d) == R.transpose() * U.inverse() * Sc(i)
  Superlattice parent;

  /// \brief PrimGrid for supercell of child structure
  /// The child lattice is recorded in its idealized state (de-rotated and
  /// un-deformed) The transformation matrices 'isometry' and 'stretch' record
  /// the idealization transformation that yielded the ideal child lattice from
  /// its initial state
  ///   Define:
  ///      - 'Sp(i)' as the ideal parent supercell
  ///      - 'Sc(i)' as the ideal child supercell
  ///      - 'Sc(d)' as the deformed child supercell
  ///      - 'R' as cartesian rotation/reflection/rotoreflection ('isometry')
  ///      - 'U' as stretch tensor ('stretch')
  ///   Then, the ideal and deformed superlattices satisfy
  ///      Sp(i) == Sc(i)
  ///   and
  ///      Sc(i) == U * R * Sc(d)
  ///   and
  ///      Sc(d) == R.transpose() * U.inverse() * Sc(i)
  Superlattice child;

  /// \brief strain_cost of the LatticeNode
  double cost;

  /// \brief Construct with ideal parent_scel and deformed child_scel, which are
  /// related by a deformation tensor
  /// @param parent_scel and @param child_scel are integer combinations of the
  /// primitive cells 'parent_prim' and 'child_prim', respectively
  /// @param child_N_atom is number of sites in the child
  /// \param _cost is used to specify mapping cost (in default case -- big_inf()
  /// -- cost will be calculated from scratch)
  LatticeNode(Lattice const &parent_prim, Lattice const &parent_scel,
              Lattice const &child_prim, Lattice const &child_scel,
              Index child_N_atom, double _cost = StrucMapping::big_inf());

  /// \brief Construct with LatticeMap, which relates a supercell of parent_prim
  /// to a supercell of child_prim
  /// @param lat_map specifies a supercell that is a supercell of parent_prim,
  /// but also an idealized supercell of child_prim
  /// @param child_N_atom is number of sites in the child
  LatticeNode(LatticeMap const &_lat_map, Lattice const &parent_prim,
              Lattice const &child_prim);
};

/// \brief Compare two LatticeMap objects, based on their mapping cost first,
/// followed by PrimGrid transformation matrices
bool less(LatticeNode const &A, LatticeNode const &B, double cost_tol);

/// \brief returns true if cost values and parent/child supercell
/// transformations are same for A and B
bool identical(LatticeNode const &A, LatticeNode const &B, double cost_tol);

/// \brief Structure to encode the solution of a constrained atomic assignmnet
/// problem This describes the permutation, translation, and time-reversal of
/// the atoms of a child structure to bring them into registration with the
/// atoms of a parent structure (assuming periodic boundary conditions). Also
/// records the constrained and unconstrained assignment costs
struct AssignmentNode {
  AssignmentNode(double _cost_tol = 1e-6)
      : time_reversal(false), cost(0), m_cost_tol(_cost_tol) {}

  /// \brief Mapping translation from child to parent
  /// Defined such that
  ///   translation=parent_coord.col(i)-child_coord.col(permutation[i])+displacement.col(i)
  ///   (for each 'i', modulo a lattice translation)
  /// This definition assumes that 'child_coord' has been de-rotated and
  /// un-deformed according to a particular lattice mapping (as defined by a
  /// LatticeNode object)
  Eigen::Vector3d translation;

  /// \brief time_reversal relationship between child and parent
  bool time_reversal;

  /// \brief parent->child site assignments that have been forced on at this
  /// node for element 'el' such that forced_on.count(el)==1, 'el' denotes that
  /// child_coord.col(el.second) maps onto parent_coord.col(el.first)
  std::set<std::pair<Index, Index>> forced_on;

  /// \brief 'real' indices of rows in the reduced 'cost_mat'
  /// When a site assignment {i,j} is added to forced_on, row 'i' and col 'j'
  /// are removed from cost_mat An element cost_mat(k,l) in the 'cost_mat'
  /// corresponds to the original element at (irow[k],icol[l]) in the original
  /// cost_mat
  std::vector<Index> irow;

  /// \brief 'real' indices of columns in the reduced 'cost_mat'
  /// When a site assignment {i,j} is added to forced_on, row 'i' and col 'j'
  /// are removed from cost_mat An element cost_mat(k,l) in the 'cost_mat'
  /// corresponds to the original element at (irow[k],icol[l]) in the original
  /// cost_mat
  std::vector<Index> icol;

  /// \brief Solution of the assignment problem for the reduced 'cost_mat'
  /// An assignment {k,l} in the reduced problem occurs when assignment[k]==l
  /// In the unreduced problem, this assignment corresponds to {irow[k],icol[l]}
  std::vector<Index> assignment;

  /// \brief Cost matrix for an assignment problem, which may be a reduced
  /// assignment problem if forced_on.size()>0 cost_mat(i,j) is the cost of
  /// mapping child site 'j' onto parent atom 'i'. If parent structure allows
  /// vacancies cost_mat may have more columns than child sites. These
  /// correspond to 'virtual vacancies', which have zero cost of mapping onto a
  /// parent site that allows vacancies and infinite cost of mapping onto a
  /// parent site that does not allow vacancies. In the case of a reduced
  /// assignment problem, cost_mat.cols()=icol.size() and
  /// cost_mat.rows()=irow.size()
  Eigen::MatrixXd cost_mat;

  /// \brief Total cost of best solution to the constrained assignment problem
  /// having some forced_on assignments
  double cost;

  double cost_tol() const { return m_cost_tol; }

  /// \brief True if cost matrix and assignment vector are uninitialized
  bool empty() const { return cost_mat.size() == 0 && assignment.empty(); }

  /// \brief Compares time_reversal and translation
  /// for time_reversal, false is less than true
  /// for translation, elements are compared lexicographically
  bool operator<(AssignmentNode const &other) const;

  /// \brief Combines constrained vector HungarianNode::assignment and
  /// HungarianNode::forced_on to obtain total permutation vector.
  /// child_struc_after.site[i] after assignment is given by
  /// child_struc_before.site[result[i]] before assignment
  std::vector<Index> permutation() const {
    std::vector<Index> result(assignment.size() + forced_on.size(), 0);
    for (auto const &pair : forced_on) {
      result[pair.first] = pair.second;
    }
    for (Index i = 0; i < assignment.size(); ++i) {
      result[irow[i]] = icol[assignment[i]];
    }
    return result;
  }

 private:
  double m_cost_tol;
};

/// \brief true if time_reversal and translation are identical
bool identical(AssignmentNode const &A, AssignmentNode const &B);

/// Data structure holding a potential mapping, which consists of deformation,
/// occupation array, and displacement field
struct MappingNode {
  // typedefs to provide flexibility if we eventually change to a
  // Eigen::Matrix<3,Eigen::Dynamic>
  typedef Eigen::MatrixXd DisplacementMatrix;

  // Can treat as a Eigen::VectorXd
  using Displacement = DisplacementMatrix::ColXpr;
  using ConstDisplacement = DisplacementMatrix::ConstColXpr;

  // Label molecules as name and occupant index (optional, default 0)
  using MoleculeLabel = std::pair<std::string, Index>;

  using AtomIndexSet = std::set<Index>;

  using MoleculeMap = std::vector<AtomIndexSet>;

  LatticeNode lattice_node;
  AssignmentNode atomic_node;
  double lattice_weight;
  double atomic_weight;

  /// \brief true if assignment problem is not yet known to be insoluable --
  /// default true
  bool is_viable;

  /// \brief true if assignment has been checked for physical validity and
  /// passes -- default false
  mutable bool is_valid;

  /// \brief true if node has been partitioned into sub-nodes for generalized
  /// k-best assignment problem -- default false
  mutable bool is_partitioned;

  /// \brief total, finalized cost, populated by a StrucMapCalculator.
  /// Not guaranteed to be a linear function of lattice_node.cost and
  /// atomic_node.cost
  double cost;

  /// \brief 3xN matrix of displacements for all sites in parent supercell (Va
  /// are included, but set to Zero)
  Eigen::MatrixXd atom_displacement;

  /// permutation lists indices of sites in input structure, as-read, so that
  /// they constitute particular mapping onto parent structure.
  /// If we define index j = atom_permutation[i], this indicates that atom 'j'
  /// of the child superstructure maps onto site 'i' of parent superstructure If
  /// parent has N sites and child has M<N atoms, vacancies are designated by
  /// values j>=M
  std::vector<Index> atom_permutation;

  /// mol_map[j] lists atom indices of parent superstructure that comprise the
  /// molecule at its j'th molecular site
  MoleculeMap mol_map;

  /// list of assigned molecule names
  std::vector<MoleculeLabel> mol_labels;

  /// \brief Static constructor to build an invalid MappingNode, can be used as
  /// return value when no valid mapping exists
  static MappingNode invalid();

  /// \brief construct with lattice node and lattice_weight. Cost is initialized
  /// assuming zero atomic_node cost
  MappingNode(LatticeNode _lattice_node, double _lattice_weight)
      : lattice_node(std::move(_lattice_node)),
        is_viable(true),
        is_valid(false),
        is_partitioned(false) {
    set_lattice_weight(_lattice_weight);
    cost = lattice_weight * lattice_node.cost;
  }

  double cost_tol() const { return atomic_node.cost_tol(); }

  /// \brief set the lattice_weight. Cost is calculated as
  /// cost = lattice_weight*lattice_node.cost + atomic_weight*atomic_node.cost
  /// lattice_weight must be on interval (0.,1.]; atomic_weight
  /// is 1.-lattice_weight
  void set_lattice_weight(double _lw) {
    lattice_weight = max(min(_lw, 1.0), 1e-9);
    atomic_weight = 1. - lattice_weight;
  }

  /// \brief Return pair of integer volumes {Vp, Vc}, where Vp is parent
  /// supercell volume and Vc is child supercell volume
  std::pair<Index, Index> vol_pair() const {
    return std::pair<Index, Index>(lattice_node.parent.size(),
                                   lattice_node.child.size());
  }

  /// \brief non-const calc method solves the assignment problem via
  /// hungarian_method sets is_viable -> false if no solution
  void calc();

  void clear() {
    throw std::runtime_error("MappingNode::clear() not implemented");
  }

  /// \brief convenience method to access MappingNode::lattice_node.isometry
  Eigen::Matrix3d const &isometry() const { return lattice_node.isometry; }

  /// \brief convenience method to access MappingNode::lattice_node.stretch
  Eigen::Matrix3d const &stretch() const { return lattice_node.stretch; }

  /// \brief convenience method to access MappingNode::atomic_node.translation
  Eigen::Vector3d const &translation() const { return atomic_node.translation; }

  /// \brief convenience method to access MappingNode::atomic_node.time_reveral
  bool time_reversal() const { return atomic_node.time_reversal; }

  /// \brief non-const access of i'th atomic displacement of mapped structure
  Displacement disp(Index i) { return atom_displacement.col(i); }

  /// \brief const access i'th atomic displacement of mapped structure
  ConstDisplacement disp(Index i) const { return atom_displacement.col(i); }

  /// Compares cost, lattice_node, atomic_node, and permutation
  /// if costs are tied, compares lattice_node.cost
  ///   if tied, uninitialized atomic_node comes before initialized atomic_node
  ///   if tied, compares lattice_node, using defined comparator
  ///   if tied, compares atomic_node, using defined comparator
  ///   if tied, does lexicographic comparison of permutation
  /// This order is essential for proper behavior of mapping algorithm
  bool operator<(MappingNode const &other) const;
};

/// \brief External accessor for isometry, to provide xtal::SymOp adaptability
inline Eigen::Matrix3d const &get_matrix(MappingNode const &_node) {
  return _node.isometry();
}

/// \brief External accessor for translation, to provide xtal::SymOp
/// adaptability
inline Eigen::Vector3d const &get_translation(MappingNode const &_node) {
  return _node.translation();
}

/// \brief External accessor for time_reversal, to provide xtal::SymOp
/// adaptability
inline bool get_time_reversal(MappingNode const &_node) {
  return _node.time_reversal();
}

/// A class for mapping an arbitrary 'child' crystal structure as a deformation
/// of a 'parent' crystal structure StrucMapper manages options for the mapping
/// algorithm and mapping cost function It also caches some information about
/// supercell lattices to improve speed of mapping multiple children crystals
/// onto a single parent structure
///
class StrucMapper {
 public:
  using LatMapType = std::map<Index, std::vector<Lattice>>;

  enum Options {
    none = 0,
    strict = (1u << 0),
    robust = (1u << 1),
    sym_strain = (1u << 2),
    sym_basis = (1u << 3),
    // soft_va_limit ensures that if no supercell volume satisfies vacancy
    // constraints, the smallest possible volume is used. Default behavior
    // results in no valid mapping
    soft_va_limit = (1u << 4)
  };

  ///\brief Construct and initialize a StrucMapper
  ///
  ///\param _calculator
  ///\parblock
  ///          specialization of base class StrucMapCalculatorInterface, which
  ///          controls the way that the cost_matrix and costs are calculated,
  ///          determines validity of proposed mappings, and finalizes the
  ///          representation of proposed mappings
  ///\endparblock
  ///
  ///\param _lattice_weight
  ///\parblock
  ///          free parameter 'w' in the cost function: total_cost =
  ///          w*lattice_deformation+(1-w)*atomic_deformation can vary between 0
  ///          (completely basis-focused) and 1 (completely lattice-focused)
  ///\endparblock
  ///
  ///\param _max_volume_change
  ///\parblock
  ///          constrains the search space by assuming a limit on allowed volume
  ///          change only taken into account when non-interstitial vacancies
  ///          are allowed in parent structure
  ///\endparblock
  ///
  ///\param _options
  ///\parblock
  ///          specify a combination of StrucMapper::Options using bitwise OR:
  ///          Ex. _options=StrucMapper::robust|StrucMapper::strict Options are:
  ///             'robust': does not assume the imported structure might be
  ///             ideal ('robust' is much slower for importing ideal structures,
  ///                       but if 'robust' is not set and a non-ideal structure
  ///                       is passed, this will be almost always be detected
  ///                       and robust methods will be used instead. Thus,
  ///                       'robust' is slightly faster if imported Structures
  ///                       are *not* ideal)
  ///
  ///             'strict': prevents transformation into canonical form. Tries
  ///             to preserve original orientation of imported structure if
  ///             possible
  ///\endparblock
  ///
  ///\param _cost_tol tolerance for mapping comparisons
  ///
  ///\param _min_va_frac minimum fraction of vacant sites, below this fraction a
  /// mapping will not be considered
  ///
  ///\param _max_va_frac maximum fraction of vacant sites, above this fraction a
  /// mapping will not be considered
  StrucMapper(StrucMapCalculatorInterface const &_calculator,
              double _lattice_weight = 0.5, double _max_volume_change = 0.5,
              int _options = 0,  // this should actually be a bitwise-OR of
                                 // StrucMapper::Options
              double _cost_tol = TOL, double _min_va_frac = 0.,
              double _max_va_frac = 1.);

  ///\brief Tolerance for determining if two mapping-cost values are identical
  double cost_tol() const { return m_cost_tol; }

  ///\brief Tolerance for initializing lattices. For now it is initialized to
  /// CASM::TOL
  double xtal_tol() const { return m_xtal_tol; }

  double lattice_weight() const { return m_lattice_weight; }

  double atomic_weight() const { return 1. - m_lattice_weight; }

  void set_lattice_weight(double _lw) {
    m_lattice_weight = max(min(_lw, 1.0), 1e-9);
  }

  /// \brief Max element considered for integer unimodular matrix
  /// transformations (which define orientation relationship of mapping)
  Index lattice_transformation_range() const {
    return m_lattice_transformation_range;
  }

  /// \brief Max element considered for integer unimodular matrix
  /// transformations (which define orientation relationship of mapping)
  void set_lattice_transformation_range(Index _new_range) {
    m_lattice_transformation_range = _new_range;
  }

  /// \brief Flag that enables the calculation of a symmetrized lattice cost
  /// when performing the lattice maps. This cost only accounts for that part of
  /// the deformation gradient that breaks the symmetry of the parent crystal
  /// structure
  void set_symmetrize_lattice_cost(bool _sym_lat_cost) {
    m_symmetrize_lattice_cost = _sym_lat_cost;
  }

  /// \brief Flag that enables the calculation of a symmetrized atomic cost
  /// when performing the atomic maps. This cost only accounts for that part of
  /// the displacement field that breaks the symmetry of the parent crystal
  /// structure
  void set_symmetrize_atomic_cost(
      bool _sym_atomic_cost, const SymOpVector &factor_group,
      const std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,
                                                 Index>> &permutation_group) {
    m_symmetrize_atomic_cost = _sym_atomic_cost;
    m_calc_ptr->set_sym_invariant_displacement_modes(
        generate_invariant_shuffle_modes(factor_group, permutation_group));
  }

  bool symmetrize_lattice_cost() const { return m_symmetrize_lattice_cost; }
  bool symmetrize_atomic_cost() const { return m_symmetrize_atomic_cost; }

  /// \brief Returns the minimum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  double min_va_frac() const { return m_min_va_frac; }

  /// \brief Sets the minimum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  void set_min_va_frac(double _min_va) { m_min_va_frac = max(_min_va, 0.); }

  /// \brief Returns the maximum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  double max_va_frac() const { return m_max_va_frac; }

  /// \brief Sets the maximum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  void set_max_va_frac(double _max_va) { m_max_va_frac = min(_max_va, 0.99); }

  /// \brief returns bit flag of selected options for this StrucMapper
  int options() const { return m_options; }

  /// \brief returns reference to parent structure
  SimpleStructure const &parent() const;

  ///\brief specify a superlattice of the parent to be searched during mapping
  void add_allowed_lattice(Lattice const &_lat) {
    m_allowed_superlat_map
        [std::abs(round(volume(_lat) / parent().lat_column_mat.determinant()))]
            .push_back(_lat);
  }

  ///\brief clear the list of allowed parent superlattices;
  /// all superlattices will be generated automatically, as needed (default)
  void clear_allowed_lattices() const { m_allowed_superlat_map.clear(); }

  ///\brief returns true if the search of parent superlattices is constrained to
  /// a pre-specified list
  bool lattices_constrained() const { return m_allowed_superlat_map.size(); }

  ///\brief specify to use filtered lattices for mapping. The filter function is
  /// of the form
  ///  bool filter(parent_lattice, proposed_lattice)
  /// where parent_lattice is the primitive lattice of the parent structure, and
  /// proposed_lattice is a proposed superlattice of the parent structure
  void set_filter(
      std::function<bool(Lattice const &, Lattice const &)> _filter_f) {
    m_filtered = true;
    m_filter_f = _filter_f;
    m_superlat_map.clear();
  }

  ///\brief specify not to use filtered lattice for mapping
  void unset_filter() {
    m_filtered = false;
    m_superlat_map.clear();
  }

  ///\brief k-best mappings of ideal child structure onto parent structure
  /// Assumes that child_struc and parent_struc have lattices related by an
  /// integer transformation so that search over lattices can be replaced with
  /// simple matrix solution (faster than typical search) Atomic positions need
  /// not be ideal
  ///
  ///\param child_struc Input structure to be mapped onto parent structure
  ///\param k Number of k-best mapping relations to return.
  ///\param max_cost Search will terminate once no mappings better than max_cost
  /// are found \param min_cost All mappings better than min_cost will be
  /// reported, without contributing to 'k'
  ///
  ///\result std::set<MappingNode> a list of valid mappings, sorted first by
  /// cost, and then other attributes
  std::set<MappingNode> map_ideal_struc(
      const SimpleStructure &child_struc, Index k = 1,
      double max_cost = StrucMapping::big_inf(), double min_cost = -TOL,
      bool keep_invalid = false) const;

  ///\brief k-best mappings of arbitrary child structure onto parent structure,
  /// without simplifying assumptions
  ///       If k and k+1, k+2, etc mappings have identical cost, k will be
  ///       increased until k+n mapping has cost greater than mapping k
  ///
  ///\param child_struc Input structure to be mapped onto parent structure
  ///\param k Number of k-best mapping relations to return
  ///\param max_cost Search will terminate once no mappings better than max_cost
  /// are found \param min_cost All mappings better than min_cost will be
  /// reported, without contributing to 'k' \param keep_invalid If true, invalid
  /// mappings are retained in output; otherwise they are discarded
  ///
  ///\result std::set<MappingNode> a list of valid mappings, sorted first by
  /// cost, and then other attributes
  std::set<MappingNode> map_deformed_struc(
      const SimpleStructure &child_struc, Index k = 1,
      double max_cost = StrucMapping::big_inf(), double min_cost = -TOL,
      bool keep_invalid = false,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  ///\brief k-best mappings of arbitrary child structure onto parent structure
  ///       imposes simplifying assumption that solution is a superstructure of
  ///       parent structure having integer volume on range [min_vol, max_vol]
  ///       If k and k+1, k+2, etc mappings have identical cost, k will be
  ///       increased until k+n mapping has cost greater than mapping k
  ///
  ///\param child_struc Input structure to be mapped onto parent structure
  ///\param k Number of k-best mapping relations to return
  ///\param max_cost Search will terminate once no mappings better than max_cost
  /// are found \param min_cost All mappings better than min_cost will be
  /// reported, without contributing to 'k' \param keep_invalid If true, invalid
  /// mappings are retained in output; otherwise they are discarded
  ///
  ///\result std::set<MappingNode> a list of valid mappings, sorted first by
  /// cost, and then other attributes
  std::set<MappingNode> map_deformed_struc_impose_lattice_vols(
      const SimpleStructure &child_struc, Index min_vol, Index max_vol,
      Index k = 1, double max_cost = StrucMapping::big_inf(),
      double min_cost = -TOL, bool keep_invalid = false,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  ///\brief k-best mappings of arbitrary child structure onto parent structure
  ///       imposes simplifying assumption that solution is THE superstructure
  ///       of parent given by
  ///       @param imposed_lat. Significantly faster than unconstrained search,
  ///       but may miss best solution If k and k+1, k+2, etc mappings have
  ///       identical cost, k will be increased until k+n mapping has cost
  ///       greater than mapping k
  ///
  ///\param child_struc Input structure to be mapped onto parent structure
  ///\param k Number of k-best mapping relations to return
  ///\param max_cost Search will terminate once no mappings better than max_cost
  /// are found \param min_cost All mappings better than min_cost will be
  /// reported, without contributing to 'k' \param keep_invalid If true, invalid
  /// mappings are retained in output; otherwise they are discarded
  ///
  ///\result std::set<MappingNode> a list of valid mappings, sorted first by
  /// cost, and then other attributes
  std::set<MappingNode> map_deformed_struc_impose_lattice(
      const SimpleStructure &child_struc, const Lattice &imposed_lat,
      Index k = 1, double max_cost = StrucMapping::big_inf(),
      double min_cost = -TOL, bool keep_invalid = false,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  ///\brief k-best mappings of arbitrary child structure onto parent structure
  ///       imposes simplifying assumption that the lattice mapping is known and
  ///       specified by
  ///       @param imposed_node. This sets the lattice and setting that defines
  ///       the mapping relation This assumption is much than unconstrained
  ///       search, but may miss best solution If k and k+1, k+2, etc mappings
  ///       have identical cost, k will be increased until k+n mapping has cost
  ///       greater than mapping k

  ///
  ///\param child_struc Input structure to be mapped onto parent structure
  ///\param k Number of k-best mapping relations to return
  ///\param max_cost Search will terminate once no mappings better than max_cost
  /// are found \param min_cost All mappings better than min_cost will be
  /// reported, without contributing to 'k' \param keep_invalid If true, invalid
  /// mappings are retained in output; otherwise they are discarded
  ///
  ///\result std::set<MappingNode> a list of valid mappings, sorted first by
  /// cost, and then other attributes
  std::set<MappingNode> map_deformed_struc_impose_lattice_node(
      const SimpleStructure &child_struc, const LatticeNode &imposed_node,
      Index k = 1, double max_cost = StrucMapping::big_inf(),
      double min_cost = -TOL, bool keep_invalid = false) const;

  ///\brief low-level function. Takes a queue of mappings and use it to seed a
  /// search for k-best
  ///       total mappings. The seed mappings may be partial (having only
  ///       LatticeNode specified, and not HungarianNode) or complete. The
  ///       result is appended to @param queue, and partial and invalid mappings
  ///       are removed from the queue. If k and k+1, k+2, etc mappings have
  ///       identical cost, k will be increased until k+n mapping has cost
  ///       greater than mapping k

  ///
  ///
  ///\param child_struc Input structure to be mapped onto parent structure
  ///\param queue
  ///\parblock
  ///  list of partial mappings to seed search. Partial mappings are removed
  ///  from list as they are searched finalized mappings are inserted into list
  ///  as they are found.
  ///\endparblock
  ///
  ///\param k Number of k-best mapping relations to return
  ///\param max_cost Search will terminate once no mappings better than max_cost
  /// are found \param min_cost All mappings better than min_cost will be
  /// reported, without contributing to 'k' \param keep_invalid If true, invalid
  /// mappings are retained in output; otherwise they are discarded \param
  /// keep_tail If true, any partial or unresolved mappings at back of queue
  /// will be retained (they are deleted otherwise) \param no_partition
  /// \parblock
  ///   If true, search over suboptimal basis permutations will be skipped (The
  ///   optimal mapping will be found, but suboptimal mappings may be excluded
  ///   in the k-best case. Degenerate permutations of atoms will not be
  ///   identified.)
  /// \endparblock
  ///
  ///\result Index specifying the number of complete, valid mappings in the
  /// solution set
  Index k_best_maps_better_than(SimpleStructure const &child_struc,
                                std::set<MappingNode> &queue, Index k = 1,
                                double max_cost = StrucMapping::big_inf(),
                                double min_cost = -TOL,
                                bool keep_invalid = false,
                                bool keep_tail = false,
                                bool no_partiton = false) const;

  StrucMapCalculatorInterface const &calculator() const { return *m_calc_ptr; }

 private:
  /// \brief generate list of partial mapping seeds (i.e., LatticeNode only)
  /// from a list of supercells of the parent structure and a list supercells of
  /// the child structure
  ///
  ///\param k Number of k-best mapping relations to return
  ///\param max_lattice_cost Search will terminate once no lattice mappings
  /// better than max_lattice_cost are found \param min_lattice_cost All lattice
  /// mappings better than min_lattice_cost will be returned, without
  /// contributing to 'k'
  std::set<MappingNode> _seed_k_best_from_super_lats(
      SimpleStructure const &child_struc,
      std::vector<Lattice> const &_parent_scels,
      std::vector<Lattice> const &_child_scels, Index k, double max_strain_cost,
      double min_strain_cost,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  /// \brief construct partial mapping nodes (with uninitialized atomic_node)
  /// based on current settings considers supercells with integer volume between
  /// min_vol and max_vol
  std::set<MappingNode> _seed_from_vol_range(
      SimpleStructure const &child_struc, Index k, Index min_vol, Index max_vol,
      double max_strain_cost, double min_strain_cost,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  ///\brief returns number of species in a SimpleStructure given the current
  /// calculator settings.
  ///       Use instead of sstruc.n_atom() for consistency
  Index _n_species(SimpleStructure const &sstruc) const;

  // Implement filter as std::function<bool(Lattice const&,Lattice const&)>
  bool _filter_lat(Lattice const &_parent_lat,
                   Lattice const &_child_lat) const {
    return m_filter_f(_parent_lat, _child_lat);
  }

  std::pair<Index, Index> _vol_range(const SimpleStructure &child_struc) const;

  notstd::cloneable_ptr<StrucMapCalculatorInterface> m_calc_ptr;

  Eigen::MatrixXd m_strain_gram_mat;

  double m_lattice_weight;
  double m_max_volume_change;
  int m_options;
  double m_cost_tol;
  double m_xtal_tol;
  double m_min_va_frac;
  double m_max_va_frac;

  Index m_lattice_transformation_range;

  bool m_symmetrize_lattice_cost;
  bool m_symmetrize_atomic_cost;

  bool m_filtered;
  std::function<bool(Lattice const &, Lattice const &)> m_filter_f;

  /// Maps the supercell volume to a vector of Lattices with that volume
  mutable LatMapType m_superlat_map;
  mutable LatMapType m_allowed_superlat_map;

  std::vector<Lattice> _lattices_of_vol(Index prim_vol) const;
};

}  // namespace xtal
}  // namespace CASM
#endif
