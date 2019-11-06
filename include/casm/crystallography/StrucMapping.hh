#ifndef CASM_StrucMapping
#define CASM_StrucMapping

#include <vector>
#include <unordered_set>
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/PrimGrid.hh"

namespace CASM {
  namespace xtal {
    class Lattice;
    class LatticeMap;
    class SimpleStructure;

    // In this file:
    struct LatticeNode;
    struct HungarianNode;
    struct MappingNode;
    class StrucMapCalculatorInterface;
    class StrucMapper;

    namespace StrucMapping {
      inline
      double big_inf() {
        return 10E20;
      }

      inline
      double small_inf() {
        return 10E10;
      }

      inline
      bool is_inf(double _val) {
        return _val > small_inf() / 2.;
      }



      /// \brief Calculate the strain cost function of a MappingNode using LatticeMap::calc_strain_cost()
      /// \param Nsites number of atoms in the relaxed structure, for proper normalization
      double strain_cost(double relaxed_lat_vol, const MappingNode &mapped_config, Index Nsites);

      /// \brief Calculate the basis cost function of a MappingNode as the mean-square displacement of its atoms
      /// \param Nsites number of atoms in the relaxed structure, for proper normalization
      double basis_cost(const MappingNode &mapped_config, Index Nsites);


    }

    /// \brief Class describing the lattice-mapping portion of a particular mapping
    /// A general map for child_struc onto parent_struc may require forming a supercell of
    /// parent_struc (most commonly) and/or of child_struc.
    /// As such, the LatticeNode is specified in terms of superlattices of both parent_struc and child_struc,
    /// as well as deformation and rotation information sufficient to fully define the lattice map
    struct LatticeNode {

      /// \brief stretch tensor that takes child superlattice from its de-rotated, deformed state to its ideal, parent-mapped state
      /// We utilize a convention in which 'stretch' is applied *after* 'isometry', so 'stretch' is the *left* stretch tensor
      Eigen::Matrix3d stretch;

      /// \brief cartesian rotation/reflection/rotoreflection that rotates the child superlattice to its de-rotated, deformed state
      Eigen::Matrix3d isometry;

      /// \brief PrimGrid for supercell of parent structure
      /// The parent structure defines the ideal strain state, and the child structure must undergo an
      /// 'idealization' deformation to perfectly tile the parent structure supercell.
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
      PrimGrid parent;

      /// \brief PrimGrid for supercell of child structure
      /// The child lattice is recorded in its idealized state (de-rotated and un-deformed)
      /// The transformation matrices 'isometry' and 'stretch' record the idealization transformation
      /// that yielded the ideal child lattice from its initial state
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
      PrimGrid child;

      /// \brief strain_cost of the LatticeNode
      double cost;

      /// \brief Construct with ideal parent_scel and deformed child_scel, which are related by a deformation tensor
      /// @param parent_scel and @param child_scel are integer combinations of the primitive cells 'parent_prim' and 'child_prim', respectively
      /// @param child_N_atom is number of sites in the child
      LatticeNode(Lattice const &parent_prim,
                  Lattice const &parent_scel,
                  Lattice const &child_prim,
                  Lattice const &child_scel,
                  Index child_N_atom,
                  double _cost = StrucMapping::big_inf());

      /// \brief Construct with LatticeMap, which relates a supercell of parent_prim to a supercell of child_prim
      /// @param lat_map specifies a supercell that is a supercell of parent_prim, but also an idealized supercell of child_prim
      /// @param child_N_atom is number of sites in the child
      LatticeNode(LatticeMap const &_lat_map,
                  Lattice const &parent_prim,
                  Lattice const &child_prim,
                  Index child_N_atom);

      /// \brief Compare two LatticeMap objects, based on their mapping cost
      bool operator<(LatticeNode const &other)const {
        return cost < other.cost - 1e-6;
      }

    };

    inline
    bool identical(LatticeNode const &A, LatticeNode const &B) {
      if(!almost_equal(A.cost, B.cost, 1e-6))
        return false;
      if(A.parent.trans_mat() != B.parent.trans_mat())
        return false;
      if(A.child.trans_mat() != B.child.trans_mat())
        return false;
      return true;
    }

    struct HungarianNode {
      HungarianNode(double _tol = 1e-6) : cost(0), cost_offset(0), m_tol(_tol) {}

      /// \brief Mapping translation from child to parent
      /// Defined such that
      ///   translation=parent_coord.col(i)-child_coord.col(permutation[i])+displacement.col(i)
      ///   (for each 'i', modulo a lattice translation)
      /// This definition assumes that 'child_coord' has been de-rotated and un-deformed according
      /// to a particular lattice mapping (as defined by a LatticeNode object)
      Eigen::Vector3d translation;

      /// \brief parent->child site assignments that have been forced on at this node
      /// for element 'el' such that forced_on.count(el)==1, 'el' denotes that
      /// child_coord.col(el.second) maps onto parent_coord.col(el.first)
      std::set<std::pair<Index, Index> > forced_on;

      /// \brief 'real' indices of rows in the reduced 'cost_mat'
      /// When a site assignment {i,j} is added to forced_on, row 'i' and col 'j' are removed from cost_mat
      /// An element cost_mat(k,l) in the 'cost_mat' corresponds to the original element at (irow[k],icol[l])
      /// in the original cost_mat
      std::vector<Index> irow;

      /// \brief 'real' indices of columns in the reduced 'cost_mat'
      /// When a site assignment {i,j} is added to forced_on, row 'i' and col 'j' are removed from cost_mat
      /// An element cost_mat(k,l) in the 'cost_mat' corresponds to the original element at (irow[k],icol[l])
      /// in the original cost_mat
      std::vector<Index> icol;

      /// \brief Solution of the assignment problem for the reduced 'cost_mat'
      /// An assignment {k,l} in the reduced problem occurs when assignment[k]==l
      /// In the unreduced problem, this assignment corresponds to {irow[k],icol[l]}
      std::vector<Index> assignment;

      /// \brief Cost matrix for an assignment problem, which may be a reduced assignment problem if forced_on.size()>0
      /// In the case of a reduced assignment problem, cost_mat.cols()=icol.size() and cost_mat.rows()=irow.size()
      Eigen::MatrixXd cost_mat;

      /// \brief Total cost of best solution to the constrained assignment problem having some forced_on assignments
      double cost;

      /// \brief Cumulative cost of all forced_on assignments
      /// This is added to best solution of reduced assignment problem to obtain cost of total constrained assignment problem
      double cost_offset;

      double tol() const {
        return m_tol;
      }

      bool empty() const {
        return cost_mat.size() == 0 && assignment.empty();
      }

      void solve() {
        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        cost = hungarian_method(cost_mat, assignment, m_tol) + cost_offset + m_tol * translation.norm() / 10.0;
      }

      bool operator<(HungarianNode const &other)const {
        return cost < other.cost - m_tol;
      }

      std::vector<Index> permutation() const {
        std::vector<Index> result(assignment.size() + forced_on.size(), 0);
        for(auto const &pair : forced_on) {
          result[pair.first] = pair.second;
        }
        for(Index i = 0; i < assignment.size(); ++i) {
          result[irow[i]] = icol[assignment[i]];
        }
        return result;
      }


    private:
      double m_tol;
    };

    inline
    bool identical(HungarianNode const &A, HungarianNode const &B) {
      if(!almost_equal(A.cost, B.cost, 1e-6))
        return false;
      if(!almost_equal(A.translation, B.translation, 1e-6))
        return false;
      if(A.forced_on != B.forced_on)
        return false;
      return true;
    }


    /// Data structure holding a potential mapping, which consists of deformation, occupation array, and displacement field
    struct MappingNode {
      // typedefs to provide flexibility if we eventually change to a Eigen::Matrix<3,Eigen::Dynamic>
      typedef Eigen::MatrixXd DisplacementMatrix;

      // Can treat as a Eigen::VectorXd
      typedef DisplacementMatrix::ColXpr Displacement;
      typedef DisplacementMatrix::ConstColXpr ConstDisplacement;

      /// \brief Static constructor to build an invalid MappingNode, can be used as return value when no valid mapping exists
      static MappingNode invalid();

      MappingNode(LatticeNode _lat_node, double _strain_weight) :
        lat_node(std::move(_lat_node)),
        is_viable(true),
        is_valid(false),
        is_partitioned(false) {
        set_strain_weight(_strain_weight);
        cost = strain_weight * lat_node.cost;
      }

      double tol() const {
        return basis_node.tol();
      }

      void set_strain_weight(double _lw) {
        strain_weight = max(min(_lw, 1.0), 1e-9);
        basis_weight = 1. - strain_weight;
      }

      std::pair<Index, Index> vol_pair() const {
        return std::pair<Index, Index>(lat_node.parent.size(), lat_node.child.size());
      }

      /// \brief non-const calc method solves the assignment problem via hungarian_method
      /// sets is_viable -> false if no solution
      void calc();

      void clear() {
        throw std::runtime_error("MappingNode::clear() not implemented");
      }

      Eigen::Matrix3d isometry() const {
        return lat_node.isometry;
      }

      Eigen::Matrix3d stretch() const {
        return lat_node.stretch;
      }

      Eigen::Vector3d translation() const {
        return basis_node.translation;
      }

      Displacement disp(Index i) {
        return displacement.col(i);
      }

      ConstDisplacement disp(Index i) const {
        return displacement.col(i);
      }


      LatticeNode lat_node;
      HungarianNode basis_node;
      double strain_weight;
      double basis_weight;

      /// \brief true if assignment problem is not yet known to be insoluable -- default true
      bool is_viable;

      /// \brief true if assignment has been checked for physical validity and passes -- default false
      mutable bool is_valid;

      /// \brief true if node has been partitioned into sub-nodes for generalized k-best assignment problem -- default false
      mutable bool is_partitioned;

      double cost;

      Eigen::MatrixXd displacement;

      /// permutation lists indices of sites in input structure, as-read, so that
      /// they constitute particular mapping onto parent structure
      std::vector<Index> permutation;

      bool operator<(MappingNode const &other) const {
        if(!almost_equal(cost, other.cost)) {
          return cost < other.cost;
        }
        if(!identical(lat_node, other.lat_node)) {
          return lat_node < other.lat_node;
        }
        if(!identical(basis_node, other.basis_node)) {
          return basis_node < other.basis_node;
        }
        return false;
      }
    };




    /// A class for mapping an arbitrary crystal structure as a configuration of a crystal template
    /// as described by a PrimClex.  StrucMapper manages options for the mapping algorithm and mapping cost function
    /// It also caches some information about supercell lattices so that batch imports are more efficient
    ///
    /// \ingroup Configuration
    class StrucMapper {
    public:
      using LatMapType = std::map<Index, std::vector<Lattice> >;

      enum Options {none = 0,
                    strict = (1u << 0),
                    robust = (1u << 1),
                    sym_strain = (1u << 2),
                    sym_basis = (1u << 3),
                   };

      ///\brief Construct and initialize a StrucMapper
      ///\param _pclex the PrimClex that describes the crystal template
      ///
      ///\param _strain_weight
      ///\parblock
      ///          free parameter 'w' in the cost function: total_cost = w*lattice_deformation+(1-w)*basis_deformation
      ///          can vary between 0 (completely basis-focused) and 1 (completely lattice-focused)
      ///\endparblock
      ///
      ///\param _max_volume_change
      ///\parblock
      ///          constrains the search space by assuming a limit on allowed volume change
      ///          only taken into account when non-interstitial vacancies are allowed
      ///\endparblock
      ///
      ///\param _options
      ///\parblock
      ///          specify a combination of StrucMapper::Options using bitwise OR: Ex. _options=StrucMapper::robust|StrucMapper::strict
      ///          Options are:
      ///             'robust': does not assume the imported structure might be ideal ('robust' is much slower for importing ideal structures,
      ///                       but if 'robust' is not set and a non-ideal structure is passed, this will be almost always be detected and
      ///                       robust methods will be used instead. Thus, 'robust' is slightly faster if imported Structures are *not* ideal)
      ///
      ///             'strict': prevents transformation into canonical form. Tries to preserve original orientation of imported structure if possible
      ///\endparblock
      ///
      ///\param _tol tolerance for mapping comparisons
      StrucMapper(StrucMapCalculatorInterface const &_calculator,
                  double _strain_weight = 0.5,
                  Index _Nbest = 1,
                  double _max_volume_change = 0.5,
                  int _options = robust, // this should actually be a bitwise-OR of StrucMapper::Options
                  double _tol = TOL,
                  double _min_va_frac = 0.,
                  double _max_va_frac = 1.);

      double tol() const {
        return m_tol;
      }

      double strain_weight() const {
        return m_strain_weight;
      }

      double basis_weight() const {
        return 1. - m_strain_weight;
      }

      void set_strain_weight(double _lw) {
        m_strain_weight = max(min(_lw, 1.0), 1e-9);
      }

      /// \brief Returns the minimum fraction of sites allowed to be vacant in the mapping relation
      /// Vacancy fraction is used to constrain the mapping supercell search, but is only used
      /// when the supercell volume cannot is not uniquely determined by the number of each
      /// species in the child structure
      double min_va_frac() const {
        return m_min_va_frac;
      }

      /// \brief Sets the minimum fraction of sites allowed to be vacant in the mapping relation
      /// Vacancy fraction is used to constrain the mapping supercell search, but is only used
      /// when the supercell volume cannot is not uniquely determined by the number of each
      /// species in the child structure
      void set_min_va_frac(double _min_va) {
        m_min_va_frac = max(_min_va, 0.);
      }

      /// \brief Returns the maximum fraction of sites allowed to be vacant in the mapping relation
      /// Vacancy fraction is used to constrain the mapping supercell search, but is only used
      /// when the supercell volume cannot is not uniquely determined by the number of each
      /// species in the child structure
      double max_va_frac() const {
        return m_max_va_frac;
      }

      /// \brief Sets the maximum fraction of sites allowed to be vacant in the mapping relation
      /// Vacancy fraction is used to constrain the mapping supercell search, but is only used
      /// when the supercell volume cannot is not uniquely determined by the number of each
      /// species in the child structure
      void set_max_va_frac(double _max_va) {
        m_max_va_frac = min(_max_va, 1.);
      }

      /// \brief returns bit flag of selected options for this StrucMapper
      int options() const {
        return m_options;
      }

      /// \brief returns reference to parent structure
      SimpleStructure const &parent() const;

      ///\brief specify which lattices should be searched when mapping configurations
      void add_allowed_lattice(Lattice const &_lat) {
        m_allowed_superlat_map[std::abs(round(volume(_lat) / parent().lat_column_mat.determinant()))].push_back(_lat);
      }

      ///\brief unset the enforcement of particular lattices (default behavior)
      void clear_allowed_lattices() const {
        m_allowed_superlat_map.clear();
      }

      ///\brief specify to use restricted hermites when mapping
      void restricted() {
        m_restricted = true;
      }

      ///\brief returns true if lattices were set to be allowed as candidates
      bool lattices_constrained() const {
        return m_allowed_superlat_map.size();
      }

      std::set<MappingNode> seed_from_vol_range(SimpleStructure const &child_struc,
                                                Index min_vol,
                                                Index max_vol) const;

      ///\brief Low-level routine to map a structure onto a ConfigDof if it is known to be ideal
      ///\param mapped_config[out] MappingNode that is result of mapping procedure
      ///\parblock
      ///                   populated by the permutation of sites in the imported structure
      ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
      ///\endparblock
      MappingNode map_ideal_struc(const SimpleStructure &child_struc) const;



      ///\brief Low-level routine to map a structure onto a ConfigDof. Does not assume structure is ideal
      ///\parblock
      ///                   populated by the permutation of sites in the imported structure
      ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
      ///\endparblock
      ///\param best_cost[in] optional parameter. Method will return false of no mapping is found better than 'best_cost'
      std::set<MappingNode> map_deformed_struc(const SimpleStructure &child_struc,
                                               double best_cost = StrucMapping::big_inf(),
                                               bool keep_invalid = false) const;

      ///\brief Low-level routine to map a structure onto a ConfigDof assuming a specific Lattice, without assuming structure is ideal
      ///       Will only identify mappings better than best_cost, and best_cost is updated to reflect cost of best mapping identified
      ///\param imposed_lat[in] Supercell Lattice onto which struc will be mapped
      ///\param best_cost Imposes an upper bound on cost of any mapping considered, and is updated to reflect best mapping encountered
      ///\parblock
      ///                   populated by the permutation of sites in the imported structure
      ///                   that maps them onto sites of the ideal crystal (excluding vacancies)
      ///\endparblock
      std::set<MappingNode> map_deformed_struc_impose_lattice(const SimpleStructure &child_struc,
                                                              const Lattice &imposed_lat,
                                                              double best_cost = StrucMapping::big_inf(),
                                                              bool keep_invalid = false) const;

      std::set<MappingNode> map_deformed_struc_impose_lattice_node(const SimpleStructure &child_struc,
                                                                   const LatticeNode &imposed_node,
                                                                   double best_cost = StrucMapping::big_inf(),
                                                                   bool keep_invalid = false) const;



      Index k_best_maps_better_than(SimpleStructure const &child_struc,
                                    std::set<MappingNode> &queue,
                                    double max_cost,
                                    bool keep_invalid,
                                    bool erase_tail = true) const;


      StrucMapCalculatorInterface const &calculator() const {
        return *m_calc_ptr;
      }

      std::set<MappingNode> seed_k_best_from_super_lats(SimpleStructure const &child_struc,
                                                        std::vector<Lattice> const &_parent_scels,
                                                        std::vector<Lattice> const &_child_scels,
                                                        Index k,
                                                        double min_cost = 1e-6,
                                                        double max_cost = StrucMapping::small_inf()) const;
    private:
      //Implement filter as std::function<bool(Lattice const&)> or as polymorphic type or as query (i.e., DataFormatter)
      bool _filter_lat(Lattice const &_lat)const {
        return true;
      }

      std::pair<Index, Index> _vol_range(const SimpleStructure &child_struc) const;

      notstd::cloneable_ptr<StrucMapCalculatorInterface> m_calc_ptr;

      Eigen::MatrixXd m_strain_gram_mat;

      double m_strain_weight;
      Index m_Nbest;
      double m_max_volume_change;
      int m_options;
      bool m_restricted = false;
      double m_tol;
      double m_min_va_frac;
      double m_max_va_frac;

      ///Maps the supercell volume to a vector of Lattices with that volume
      mutable LatMapType m_superlat_map;
      mutable LatMapType m_allowed_superlat_map;

      std::vector<Lattice> _lattices_of_vol(Index prim_vol) const;


    };



    namespace StrucMapping_impl {

      std::set<MappingNode> k_best_lat_maps_better_than(Lattice const &_parent_prim,
                                                        std::vector<Lattice> const &_parent_scels,
                                                        Lattice const &_child_prim,
                                                        std::vector<Lattice> const &_child_scels,
                                                        Index _num_atoms,
                                                        std::vector<SymOp> const &_parent_pg,
                                                        Index k,
                                                        double strain_weight,
                                                        double basis_weight,
                                                        double _tol,
                                                        double min_cost = 1e-6,
                                                        double max_cost = StrucMapping::small_inf());

    }


  }
}
#endif
