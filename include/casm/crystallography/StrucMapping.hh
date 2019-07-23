#ifndef CASM_StrucMapping
#define CASM_StrucMapping

#include <vector>
#include <unordered_set>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/SupercellEnumerator.hh"

namespace CASM {
  class Lattice;
  class LatticeMap;
  class SimpleStructure;

  // In this file:
  class LatticeNode;
  class HungarianNode;
  class MappingNode;
  class StrucMapCalculatorInterface;
  class SimpleStrucMapCalculator;
  class StrucMapper;

  /// \brief Class describing the lattice-mapping portion of a particular mapping
  /// A general map for child_struc onto parent_struc may require forming a supercell of
  /// parent_struc (most commonly) and/or of child_struc.
  /// As such, the LatticeNode is specified in terms of superlattices of both parent_struc and child_struc,
  /// as well as deformation and rotation information sufficient to fully define the lattice map
  struct LatticeNode {

    /// \brief Right stretch tensor that takes child superlattice to its deformed state
    Eigen::Matrix3d rstretch;

    /// \brief cartesian rotation/reflection/rotoreflection that rotates the deformed child superlattice to
    Eigen::Matrix3d isometry;

    /// \brief PrimGrid for supercell of parent structure
    /// The parent structure defines the ideal strain state, and the child structure must undergo an
    /// 'idealization' deformation to perfectly tile the parent structure supercell.
    ///   Define:
    ///      - 'Sp(i)' as the ideal parent supercell
    ///      - 'Sc(i)' as the ideal child supercell
    ///      - 'Sc(d)' as the deformed child supercell
    ///      - 'R' as cartesian rotation/reflection/rotoreflection ('isometry')
    ///      - 'U' as right stretch tensor
    ///   Then, the ideal and deformed superlattices satisfy
    ///      Sp(i) == Sc(i)
    ///   and
    ///      Sc(d) == R * U * Sc(i)
    PrimGrid parent;

    /// \brief PrimGrid for supercell of child structure
    /// LatticeNode is constructed using the 'idealized' child superlattice. For a particular superlattice
    /// of the child structure, the superlattice is rotated and deformed to be exactly equal to the
    /// superlattice of the parent structure.
    ///   Define:
    ///      - 'Sp(i)' as the ideal parent supercell
    ///      - 'Sc(i)' as the ideal child supercell
    ///      - 'Sc(d)' as the deformed child supercell
    ///      - 'R' as cartesian rotation/reflection/rotoreflection ('isometry')
    ///      - 'U' as right stretch tensor
    ///   Then, the ideal and deformed superlattices satisfy
    ///      Sp(i) == Sc(i)
    ///   and
    ///      Sc(d) == R * U * Sc(i)
    PrimGrid child;

    /// \brief strain_cost of the LatticeNode
    double cost;

    LatticeNode();

    /// \brief Construct with ideal parent_scel and deformed child_scel, which are related by a deformation tensor
    /// @param parent_scel and @param child_scel are integer combinations of the primitive cells 'parent_prim' and 'child_prim', respectively
    /// @param child_N_atom is number of sites in the child
    LatticeNode(Lattice const &parent_prim,
                Lattice const &parent_scel,
                Lattice const &child_prim,
                Lattice const &child_scel,
                Index child_N_atom);

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

    Eigen::Vector3d translation;

    std::set<std::pair<Index, Index> > forced_on;
    //std::set<std::pair<Index,Index> > forced_off;
    std::vector<Index> irow;
    std::vector<Index> icol;
    std::vector<Index> assignment;
    Eigen::MatrixXd cost_mat;
    double cost;
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

    MappingNode() :
      is_viable(true),
      is_valid(false),
      is_partitioned(false) {
    }

    MappingNode(LatticeNode _lat_node, double _lat_weight, double _basis_weight) :
      lat_node(std::move(_lat_node)),
      lat_weight(_lat_weight),
      basis_weight(_basis_weight),
      is_viable(true),
      is_valid(false),
      is_partitioned(false),
      cost(lat_weight * lat_node.cost) {

    }

    double tol() const {
      return basis_node.tol();
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

    Displacement disp(Index i) {
      return displacement.col(i);
    }

    ConstDisplacement disp(Index i) const {
      return displacement.col(i);
    }


    LatticeNode lat_node;
    HungarianNode basis_node;
    double lat_weight;
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

  namespace StrucMapping {

    // Denotes (name, number composition) of any species whose number-composition is fixed for the parent primitive cell
    using FixedSpecies = std::map<std::string, Index>;

    // List of species allowed at each site of primitive
    using AllowedSpecies = std::vector<std::unordered_set<std::string>>;

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

  class StrucMapCalculatorInterface {
  public:

    StrucMapCalculatorInterface(SimpleStructure _parent,
                                std::vector<SymOp> _point_group = {SymOp()},
                                SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                StrucMapping::AllowedSpecies allowed_species = {}) :
      m_parent(std::move(_parent)),
      m_point_group(std::move(_point_group)),
      m_species_mode(_species_mode),
      m_allowed_species(std::move(allowed_species)) {

      if(m_allowed_species.empty()) {
        auto const &p_info(info(parent()));
        m_allowed_species.resize(p_info.size());
        for(Index i = 0; i < m_allowed_species.size(); ++i)
          m_allowed_species[i].insert(p_info.names[i]);
      }


      for(auto const &slist : m_allowed_species) {
        if(slist.size() != 1) {
          for(auto const &sp : slist) {
            _fixed_species()[sp] = 0;
          }
        }
        auto it = _fixed_species().find(*slist.begin());
        if(it == _fixed_species().end())
          _fixed_species()[*slist.begin()] = 1;
        else if((it->second) > 0)
          ++(it->second);
      }
      for(auto it = _fixed_species().begin(); it != _fixed_species().end(); ++it) {
        if((it->second) == 0)
          _fixed_species().erase(it);
      }
    }

    virtual ~StrucMapCalculatorInterface() {}

    SimpleStructure::Info const &info(SimpleStructure const &_struc) const {
      return _struc.info(m_species_mode);
    }

    SimpleStructure const &parent() const {
      return m_parent;
    }

    std::vector<SymOp> const &point_group() const {
      return m_point_group;
    }

    void set_point_group(std::vector<SymOp> _point_group) {
      m_point_group = std::move(_point_group);
    }

    StrucMapping::FixedSpecies const &fixed_species() const {
      return m_fixed_species;
    }

    virtual Index max_n_va() const = 0;

    virtual std::vector<Eigen::Vector3d> translations(SimpleStructure const &child_struc,
                                                      LatticeNode const &lat_node) const = 0;


    virtual void finalize(MappingNode &_node,
                          SimpleStructure const &child_struc) const = 0;

    virtual bool populate_cost_mat(MappingNode &_node,
                                   SimpleStructure const &child_struc) const = 0;

    virtual void populate_properties(MappingNode &_node,
                                     const SimpleStructure &child_struc) const = 0;

    virtual bool validate(MappingNode const &_node) const = 0;

    /// \brief Make an exact copy of the calculator (including any initialized members)
    std::unique_ptr<StrucMapCalculatorInterface> clone() const {
      return std::unique_ptr<StrucMapCalculatorInterface>(this->_clone());
    }

    /// \brief Make an exact copy of the calculator (including any initialized members)
    std::unique_ptr<StrucMapCalculatorInterface> quasi_clone(SimpleStructure _parent,
                                                             std::vector<SymOp> _point_group = {SymOp()},
                                                             SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                             StrucMapping::AllowedSpecies _allowed_species = {}) const {
      return std::unique_ptr<StrucMapCalculatorInterface>(this->_quasi_clone(std::move(_parent),
                                                                             std::move(_point_group),
                                                                             _species_mode,
                                                                             std::move(_allowed_species)));
    }

  protected:
    StrucMapping::AllowedSpecies &_allowed_species() {
      return m_allowed_species;
    }

    StrucMapping::AllowedSpecies const &_allowed_species() const {
      return m_allowed_species;
    }

    StrucMapping::FixedSpecies &_fixed_species() {
      return m_fixed_species;
    }

  private:
    SimpleStructure m_parent;

    std::vector<SymOp> m_point_group;

    SimpleStructure::SpeciesMode m_species_mode;

    StrucMapping::AllowedSpecies m_allowed_species;

    StrucMapping::FixedSpecies m_fixed_species;

    /// \brief Make an exact copy of the calculator (including any initialized members)
    virtual StrucMapCalculatorInterface *_clone() const = 0;

    /// \brief Make an exact copy of the calculator (including any initialized members)
    virtual StrucMapCalculatorInterface *_quasi_clone(SimpleStructure _parent,
                                                      std::vector<SymOp> _point_group = {SymOp()},
                                                      SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                      StrucMapping::AllowedSpecies _allowed_species = {}) const = 0;

  };

  class SimpleStrucMapCalculator : public StrucMapCalculatorInterface {
  public:
    SimpleStrucMapCalculator(SimpleStructure _parent,
                             std::vector<SymOp> _point_group = {SymOp()},
                             SimpleStructure::SpeciesMode species_mode = SimpleStructure::SpeciesMode::ATOM,
                             StrucMapping::AllowedSpecies allowed_species = {}) :
      StrucMapCalculatorInterface(std::move(_parent),
                                  std::move(_point_group),
                                  species_mode,
                                  std::move(allowed_species)) {

      for(Index i = 0; i < _allowed_species().size(); ++i) {
        if(_allowed_species()[i].count("Va") || _allowed_species()[i].count("VA") || _allowed_species()[i].count("va"))
          m_va_allowed.insert(i);
      }
    }

    Index max_n_va() const override {
      return m_va_allowed.size();
    }

    std::vector<Eigen::Vector3d> translations(SimpleStructure const &child_struc,
                                              LatticeNode const &lat_node) const override;
    void finalize(MappingNode &_node,
                  SimpleStructure const &child_struc) const override ;//{TODO();};


    bool populate_cost_mat(MappingNode &_node,
                           SimpleStructure const &child_struc) const override;

    void populate_displacement(MappingNode &_node,
                               SimpleStructure const &child_struc) const;

    void populate_properties(MappingNode &_node,
                             const SimpleStructure &child_struc) const override;

    bool validate(MappingNode const &_node) const override {
      return true;
    }

  private:
    /// \brief Make an exact copy of the calculator (including any initialized members)
    virtual StrucMapCalculatorInterface *_clone() const override {
      return new SimpleStrucMapCalculator(*this);
    }

    /// \brief Make an exact copy of the calculator (including any initialized members)
    virtual StrucMapCalculatorInterface *_quasi_clone(SimpleStructure _parent,
                                                      std::vector<SymOp> _point_group = {SymOp()},
                                                      SimpleStructure::SpeciesMode _species_mode = SimpleStructure::SpeciesMode::ATOM,
                                                      StrucMapping::AllowedSpecies _allowed_species = {}) const override {
      return new SimpleStrucMapCalculator(std::move(_parent), std::move(_point_group), _species_mode, std::move(_allowed_species));
    }

    std::unordered_set<Index> m_va_allowed;
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
                  robust = (1u << 1)
                 };

    ///\brief Construct and initialize a StrucMapper
    ///\param _pclex the PrimClex that describes the crystal template
    ///
    ///\param _lattice_weight
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
                double _lattice_weight = 0.5,
                Index _Nbest = 1,
                double _max_volume_change = 0.5,
                int _options = robust, // this should actually be a bitwise-OR of StrucMapper::Options
                double _tol = TOL,
                double _min_va_frac = 0.,
                double _max_va_frac = 1.);

    double tol() const {
      return m_tol;
    }

    double lattice_weight() const {
      return m_lattice_weight;
    }

    double basis_weight() const {
      return 1. - m_lattice_weight;
    }

    void set_lattice_weight(double _lw) {
      m_lattice_weight = max(min(_lw, 1.0), 1e-9);
    }

    double min_va_frac() const {
      return m_min_va_frac;
    }

    void set_min_va_frac(double _min_va) {
      m_min_va_frac = max(_min_va, 0.);
    }

    double max_va_frac() const {
      return m_max_va_frac;
    }

    void set_max_va_frac(double _max_va) {
      m_max_va_frac = min(_max_va, 1.);
    }

    int options() const {
      return m_options;
    }

    SimpleStructure const &parent() const {
      return _calculator().parent();
    }

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

    template<typename OutputIterator>
    bool insert_at_most_k_maps(SimpleStructure child_struc,
                               MappingNode const &seed,
                               double max_cost,
                               Index k,
                               OutputIterator it) const {

      //derotate first
      child_struc.rotate(seed.lat_node.isometry.transpose());

      //Then undeform by inverse of right stretch
      child_struc.deform(seed.lat_node.rstretch.inverse());

      // We want to get rid of translations.
      // define translation such that:
      //    IDEAL = RELAXED + translation
      // and use it when calculating cost matrix

      for(Eigen::Vector3d const &translation : _calculator().translations(child_struc,
                                                                          seed.lat_node)) {
        MappingNode node = seed;
        node.basis_node.translation = translation;
        if(!_calculator().populate_cost_mat(node,
                                            child_struc)) {
          // Indicates that structure is incompatible with supercell, regardless of translation so return false
          return false;
        }

        // The mapping routine is called here
        node.calc();

        // if assignment is smaller than child_struc.basis().size(), then child_struc is incompattible with supercell
        // (assignment.size()==0 if the hungarian routine detects an incompatibility, regardless of translation)
        if(!node.is_viable) {
          return false;
        }

        // Now we are filling up displacements
        _calculator().finalize(node, child_struc);

        // add small penalty (~_tol) for larger translation distances, so that shortest equivalent translation is used
        //node.cost = lattice_weight() * node.lat_node.cost + basis_weight() * (StrucMapping::basis_cost(node, c_info.size() * node.lat_node.child.size()) + m_tol * node.basis_node.translation.norm() / 10.0);
        if(node.cost < max_cost) {
          *it = node;
        }

      }
      return true;
    }

    StrucMapCalculatorInterface const &calculator() const {
      return *m_calc_ptr;
    }

  private:
    //Implement filter as std::function<bool(Lattice const&)> or as polymorphic type or as query (i.e., DataFormatter)
    bool _filter_lat(Lattice const &_lat)const {
      return true;
    }

    std::pair<Index, Index> _vol_range(const SimpleStructure &child_struc) const;

    StrucMapCalculatorInterface const &_calculator() const {
      return *m_calc_ptr;
    }


    notstd::cloneable_ptr<StrucMapCalculatorInterface> m_calc_ptr;

    double m_lattice_weight;
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

  template<typename OutputIterator>
  void partition_node(MappingNode const &_node,
                      StrucMapCalculatorInterface const &_calculator,
                      SimpleStructure const &_child_struc,
                      OutputIterator it) {
    Index j, jj, currj, tj;
    Index n = _node.basis_node.assignment.size();
    MappingNode t1(_node),
                t2(_node.lat_node, _node.lat_weight, _node.basis_weight);

    MappingNode *p1 = &t1;
    MappingNode *p2 = &t2;
    for(Index m = 0; m < (n - 1); ++m) {
      HungarianNode &n1(p1->basis_node);
      HungarianNode &n2(p2->basis_node);
      //clear assignment and cost_mat
      n2.assignment.clear();
      n2.cost_mat.resize(n - m - 1, n - m - 1);

      // We are forcing on the first site assignment in t1: i.e.,  [0,currj]
      // this involves striking row 0 and col currj from t2's cost_mat
      // and augmenting the cost_offset of t2 by the cost of [0,currj]
      currj = n1.assignment[0];
      n2.cost_offset = n1.cost_offset + n1.cost_mat(0, currj);

      // [0,currj] have local context only. We store the forced assignment in forced_on
      // using the original indexing
      n2.forced_on.emplace(n2.irow[0], n2.icol[currj]);

      // Strike row 0 and col currj to form new HungarianNode for t2
      // (i,j) indexes the starting cost_mat, (i-1, jj) indexes the resulting cost_mat
      n2.irow = std::vector<Index>(++n1.irow.begin(), n1.irow.end());
      for(j = 0, jj = 0; j < (n - m); ++j, ++jj) {
        if(j == currj) {
          --jj;
          continue;
        }
        n2.icol.push_back(n1.icol[j]);

        // We will also store an updated assignment vector in t2, which will be
        // used to construct next node of partition
        tj = n1.assignment[j];
        if(tj > currj)
          tj--;
        n2.assignment.push_back(tj);

        // Fill col jj of t2's cost mat
        for(Index i = 1; i < n; ++i)
          n2.cost_mat(i - 1, jj) = n1.cost_mat(i, j);
      }
      //t2 properly initialized; we can now force OFF [0,currj] in t1, and add it to node list
      n1.cost_mat(0, currj) = StrucMapping::big_inf();
      n1.assignment.clear();
      p1->calc();
      _calculator.finalize(*p1, _child_struc);
      it = *p1;
      std::swap(p1, p2);
    }
  }


  namespace StrucMapping_impl {

    std::set<MappingNode> k_best_nodes_better_than(Lattice const &_parent_prim,
                                                   std::vector<Lattice> const &_parent_scels,
                                                   Lattice const &_child_prim,
                                                   std::vector<Lattice> const &_child_scels,
                                                   Index _num_atoms,
                                                   std::vector<SymOp> const &_parent_pg,
                                                   Index k,
                                                   double lattice_weight,
                                                   double basis_weight,
                                                   double _tol,
                                                   double min_cost = 1e-6,
                                                   double max_cost = StrucMapping::small_inf());

  }


}

#endif
