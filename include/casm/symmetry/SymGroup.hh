#ifndef SYMGROUP_HH
#define SYMGROUP_HH

#include <iostream>
#include <string>
#include <iomanip>

#include "casm/symmetry/SymOp.hh"

namespace CASM {

  class Lattice;

  class SymGroup;
  class MasterSymGroup;
  class SymGroupRep;
  struct SymInfo;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /** \defgroup SymGroup
   *  \ingroup Symmetry
   *  \brief Relates to symmetry groups
   *  @{
   */

  ///\brief SymGroup is a collection of symmetry operations that satisfy the group property
  /// The symmetry operations are stored as their coordinate representation, as described by the SymOp class
  /// i.e., if SymOps 'A' and 'B' are in SymGroup, C=A*B is also in SymGroup
  ///       if 'A' is in SymGroup, then A.inverse() is in SymGroup
  ///       SymGroup always contains an identity operation
  class SymGroup : public Array<SymOp> {
  public:
    typedef SymOp::vector_type vector_type;
    typedef SymOp::matrix_type matrix_type;
    /// Initialize by setting periodicity mode (default mode is PERIODIC)
    SymGroup(PERIODICITY_TYPE init_type = PERIODIC) :
      m_lat_ptr(nullptr),
      m_group_periodicity(init_type),
      m_max_error(-1) {
      name.clear();
      latex_name.clear();
    }

    SymGroup(const Array<SymOp> &from_array, PERIODICITY_TYPE init_type = PERIODIC);

    virtual void push_back(const SymOp &new_op);
    virtual void clear();
    virtual void clear_tables();

    void set_lattice(const Lattice &new_lat);

    const Lattice &lattice() const;

    const MasterSymGroup &master_group() const {
      assert(size() && at(0).has_valid_master());
      return at(0).master_group();
    }

    ReturnArray<Index> op_indices() const;

    ///Check to see if a SymOp is contained in in SymGroup
    //maybe contains and find should account for group_periodicity
    //bool contains(const SymOp &test_op) const; //Donghee

    bool contains_periodic(const SymOp &test_op, double tol = TOL) const;


    ///Check to see if a SymOp is contained in in SymGroup and return its index
    //Index find(const SymOp &test_op) const;   //Donghee

    /// Check to see if a SymOp matrix ONLY is contained in SymGroup and return the index of this operation.
    /// This was originally written for pruning the factor groups of primitive structures to construct the
    /// factor groups of their superstructures to be consistent with the supercell lattice point groups.
    Index find_no_trans(const SymOp &test_op) const;

    /// This is meant for factor groups. It will compare the Cartesian matrix of the test_op with those
    /// of the SymOps in the group. Upon a successful matrix match, it will attempt to match the shift
    /// shift vector with min_dist.
    Index find_periodic(const SymOp &test_op, double tol = TOL) const;
    ReturnArray<Index> find_all_periodic(const Array<SymOp> &subgroup, double tol = TOL) const;

    /// \brief Sort SymOp in the SymGroup
    virtual void sort();

    /// Adds SymOps from 'other_group' and enforces the group property
    // Should we change the name to something else?
    SymGroup get_union(const SymGroup &other_group) const;  //Ivy do this

    /// Calls 'apply_sym' on all SymOps in the group
    SymGroup &apply_sym(const SymOp &op);

    /// Get index of operation that is result of multiplication of at(i)*at(j)
    Index ind_prod(Index i, Index j) const;

    /// Get index of operation that is inverse of operation at(i)
    Index ind_inverse(Index i) const;

    /// Get conjugacy class index of operation at(i)
    Index class_of_op(Index i) const;

    /// set symrep ID of a particular irrep
    void set_irrep_ID(Index i, SymGroupRepID ID) const;

    /// Get symrep ID of a particular irrep
    SymGroupRepID get_irrep_ID(Index i) const;

    /// Get symrep ID of the representation that stores the Cartesian symop matrices
    SymGroupRepID coord_rep_ID() const;

    /// Get symrep for a particular irrep
    SymGroupRep const &get_irrep(Index i) const;

    /// Add a new empty representation
    SymGroupRepID add_empty_representation() const;


    /// Gets all the space group operations in unit cell and stores them in space_group
    /// assuming that this SymGroup contains the factor group
    void calc_space_group_in_cell(SymGroup &space_group, const Lattice &_cell) const;

    /// gets all teh space group operations corresponding to translations in the specified range
    /// max_trans sets boundary of parillellipiped centered at origin.
    void calc_space_group_in_range(SymGroup &space_group, const Lattice &_cell, Eigen::Vector3i min_trans,  Eigen::Vector3i max_trans) const;

    /// Check to see if SymGroup satisfies the group property
    bool is_group(double tol = TOL) const;
    /// Enforce group property by adding products of operations to the group
    void enforce_group(double tol = TOL, Index max_size = 200);  //AAB

    /// print locations of the symmetry-generating element of each SymOp
    void print_locations(std::ostream &stream) const;

    /// Write the SymGroup to a file
    void write(std::string filename, COORD_TYPE mode) const;

    /// Print the SymGroup to a stream
    void print(std::ostream &out, COORD_TYPE mode) const;

    /// Cartesian translation of SymGroup origin by vector 'shift'
    SymGroup &operator+=(const Eigen::Ref<const SymOp::vector_type> &shift);
    SymGroup &operator-=(const Eigen::Ref<const SymOp::vector_type> &shift);

    /// This returns the group's max_error
    double max_error();

    ReturnArray<Array<Index> > left_cosets(const Array<SymOp> &subgroup, double tol = TOL) const;
    ReturnArray<Array<Index> > left_cosets(const Array<Index> &subgroup_inds) const;

    const Array<Index>::X2 &get_multi_table() const;
    const Array<Index>::X2 &get_alt_multi_table() const;
    void invalidate_multi_tables() const;
    const Array<Index>::X2 &get_conjugacy_classes() const;
    const Array<std::complex<double> >::X2 &character_table() const;
    const Array<bool> &get_complex_irrep_list() const;
    const std::string &get_name() const;
    const std::string &get_latex_name() const;

    PERIODICITY_TYPE periodicity() const {
      return m_group_periodicity;
    }

    std::string possible_space_groups() const {
      return comment;
    }


    const Array<Index>::X3 &subgroups() const;

    void print_character_table(std::ostream &stream);
    ReturnArray<Index> get_irrep_decomposition() const;
    bool is_irreducible() const;

    std::vector<SymGroup> unique_subgroups() const;

    ///Space group (added by Donghee );
    void get_rotation_groups()const;
    void get_point_group_type()const;
    void print_space_group_info(std::ostream &out) const;

    ///Fill up a SymGroup with *this minus the shifts
    void copy_no_trans(SymGroup &shiftless, bool keep_repeated = false) const;

    jsonParser &to_json(jsonParser &json) const;

    void from_json(const jsonParser &json);

    SymInfo info(Index i) const;

  protected:
    void _generate_conjugacy_classes() const;
    void _generate_character_table() const;
    void _generate_centralizers() const;
    void _generate_elem_order_table() const;
    void _generate_class_names() const;
    void _generate_irrep_names() const;
    bool _generate_multi_table() const;
    void _generate_alt_multi_table() const;
    void _generate_subgroups() const;

    // small_groups are cyclic subgroups.  Found by taking a group
    // element and multiplying it by itself until a group is generated
    // small_groups[i][j][k] is index of symop 'k' in subgroup (i,j) -- the equivalent subroup 'j' of an orbit 'i' of equivalent subgroups
    Array<Index>::X3 _small_subgroups() const;


    /// Pointer to a lattice for doing periodic comparisons
    Lattice const *m_lat_ptr;

    /// Specifies whether to use lattice periodicity when testing for equivalence
    PERIODICITY_TYPE m_group_periodicity;

    /// multi_table[i][j] gives index of operation that is result of at(i)*at(j)
    mutable Array<Array<Index> > multi_table;

    /// alt_multi_table[i][j] gives index of operation that is result of at(i).inverse()*at(j)
    mutable Array<Array<Index> > alt_multi_table;

    // information about conjugacy classes
    // conjugacy_classes[i][j] gives index of SymOp 'j' in class 'i'
    mutable Array<Array<Index> > conjugacy_classes;
    mutable Array<std::string> class_names;
    mutable Array<Index> index2conjugacy_class;

    // Information about irreducible representations
    // m_character_table[i][j] is character of conjugacy class 'j' in irrep 'i'
    mutable Array<Array<std::complex<double> > > m_character_table;
    mutable Array<SymGroupRepID> irrep_IDs;
    mutable Array<bool> complex_irrep;
    mutable Array<std::string> irrep_names;

    // subgroups are found by finding the closure for each possible union of small_subgroups
    // organized the same way as small_subgroups
    mutable Array<Array<Array<Index> > > m_subgroups;


    mutable Array<Array<Index> > centralizer_table;
    mutable Array<Array<Index> > elem_order_table;

    mutable std::string name;
    mutable std::string latex_name;
    mutable std::string comment;

    mutable double m_max_error;

    ///Space group (added by Donghee );
    mutable Array<Array<SymOp> > rotation_groups;
    mutable std::string crystal_system;
    mutable bool centric; // if it is centric, special point is same in reciprocal space
    mutable Array<int> group_number; // space group number (min and max)
    mutable Array<std::string> group_name; // 0: International 1: Schonflies


  };

  jsonParser &to_json(const SymGroup &group, jsonParser &json);

  // Note: as a hack this expects group[0] to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the Molecules
  void from_json(SymGroup &group, const jsonParser &json);

  // return SymGroup with all molecular point group sym ops
  // I will centerize your coord_map, fyi.
  SymGroup molecular_point_group(std::map<int, std::vector<Eigen::Vector3d> > coord_map);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class MasterSymGroup : public SymGroup {
    //NOTE: It may be useful to store a user-specified set of "favored directions" in MasterSymGroup
    //      e.g., {(0,0,1), (1,0,0), (0,1,0)}, which would be used to arbitrate decisions in symmetry-related algorithms
    //            Relevant applications include
    //                - constucting orbit equivalence map
    //                - optimizing the coordinate system of irreps
    //                - deciding prototype cluster

  public:
    MasterSymGroup(PERIODICITY_TYPE init_type = PERIODIC) :
      SymGroup(init_type),
      m_group_index(GROUP_COUNT++) {
    }

    MasterSymGroup(const MasterSymGroup &RHS);
    ~MasterSymGroup();

    Index group_index() const {
      return m_group_index;
    }

    MasterSymGroup &operator=(const MasterSymGroup &RHS);

    /// push_back sets home_group and op_index of added SymOp;
    /// virtual in SymGroup, so this overrides
    void push_back(const SymOp &op);

    /// Reset everything and delete all representations
    /// virtual in SymGroup, so this overrides
    void clear();

    void sort();
    //void sort_by_class();

    const SymGroup &point_group() const;

    /// Const access of alternate Representations of a SymGroup
    SymGroupRep const &representation(SymGroupRepID i) const;

    SymGroupRep const &reg_rep() const;

    SymGroupRep const &coord_rep() const;

    /// Add a new representation by passing a reference.  SymGroup will store a copy
    SymGroupRepID add_representation(const SymGroupRep &new_rep) const;

    /// Add a new empty representation
    SymGroupRepID add_empty_representation() const;

    SymGroupRepID reg_rep_ID() const;
    SymGroupRepID coord_rep_ID() const;
    SymGroupRepID identity_rep_ID(Index dim) const;
    SymGroupRepID add_kronecker_rep(SymGroupRepID ID1, SymGroupRepID ID2) const;
    SymGroupRepID add_direct_sum_rep(const Array<SymGroupRepID> &rep_IDs) const;
    SymGroupRepID add_transformed_rep(SymGroupRepID orig_ID, const Eigen::MatrixXd &trans_mat) const;
    SymGroupRepID add_rotation_rep() const;

    jsonParser &to_json(jsonParser &json) const;

    void from_json(const jsonParser &json);

  private:

    SymGroupRep *_representation_ptr(SymGroupRepID _id) const;

    SymGroupRepID _add_reg_rep() const;
    SymGroupRepID _add_coord_rep() const;

    /// Add a new representation by passing a pointer to SymGroupRep allocated via 'new'.
    /// MasterSymGroup will store the pointer, so don't delete it after calling
    SymGroupRepID _add_representation(SymGroupRep *_rep_ptr) const;

    /// Counts number of instantiated MasterSymGroups, excluding ones created via copy
    static Index GROUP_COUNT;

    /// Index of this group, initialized from value of GROUP_COUNT upon construction
    Index m_group_index;

    /// Collection of alternate representations of this symmetry group
    /// Stored as pointers to avoid weird behavior with resizing
    mutable Array<SymGroupRep *> m_rep_array;

    /// ID of Cartesian representation
    mutable SymGroupRepID m_coord_rep_ID;

    /// ID of 'regular representation', which is (size() X size()) representation constructed from alt_multi_table()
    mutable SymGroupRepID m_reg_rep_ID;

    /// identity representations: m_identity_rep_IDs[dim] refers to the Identity representation of dimention 'dim'
    mutable Array<SymGroupRepID> m_identity_rep_IDs;

    /// Copy of *this with translations removed
    mutable SymGroup m_point_group;
  };

  jsonParser &to_json(const SymGroup &group, jsonParser &json);

  // Note: as a hack this expects group[0] to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the Molecules
  void from_json(SymGroup &group, const jsonParser &json);

  bool compare_periodic(const SymOp &a,
                        const SymOp &b,
                        const Lattice &lat,
                        PERIODICITY_TYPE periodicity,
                        double _tol);

  SymOp within_cell(const SymOp &a, const Lattice &lat, PERIODICITY_TYPE periodicity);

  /** @} */
}

#endif
