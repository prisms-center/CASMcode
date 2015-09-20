#ifndef SYMGROUP_HH
#define SYMGROUP_HH

#include <iostream>
#include <string>
#include <iomanip>

#include "casm/symmetry/SymOp.hh"

//#include "casm/../CASM_global_definitions.cc"
//#include "casm/../crystallography/CoordinateSystems.cc"
//#include "casm/../container/Array.cc"
//#include "casm/../symmetry/SymOp.hh"
//#include "casm/../symmetry/SymGroup.hh"
//#include "casm/../container/Counter.cc"

namespace CASM {
  class SymGroup;
  class MasterSymGroup;
  class SymGroupRep;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ///\brief SymGroup is a collection of symmetry operations that satisfy the group property
  /// The symmetry operations are stored as their coordinate representation, as described by the SymOp class
  /// i.e., if SymOps 'A' and 'B' are in SymGroup, C=A*B is also in SymGroup
  ///       if 'A' is in SymGroup, then A.inverse() is in SymGroup
  ///       SymGroup always contains an identity operation
  class SymGroup : public Array<SymOp> {
  protected:

    /// Specifies whether to use lattice periodicity when testing for equivalence
    PERIODICITY_TYPE group_periodicity;

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
    // character_table[i][j] is character of conjugacy class 'j' in irrep 'i'
    mutable Array<Array<std::complex<double> > > character_table;
    mutable Array<Index> irrep_IDs;
    mutable Array<bool> complex_irrep;
    mutable Array<std::string> irrep_names;

    Array<SymGroup> unique_subgroups;

    // small_groups are cyclic subgroups.  Found by taking a group
    // element and multiplying it by itself until a group is generated
    // small_groups[i][j][k] is index of symop 'k' in subgroup (i,j) -- the equivalent subroup 'j' of an orbit 'i' of equivalent subgroups
    mutable Array<Array<Array<Index> > > small_groups;

    // large_groups are found by finding the closure for each possible union of small_groups
    // organized the same way as small_groups
    mutable Array<Array<Array<Index> > > large_groups;


    mutable Array<Array<Index> > centralizer_table;
    mutable Array<Array<Index> > elem_order_table;

    mutable std::string name;
    mutable std::string latex_name;
    mutable std::string comment;

    mutable double max_error;

    ///Space group (added by Donghee );
    mutable Array<Array<SymOp> > rotation_groups;
    mutable std::string crystal_system;
    mutable bool centric; // if it is centric, special point is same in reciprocal space
    mutable Array<int> group_number; // space group number (min and max)
    mutable Array<std::string> group_name; // 0: International 1: Schonflies

  protected:
    void calc_conjugacy_classes() const;
    void calc_character_table() const;
    void calc_centralizers() const;
    void calc_elem_order_table() const;
    void generate_class_names() const;
    void generate_irrep_names() const;
    bool calc_multi_table() const;
    void calc_alt_multi_table() const;
    void calc_small_subgroups() const;
    void calc_large_subgroups() const;


  public:
    /// Initialize by setting periodicity mode (default mode is PERIODIC)
    SymGroup(PERIODICITY_TYPE init_type = PERIODIC) : group_periodicity(init_type), max_error(-1) {
      name.clear();
      latex_name.clear();
    };

    virtual void push_back(const SymOp &new_op);
    virtual void clear();
    virtual void clear_tables();

    ///Check to see if a SymOp is contained in in SymGroup
    //maybe contains and find should account for group_periodicity
    bool contains(const SymOp &test_op) const; //Donghee

    bool contains_periodic(const SymOp &test_op, double tol = TOL) const;


    ///Check to see if a SymOp is contained in in SymGroup and return its index
    Index find(const SymOp &test_op) const;   //Donghee

    /// Check to see if a SymOp matrix ONLY is contained in SymGroup and return the index of this operation.
    /// This was originally written for pruning the factor groups of primitive structures to construct the
    /// factor groups of their superstructures to be consistent with the supercell lattice point groups.
    Index find_no_trans(const SymOp &test_op) const;

    /// This is meant for factor groups. It will compare the Cartesian matrix of the test_op with those
    /// of the SymOps in the group. Upon a successful matrix match, it will attempt to match the shift
    /// shift vector with min_dist.
    Index find_periodic(const SymOp &test_op, double tol = TOL) const;

    /// Sort SymOps in SymGroup
    /// Positive determinant SymOps come before negative determinant
    /// Sorted in order of decreasing character (trace)
    virtual void sort();
    virtual void sort_by_class(); //AAB

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
    void set_irrep_ID(Index i, Index ID) const;

    /// Get symrep ID of a particular irrep
    Index get_irrep_ID(Index i) const;

    /// Get symrep for a particular irrep
    SymGroupRep const *get_irrep(Index i) const;

    Index make_empty_representation() const;

    void set_lattice(const Lattice &new_lat, COORD_TYPE mode);

    /// Gets all the space group operations in unit cell and stores them in space_group
    /// assuming that this SymGroup contains the factor group
    void calc_space_group_in_cell(SymGroup &space_group) const;

    /// gets all teh space group operations corresponding to translations in the specified range
    /// max_trans sets boundary of parillellipiped centered at origin.
    void calc_space_group_in_range(SymGroup &space_group, Vector3<int> min_trans, Vector3<int> max_trans) const;

    /// Check to see if SymGroup satisfies the group property
    bool is_group(double tol = TOL) const;
    /// Enforce group property by adding products of operations to the group
    void enforce_group(double tol = TOL, Index max_size = 1000);  //AAB

    /// print locations of the symmetry-generating element of each SymOp
    void print_locations(std::ostream &stream);

    /// Write the SymGroup to a file
    void write(std::string filename, COORD_TYPE mode) const;

    /// Print the SymGroup to a stream
    void print(std::ostream &out, COORD_TYPE mode) const;

    /// Cartesian translation of SymGroup origin Coordinate shift
    SymGroup &operator+=(const Coordinate &shift);
    SymGroup &operator-=(const Coordinate &shift);

    Eigen::MatrixXd const *get_MatrixXd(Index i) const;

    /// This returns the group's max_error
    double get_max_error();

    const Array<Array<Index> > &get_multi_table() const;
    const Array<Array<Index> > &get_alt_multi_table() const;
    void invalidate_multi_tables() const;
    const Array<Array<Index> > &get_conjugacy_classes() const;
    const Array<Array<std::complex<double> > > &get_character_table() const;
    const Array<bool> &get_complex_irrep_list() const;
    const std::string &get_name() const;
    const std::string &get_latex_name() const;
    PERIODICITY_TYPE get_periodicity() const {
      return group_periodicity;
    }
    std::string possible_space_groups() const {
      return comment;
    }

    const Array<Array<Array<Index> > > &get_large_subgroups() const;
    const Array<Array<Array<Index> > > &get_small_subgroups() const;

    void print_character_table(std::ostream &stream);
    ReturnArray<Index> get_irrep_decomposition() const;
    bool is_irreducible() const;

    ///Space group (added by Donghee );
    void get_rotation_groups()const;
    void get_point_group_type()const;
    void print_space_group_info(std::ostream &out) const;

    void make_unique_subgroups();

    ///Fill up a SymGroup with *this minus the shifts
    void copy_no_trans(SymGroup &shiftless, bool keep_repeated = false) const;


    jsonParser &to_json(jsonParser &json) const;

    // Note: as a hack this expects at(0) to be present and have the right lattice!!!
    //   it's just used to set the lattice for all the SymOp
    void from_json(const jsonParser &json);

  };

  jsonParser &to_json(const SymGroup &group, jsonParser &json);

  // Note: as a hack this expects group[0] to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the Molecules
  void from_json(SymGroup &group, const jsonParser &json);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class MasterSymGroup : public SymGroup {
    ///  Collection of alternate representations of this symmetry group
    /// (default is the coordinate representation)
    mutable Array<SymGroupRep *> rep_array;
    mutable Index coord_rep_ID, reg_rep_ID;

    mutable SymGroup point_group_internal;

  public:
    MasterSymGroup(PERIODICITY_TYPE init_type = PERIODIC) : SymGroup(init_type), coord_rep_ID(-1), reg_rep_ID(-1) {};
    MasterSymGroup(const MasterSymGroup &RHS);
    ~MasterSymGroup();

    MasterSymGroup &operator=(const MasterSymGroup &RHS);

    /// push_back sets home_group and op_index of added SymOp;
    /// virtual in SymGroup, so this overrides
    void push_back(const SymOp &op);

    /// Reset everything and delete all representations
    /// virtual in SymGroup, so this overrides
    void clear();

    void sort();
    void sort_by_class();

    const SymGroup &point_group() const;

    /// Const access of alternate Representations of a SymGroup
    SymGroupRep const *representation(Index i) const;

    /// Add a new representation by passing a pointer;
    /// SymGroup assumes ownership (user promises that no other object will delete new_ptr)
    Index add_representation(SymGroupRep *new_ptr) const;

    /// Add a new empty representation
    Index make_empty_representation() const;

    /// Add a new representation by passing a reference.  SymGroup will store a copy
    Index add_representation(const SymGroupRep &new_rep) const;
    Index get_reg_rep_ID() const;
    Index get_coord_rep_ID() const;
    SymGroupRep const *get_reg_rep() const;
    SymGroupRep const *get_coord_rep() const;
    Index add_reg_rep() const;
    Index add_coord_rep() const;
    Index add_kronecker_rep(Index ID1, Index ID2) const;
    Index add_direct_sum_rep(const Array<Index> &rep_IDs) const;
    //Index add_transformed_rep(Index orig_ID, const Eigen::MatrixXd &trans_mat) const;

    //Eigen::MatrixXd const* get_MatrixXd(Index i) const;

    jsonParser &to_json(jsonParser &json) const;

    // Note: as a hack this expects at(0) to be present and have the right lattice!!!
    //   it's just used to set the lattice for all the SymOp
    void from_json(const jsonParser &json);

  };

  jsonParser &to_json(const SymGroup &group, jsonParser &json);

  // Note: as a hack this expects group[0] to be present and have the right lattice!!!
  //   it's just used to set the lattice for all the Molecules
  void from_json(SymGroup &group, const jsonParser &json);


};

#endif
