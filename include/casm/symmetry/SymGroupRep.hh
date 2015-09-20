#ifndef SYMGROUPREP_HH
#define SYMGROUPREP_HH

#include <iostream>
#include <string>
#include <iomanip>

#include "casm/container/Array.hh"
#include "casm/symmetry/SymOp.hh"
#include "casm/symmetry/SymGroup.hh"

namespace CASM {
  class SymGroupRep;
  class Lattice;
  class MasterSymGroup;

  class SymOpRepresentation;
  class SymGroup;
  class Permutation;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class SymGroupRepHandle;

  ///\brief SymGroupRep is an alternative representation of a SymGroup for something other than real space.
  /// There is a one-to-one correspondence of SymOps in some SymGroup with the SymOpRepresentations in SymGroupRep
  /// SymGroupRep does not know or care about the specifics of what the SymOpRepresentations describe
  /// or how they are implemented
  class SymGroupRep : public Array<SymOpRepresentation *> {
  private:
    /// REP_COUNT keeps track of the rep_ID that will be given to the next SymGroupRep to be created
    static Index REP_COUNT;

    /// rep_ID is unique identifier of a specific SymGroupRep instantiation
    mutable Index m_rep_ID;

    /// Pointer to the home_group that generated this SymGroupRep
    MasterSymGroup const *m_home_group;
    // don't use default constructor
    SymGroupRep() {
      std::cerr << "Cannot perform default construction of SymGroupRep.\nExiting...\n";
      exit(1);
    };
    using Array<SymOpRepresentation *>::push_back;
  public:
    typedef SymGroupRepHandle RemoteHandle;
    enum NullInitializer {NO_HOME};
    static Index CURR_REP_COUNT() {
      return REP_COUNT;
    };

    /// Use this constructor when MasterSymGroup is unknown or doesn't exist
    /// You must promise that you know what you're doing
    SymGroupRep(SymGroupRep::NullInitializer init) : m_rep_ID(REP_COUNT++), m_home_group(NULL) { };

    SymGroupRep(const MasterSymGroup &_home) : m_rep_ID(REP_COUNT++), m_home_group(&_home) { };
    SymGroupRep(const SymGroupRep &RHS);

    ~SymGroupRep();

    void push_back(const SymOpRepresentation &new_op);

    /// SPECIAL METHOD: const, but it modifies symop representation 'i',
    /// after performing appropriate checks. Do not refer to this as an example
    /// of how to write mutators.
    /// Critical Failure if symop representation 'i' is already initialized
    void set_rep(Index i, const SymOpRepresentation &new_op) const;

    void set_master_group(const MasterSymGroup &master);

    /// Adds this representation its home_group, if it hasn't been added yet
    /// returns its rep_ID
    Index add_self_to_home();

    SymGroupRep &operator=(const SymGroupRep &RHS);

    void clear();

    SymGroupRep *copy()const {
      return new SymGroupRep(*this);
    };

    //John G 050513

    void print_permutation(std::ostream &stream) const;
    void print_MatrixXd(std::ostream &stream) const;
    void print_MatrixXd(std::ostream &stream, const SymGroup &subgroup) const;

    //\John G

    // block_shape_matrix is sum of squares of each (i,j) matrix element over all operations in SymGroupRep
    // It reveals the block_diagonalization of the symgrouprep
    Eigen::MatrixXd block_shape_matrix() const;
    Eigen::MatrixXd block_shape_matrix(const SymGroup &subgroup) const;

    // counts number of blocks in the block_shape_matrix
    Index num_blocks() const;
    Index num_blocks(const SymGroup &subgroup) const;

    SymOpRepresentation const *get_representation(const SymOpRepresentation &) const;

    Eigen::MatrixXd const *get_MatrixXd(Index i) const;
    Eigen::MatrixXd const *get_MatrixXd(const SymOpRepresentation &) const;
    Permutation const *get_permutation(Index i) const;
    Permutation const *get_permutation(const SymOpRepresentation &) const;
    Index get_ID() const {
      return m_rep_ID;
    };

    Index refresh_ID() {
      return m_rep_ID = REP_COUNT++;
    };

    /// Make a copy of representation on vector space 'V' that is transformed into a representation on vector space 'W'
    /// 'trans_mat' is the unitary matrix that isomorphically maps 'V'->'W' (i.e., [w = trans_mat * v] and [v = trans_mat.transpose() * w] )
    /// If the original representation to be transformed is just a temporary standalone SymGroupRep, be sure to delete it before falling out of scope
    //SymGroupRep *coord_transformed_copy(const Eigen::MatrixXd &trans_mat) const;

    /// Find a new coordinate system oriented along high-symmetry directions in vector space 'V' as determined by
    /// the subset of SymOpRepresentations specified by 'subgroup'. Then makes a copy of SymGroupRep in the new coordinate system
    /// If the original representation to be transformed is just a temporary standalone SymGroupRep, be sure to delete it before falling out of scope
    //SymGroupRep *coord_symmetrized_copy(const SymGroup &subgroup) const;

    //ReturnArray<Array<Eigen::VectorXd> > calc_special_directions(const SymGroup &subgroup)const;
    //ReturnArray<Array< Eigen::MatrixXd> > calc_special_subspaces(const SymGroup &subgroup)const;



    //ReturnArray<Index> num_each_irrep() const;
    //ReturnArray<Index> num_each_irrep(const SymGroup &sub_group) const;
    //ReturnArray<Index> num_each_real_irrep(const SymGroup &subgroup) const;

    //ReturnArray<Index> get_irrep_IDs(const SymGroup &subgroup) const;

    //    void calc_new_irreps(int max_iter = 1000) const;
    //  void calc_new_irreps(const SymGroup &sub_group, int max_iter = 1000) const;
    //bool is_irrep() const;
    //bool is_irrep(const SymGroup &head_group) const;

    //Eigen::MatrixXd get_irrep_trans_mat(const SymGroup &head_group) const;
    //ReturnArray<Eigen::MatrixXd> get_projection_operators() const;

    jsonParser &to_json(jsonParser &json) const;

    // If 'm_home_group' is not NULL, should be initialized accordingly
    // Lattice may be necessary for constructing SymOps
    void from_json(const jsonParser &json, const Lattice &lat);
  };

  jsonParser &to_json(const SymGroupRep &rep, jsonParser &json);

  // If 'm_home_group' is not NULL, should be initialized accordingly
  // Lattice may be necessary for constructing SymOps
  void from_json(SymGroupRep &rep, const jsonParser &json, const Lattice &lat);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /// SymGroupRepHandle is used to provide easy remote access to a SymGroup representation
  /// The head group may be a subgroup of the MasterSymGroup where the SymGroupRep is stored.
  class SymGroupRepHandle { // <-- typedefed as SymGroupRep::RemoteHandle
    SymGroup const *m_head_group;
    SymGroupRep const *m_group_rep;
  public:
    SymGroupRepHandle():
      m_head_group(NULL), m_group_rep(NULL) {}

    SymGroupRepHandle(const SymGroup &head_group, Index symrep_ID) :
      m_head_group(&head_group), m_group_rep(NULL) {
      assert(m_head_group->size() && valid_index(symrep_ID) && (*m_head_group)[0].has_valid_master());
      m_group_rep = (*m_head_group)[0].master_group().representation(symrep_ID);
      assert(m_group_rep);
    }

    /// Size of the associate SymGroup
    Index size() const {
      return m_head_group->size();
    }

    SymOpRepresentation const *operator[](Index i) const {
      return m_group_rep->at((*m_head_group)[i].index());
    }

    const SymOp &sym_op(Index i) const {
      return (*m_head_group)[i];
    }

    Index ind_inverse(Index i) const {
      return m_head_group->ind_inverse(i);
    }

    Index ind_prod(Index i, Index j) const {
      return m_head_group->ind_prod(i, j);
    }

    bool operator==(const SymGroupRepHandle &RHS) {
      return m_head_group == RHS.m_head_group && m_group_rep == RHS.m_group_rep;
    }
  };
};

#endif
