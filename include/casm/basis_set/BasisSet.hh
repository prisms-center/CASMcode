#ifndef BASISSET_HH
#define BASISSET_HH

#include <iostream>
#include <sstream>  //John G 010413
#include <map>

#include "casm/basis_set/BasisFunction.hh"

namespace CASM {

  class DiscreteDoF;
  class SymGroup;
  class SymGroupRep;



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class BasisSet: public Array<Function *> {
    mutable Index m_basis_symrep_ID;

    // When forming tensor products with other basis sets,
    // what is the minimum and maximum order of polynomial function
    // that elements of this basis set can be used in
    // For example, if basis_set1 = {u,v,w} with
    // min_poly_order 0 and max_poly_order 2
    // and basis_set2 = {x, y} with min_poly_order 1 and max_poly_order 1
    // then the product space of basis_set1 and basis_set2 yields
    // {x, y, x*u, x*v, x*w, y*u, y*v,y*w, x*u*u, x*u*v, x*u*w, ..., y*w*w}
    Index m_min_poly_order, m_max_poly_order;

    // public push_back is unsafe. Use BasisSet::append() instead
    //using Array<Function *> :: push_back;
    Function *linear_combination(const Eigen::VectorXd &coeffs) const;
  public:
    using Array<Function *> :: back;
    using Array<Function *> :: swap_elem;

    Array<BasisSet> subspaces;

    BasisSet() : m_basis_symrep_ID(-1), m_min_poly_order(0), m_max_poly_order(-1) {};

    BasisSet(const BasisSet &init_basis);
    const BasisSet &operator=(const BasisSet &RHS);

    //Delete Functions upon destruction
    ~BasisSet();

    void clear();

    Index basis_symrep_ID()const {
      return m_basis_symrep_ID;
    };

    Index max_poly_order()const {
      return m_max_poly_order;
    };
    Index min_poly_order()const {
      return m_min_poly_order;
    };

    BasisSet poly_quotient_set(const Function *divisor) const;

    void accept(const FunctionVisitor &visitor);

    /// Remotely evaluate each basis function and add it to the respective value in cumulant
    void remote_eval_and_add_to(Array<double> &cumulant)const;

    template<typename IteratorType>
    void remote_eval_and_add_to(IteratorType begin, IteratorType end)const;

    void update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs);

    template<typename T>
    int register_remotes(const std::string &dof_name, const Array<T> &remote_vals);

    void append(const BasisSet &RHS);

    /// Construct a polynomial basis set that contains all allowed polynomials of functions specified by tsubs
    /// if tsubs specifies, e.g., {{x,y,z}, {x,y,z}, {w,v}}, the resulting basis set will be {x*x*w, x*x*v, x*y*w+y*x*w,..., z*z*v}
    void construct_discrete_occupations(const DiscreteDoF &allowed_occs, Index basis_ind, Index sym_rep_ind = -2);

    void construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs, const Eigen::MatrixXd &gram_mat, Index basis_ind, Index sym_rep_ind = -2);

    void construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs, const Array<double> &occ_probs, Index basis_ind, Index sym_rep_ind = -2);

    //void construct_invariant_new_polynomials(const BasisSet &args, Index order, const SymGroup &head_sym_group);

    void construct_invariant_cluster_polynomials(const Array<Array<BasisSet const *> > &site_args, const Array<BasisSet const *> &global_args, const SymGroup &head_group,
                                                 const SymGroupRep &permute_group, const Array<Index> &min_site_order, Index max_poly_order);
    void construct_invariant_cluster_polynomials(const Array<Array<BasisSet const *> > &site_args, const Array<BasisSet const *> &global_args, const SymGroup &head_group,
                                                 const SymGroupRep &permute_group, Index max_poly_order);

    void calc_invariant_functions(const SymGroup &head_sym_group);

    bool is_normal_basis_for(const SymGroup &head_sym_group);
    bool is_normal_basis_for_old(const SymGroup &head_sym_group);

    BasisSet calc_normal_basis(const SymGroup &head_sym_group, Eigen::MatrixXd &trans_mat)const;
    BasisSet calc_normal_basis(const SymGroup &head_sym_group, Eigen::MatrixXd &trans_mat, int principle_axis = -1, int nrecursions = 0);


    BasisSet transform_copy(const Eigen::MatrixXd &trans_mat) const;

    BasisSet &apply_sym(const SymOp &op);

    //returns true if Gram_Schmidt leave BasisSet unchanged.
    bool Gram_Schmidt();

    bool Gaussian_Elim();

    void get_symmetry_representation(const SymGroup &head_sym_group) const;

    //return true if basis_set was already orthogonal to argument
    bool make_orthogonal_to(Function const *ortho_func);
    bool make_orthogonal_to(const BasisSet &ortho_basis);

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  };

  jsonParser &to_json(const BasisSet &bset, jsonParser &json);
  void from_json(BasisSet &bset, const jsonParser &json);

  //***************************************************************************

  template<typename T>
  int BasisSet::register_remotes(const std::string &dof_name, const Array<T> &remote_vals) {
    int sum(0);
    for(int i = 0; i < size(); i++)
      sum += at(i)->register_remotes(dof_name, remote_vals);

    return sum;
  }

  template<typename IteratorType>
  void BasisSet::remote_eval_and_add_to(IteratorType it_begin, IteratorType it_end)const {
    Index i = 0;
    for(; it_begin != it_end; ++it_begin) {
      if(!at(i)) continue;
      (*it_begin) += at(i++)->remote_eval();
    }
  }

}
#endif
