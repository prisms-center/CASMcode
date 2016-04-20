#ifndef BASISSET_HH
#define BASISSET_HH

#include <iostream>
#include <sstream>  //John G 010413
#include <map>
#include <memory>
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/BasisFunction.hh"

//#include "../container/SparseTensor.cc"
//#include "../container/Tensor.cc"
//#include "../basis_set/DoF.cc"
//#include "../symmetry/SymOp.hh"
//#include "../symmetry/SymGroup.hh"
//#include "../misc/HierarchyID.cc"

namespace CASM {

  class FunctionVisitor;
  class SymOp;
  class SymGroup;
  class SymGroupRep;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class BasisSet: private Array<Function *>, public std::enable_shared_from_this<BasisSet> {
  public:
    typedef Array<Index> SubBasis;
    typedef std::pair<SubBasis, Index> PolyConstraint;

    using Array<Function *> :: size;
    using Array<Function *> :: find;
    using Array<Function *> :: swap_elem;

    BasisSet(const std::string &name = "") :
      m_name(name), m_basis_ID(_new_ID()), m_min_poly_order(0), m_max_poly_order(-1) {}
    BasisSet(const BasisSet &init_basis);
    const BasisSet &operator=(const BasisSet &RHS);

    //Delete Functions upon destruction
    ~BasisSet();

    void clear();

    const std::string &name() const {
      return m_name;
    }

    SymGroupRepID basis_symrep_ID()const {
      return m_basis_symrep_ID;
    }

    Index max_poly_order()const {
      return m_max_poly_order;
    }
    Index min_poly_order()const {
      return m_min_poly_order;
    }

    const Array<PolyConstraint> &min_poly_constraints() const {
      return m_min_poly_constraints;
    }
    const Array<PolyConstraint> &max_poly_constraints() const {
      return m_max_poly_constraints;
    }
    BasisSet poly_quotient_set(const Function *divisor) const;

    std::shared_ptr<BasisSet> shared_copy() const {
      return std::make_shared<BasisSet>(*this);
    }

    Function const *operator[](Index i) const {
      return Array<Function *> :: operator[](i);
    }

    Function const *back() const {
      return Array<Function *> :: back();
    }

    const double &eval_cache(Index i) const {
      assert(i < size() && i < m_eval_cache.size());
      return m_eval_cache[i];
    }

    const double &deval_cache(Index i) const {
      assert(i < size() && i < m_deval_cache.size());
      return m_deval_cache[i];
    }

    bool compare(const BasisSet &RHS) const;

    int dependency_layer() const;

    void clear_formulae() {
      for(Index i = 0; i < size(); i++) {
        at(i)->clear_formula();
      }
    }

    void set_name(const std::string &new_name) {
      m_name = new_name;
    }

    void set_basis_symrep_ID(SymGroupRepID new_ID) {
      m_basis_symrep_ID = new_ID;
    }

    void add_min_poly_constraint(const Array<Index> &expons, Index expon_sum) {
      m_min_poly_constraints.push_back(PolyConstraint(expons, expon_sum));
    }

    void add_max_poly_constraint(const Array<Index> &expons, Index expon_sum) {
      m_max_poly_constraints.push_back(PolyConstraint(expons, expon_sum));
    }

    bool satisfies_exponent_constraints(const Array<Index> &expons)const;

    bool accept(const FunctionVisitor &visitor);

    /// Remotely evaluate each basis function and add it to the respective value in cumulant
    void remote_eval_and_add_to(Array<double> &cumulant)const;

    /// Remotely evaluate derivative of each basis function (w.r.t. dvar) and add it to the respective value in cumulant
    void remote_deval_and_add_to(Array<double> &cumulant, const DoF::RemoteHandle &dvar)const;

    /// Remotely evaluate each basis function and add it to the respective value in cumulant
    template<typename IteratorType>
    void remote_eval_to(IteratorType result_begin, IteratorType result_end)const;

    /// Remotely evaluate derivative of each basis function (w.r.t. dvar) and add it to the respective value in cumulant
    template<typename IteratorType>
    void remote_deval_to(IteratorType result_begin, IteratorType result_end, const DoF::RemoteHandle &dvar)const;

    /// Define the basis set to contain only variables (e.g., x,y,z)
    void set_variable_basis(const Array<ContinuousDoF> &tvar_compon, SymGroupRepID _sym_rep_ID);

    void set_dof_IDs(const Array<Index> &new_IDs);

    const Array<Index> &dof_IDs() const {
      return m_dof_IDs;
    }

    const SubBasis &dof_sub_basis(Index i) const {
      return m_dof_subbases[i];
    }

    const Array<SubBasis> &dof_sub_bases() const {
      return m_dof_subbases;
    }

    std::vector<std::set<Index> > independent_sub_bases() const;

    int register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles);

    void append(const BasisSet &RHS);

    /// Construct a polynomial basis set that contains all allowed polynomials of functions specified by tsubs
    /// if tsubs specifies, e.g., {{x,y,z}, {x,y,z}, {w,v}}, the resulting basis set will be {x*x*w, x*x*v, x*y*w+y*x*w,..., z*z*v}
    void construct_polynomials_by_order(const Array<BasisSet const *> &tsubs, Index order);

    // expects single basis set of polynomials with the same order( usually just order 1).
    // makes a basis set of all of the combinations of polynomials with order
    // spcecified by order.
    //void construct_polynomials_by_order(const Array<BasisSet const * > &tsubs, Index order);

    void construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs,
                                                  const Eigen::MatrixXd &gram_mat,
                                                  Index basis_ind,
                                                  const SymGroup &symgroup);

    void construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs,
                                                  const Array<double> &occ_probs,
                                                  Index basis_ind,
                                                  const SymGroup &symgroup);

    void construct_invariant_polynomials(const Array<BasisSet const *> &tsubs,
                                         const SymGroup &head_sym_group,
                                         Index order,
                                         Index min_dof_order = 1);

    //void construct_green_lagrange_dot_prods(const Array<BasisSet const *> &site_disp_dofs, BasisSet const *LG_strain_dofs, const Eigen::MatrixXd &ref_clust);
    //void construct_disp_grad_dot_prods(const Array<BasisSet const *> &site_disp_dofs, BasisSet const *F_strain_dofs, const Eigen::MatrixXd &ref_clust);
    //void construct_quadratic_ccds(const Array<BasisSet const *> &site_disp_dofs, BasisSet const *strain_dofs, const Eigen::MatrixXd &ref_clust,
    //                              const Array<Array<SymOp> > &equiv_map, const SymGroupRep &permute_group, double sigma);

    //void construct_harmonic_polynomials(const Array<BasisSet const *> &rotation_dofs, Index order, Index min_order, bool even);

    void calc_invariant_functions(const SymGroup &head_sym_group);

    bool is_normal_basis_for(const SymGroup &head_sym_group);

    BasisSet calc_normal_basis(const SymGroup &head_sym_group, Eigen::MatrixXd &trans_mat)const;

    BasisSet transform_copy(const Eigen::MatrixXd &trans_mat) const;

    BasisSet &apply_sym(const SymOp &op, int dependency_layer = 1);

    //returns true if Gram_Schmidt leave BasisSet unchanged.
    bool Gram_Schmidt();

    bool Gaussian_Elim();

    void get_symmetry_representation(const SymGroup &head_sym_group) const;

    //return true if basis_set was already orthogonal to argument
    bool make_orthogonal_to(Function const *ortho_func);
    bool make_orthogonal_to(const BasisSet &ortho_basis);

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);

  private:
    mutable SymGroupRepID m_basis_symrep_ID;

    std::string m_name;
    Index m_basis_ID;

    std::vector<std::shared_ptr<BasisSet> > m_argument;

    // When forming tensor products with other basis sets,
    // what is the minimum and maximum order of polynomial function
    // that elements of this basis set can be used in
    // For example, if basis_set1 = {u,v,w} with
    // min_poly_order 0 and max_poly_order 2
    // and basis_set2 = {x, y} with min_poly_order 1 and max_poly_order 1
    // then the product space of basis_set1 and basis_set2 yields
    // {x, y, x*u, x*v, x*w, y*u, y*v,y*w, x*u*u, x*u*v, x*u*w, ..., y*w*w}
    Index m_min_poly_order, m_max_poly_order;

    Array<Index> m_dof_IDs;
    Array<SubBasis> m_dof_subbases;
    Array<PolyConstraint> m_min_poly_constraints, m_max_poly_constraints;

    mutable std::vector<double> m_eval_cache;
    mutable std::vector<double> m_deval_cache;

    // public push_back is unsafe. Use BasisSet::append() instead,
    // which calls this private push_back method
    void push_back(Function *new_func);

    // non-const access to last function
    Function *_back() {
      return Array<Function *> :: back();
    }

    // non-const access to function at Index (i)
    Function *&_at(Index i) {
      return Array<Function *> :: at(i);
    }

    static Index _new_ID() {
      static Index ID_COUNT(0);
      return ID_COUNT++;
    }

    void _refresh_ID() {
      m_basis_ID = _new_ID();
    }

    void _eval_to_cache() const;
    void _deval_to_cache(const DoF::RemoteHandle &_dvar) const;

    Function *_linear_combination(const Eigen::VectorXd &coeffs) const;

    void _set_arguments(const Array<BasisSet const *> &new_args);
    void _set_arguments(const std::vector<std::shared_ptr<BasisSet> > &new_args) {
      m_argument = new_args;
    }

    void _update_dof_IDs(const Array<Index> before_IDs, const Array<Index> &after_IDs);

    //non-const access to constraints
    Array<PolyConstraint> &_min_poly_constraints() {
      return m_min_poly_constraints;
    }
    Array<PolyConstraint> &_max_poly_constraints() {
      return m_max_poly_constraints;
    }



  };

  jsonParser &to_json(const BasisSet &bset, jsonParser &json);
  void from_json(BasisSet &bset, const jsonParser &json);

  BasisSet operator*(const SymOp &LHS, const BasisSet &RHS);


  //*******************************************************************************************
  /// Remotely evaluate each basis function and add it to the respective value in cumulant
  template<typename IteratorType>
  void BasisSet::remote_eval_to(IteratorType result_begin, IteratorType result_end)const {
    for(Index i = 0; i < m_argument.size(); i++)
      m_argument[i]->_eval_to_cache();

    Index i(0);
    for(; result_begin != result_end; ++result_begin) {
      assert(i < size());
      if(at(i))
        *result_begin = at(i)->cache_eval();
      i++;
    }
  }

  //*******************************************************************************************
  /// Remotely evaluate derivative of each basis function (w.r.t. dvar) and add it to the respective value in result
  template<typename IteratorType>
  void BasisSet::remote_deval_to(IteratorType result_begin, IteratorType result_end, const DoF::RemoteHandle &dvar)const {
    for(Index i = 0; i < m_argument.size(); i++)
      m_argument[i]->_deval_to_cache(dvar);

    Index i(0);
    for(; result_begin != result_end; ++result_begin) {
      assert(i < size());
      if(at(i))
        *result_begin = at(i)->cache_deval(dvar);
      i++;
    }
  }

}
#endif
