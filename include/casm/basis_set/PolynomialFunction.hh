#ifndef POLYNOMIALFUNCTION_HH
#define POLYNOMIALFUNCTION_HH

#include <iostream>
#include <sstream>
#include <map>
#include <memory>
#include <vector>

#include "casm/basis_set/BasisFunction.hh"
#include "casm/container/PolyTrie.hh"

namespace CASM {

  class Function;
  class InnerProduct;
  class FunctionOperation;
  class PolynomialFunction;
  class OccupantFunction;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class PolynomialFunction :
    public Function, public DerivedID<PolynomialFunction, Function> {
  public:
    PolynomialFunction() : m_coeffs(0) {  };
    PolynomialFunction(const std::vector<std::shared_ptr<BasisSet> > &_args);
    PolynomialFunction(const std::vector<std::shared_ptr<BasisSet> > &_args, const PolyTrie<double> &_coeffs);

    PolynomialFunction(const PolynomialFunction &RHS);
    PolynomialFunction(const PolynomialFunction &RHS, const PolyTrie<double> &_coeffs);

    //Create new polynomial function that is product of two others.
    PolynomialFunction(const PolynomialFunction &LHS, const PolynomialFunction &RHS);


    static void fill_dispatch_table();
    SparseTensor<double> const *get_coeffs()const override;

    const PolyTrie<double> &poly_coeffs()const {
      return m_coeffs;
    };

    std::string type_name()const override {
      return "PolynomialFunction";
    };

    Function *copy()const override;
    //PolynomialFunction *copy(const PolyTrie<double> &new_coeffs)const;

    Function *copy(const PolyTrie<double> &new_coeffs)const;

    bool depends_on(const Function *test_func) const override;
    bool is_zero() const override;
    void small_to_zero(double tol = TOL) override;
    Index num_terms() const override;

    double leading_coefficient() const override;
    double leading_coefficient(Index &index) const override;
    double get_coefficient(Index i) const override;

    void make_formula() const override;
    void make_formula(double prefactor) const;

    std::set<Index> dof_IDs() const override;
    int class_ID() const override {
      return DerivedID<PolynomialFunction, Function>::get_class_ID();
    };
    static int sclass_ID() {
      return DerivedID<PolynomialFunction, Function>::get_class_ID();
    };

    Function *transform_monomial_and_add_new(double prefactor, const Array<Index> &ind, const SymOp &op, const std::vector<bool> &transform_flags);
    Function *transform_monomial_and_add(double prefactor, const Array<Index> &ind, const SymOp &op);
    void scale(double scale_factor) override;

    double frobenius_scalar_prod(const PolynomialFunction &RHS)const;
    double box_integral_scalar_prod(const PolynomialFunction &RHS, double edge_length)const;
    double gaussian_integral_scalar_prod(const PolynomialFunction &RHS, double std_dev)const;

    bool compare(const PolynomialFunction *RHS)const;
    bool prune_zeros();

    Function *plus_equals(const PolynomialFunction *RHS);
    Function *minus_equals(const PolynomialFunction *RHS);

    double remote_eval() const override;
    double remote_deval(const DoF::RemoteHandle &dvar) const override;

    double cache_eval() const override;
    double cache_deval(const DoF::RemoteHandle &dvar) const override;

    Function *apply_sym_coeffs(const SymOp &op, int dependency_layer) override;

    Function *poly_quotient(const Variable *RHS) const;
    Function *poly_remainder(const Variable *RHS) const;


    Function *poly_quotient(const OccupantFunction *RHS) const;
    Function *poly_remainder(const OccupantFunction *RHS) const;

  protected:

    //**Inherited from Function:**
    //  Array<std::shared_ptr<BasisSet> > m_argument;
    //  mutable std::string m_formula, m_tex_formula;
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // sub_sym_reps has sym_rep index of each subspace
    //Array<Index> m_sub_sym_reps;

    // subspace_array specifies how Functions in 'argument'
    // are divided into subspaces that transform distinctly
    // under symmetry.  There is one interior array for each subspace
    // e.g.: if argument contains {x, y, v, w}
    //       and {x,y} transforms separately from {v,w}, then
    //       dim2arg_array is {{0,1},{2,3}}
    //Array<Array<Index> > m_subspaces;

    // arg2sub tells which subspace each argument belongs to
    // in example above for polynomials of {x, y, v, w}
    // arg2sub will be {0,0,1,1}
    //Array<Index> m_arg2sub;
    //Array<Index> m_arg2fun;

    PolyTrie<double> m_coeffs;

    void _set_arguments(const ArgumentContainer &new_arg, std::vector<Index> const &compatibility_map) override;

    Function *_apply_sym(const SymOp &op) override;
    Function *_apply_sym(const SymOp &op, const std::vector<bool> &transform_flags);

    bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL) override;

    bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL)const override;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class BasicPolyPolyScalarProd : public InnerProduct {
  public:
    double dot(Function const *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class PolyPolyOperation : public FunctionOperation {
  public:
    bool compare(Function const *LHS, Function const *RHS) const;

    Function *multiply(Function const *LHS, Function const *RHS) const;
    Function *multiply_by(Function *LHS, Function const *RHS) const;

    Function *add(Function const *LHS, Function const *RHS) const;
    Function *add_to(Function *LHS, Function const *RHS) const;

    Function *subtract(Function const *LHS, Function const *RHS) const;
    Function *subtract_from(Function *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class PolyVarOperation : public FunctionOperation {
    Function *poly_quotient(Function const *LHS, Function const *RHS) const;
    Function *poly_remainder(Function const *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class PolyOccOperation : public FunctionOperation {
    Function *poly_quotient(Function const *LHS, Function const *RHS) const;
    Function *poly_remainder(Function const *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


}
#endif
