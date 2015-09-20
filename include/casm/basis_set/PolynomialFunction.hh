#ifndef POLYNOMIALFUNCTION_HH
#define POLYNOMIALFUNCTION_HH

#include <iostream>
#include <sstream>  //John G 010413
#include <map>

#include "casm/container/PolyTrie.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {

  class OccupantFunction;
  class InnerProduct;
  class FunctionOperation;
  class PolynomialFunction;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class PolynomialFunction :
    public Function, public DerivedID<PolynomialFunction, Function> {

    //**Inherited from Function:**
    //  Array<Function*> m_argument;
    //  mutable std::string m_formula, m_tex_formula;
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // sub_sym_reps has sym_rep index of each subspace
    Array<Index> m_sub_sym_reps;

    // subspace_array specifies how Functions in 'argument'
    // are divided into subspaces that transform distinctly
    // under symmetry.  There is one interior array for each subspace
    // e.g.: if argument contains {x, y, v, w}
    //       and {x,y} transforms separately from {v,w}, then
    //       dim2arg_array is {{0,1},{2,3}}
    Array<Array<Index> > m_subspaces;

    // arg2sub tells which subspace each argument belongs to
    // in example above for polynomialcs of {x, y, v, w}
    // arg2sub will be {0,0,1,1}
    Array<Index> m_arg2sub;

    PolyTrie<double> m_coeffs;

  public:
    PolynomialFunction() : m_coeffs(0) {  };
    PolynomialFunction(const Array<BasisSet > &init_args);
    PolynomialFunction(const Array<BasisSet const *> &init_args);
    PolynomialFunction(const PolynomialFunction &RHS);
    PolynomialFunction(const PolynomialFunction &RHS, const PolyTrie<double> &_coeffs);

    //Create new polynomial function that is product of two others.
    PolynomialFunction(const PolynomialFunction &LHS, const PolynomialFunction &RHS);


    static void fill_dispatch_table();
    //SparseTensor<double> const *get_coeffs()const;

    std::string type_name()const {
      return "PolynomialFunction";
    };

    Function *copy()const;

    bool accept(const FunctionVisitor &visitor);

    bool depends_on(const Function *test_func) const;
    bool is_zero() const;
    void small_to_zero(double tol = TOL);
    Index num_terms() const;

    double leading_coefficient() const;
    double leading_coefficient(Index &index) const;
    double get_coefficient(Index i) const;

    void make_formula() const;
    void make_formula(double prefactor) const;
    int class_ID() const {
      return DerivedID<PolynomialFunction, Function>::get_class_ID();
    };
    static int sclass_ID() {
      return DerivedID<PolynomialFunction, Function>::get_class_ID();
    };

    //double dot(const PolynomialFunction &RHS) const;
    Function *apply_sym(const SymOp &op);

    Function *transform_monomial_and_add_new(double prefactor, const Array<Index> &ind, const SymOp &op);
    Function *transform_monomial_and_add(double prefactor, const Array<Index> &ind, const SymOp &op);
    void scale(double scale_factor);

    double frobenius_scalar_prod(const PolynomialFunction &RHS)const;
    double box_integral_scalar_prod(const PolynomialFunction &RHS, double edge_length)const;
    double gaussian_integral_scalar_prod(const PolynomialFunction &RHS, double std_dev)const;

    Function *plus_equals(const PolynomialFunction *RHS);
    Function *minus_equals(const PolynomialFunction *RHS);

    double remote_eval() const;

    double eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const;
    double eval(const Array<Index> &dof_IDs, const Array<double> &arg_states) const;
    double poly_eval(const Array<double> &arg_states) const;


    Function *poly_quotient(const OccupantFunction *RHS) const;
    Function *poly_remainder(const OccupantFunction *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class BasicPolyPolyScalarProd : public InnerProduct {
  public:
    double dot(Function const *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class PolyPolyOperation : public FunctionOperation {
  public:
    //bool compare(Function const *LHS, Function const *RHS) const;

    Function *multiply(Function const *LHS, Function const *RHS) const;
    Function *multiply_by(Function *LHS, Function const *RHS) const;

    Function *add(Function const *LHS, Function const *RHS) const;
    Function *add_to(Function *LHS, Function const *RHS) const;

    Function *subtract(Function const *LHS, Function const *RHS) const;
    Function *subtract_from(Function *LHS, Function const *RHS) const;
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
