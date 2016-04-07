#ifndef VARIABLE_HH
#define VARIABLE_HH

#include <iostream>
#include <sstream>
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/basis_set/BasisFunction.hh"


namespace CASM {

  class Function;
  class FunctionVisitor;
  class FunctionOperation;
  class InnerProduct;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //Function is a virtual class, which means you use it like this:
  //
  //  Function* my_funct_ptr=new VariableArgument("x");
  //  my_funct_ptr->make_formula();
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class Variable :
    public Function, public DerivedID<Variable, Function> {
    //**Inherited from Function:**
    //  int func_ID;
    //  Array<Function*> m_argument;
    //  mutable std::string m_formula, m_tex_formula;
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // m_var_compon defines coordinate system of variables
    Array<ContinuousDoF> m_var_compon;

    //Symmetry representation for variable space (e.g., {x,y,z})
    SymGroupRepID m_sym_rep_ID;

    //coeffs defines linear combination (i.e., vector) of m_var_compon
    Eigen::VectorXd m_coeffs;

    //Default construction not allowed
    //Variable(){}

    // Copy construction is private.
    // Copy construction should only occur in Variable::copy()
    Variable(const Variable &old_var);

  public:

    Variable(const Array<ContinuousDoF> &tvar, int var_ind, SymGroupRepID rep_ID);
    Variable(const Array<ContinuousDoF> &tvar, const Eigen::VectorXd &init_coeffs, SymGroupRepID rep_ID);

    static void fill_dispatch_table();

    std::string type_name() const {
      return "Variable";
    }

    Function *copy()const {
      return new Variable(*this);
    }

    bool is_zero() const;
    void small_to_zero(double tol = TOL);
    Index num_terms() const;

    double leading_coefficient() const;
    double leading_coefficient(Index &index) const;
    double get_coefficient(Index i) const;

    const Array<ContinuousDoF> &var_compon() const {
      return m_var_compon;
    }
    SymGroupRepID sym_rep_ID() const {
      return m_sym_rep_ID;
    }

    const Eigen::VectorXd &coeffs() const {
      return m_coeffs;
    }

    void make_formula()const;
    void make_formula(double prefactor)const;

    int register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles);

    bool compare(const Variable *RHS) const;

    int class_ID() const {
      return DerivedID<Variable, Function>::get_class_ID();
    }
    static int sclass_ID() {
      return DerivedID<Variable, Function>::get_class_ID();
    }

    double dot(Function const *RHS) const;
    void scale(double scale_factor);

    double remote_eval() const;

    double remote_deval(const DoF::RemoteHandle &dvar) const;

    double cache_eval() const {
      return remote_eval();
    }

    double cache_deval(const DoF::RemoteHandle &dvar) const {
      return remote_deval(dvar);
    }

    Function *minus_equals(const Variable *RHS);
    Function *plus_equals(const Variable *RHS);

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  protected:
    Function *_apply_sym(const SymOp &op);
    bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL);

    bool _update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs);

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class BasicVarVarScalarProd : public InnerProduct {
  public:
    double dot(Function const *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class VarVarOperation : public FunctionOperation {
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

  jsonParser &to_json(const Variable &var, jsonParser &json);
  void from_json(Variable &var, const jsonParser &json);

}
#endif
