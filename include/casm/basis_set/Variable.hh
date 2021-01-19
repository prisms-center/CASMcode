#ifndef VARIABLE_HH
#define VARIABLE_HH

#include <iostream>
#include <sstream>

#include "casm/basis_set/BasisFunction.hh"
#include "casm/basis_set/DoFSet.hh"
#include "casm/symmetry/SymGroupRepID.hh"

namespace CASM {

class Function;
class FunctionVisitor;
class FunctionOperation;
class InnerProduct;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Function is a virtual class, which means you use it like this:
//
//  Function* my_funct_ptr=new VariableArgument("x");
//  my_funct_ptr->make_formula();
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Variable : public Function, public DerivedID<Variable, Function> {
  //**Inherited from Function:**
  //  int func_ID;
  //  Array<Function*> m_argument;
  //  mutable std::string m_formula, m_tex_formula;
  //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // m_dof_set defines coordinate system of variables
  DoFSet m_dof_set;

  // coeffs defines linear combination (i.e., vector) of m_dof_set
  Eigen::VectorXd m_coeffs;

  // Default construction not allowed
  // Variable(){}

  // Copy construction is private.
  // Copy construction should only occur in Variable::copy()
  Variable(const Variable &old_var) = default;

 public:
  Variable(const DoFSet &tvar, int var_ind);
  Variable(const DoFSet &tvar, const Eigen::VectorXd &init_coeffs);

  static void fill_dispatch_table();

  std::string type_name() const override { return "Variable"; }

  Function *copy() const override { return new Variable(*this); }

  bool is_zero() const override;
  void small_to_zero(double tol = TOL) override;
  Index num_terms() const override;

  double leading_coefficient() const override;
  double leading_coefficient(Index &index) const override;
  double get_coefficient(Index i) const override;

  const DoFSet &dof_set() const { return m_dof_set; }

  SymGroupRepID symrep_ID() const { return dof_set().symrep_ID(); }

  const Eigen::VectorXd &coeffs() const { return m_coeffs; }

  void make_formula() const override;

  std::set<Index> dof_IDs() const override;

  int register_remotes(
      const std::vector<DoF::RemoteHandle> &remote_handles) override;

  bool compare(const Variable *RHS) const;

  int class_ID() const override {
    return DerivedID<Variable, Function>::get_class_ID();
  }
  static int sclass_ID() {
    return DerivedID<Variable, Function>::get_class_ID();
  }

  double dot(Function const *RHS) const;
  void scale(double scale_factor) override;

  double remote_eval() const override;

  double remote_deval(const DoF::RemoteHandle &dvar) const override;

  double cache_eval() const override { return remote_eval(); }

  double cache_deval(const DoF::RemoteHandle &dvar) const override {
    return remote_deval(dvar);
  }

  Function *minus_equals(const Variable *RHS);
  Function *plus_equals(const Variable *RHS);

  jsonParser &to_json(jsonParser &json) const override;
  void from_json(const jsonParser &json);

 protected:
  Function *_apply_sym(const SymOp &op) override;
  bool _accept(const FunctionVisitor &visitor,
               BasisSet const *home_basis_ptr = NULL) override;
  bool _accept(const FunctionVisitor &visitor,
               BasisSet const *home_basis_ptr = NULL) const override;

  bool _update_dof_IDs(const std::vector<Index> &before_IDs,
                       const std::vector<Index> &after_IDs) override;
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

}  // namespace CASM
#endif
