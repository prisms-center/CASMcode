#ifndef OCCUPANTFUNCTION_HH
#define OCCUPANTFUNCTION_HH

#include <iostream>
#include <sstream>
#include <map>
#include "casm/basis_set/BasisFunction.hh"
#include "casm/basis_set/DoF.hh"

namespace CASM {

  class Function;
  class FunctionOperation;
  class InnerProduct;
  class OccupantFunction;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //
  //
  //
  //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  class OccupantFunction :
    public Function, public DerivedID<OccupantFunction, Function> {
  public:
    OccupantFunction(const DiscreteDoF &init_var, const Eigen::VectorXd &init_eval, int _occ_func_ind, int _basis_ind, SymGroupRepID _sym_rep_ID):
      m_var(init_var.copy()), m_eval_table(init_eval), m_sym_rep_ID(_sym_rep_ID), m_occ_func_ind(_occ_func_ind), m_basis_ind(_basis_ind) { }

    OccupantFunction(const OccupantFunction &RHS) : Function(RHS), m_var(RHS.m_var->copy()), m_eval_table(RHS.m_eval_table),
      m_sym_rep_ID(RHS.m_sym_rep_ID), m_occ_func_ind(RHS.occ_func_ind()), m_basis_ind(RHS.basis_ind()) {}

    ~OccupantFunction() {
      if(m_var)
        delete m_var;
    }

    static int sclass_ID();
    int class_ID() const;

    std::string type_name() const {
      return "OccupantFunction";
    }

    Index occ_func_ind()const {
      return m_occ_func_ind;
    }

    Index basis_ind()const {
      return m_basis_ind;
    }

    void set_basis_ind(int new_ind) {
      m_basis_ind = new_ind;
    }

    const DiscreteDoF &dof() const {
      return *m_var;
    }


    Function *copy() const;

    bool is_zero() const;
    Index num_terms() const;

    const Eigen::VectorXd &eval_table() const {
      return m_eval_table;
    }

    double leading_coefficient() const;
    double leading_coefficient(Index &index) const;
    double get_coefficient(Index i) const;


    void small_to_zero(double tol = TOL);
    void scale(double scale_factor);
    void make_formula() const;

    int register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles);

    bool compare(const OccupantFunction *RHS) const;

    static void fill_dispatch_table();
    Eigen::VectorXd const *get_eigen_coeffs() const;

    double remote_eval() const;

    double remote_deval(const DoF::RemoteHandle &dvar) const;

    double cache_eval() const {
      return remote_eval();
    }

    double cache_deval(const DoF::RemoteHandle &dvar)const {
      return remote_deval(dvar);
    }

    double eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const;

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);
  protected:
    Function *_apply_sym(const SymOp &op);

    bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL);

    bool _update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs);

  private:
    //**Inherited from Function:**
    //  Index func_ID;
    //  Array<Function*> m_argument;
    //  mutable std::string m_formula, m_tex_formula;
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //Array<std::string> m_formula_bits;     //mutable?

    DiscreteDoF *m_var;
    Eigen::VectorXd m_eval_table;
    SymGroupRepID m_sym_rep_ID;
    Index m_occ_func_ind, m_basis_ind;

    OccupantFunction() : m_var(nullptr) {} // no default construction

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class OccOccOperation : public FunctionOperation {
  public:
    bool compare(Function const *LHS, Function const *RHS) const;

    //Function *multiply(Function const *LHS, Function const *RHS) const;
    //Function *multiply_by(Function *LHS, Function const *RHS) const;

    //Function *add(Function const *LHS, Function const *RHS) const;
    //Function *add_to(Function *LHS, Function const *RHS) const;

    //Function *subtract(Function const *LHS, Function const *RHS) const;
    //Function *subtract_from(Function *LHS, Function const *RHS) const;
  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  jsonParser &to_json(const OccupantFunction &func, jsonParser &json);
  void from_json(OccupantFunction &func, const jsonParser &json);


}
#endif
