#ifndef OCCUPANTFUNCTION_HH
#define OCCUPANTFUNCTION_HH

#include <iostream>
#include <sstream>
#include <map>
#include "casm/basis_set/BasisFunction.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/misc/cloneable_ptr.hh"

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
      m_var(init_var.clone()), m_eval_table(init_eval), m_sym_rep_ID(_sym_rep_ID), m_occ_func_ind(_occ_func_ind), m_basis_ind(_basis_ind) { }

    static int sclass_ID();
    int class_ID() const override;

    std::string type_name() const override {
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


    Function *copy() const override;

    bool is_zero() const override;
    Index num_terms() const override;

    const Eigen::VectorXd &eval_table() const {
      return m_eval_table;
    }

    double leading_coefficient() const override;
    double leading_coefficient(Index &index) const override;
    double get_coefficient(Index i) const override;


    void small_to_zero(double tol = TOL) override;
    void scale(double scale_factor) override;
    void make_formula() const override;

    int register_remotes(const std::vector<DoF::RemoteHandle> &remote_handles) override;


    std::set<Index> dof_IDs() const override {
      return std::set<Index>({dof().ID()});
    }
    bool compare(const OccupantFunction *RHS) const;

    static void fill_dispatch_table();
    Eigen::VectorXd const *get_eigen_coeffs() const override;

    double discrete_eval(int state) const;

    double remote_eval() const override;

    double remote_deval(const DoF::RemoteHandle &dvar) const override;

    double cache_eval() const override {
      return remote_eval();
    }

    double cache_deval(const DoF::RemoteHandle &dvar)const override {
      return remote_deval(dvar);
    }

    jsonParser &to_json(jsonParser &json) const override;
    void from_json(const jsonParser &json);
  protected:
    Function *_apply_sym(const SymOp &op) override;

    bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL) override;

    bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL) const override;

    bool _update_dof_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs) override;

  private:
    //**Inherited from Function:**
    //  Index func_ID;
    //  Array<Function*> m_argument;
    //  mutable std::string m_formula, m_tex_formula;
    //  ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //Array<std::string> m_formula_bits;     //mutable?

    notstd::cloneable_ptr<DiscreteDoF> m_var;
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
