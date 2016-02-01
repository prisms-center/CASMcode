#ifndef BASISFUNCTION_HH
#define BASISFUNCTION_HH

#include <iostream>
#include <sstream>  //John G 010413
#include <map>

#include "casm/misc/HierarchyID.hh"
#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"

namespace CASM {

  class Function;
  class FunctionVisitor;
  class FunctionOperation;
  class InnerProduct;
  class MonomialFunction;
  class BasisSet;
  class ContinuousDoF;
  class DiscreteDoF;
  class SymOp;


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /**
     Function is a virtual class from which all basis function types
     (e.g., Variable, MonomialFunction, etc.) are derived
   **/

  class Function : public HierarchyID<Function> {
    friend class HierarchyID<Function>;
    static Index ID_count;
  protected:
    Index func_ID;

    //    double scale_val;

    // static members hold tables of objects that define DerivedA--DerivedB operations
    static Array< Array< InnerProduct * > > inner_prod_table;
    static Array< Array< FunctionOperation * > > operation_table;

    // Function::extend_hierarchy() is called at first object initialization of a new derived type
    static void extend_hierarchy() {
      for(Index i = 0; i < inner_prod_table.size(); i++) {
        operation_table[i].push_back(nullptr);
        inner_prod_table[i].push_back(nullptr);
      }

      inner_prod_table.push_back(Array<InnerProduct * > (inner_prod_table.size() + 1, nullptr));
      operation_table.push_back(Array<FunctionOperation * > (operation_table.size() + 1, nullptr));
    };

    //Function 'owns' the Function pointers in argument.
    //The Function objects they point to are deleted when the Function is destroyed.
    Array<Function *> m_argument;

    mutable std::string m_formula, m_tex_formula;

  public:

    Function(std::string &init_formula) : func_ID(ID_count++), m_formula(init_formula), m_tex_formula(init_formula) {};
    Function(const Function &RHS) : func_ID(RHS.func_ID),
      m_formula(RHS.m_formula), m_tex_formula(RHS.m_tex_formula) { };
    Function() : func_ID(ID_count++) {};

    //WARNING: If you write a destructor for a Function-derived class,
    //         it should internally call Function::~Function()
    virtual ~Function();
    void refresh_ID();
    Index ID() const {
      return func_ID;
    };


    std::string formula() const;
    std::string tex_formula() const;

    void print(std::ostream &stream) const;
    void print_tex(std::ostream &stream) const;
    void set_formula(const std::string &new_formula) {
      m_formula = new_formula;
      m_tex_formula = new_formula;
    };
    void set_tex_formula(const std::string &new_formula) {
      m_tex_formula = new_formula;
    };
    void clear_formula() {
      m_formula.clear();
      m_tex_formula.clear();
      for(Index i = 0; i < m_argument.size(); i++)
        m_argument[i]->clear_formula();
    }

    virtual std::string type_name() const = 0;

    //virtual void make_formula(double prefactor)const =0;
    virtual Function *copy() const = 0;
    virtual void make_formula()const = 0;
    virtual bool is_zero() const = 0;

    // very basic dependency check -- only compares pointer values for now and checks that coefficient is nonzero
    virtual bool depends_on(const Function *test_func) const {
      return this == test_func;
    };

    //for derived accept method, always do:
    //accept(const FunctionVisitor &visitor){Function::accept(visitor); visitor->visit(*this);};
    virtual bool accept(const FunctionVisitor &visitor);
    virtual void small_to_zero(double tol = TOL) = 0;
    virtual Index num_terms() const = 0;
    virtual double leading_coefficient() const = 0;
    virtual double leading_coefficient(Index &index) const = 0;
    virtual double get_coefficient(Index i) const = 0;
    virtual int class_ID() const = 0;
    virtual Function *apply_sym(const SymOp &op) = 0;
    virtual void scale(double scale_factor) = 0;

    virtual Eigen::VectorXd const *get_eigen_coeffs() const {
      return nullptr;
    };

    virtual double remote_eval() const = 0;

    //virtual double eval(int var_state) const;
    virtual double eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const;
    virtual double eval(const Array<Index> &dof_IDs, const Array<double> &arg_states) const;

    virtual int register_remotes(const std::string &dof_name, const Array<int> &discrete_remotes);
    virtual int register_remotes(const std::string &dof_name, const Array<double> &continuous_remotes);

    virtual bool update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs);

    Function *sym_copy(const SymOp &op) const;

    void normalize();

    double dot(Function const *RHS) const;
    bool compare(Function const *RHS) const;
    Function *minus(Function const *RHS) const;
    Function *plus(Function const *RHS) const;
    Function *multiply(Function const *RHS) const;
    Function *poly_quotient(Function const *RHS) const;
    Function *poly_remainder(Function const *RHS) const;
    Function *minus_in_place(Function const *RHS);
    Function *plus_in_place(Function const *RHS);


    static void print_table() {
      for(Index i = 0; i < inner_prod_table.size(); i++) {
        for(Index j = 0; j < inner_prod_table[i].size(); j++)
          std::cout << inner_prod_table[i][j] << "  ";
        std::cout << '\n';
      }
    };

    virtual jsonParser &to_json(jsonParser &json) const;
  };

  jsonParser &to_json(const Function *func, jsonParser &json);

  // This does not exist: void from_json(Function *func, const jsonParser &json);
  // Use the json (i.e. json["Function_type"] = "MonomialFunction")
  //   to know which Function's from_json to call


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  // InnerProduct is separate from FunctionOperation (for now)
  // Because InnerProduct behavior may be different for two Functions of the same Derived type
  class InnerProduct {
  public:
    virtual double dot(Function const *LHS, Function const *RHS) const = 0;
  };


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class FunctionOperation {
  public:

    // Comparison
    virtual bool compare(Function const *LHS, Function const *RHS) const;

    // Multiplication (tensor product)
    virtual Function *multiply(Function const *LHS, Function const *RHS) const;
    virtual Function *multiply_by(Function *LHS, Function const *RHS) const;

    // Addition
    virtual Function *add(Function const *LHS, Function const *RHS) const;
    virtual Function *add_to(Function *LHS, Function const *RHS) const;

    // Subtraction
    virtual Function *subtract(Function const *LHS, Function const *RHS) const;
    virtual Function *subtract_from(Function *LHS, Function const *RHS) const;

    // Polynomial division
    virtual Function *poly_quotient(Function const *LHS, Function const *RHS) const;
    virtual Function *poly_remainder(Function const *LHS, Function const *RHS) const;

  };

}
#endif
