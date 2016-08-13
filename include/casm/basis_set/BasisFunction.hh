#ifndef BASISFUNCTION_HH
#define BASISFUNCTION_HH

#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <memory>

#include "casm/external/Eigen/Dense"
#include "casm/CASM_global_definitions.hh"
#include "casm/container/Array.hh"
#include "casm/misc/HierarchyID.hh"
#include "casm/basis_set/DoF.hh"

//#include "../container/SparseTensor.cc"
//#include "../container/Tensor.cc"
//#include "../basis_set/DoF.cc"
//#include "../symmetry/SymOp.hh"
//#include "../symmetry/SymGroup.hh"
//#include "../misc/HierarchyID.cc"

namespace CASM {

  class Function;
  class FunctionVisitor;
  class FunctionOperation;
  class InnerProduct;
  class Variable;
  class BasisSet;
  class DoF;
  class SymOp;
  template <typename T>
  class SparseTensor;



  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /**
     Function is a virtual class from which all basis function types
     (e.g., Variable, PolynomialFunction, etc.) are derived
   **/

  class Function : public HierarchyID<Function> {
  protected:
    typedef std::vector<std::shared_ptr<BasisSet> > ArgumentContainer;
  public:
    //Function(std::string &init_formula) : func_ID(ID_count++), m_formula(init_formula), m_tex_formula(init_formula) {}
    Function(const Function &RHS) : func_ID(RHS.func_ID), m_argument(RHS.m_argument), m_arg2sub(RHS.m_arg2sub), m_arg2fun(RHS.m_arg2fun),
      m_formula(RHS.m_formula), m_tex_formula(RHS.m_tex_formula) { }
    Function(const ArgumentContainer &_args);

    Function() : func_ID(ID_count++) {}

    //WARNING: If you write a destructor for a Function-derived class,
    //         it should internally call Function::~Function()
    virtual ~Function() {};
    void refresh_ID();
    Index ID() const {
      return func_ID;
    }

    Index num_args() const {
      return m_arg2sub.size();
    }


    std::string formula() const;
    std::string tex_formula() const;

    void print(std::ostream &stream) const;
    void print_tex(std::ostream &stream) const;
    void set_label_format(const std::string &format) {
      m_label_format = format;
    }
    const std::string &label_format()const {
      return m_label_format;
    }
    void set_formula(const std::string &new_formula) {
      m_formula = new_formula;
      m_tex_formula = new_formula;
    }
    void set_tex_formula(const std::string &new_formula) {
      m_tex_formula = new_formula;
    }
    void clear_formula() {
      m_formula.clear();
      m_tex_formula.clear();
      //for(Index i = 0; i < m_argument.size(); i++)
      //m_argument[i]->clear_formulae();
    }

    virtual std::string type_name() const = 0;

    //virtual void make_formula(double prefactor)const =0;
    virtual Function *copy() const = 0;
    virtual void make_formula()const = 0;
    virtual bool is_zero() const = 0;

    // very basic dependency check -- only compares pointer values for now and checks that coefficient is nonzero
    virtual bool depends_on(const Function *test_func) const {
      return this == test_func;
    }

    //for derived accept method, always do:
    //accept(const FunctionVisitor &visitor){Function::accept(visitor); visitor->visit(*this);}
    bool accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL);
    virtual void small_to_zero(double tol = TOL) = 0;
    virtual Index num_terms() const = 0;
    virtual double leading_coefficient() const = 0;
    virtual double leading_coefficient(Index &index) const = 0;
    virtual double get_coefficient(Index i) const = 0;
    virtual int class_ID() const = 0;
    virtual void scale(double scale_factor) = 0;
    virtual SparseTensor<double> const *get_coeffs()const {
      return NULL;
    }

    virtual Eigen::VectorXd const *get_eigen_coeffs() const {
      return nullptr;
    }

    virtual double remote_eval() const = 0;
    virtual double remote_deval(const DoF::RemoteHandle &dvar) const = 0;

    virtual double cache_eval() const = 0;
    virtual double cache_deval(const DoF::RemoteHandle &dvar) const = 0;

    //virtual double eval(int var_state) const;
    virtual double eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const;
    virtual double eval(const Array<Index> &dof_IDs, const Array<double> &arg_states) const;

    virtual int register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles);

    bool update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs);


    virtual Function *apply_sym_coeffs(const SymOp &op, int dependency_layer = 1) {
      if(_dependency_layer() == dependency_layer)
        return _apply_sym(op);

      //for(Index i = 0; i < m_argument.size(); ++i)
      //m_argument[i]->apply_sym(op, dependency_layer);

      return this;
    }

    Function *sym_copy_coeffs(const SymOp &op, int dependency_layer = 1) const;

    void normalize();

    double dot(Function const *RHS) const;

    // shallow comparison assumes that m_argument is equivalent to RHS.m_argument
    bool shallow_compare(Function const *RHS) const;

    // 'deep' comparison also checks that m_argument is equivalent to RHS.m_argument
    bool compare(Function const *RHS) const;

    Function *minus(Function const *RHS) const;
    Function *plus(Function const *RHS) const;
    Function *multiply(Function const *RHS) const;
    Function *poly_quotient(Function const *RHS) const;
    Function *poly_remainder(Function const *RHS) const;
    Function *minus_in_place(Function const *RHS);
    Function *plus_in_place(Function const *RHS);


    void set_arguments(const ArgumentContainer &new_arg) {
      m_argument = new_arg;
    }


    const ArgumentContainer &argument_bases() const {
      return m_argument;
    }

    static void print_table() {
      for(Index i = 0; i < inner_prod_table.size(); i++) {
        for(Index j = 0; j < inner_prod_table[i].size(); j++)
          std::cout << inner_prod_table[i][j] << "  ";
        std::cout << '\n';
      }
    }

    virtual jsonParser &to_json(jsonParser &json) const;

  protected:
    // static members hold tables of objects that define DerivedA--DerivedB operations
    static Array< Array< InnerProduct * > > inner_prod_table;
    static Array< Array< FunctionOperation * > > operation_table;

    // Function::extend_hierarchy() is called at first object initialization of a new derived type
    static void extend_hierarchy() {
      for(Index i = 0; i < inner_prod_table.size(); i++) {
        operation_table[i].push_back(NULL);
        inner_prod_table[i].push_back(NULL);
      }

      inner_prod_table.push_back(Array<InnerProduct * > (inner_prod_table.size() + 1, NULL));
      operation_table.push_back(Array<FunctionOperation * > (operation_table.size() + 1, NULL));
    }

    Index func_ID;

    ArgumentContainer m_argument;

    /// m_label_format sets the label format used to generate a label string for a Function object.
    /// It is specified as a string of the form (substr1 + "%a" + substr2 + "%b" + substr3 + ... ),
    /// where "%a" and "%b" are flags that specify object-specific values.
    /// The following flags are allowed:
    ///  - %f : function index (only available for some derived types)
    ///  - %b : basis index (only available for OccupantFunction and Variable)
    ///  - %n : neighbor list index -- the DoF ID if dof.is_locked()==false
    ///  - %g : global DoF index -- the DoF ID if dof.is_locked()==true
    /// For objects where the flag does not uniquely specify a single value, it evaluates to a substring
    /// that concatenates multiple values in ascending order.
    /// Example: For a polynomial function that combines DoFs from the sites {8, 2, 4} of the neighborlist,
    ///          %n will evaluate to the substring expression "2_4_8"
    std::string m_label_format;


    //    double scale_val;

    // arg2sub tells which subspace each argument belongs to
    // in example above for polynomials of {x, y, v, w}
    // arg2sub will be {0,0,1,1}
    Array<Index> m_arg2sub;
    Array<Index> m_arg2fun;


    mutable std::string m_formula, m_tex_formula;

    ReturnArray<SymGroupRepID> _sub_sym_reps() const;

    Function const *_argument(Index i) const;

    double _arg_eval_cache(Index i) const;

    double _arg_deval_cache(Index i) const;

    //Function *_argument(Index i);

    // dependency layer is number of functional layers between the degree of freedom and the current layer
    // Consider the function H(G(F(x))), where 'x' is a degree of freedom.
    //            'x' has dependency layer = 0
    //            F(x) has dependency layer = 1
    //            G(F) has dependency layer = 2
    //            H(G) has dependency layer = 3
    //            G(F) * F(x) has dependency layer = 3  // because it is a function K(G,F)
    int _dependency_layer() const;

    virtual Function *_apply_sym(const SymOp &op) = 0;

    virtual bool _accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr = NULL) = 0;

    virtual bool _update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs) {
      // default action: do nothing, report that function does not change (via 'return false');
      return false;
    }
  private:
    friend class HierarchyID<Function>;
    static Index ID_count;
  };

  jsonParser &to_json(const Function *func, jsonParser &json);

  // This does not exist: void from_json(Function *func, const jsonParser &json);
  // Use the json (i.e. json["Function_type"] = "PolynomialFunction")
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
