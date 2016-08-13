#include "casm/basis_set/BasisFunction.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {

  Index Function::ID_count(0);
  Array<Array< InnerProduct * > > Function::inner_prod_table = Array<Array< InnerProduct * > > ();
  Array<Array< FunctionOperation * > > Function::operation_table = Array<Array< FunctionOperation * > > ();


  //*******************************************************************************************

  Function::Function(const std::vector<std::shared_ptr<BasisSet> > &_args) : func_ID(ID_count++), m_argument(_args) {
    //Index linear_ind(0);
    for(Index i = 0; i < m_argument.size(); i++) {
      for(Index j = 0; j < m_argument[i]->size(); j++) {
        m_arg2sub.push_back(i);
        m_arg2fun.push_back(j);
      }
    }
  }
  //*******************************************************************************************

  double Function::dot(Function const *RHS) const {
    return inner_prod_table[this->class_ID()][RHS->class_ID()]->dot(this, RHS);
  }

  //*******************************************************************************************
  void Function::normalize() {
    double mag(sqrt(this->dot(this)));
    if(almost_zero(mag)) {
      return;
    }
    m_formula.clear();
    m_tex_formula.clear();
    this->scale(1.0 / mag);
  }

  //*******************************************************************************************

  bool Function::shallow_compare(Function const *RHS) const {
    return (operation_table[this->class_ID()][RHS->class_ID()])->compare(this, RHS);
  }

  //*******************************************************************************************

  bool Function::compare(Function const *RHS) const {
    if(m_argument.size() != (RHS->m_argument).size())
      return false;

    //if operations between this and RHS are undefined, return false
    if(!operation_table[this->class_ID()][RHS->class_ID()])
      return false;

    for(Index i = 0; i < m_argument.size(); i++) {
      if(!(m_argument[i]->compare(*(RHS->m_argument)[i])))
        return false;
    }

    return shallow_compare(RHS);
  }

  //*******************************************************************************************

  Function *Function::plus(Function const *RHS) const {
    return operation_table[this->class_ID()][RHS->class_ID()]->add(this, RHS);
  }
  //*******************************************************************************************

  Function *Function::minus(Function const *RHS) const {
    return operation_table[this->class_ID()][RHS->class_ID()]->subtract(this, RHS);
  }
  //*******************************************************************************************

  Function *Function::poly_quotient(Function const *RHS) const {
    return operation_table[this->class_ID()][RHS->class_ID()]->poly_quotient(this, RHS);
  }
  //*******************************************************************************************

  Function *Function::poly_remainder(Function const *RHS) const {
    return operation_table[this->class_ID()][RHS->class_ID()]->poly_remainder(this, RHS);
  }
  //*******************************************************************************************

  Function *Function::multiply(Function const *RHS) const {
    return (operation_table[this->class_ID()][RHS->class_ID()])->multiply(this, RHS);
  }
  //*******************************************************************************************

  Function *Function::plus_in_place(Function const *RHS) {
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    return operation_table[this->class_ID()][RHS->class_ID()]->add_to(this, RHS);
  }
  //*******************************************************************************************

  Function *Function::minus_in_place(Function const *RHS) {
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    return operation_table[this->class_ID()][RHS->class_ID()]->subtract_from(this, RHS);
  }

  //*******************************************************************************************

  ReturnArray<SymGroupRepID> Function::_sub_sym_reps() const {
    Array<SymGroupRepID> t_result(m_argument.size());
    for(Index i = 0; i < m_argument.size(); i++) {
      t_result[i] = m_argument[i]->basis_symrep_ID();
    }
    return t_result;
  }

  //*******************************************************************************************

  Function const *Function::_argument(Index i) const {
    return (*m_argument[m_arg2sub[i]])[m_arg2fun[i]];
  }

  //*******************************************************************************************

  double Function::_arg_eval_cache(Index i) const {
    return m_argument[m_arg2sub[i]]->eval_cache(m_arg2fun[i]);
  }

  //*******************************************************************************************

  double Function::_arg_deval_cache(Index i) const {
    return m_argument[m_arg2sub[i]]->deval_cache(m_arg2fun[i]);
  }

  //*******************************************************************************************

  int Function::_dependency_layer() const {
    int tdep = -1;
    for(Index i = 0; i < m_argument.size(); i++) {
      tdep = max(tdep, m_argument[i]->dependency_layer());
    }
    return ++tdep;
  }
  //*******************************************************************************************
  /*
    Function *Function::_argument(Index i) {
    return (*m_argument[m_arg2sub[i]])[m_arg2fun[i]];
    }
  */
  //*******************************************************************************************

  bool Function::accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr) {
    //bool is_updated(false);
    //std::cout << "INSIDE BASE Function::accept\n";
    // shared_basis
    //for(Index i = 0; i < m_argument.size(); i++)
    //is_updated = (m_argument[i]->accept(visitor)) || is_updated; // should we add && depends_on to first part?

    //if(is_updated) {
    //m_formula.clear();
    //m_tex_formula.clear();
    //}

    return _accept(visitor, home_basis_ptr);// || is_updated;
  }

  //*******************************************************************************************

  void Function::refresh_ID() {
    func_ID = ID_count++;
  }

  //*******************************************************************************************

  void Function::print(std::ostream &stream)const {
    if(!m_formula.size()) {
      this->make_formula();
    }
    stream << m_formula;
    return;
  }
  //*******************************************************************************************

  void Function::print_tex(std::ostream &stream)const {
    if(!m_formula.size()) {
      this->make_formula();
    }
    stream << m_tex_formula;
    return;
  }

  //*******************************************************************************************

  Function *Function::sym_copy_coeffs(const SymOp &op, int dependency_level) const {
    Function *tptr(this->copy());
    return tptr->apply_sym_coeffs(op, dependency_level);
  }

  //*******************************************************************************************

  int Function::register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles) {
    int t_tot(0);
    // From now on, this will be handled at BasisSet level
    //for(Index i = 0; i < m_argument.size(); i++) {
    //t_tot += m_argument[i]->register_remotes(dof_name, remote_handles);
    //}
    return t_tot;
  }

  //*******************************************************************************************

  bool Function::update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs) {
    //JCT shared_basis
    /*bool is_updated(false);
    for(Index i = 0; i < m_argument.size(); i++) {
      is_updated = ((m_argument[i]->update_dof_IDs(before_IDs, after_IDs)) && depends_on(m_argument[i])) || is_updated;
    }

    if(is_updated) {
      m_formula.clear();
      m_tex_formula.clear();
    }
    */
    return  _update_dof_IDs(before_IDs, after_IDs);
  }

  //*******************************************************************************************

  std::string Function::formula() const {
    if(!m_formula.size()) {
      this->make_formula();
    }

    return m_formula;
  }

  //*******************************************************************************************

  std::string Function::tex_formula() const {
    if(!m_formula.size()) {
      this->make_formula();
    }

    return m_tex_formula;
  }

  //*******************************************************************************************
  /*
  double Function::eval(int var_state) const {

    std::cerr << "WARNING: You are trying to evaluate a function using integer argument,\n"
              << "         but that option is not available for the function:\n"
              << "         " << (this->formula()) << "\n";
    return NAN;
  }
  */
  //*******************************************************************************************

  double Function::eval(const Array<Index> &dof_IDs, const Array<double> &arg_state) const {

    std::cerr << "WARNING: You are trying to evaluate a function using an Array<double> as argument,"
              << "         but that option is not available for the function:\n"
              << "         " << (this->formula()) << "\n";
    return NAN;
  }

  //*******************************************************************************************

  double Function::eval(const Array<Index> &dof_IDs, const Array<Index> &var_state) const {

    std::cerr << "WARNING: You are trying to evaluate a function using an Array<int> as argument,\n"
              << "         but that option is not available for the function:\n"
              << "         " << (this->formula()) << "\n";
    return NAN;
  }

  //*******************************************************************************************

  // Comparison
  bool FunctionOperation::compare(Function const *LHS, Function const *RHS) const {
    // default behavior if Functions are of different types
    return false;
  }

  //*******************************************************************************************

  // Multiplication (tensor product)
  Function *FunctionOperation::multiply(Function const *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Multiplication of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  Function *FunctionOperation::multiply_by(Function *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Multiplication of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  // Addition
  Function *FunctionOperation::add(Function const *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Addition of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  Function *FunctionOperation::add_to(Function *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Addition of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  // Subtraction
  Function *FunctionOperation::subtract(Function const *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Subtraction of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  Function *FunctionOperation::subtract_from(Function *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Subtraction of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  // Polynomial division
  Function *FunctionOperation::poly_quotient(Function const *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Polynomial division of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************

  Function *FunctionOperation::poly_remainder(Function const *LHS, Function const *RHS) const {
    std::cerr << "WARNING: Polynomial division of type " << LHS->type_name() << " with type " << RHS->type_name() << " is undefined!!!\n";
    return nullptr;
  }

  //*******************************************************************************************
  //** jsonParser stuff - Function
  //*******************************************************************************************

  jsonParser &Function::to_json(jsonParser &json) const {

    // Not quite sure what to do with all this.  For now, just print 'formula'.

    //  static members:
    //
    // class Function : public HierarchyID<Function>
    //   static int Nclass;
    // static int ID_count;
    // static Array< Array< InnerProduct * > > inner_prod_table;
    // static Array< Array< TensorProduct * > > tensor_prod_table;
    // static Array< Array< FunctionDifference * > > difference_table;
    // static Array< Array< FunctionSum * > > sum_table;


    // int func_ID;
    // Array<Function *> m_argument;
    // mutable std::string m_formula, m_tex_formula;

    json.put_obj();

    json["m_formula"] = m_formula;

    return json;
  }

  //*******************************************************************************************

  /*
  void Function::from_json(const jsonParser &json) {

    // no reading functions for now

  }
  */

  //*******************************************************************************************

  jsonParser &to_json(const Function *func, jsonParser &json) {
    return func->to_json(json);
  }

  // This does not exist: void from_json(Function *func, const jsonParser &json);
  // Use the json (i.e. json["Function_type"] = "MonomialFunction")
  //   to know which Function's from_json to call


}

