#include "casm/basis_set/Variable.hh"

#include <memory>

#include "casm/container/PolyTrie.hh"

#include "casm/symmetry/SymOp.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/PolynomialFunction.hh"
#include "casm/basis_set/BasisSet.hh"

namespace CASM {

  //*******************************************************************************************

  void Variable::fill_dispatch_table() {
    Function::inner_prod_table[sclass_ID()][sclass_ID()] = new BasicVarVarScalarProd();
    Function::operation_table[sclass_ID()][sclass_ID()] = new VarVarOperation();

  }
  /*
  template< class Base >
  template< class Derived >
  int HierarchyID<Derived>::new_class_ID(const DerivedID<Derived, Base> &_ID){
    static int Nclass=0;
    Base::extend_hierarchy();
    Derived::fill_dispatch_table();
    return Nclass++;
  }
  */


  //*******************************************************************************************

  Variable::Variable(const Variable &old_var) :
    m_var_compon(old_var.var_compon()),
    m_sym_rep_ID(old_var.sym_rep_ID()),
    m_coeffs(old_var.coeffs()) {

  }
  //*******************************************************************************************

  Variable::Variable(const Array<ContinuousDoF> &_var_compon, int var_ind, SymGroupRepID _sym_rep_ID) :
    m_var_compon(_var_compon),
    m_sym_rep_ID(_sym_rep_ID),
    m_coeffs(_var_compon.size()) {

    m_coeffs.setZero();
    m_coeffs[var_ind] = 1.0;
  }

  //*******************************************************************************************

  Variable::Variable(const Array<ContinuousDoF> &_var_compon, const Eigen::VectorXd &_coeffs, SymGroupRepID _sym_rep_ID) :
    m_var_compon(_var_compon),
    m_sym_rep_ID(_sym_rep_ID),
    m_coeffs(_coeffs) {

  }

  //*******************************************************************************************

  bool Variable::_accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr/*=NULL*/) {
    return visitor.visit(*this, home_basis_ptr);
  }

  //*******************************************************************************************
  bool Variable::is_zero() const {
    //Could check to see if arguments are zero.  For now, assume this is done at time of construction
    for(EigenIndex i = 0; i < m_coeffs.size(); i++) {
      if(!almost_zero(m_coeffs[i])) {
        return false;
      }
    }
    return true;
  }

  //*******************************************************************************************

  void Variable::small_to_zero(double tol) {
    for(EigenIndex i = 0; i < m_coeffs.size(); i++) {
      if(almost_zero(m_coeffs[i], tol)) {
        m_coeffs[i] = 0.0;
      }
    }
  }

  //*******************************************************************************************

  Index Variable::num_terms() const {
    Index np(0);

    for(Index i = 0; i < m_var_compon.size(); i++) {
      if(!almost_zero(m_coeffs[i])) {
        np++;
      }
    }

    return np;

  }

  //*******************************************************************************************

  double Variable::leading_coefficient() const {
    for(Index i = 0; i < m_var_compon.size(); i++) {
      if(!almost_zero(m_coeffs[i])) {
        return m_coeffs[i];
      }
    }
    return 0.0;
  }

  //*******************************************************************************************

  double Variable::leading_coefficient(Index &index) const {
    for(index = 0; index < m_var_compon.size(); index++) {
      if(!almost_zero(m_coeffs[index])) {
        return m_coeffs[index];
      }
    }
    return 0.0;
  }

  //*******************************************************************************************

  double Variable::get_coefficient(Index i) const {
    if(EigenIndex(i) < m_coeffs.size())
      return m_coeffs[i];
    return 0.0;
  }

  //*******************************************************************************************

  void Variable::make_formula() const {
    m_formula.clear();
    m_tex_formula.clear();

    std::stringstream tformula, ttex;
    Array<int> var_ind;
    for(Index i = 0; i < m_var_compon.size(); i++) {
      if(!almost_zero(m_coeffs[i])) {
        var_ind.push_back(i);
      }
    }

    if(!var_ind.size()) {
      m_formula = "0";
      m_tex_formula = "0";
      return;
    }

    double var_scale(m_coeffs[var_ind[0]]);
    if(almost_zero(var_scale + 1)) {
      ttex << '-';
    }
    else if(!almost_zero(var_scale - 1)) {
      ttex << irrational_to_tex_string(var_scale, 10 * m_coeffs.size() * m_coeffs.size());
    }

    if(var_ind.size() > 1) {
      tformula << "(";
      ttex << "(";
    }

    for(Index i = 0; i < var_ind.size(); i++) {
      if(i > 0 && m_coeffs[var_ind[i]] > 0) {
        tformula << '+';
      }
      if(almost_zero(m_coeffs[var_ind[i]] + 1)) {
        tformula << '-';
      }
      if(!almost_zero(std::abs(m_coeffs[var_ind[i]]) - 1)) {
        tformula << m_coeffs[var_ind[i]];
        tformula << '*';
      }

      if(i > 0 && m_coeffs[var_ind[i]] / var_scale > 0) {
        ttex << '+';
      }
      if(almost_zero(m_coeffs[var_ind[i]] / var_scale + 1)) {
        ttex << '-';
      }
      if(!almost_zero(std::abs(m_coeffs[var_ind[i]] / var_scale) - 1)) {
        ttex << irrational_to_tex_string(m_coeffs[var_ind[i]] / var_scale, 10 * m_coeffs.size() * m_coeffs.size()) << ' ';
      }

      tformula << m_var_compon[var_ind[i]].type_name() << '[' << m_var_compon[var_ind[i]].ID() << ']';
      ttex << m_var_compon[var_ind[i]].type_name() << '_' << m_var_compon[var_ind[i]].ID();
    }
    if(var_ind.size() > 1) {
      tformula << ")";
      ttex << ')';
    }
    m_tex_formula = ttex.str();
    m_formula = tformula.str();
    return;
  }

  //*******************************************************************************************
  Function *Variable::_apply_sym(const SymOp &op) {
    m_formula.clear();
    refresh_ID();
    Eigen::MatrixXd const *tmat;
    tmat = op.get_matrix_rep(m_sym_rep_ID);
    if(tmat) {
      m_coeffs = (*tmat) * m_coeffs;
    }
    else {
      std::cerr << "WARNING: Attempting to reference invalid symmetry matrix (From rep_ID " << m_sym_rep_ID << ") in Variable::apply_sym!  Continuing...\n";
    }
    return this;
  }

  //*******************************************************************************************
  int Variable::register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles) {
    int t_tot(0);
    for(Index i = 0; i < m_var_compon.size(); i++) {
      if(m_var_compon[i].type_name() == dof_name) {
        if(!valid_index(m_var_compon[i].ID()) || m_var_compon[i].ID() >= remote_handles.size()) {
          std::cerr << "CRITICAL ERROR: In Variable::register_remotes(), m_var_compon[" << i << "].ID() = " << m_var_compon[i].ID() << " is out of bounds.\n"
                    << "                Exiting...\n";
          exit(1);
        }
        m_var_compon[i].register_remote(remote_handles[m_var_compon[i].ID()]);
        t_tot++;
      }
    }
    return t_tot;
  }

  //*******************************************************************************************

  bool Variable::_update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs) {

    Index ID_ind;
    bool is_updated(false);

    for(Index i = 0; i < m_var_compon.size(); i++) {
      // IMPORTANT: Do before_IDs.find(), NOT m_var_compon.find() (if such a thing existed)
      ID_ind = before_IDs.find(m_var_compon[i].ID());
      // Only set ID if DoF doesn't have an ID lock
      if(ID_ind < after_IDs.size() && !m_var_compon[i].is_locked()) {
        m_var_compon[i].set_ID(after_IDs[ID_ind]);
        // The new ID only changes the formula if the corresponding coeff is nonzero
        if(!almost_zero(m_coeffs[i]))
          is_updated = true;
      }
    }
    if(is_updated) {
      m_formula.clear();
      m_tex_formula.clear();
    }
    return is_updated;
  }

  //*******************************************************************************************

  bool Variable::compare(const Variable *RHS)const {

    if(m_var_compon.size() != RHS->m_var_compon.size())
      return false;

    for(Index i = 0; i < m_var_compon.size(); i++) {
      if(!m_var_compon[i].compare((RHS->m_var_compon)[i], false))
        return false;
    }

    return almost_equal(m_coeffs, RHS->m_coeffs);

  }

  //*******************************************************************************************

  void Variable::scale(double scale_factor) {
    m_tex_formula.clear();
    m_formula.clear();
    refresh_ID();
    m_coeffs *= scale_factor;
  }

  //*******************************************************************************************

  double Variable::remote_eval() const {
    double t_eval(0.0);
    for(Index i = 0; i < m_var_compon.size(); i++)
      t_eval += m_coeffs[i] * m_var_compon[i].remote_value();

    return t_eval;
  }

  //*******************************************************************************************

  double Variable::remote_deval(const DoF::RemoteHandle &dvar) const {
    for(Index i = 0; i < m_var_compon.size(); i++) {
      if(dvar.d_ptr() && dvar.d_ptr() == m_var_compon[i].remote_ptr()) {
        //std::cout << "derivative hit: " << m_coeffs[i] << '*' <<  m_var_compon[i].type_name() << '[' << m_var_compon[i].ID() << "]\n";
        return m_coeffs[i];
      }
    }
    return 0.0;
  }

  //*******************************************************************************************

  Function *Variable::minus_equals(const Variable *RHS) {
    m_coeffs -= RHS->coeffs();
    return this;
  }

  //*******************************************************************************************

  Function *Variable::plus_equals(const Variable *RHS) {
    m_coeffs += RHS->coeffs();
    return this;
  }

  //*******************************************************************************************
  double BasicVarVarScalarProd::dot(Function const *LHS, Function const *RHS) const {
    Variable const *VLHS(static_cast<Variable const *>(LHS)),
             *VRHS(static_cast<Variable const *>(RHS));
    return (VLHS->coeffs()).dot(VRHS->coeffs());
  }

  //*******************************************************************************************
  bool VarVarOperation::compare(const Function *LHS, const Function *RHS) const {
    return static_cast<const Variable *>(LHS)->compare(static_cast<const Variable *>(RHS));
  }

  //*******************************************************************************************
  Function *VarVarOperation::multiply_by(Function *LHS, Function const *RHS) const {
    std::cerr << "WARNING: In-place multiplication of one Variable with another is not defined!!\n"
              << "         Exiting...\n";
    exit(1);
    return NULL;

  }

  //*******************************************************************************************

  Function *VarVarOperation::multiply(Function const *LHS, Function const *RHS) const {

    Variable const *VLHS(static_cast<Variable const *>(LHS));
    Variable const *VRHS(static_cast<Variable const *>(RHS));
    BasisSet LSet, RSet;
    LSet.set_variable_basis(VLHS->var_compon(), VLHS->sym_rep_ID());
    std::vector<std::shared_ptr<BasisSet> > args;

    Index full_dim(LSet.size());
    Index LOffset(0);
    args.push_back(LSet.shared_copy());

    if((VLHS->var_compon()) != (VRHS->var_compon())) {
      RSet.set_variable_basis(VRHS->var_compon(), VRHS->sym_rep_ID());
      args.push_back(RSet.shared_copy());
      LOffset = LSet.size();
      full_dim += RSet.size();
    }
    PolyTrie<double> tcoeffs(full_dim);
    Array<Index> expon(full_dim, 0);
    for(Index i = 0; i < (VLHS->coeffs()).size(); i++) {
      if(almost_zero((VLHS->coeffs())[i]))
        continue;
      expon[i]++;
      for(Index j = 0; j < (VLHS->coeffs()).size(); j++) {
        if(almost_zero((VRHS->coeffs())[j]))
          continue;
        expon[j + LOffset]++;
        tcoeffs(expon) = (VLHS->coeffs())[i] * (VRHS->coeffs())[j];
        expon[j + LOffset]--;
      }
      expon[i]--;
    }

    return new PolynomialFunction(args, tcoeffs);

  }

  //*******************************************************************************************
  Function *VarVarOperation::subtract(Function const *LHS, Function const *RHS) const {

    Variable *VLHS_copy(static_cast<Variable *>(LHS->copy()));
    Variable const *VRHS(static_cast<Variable const *>(RHS));

    VLHS_copy->minus_equals(VRHS);
    VLHS_copy->refresh_ID();
    return VLHS_copy;
  }

  //*******************************************************************************************
  Function *VarVarOperation::subtract_from(Function *LHS, Function const *RHS) const {
    Variable *VLHS(static_cast<Variable *>(LHS));
    Variable const *VRHS(static_cast<Variable const *>(RHS));

    VLHS->minus_equals(VRHS);
    VLHS->refresh_ID();

    return VLHS;

  }

  //*******************************************************************************************
  Function *VarVarOperation::add(Function const *LHS, Function const *RHS) const {
    Variable *VLHS_copy(static_cast<Variable *>(LHS->copy()));
    Variable const *VRHS(static_cast<Variable const *>(RHS));

    VLHS_copy->plus_equals(VRHS);
    VLHS_copy->refresh_ID();

    return VLHS_copy;
  }
  //*******************************************************************************************
  Function *VarVarOperation::add_to(Function *LHS, Function const *RHS) const {
    Variable *VLHS(static_cast<Variable *>(LHS));
    Variable const *VRHS(static_cast<Variable const *>(RHS));

    VLHS->plus_equals(VRHS);

    VLHS->refresh_ID();

    return VLHS;
  }


  //*******************************************************************************************
  //** jsonParser stuff - Variable
  //*******************************************************************************************

  jsonParser &Variable::to_json(jsonParser &json) const {

    Function::to_json(json);

    // Just write Function::formula for now

    return json;
  }

  //*******************************************************************************************

  /*
  void Variable::from_json(const jsonParser &json) {

    // no reading functions for now

  }
  */

  //*******************************************************************************************

  jsonParser &to_json(const Variable &var, jsonParser &json) {
    return var.to_json(json);
  }

  /*
  // No reading functions for now
  void from_json(Variable &var, const jsonParser &json) {
    return var.from_json(json);
  }
  */

}
