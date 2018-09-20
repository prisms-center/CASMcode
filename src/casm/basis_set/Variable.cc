#include "casm/basis_set/Variable.hh"

#include <memory>

#include "casm/misc/CASM_Eigen_math.hh"

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

  //*******************************************************************************************

  Variable::Variable(const DoFSet &_dof_set, int var_ind) :
    m_dof_set(_dof_set),
    m_coeffs(_dof_set.size()) {

    m_coeffs.setZero();
    m_coeffs[var_ind] = 1.0;
  }

  //*******************************************************************************************

  Variable::Variable(const DoFSet &_dof_set, const Eigen::VectorXd &_coeffs) :
    m_dof_set(_dof_set),
    m_coeffs(_coeffs) {

  }

  //*******************************************************************************************

  bool Variable::_accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr/*=NULL*/) {
    return visitor.visit(*this, home_basis_ptr);
  }

  //*******************************************************************************************

  bool Variable::_accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr/*=NULL*/) const {
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

    for(Index i = 0; i < dof_set().size(); i++) {
      if(!almost_zero(m_coeffs[i])) {
        np++;
      }
    }

    return np;

  }

  //*******************************************************************************************

  double Variable::leading_coefficient() const {
    for(Index i = 0; i < dof_set().size(); i++) {
      if(!almost_zero(m_coeffs[i])) {
        return m_coeffs[i];
      }
    }
    return 0.0;
  }

  //*******************************************************************************************

  double Variable::leading_coefficient(Index &index) const {
    for(index = 0; index < dof_set().size(); index++) {
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
    for(Index i = 0; i < dof_set().size(); i++) {
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

      tformula << dof_set()[var_ind[i]].var_name() << '[' << dof_set()[var_ind[i]].ID() << ']';
      ttex << dof_set()[var_ind[i]].var_name() << '_' << dof_set()[var_ind[i]].ID();
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
    tmat = op.get_matrix_rep(symrep_ID());
    if(tmat) {
      m_coeffs = (*tmat) * m_coeffs;
    }
    else {
      std::cerr << "WARNING: Attempting to reference invalid symmetry matrix (From rep_ID " << sym_rep_ID() << ") in Variable::apply_sym!  Continuing...\n";
    }
    return this;
  }

  //*******************************************************************************************
  int Variable::register_remotes(const std::vector<DoF::RemoteHandle> &remote_handles) {
    int t_tot(0);
    for(Index i = 0; i < dof_set().size(); i++) {
      std::vector<DoF::RemoteHandle>::const_iterator it = find(remote_handles.begin(),
                                                               remote_handles.end(),
                                                               dof_set()[i].handle());
      if(it != remote_handles.end()) {
        if(!valid_index(dof_set()[i].ID())) {
          throw std::runtime_error("In Variable::register_remotes(), attempting to register dof with ID = "
                                   + std::to_string(dof_set()[i].ID()) + ", which is out of bounds.\n");
        }
        dof_set()[i].register_remote(*it);
        t_tot++;
      }
    }
    return t_tot;
  }

  //*******************************************************************************************

  bool Variable::_update_dof_IDs(const std::vector<Index> &before_IDs, const std::vector<Index> &after_IDs) {

    bool is_updated = m_dof_set.update_IDs(before_IDs, after_IDs);

    if(is_updated) {
      m_formula.clear();
      m_tex_formula.clear();
    }
    return is_updated;
  }

  //*******************************************************************************************

  bool Variable::compare(const Variable *RHS)const {

    if(dof_set().size() != RHS->dof_set().size())
      return false;

    for(Index i = 0; i < dof_set().size(); i++) {
      if(!dof_set()[i].compare((RHS->dof_set())[i], false))
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
    for(Index i = 0; i < dof_set().size(); i++)
      t_eval += m_coeffs[i] * dof_set()[i].remote_value();

    return t_eval;
  }

  //*******************************************************************************************

  double Variable::remote_deval(const DoF::RemoteHandle &dvar) const {
    for(Index i = 0; i < dof_set().size(); i++) {
      if(dvar.d_ptr() && dvar.d_ptr() == dof_set()[i].remote_ptr()) {
        //std::cout << "derivative hit: " << m_coeffs[i] << '*' <<  dof_set()[i].var_name() << '[' << dof_set()[i].ID() << "]\n";
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
    LSet.set_variable_basis(VLHS->dof_set());
    std::vector<std::shared_ptr<BasisSet> > args;

    Index full_dim(LSet.size());
    Index LOffset(0);
    args.push_back(LSet.shared_copy());

    if((VLHS->dof_set()) != (VRHS->dof_set())) {
      RSet.set_variable_basis(VRHS->dof_set());
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
