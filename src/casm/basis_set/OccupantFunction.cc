
#include "casm/basis_set/OccupantFunction.hh"

#include "casm/symmetry/SymOp.hh"
#include "casm/basis_set/FunctionVisitor.hh"

namespace CASM {


  //*******************************************************************************************


  void OccupantFunction::fill_dispatch_table() {
    //Function::inner_prod_table[sclass_ID()][sclass_ID()] = new BasicVarVarScalarProd();
    Function::operation_table[sclass_ID()][sclass_ID()] = new OccOccOperation();

    //BasicVarPolyProduct does not yet exist
    //Function::inner_prod_table[sclass_ID()][PolynomialFunction::sclass_ID()]=new BasicVarMonoProduct();
  }

  //*******************************************************************************************


  int OccupantFunction::class_ID() const {
    return DerivedID<OccupantFunction, Function>::get_class_ID();
  }

  //*******************************************************************************************


  int OccupantFunction::sclass_ID() {
    return DerivedID<OccupantFunction, Function>::get_class_ID();
  }

  //*******************************************************************************************


  Function *OccupantFunction::copy() const {
    return new OccupantFunction(*this);
  }

  //*******************************************************************************************


  //*******************************************************************************************

  bool OccupantFunction::_accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr/*=NULL*/) {
    return visitor.visit(*this, home_basis_ptr);
  }

  //*******************************************************************************************
  //John G 010413


  void OccupantFunction::make_formula() const {
    m_formula.clear();
    m_tex_formula.clear();

    std::stringstream tformula, ttex;
    Array<int> var_ind;

    // count the number of ones in var_ind;
    Index one_count(0);
    for(Index i = 0; i < dof().size(); i++) {
      if(!almost_zero(m_eval_table[i]))
        var_ind.push_back(i);

      if(almost_equal(m_eval_table[i], 1.0))
        one_count++;
    }
    if(!var_ind.size()) {
      m_formula = "0";
      m_tex_formula = "0";
      return;
    }

    if(one_count == dof().size()) {
      m_formula = "1";
      m_tex_formula = "1";
      return;
    }

    double var_scale(m_eval_table[var_ind[0]]);
    if(almost_equal(var_scale, -1.0)) {
      ttex << '-';
    }
    else if(!almost_equal(var_scale, 1.0)) {
      ttex << irrational_to_tex_string(var_scale, 2 * m_eval_table.size());
    }

    if(var_ind.size() > 1) {
      ttex << "(";
    }

    for(Index i = 0; i < var_ind.size(); i++) {
      if(i > 0 && m_eval_table[var_ind[i]] > 0) {
        tformula << '+';
      }
      if(almost_zero(m_eval_table[var_ind[i]] + 1)) {
        tformula << '-';
      }
      if(!almost_zero(std::abs(m_eval_table[var_ind[i]]) - 1)) {
        tformula << m_eval_table[var_ind[i]];
        tformula << '*';
      }

      if(i > 0 && m_eval_table[var_ind[i]] / var_scale > 0) {
        ttex << '+';
      }
      else if(almost_equal(m_eval_table[var_ind[i]] / var_scale, -1.0)) {
        ttex << '-';
      }

      if(!almost_zero(std::abs(m_eval_table[var_ind[i]] / var_scale) - 1)) {
        ttex << irrational_to_tex_string(m_eval_table[var_ind[i]] / var_scale, 2 * m_eval_table.size());
      }

      //tformula << "p_" << m_var[var_ind[i]].name;
      //tformula << "p_" << m_var[var_ind[i]].name;
    }
    if(var_ind.size() > 1) {
      ttex << ')';
    }
    m_tex_formula = ttex.str();
    m_formula = tformula.str();
    return;
  }

  //*******************************************************************************************

  int OccupantFunction::register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles) {

    if(dof().type_name() == dof_name) {
      if(!valid_index(dof().ID()) || dof().ID() >= remote_handles.size()) {
        std::cerr << "CRITICAL ERROR: In OccupantFunction::register_remotes(), dof().ID() = " << dof().ID() << " is out of bounds.\n"
                  << "                Exiting...\n";
        exit(1);
      }
      //std::cout << "Setting remote at Occ DoF " << dof().ID() << " of " << discrete_remotes.size() << "\n";
      m_var->register_remote(remote_handles[dof().ID()]);
      return 1;
    }

    return 0;
  }

  //*******************************************************************************************


  bool OccupantFunction::_update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs) {
    if(dof().is_locked()) return false;

    Index ID_ind = before_IDs.find(dof().ID());

    if(ID_ind < after_IDs.size()) {
      m_var->set_ID(after_IDs[ID_ind]);
      m_formula.clear();
      m_tex_formula.clear();
      return true;
    }
    //std::cout << "*****COULD NOT UPDATE ID " << dof().ID() << "!\n";
    return false;
  }

  //*******************************************************************************************


  bool OccupantFunction::compare(const OccupantFunction *RHS)const {
    return dof().ID() == RHS->dof().ID()
           && dof().type_name() == RHS->dof().type_name()
           && almost_equal(m_eval_table, RHS->m_eval_table);
  }

  //*******************************************************************************************


  Eigen::VectorXd const *OccupantFunction::get_eigen_coeffs() const {
    return &m_eval_table;
  }

  //\John G 010413

  //*******************************************************************************************


  bool OccupantFunction::is_zero() const {
    if(m_eval_table.size() == 0) return true;

    for(EigenIndex i = 0; i < m_eval_table.size(); i++) {
      if(!almost_zero(m_eval_table[i]))
        return false;
    }
    return true;
  }

  //*******************************************************************************************


  void OccupantFunction::small_to_zero(double tol) {

    for(EigenIndex i = 0; i < m_eval_table.size(); i++) {
      if(almost_zero(m_eval_table[i], tol))
        m_eval_table[i] = 0;
    }

  }

  //*******************************************************************************************


  Index OccupantFunction::num_terms() const {
    Index tnum(0);

    for(EigenIndex i = 0; i < m_eval_table.size(); i++) {
      if(!almost_zero(m_eval_table[i]))
        tnum++;
    }
    return tnum;
  }

  //*******************************************************************************************


  double OccupantFunction::leading_coefficient() const {

    for(EigenIndex i = 0; i < m_eval_table.size(); i++) {
      if(!almost_zero(m_eval_table[i]))
        return m_eval_table[i];
    }
    return 0.0;
  }

  //*******************************************************************************************


  double OccupantFunction::leading_coefficient(Index &index) const {

    for(index = 0; EigenIndex(index) < m_eval_table.size(); index++) {
      if(!almost_zero(m_eval_table[index]))
        return m_eval_table[index];
    }
    return 0.0;
  }

  //*******************************************************************************************


  double OccupantFunction::get_coefficient(Index i) const {
    if(EigenIndex(i) < m_eval_table.size())
      return m_eval_table[i];

    return 0.0;
  }

  //*******************************************************************************************


  Function *OccupantFunction::_apply_sym(const SymOp &op) {
    if(m_sym_rep_ID.is_identity()) return this;

    m_formula.clear();
    refresh_ID();
    Eigen::MatrixXd const *tmat;
    tmat = op.get_matrix_rep(m_sym_rep_ID);

    //Eigen matrix multiplication is alias-safe
    if(tmat)
      m_eval_table = (*tmat) * m_eval_table;
    //m_eval_table.transform(*tmat);
    else
      std::cerr << "WARNING: Attempting to reference invalid symmetry matrix in OccupantFunction::apply_sym!  Continuing...\n";
    return this;

  }

  //*******************************************************************************************


  void OccupantFunction::scale(double scale_factor) {
    m_formula.clear();
    refresh_ID();
    m_eval_table *= scale_factor;
  }

  //*******************************************************************************************


  double OccupantFunction::remote_eval() const {
    return m_eval_table[dof().remote_value()];
  }

  //*******************************************************************************************

  double OccupantFunction::remote_deval(const DoF::RemoteHandle &dvar) const {
    if(dvar.i_ptr() && dvar.i_ptr() == dof().remote_ptr()) {
      std::cerr << "CRITICAL ERROR: OccupantFunction::remote_deval() is not implemented! Exiting...\n";
      exit(1);
    }
    return 0.0;
  }

  //*******************************************************************************************


  double OccupantFunction::eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const {
    Index which_state = dof_IDs.find(dof().ID());
    if(which_state < var_states.size())
      return m_eval_table[var_states[which_state]];
    //else
    std::cerr << "CRITICAL ERROR:  In OccupantFunction::eval(), DoF ID of this OccupantFunction (" << dof().ID() << ") is not among specified IDs: " << dof_IDs << "\n"
              << "                 Exiting...\n";
    assert(0);
    exit(1);
  }

  //*******************************************************************************************


  bool OccOccOperation::compare(const Function *LHS, const Function *RHS) const {
    return static_cast<OccupantFunction const *>(LHS)->compare(static_cast<OccupantFunction const *>(RHS));

  }

  //*******************************************************************************************
  //** jsonParser stuff - OccupantFunction
  //*******************************************************************************************

  jsonParser &OccupantFunction::to_json(jsonParser &json) const {

    Function::to_json(json);

    // Just write Function::formula for now

    return json;
  }

  //*******************************************************************************************

  /*
  void OccupantFunction::from_json(const jsonParser &json) {

    // no reading functions for now

  }
  */

  //*******************************************************************************************

  jsonParser &to_json(const OccupantFunction &func, jsonParser &json) {
    return func.to_json(json);
  }

  /*
  // No reading functions for now
  void from_json(OccupantFunction &func, const jsonParser &json) {
    return func.from_json(json);
  }
  */

}
