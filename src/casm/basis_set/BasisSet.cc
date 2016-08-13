#include "casm/basis_set/BasisSet.hh"

#include <algorithm>

#include "casm/misc/CASM_math.hh"

#include "casm/container/Permutation.hh"
#include "casm/container/Counter.hh"
#include "casm/container/IsoCounter.hh"
#include "casm/container/MultiCounter.hh"

#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymGroupRep.hh"

#include "casm/basis_set/Variable.hh"
#include "casm/basis_set/PolynomialFunction.hh"
#include "casm/basis_set/OccupantFunction.hh"
#include "casm/basis_set/FunctionVisitor.hh"

namespace CASM {

  //*******************************************************************************************

  BasisSet::BasisSet(const BasisSet &init_basis) :
    Array<Function * >(0),
    m_basis_symrep_ID(init_basis.basis_symrep_ID()),
    m_name(init_basis.name()),
    m_basis_ID(init_basis.m_basis_ID),
    m_min_poly_order(init_basis.m_min_poly_order),
    m_max_poly_order(init_basis.m_max_poly_order),
    //m_eval_cache & m_deval_cache are taken care of by BasisSet::push_back
    //m_eval_cache(init_basis.m_eval_cache),
    //m_deval_cache(init_basis.m_deval_cache),
    m_dof_IDs(init_basis.m_dof_IDs),
    m_dof_subbases(init_basis.m_dof_subbases),
    m_min_poly_constraints(init_basis.min_poly_constraints()),
    m_max_poly_constraints(init_basis.max_poly_constraints()) {

    for(Index i = 0; i < init_basis.m_argument.size(); i++) {
      m_argument.push_back(init_basis.m_argument[i]->shared_copy());
    }

    for(Index i = 0; i < init_basis.size(); i++) {
      if(!init_basis[i])
        push_back(nullptr);
      else {
        push_back(init_basis[i]->copy());
        _back()->set_arguments(m_argument);
      }
    }

  }

  //*******************************************************************************************

  const BasisSet &BasisSet::operator=(const BasisSet &RHS) {
    if(this == &RHS) {
      return *this;
    }
    clear();
    m_basis_symrep_ID = RHS.basis_symrep_ID();
    m_min_poly_order = RHS.min_poly_order();
    m_max_poly_order = RHS.max_poly_order();
    //m_eval_cache & m_deval_cache are taken care of by BasisSet::push_back
    //m_eval_cache = RHS.m_eval_cache;
    //m_deval_cache = RHS.m_deval_cache;
    m_eval_cache.clear();
    m_deval_cache.clear();
    m_name = RHS.name();
    m_dof_IDs = RHS.m_dof_IDs;
    m_dof_subbases = RHS.m_dof_subbases;
    _min_poly_constraints() = RHS.min_poly_constraints();
    _max_poly_constraints() = RHS.max_poly_constraints();

    m_argument.clear();
    for(Index i = 0; i < RHS.m_argument.size(); i++) {
      m_argument.push_back(RHS.m_argument[i]->shared_copy());
    }
    for(Index i = 0; i < RHS.size(); i++) {
      if(!RHS[i])
        push_back(nullptr);
      else {
        push_back(RHS[i]->copy());
        _back()->set_arguments(m_argument);
      }
    }
    return *this;
  }


  //*******************************************************************************************

  BasisSet::~BasisSet() {
    clear();
  }

  //*******************************************************************************************

  void BasisSet::clear() {
    for(Index i = 0; i < size(); i++) {
      if(at(i))
        delete at(i);
    }
    Array<Function *>::clear();
    _refresh_ID();
  }

  //*******************************************************************************************

  void BasisSet::append(const BasisSet &RHS) {
    //Before appending functions, copy over  DoF IDs and subbasis info
    for(Index i = 0; i < RHS.m_dof_IDs.size(); i++) {
      Index ID_ind = m_dof_IDs.find(RHS.m_dof_IDs[i]);
      if(ID_ind == m_dof_IDs.size()) {
        assert(0 && "In BasisSet::append(), it is unsafe to append a BasisSet whose dependencies differ from (this).");
        m_dof_IDs.push_back(RHS.m_dof_IDs[i]);
        m_dof_subbases.push_back(SubBasis());
      }

      for(Index j = 0; j < RHS.m_dof_subbases[i].size(); j++) {
        //std::cout << "*** Push back " << j << " to subbases " << ID_ind << " value " << RHS.m_dof_subbases[i][j] + size() << "\n\n";
        m_dof_subbases[ID_ind].push_back(RHS.m_dof_subbases[i][j] + size());
      }
    }

    //Convert polynomial constraints
    _min_poly_constraints().reserve(min_poly_constraints().size() + RHS.min_poly_constraints().size());
    for(Index i = 0; i < RHS.min_poly_constraints().size(); i++) {
      _min_poly_constraints().push_back(RHS.min_poly_constraints()[i]);
      for(Index j = 0; j < min_poly_constraints().back().first.size(); j++)
        _min_poly_constraints().back().first[i] += size();
    }

    _max_poly_constraints().reserve(max_poly_constraints().size() + RHS.max_poly_constraints().size());
    for(Index i = 0; i < RHS.max_poly_constraints().size(); i++) {
      _max_poly_constraints().push_back(RHS.max_poly_constraints()[i]);
      for(Index j = 0; j < max_poly_constraints().back().first.size(); j++)
        _max_poly_constraints().back().first[i] += size();
    }

    //Convert RHS.min_poly_order() and RHS.max_poly_order into polynomial constraints
    if(RHS.min_poly_order() > 0) {
      _min_poly_constraints().push_back(PolyConstraint(Array<Index>::sequence(size(), size() + RHS.size() - 1), RHS.min_poly_order()));
    }
    if(valid_index(RHS.max_poly_order())) {
      _max_poly_constraints().push_back(PolyConstraint(Array<Index>::sequence(size(), size() + RHS.size() - 1), RHS.max_poly_order()));
    }
    for(Index i = 0; i < RHS.size(); i++) {
      if(!RHS[i])
        push_back(nullptr);
      else {
        push_back(RHS[i]->copy());
        _back()->set_arguments(m_argument);
      }
    }
    _refresh_ID();
  }

  //*******************************************************************************************

  BasisSet BasisSet::poly_quotient_set(const Function *divisor) const {
    BasisSet new_set;
    for(Index i = 0; i < size(); i++) {
      if(!at(i))
        new_set.push_back(nullptr);
      else
        new_set.push_back(at(i)->poly_quotient(divisor));
    }
    return new_set;
  }

  //*******************************************************************************************

  bool BasisSet::satisfies_exponent_constraints(const Array<Index> &expons)const {
    Index exp_sum;
    for(Index i = 0; i < max_poly_constraints().size(); i++) {
      exp_sum = 0;
      for(Index j = 0; j < max_poly_constraints()[i].first.size(); j++) {
        exp_sum += expons[max_poly_constraints()[i].first[j]];
      }
      if(max_poly_constraints()[i].second < exp_sum)
        return false;
    }
    for(Index i = 0; i < min_poly_constraints().size(); i++) {
      exp_sum = 0;
      for(Index j = 0; j < min_poly_constraints()[i].first.size(); j++) {
        exp_sum += expons[min_poly_constraints()[i].first[j]];
      }
      if(exp_sum < min_poly_constraints()[i].second)
        return false;
    }
    exp_sum = expons.sum();
    if(valid_index(min_poly_order()) && exp_sum < min_poly_order())
      return false;
    if(valid_index(max_poly_order()) && max_poly_order() < exp_sum)
      return false;
    return true;
  }
  //*******************************************************************************************
  bool BasisSet::accept(const FunctionVisitor &visitor) {
    bool is_changed(false), tchanged;
    //std::cout << "Accepting visitor of type " << visitor.type_name() << "... \n";
    for(Index i = 0; i < m_argument.size(); i++) {
      tchanged = m_argument[i]->accept(visitor);
      is_changed = tchanged || is_changed;
    }
    //std::cout << "BasisSet " << name() << " (" << this << ") and is_changed is " << is_changed << "\n";

    for(Index i = 0; i < size(); i++) {
      if(at(i)) {
        if(is_changed)
          at(i)->clear_formula();
        tchanged = at(i)->accept(visitor, this);
        is_changed = tchanged || is_changed;
      }
    }
    //if(is_changed)_refresh_ID(); //should we do this here?
    return is_changed;
  }

  //*******************************************************************************************
  /// Remotely evaluate each basis function and add it to the respective value in cumulant
  void BasisSet::remote_eval_and_add_to(Array<double> &cumulant)const {
    assert(size() == cumulant.size());
    for(Index i = 0; i < m_argument.size(); i++)
      m_argument[i]->_eval_to_cache();

    for(Index i = 0; i < size(); i++) {
      if(at(i))
        cumulant[i] += at(i)->cache_eval();
    }
  }

  //*******************************************************************************************
  /// Remotely evaluate derivative of each basis function (w.r.t. dvar) and add it to the respective value in cumulant
  void BasisSet::remote_deval_and_add_to(Array<double> &cumulant, const DoF::RemoteHandle &dvar)const {
    assert(size() == cumulant.size());
    for(Index i = 0; i < m_argument.size(); i++)
      m_argument[i]->_deval_to_cache(dvar);

    for(Index i = 0; i < size(); i++) {
      if(at(i))
        cumulant[i] += at(i)->cache_deval(dvar);
    }
  }

  //*******************************************************************************************

  bool BasisSet::compare(const BasisSet &RHS)const {
    if(m_argument.size() != RHS.m_argument.size() || size() != RHS.size())
      return false;

    for(Index i = 0; i < m_argument.size(); i++) {
      if(!m_argument[i]->compare(*(RHS.m_argument[i])))
        return false;
    }

    for(Index i = 0; i < size(); i++) {
      if(!(at(i)->shallow_compare(RHS[i])))
        return false;
    }
    return true;
  }

  //*******************************************************************************************

  int BasisSet::dependency_layer() const {
    int result = -1;
    for(Index i = 0; i < m_argument.size(); i++)
      result = CASM::max(result, m_argument[i]->dependency_layer());
    return ++result;
  }
  //*******************************************************************************************
  /// Remotely evaluate each basis function and add it to the respective value in cumulant
  void BasisSet::_eval_to_cache()const {
    remote_eval_to(m_eval_cache.begin(), m_eval_cache.end());
  }
  //*******************************************************************************************
  /// Remotely evaluate each basis function and add it to the respective value in cumulant
  void BasisSet::_deval_to_cache(const DoF::RemoteHandle &dvar)const {
    remote_eval_to(m_eval_cache.begin(), m_eval_cache.end());
    remote_deval_to(m_deval_cache.begin(), m_deval_cache.end(), dvar);
    //std::cout << "eval_cache: " << m_eval_cache << "\ndeval_cache: " << m_deval_cache << "\n\n";

  }
  //*******************************************************************************************
  void BasisSet::set_variable_basis(const Array<ContinuousDoF> &tvar_compon, SymGroupRepID _var_sym_rep_ID) {
    m_argument.clear();
    m_basis_symrep_ID = _var_sym_rep_ID;
    Array<Index> tdof_IDs;
    for(Index i = 0; i < tvar_compon.size(); i++) {
      if(!tvar_compon[i].is_locked() && !tdof_IDs.contains(tvar_compon[i].ID()))
        tdof_IDs.push_back(tvar_compon[i].ID());
    }
    set_dof_IDs(tdof_IDs);
    for(Index i = 0; i < tvar_compon.size(); i++) {
      push_back(new Variable(tvar_compon, i, _var_sym_rep_ID));
    }
    _refresh_ID();
  }
  //*******************************************************************************************
  void BasisSet::set_dof_IDs(const Array<Index> &new_IDs) {
    _update_dof_IDs(m_dof_IDs, new_IDs);
  }
  //*******************************************************************************************
  // Pass before_IDs by value to avoid aliasing issues when called from BasisSet::set_dof_IDs()
  void BasisSet::_update_dof_IDs(const Array<Index> before_IDs, const Array<Index> &after_IDs) {
    //std::cout << "BasisSet " << this << " --> m_dof_IDs is " << m_dof_IDs << "; before_IDs: " << before_IDs << "; after_ID: " << after_IDs << "\n" << std::endl;
    if(before_IDs.size() != after_IDs.size() && size() > 0) {
      std::cerr << "CRITICAL ERROR: In BasisSet::update_dof_IDs(), new IDs are incompatible with current IDs.\n"
                << "                Exiting...\n";
      assert(0);
      exit(1);
    }

    //update m_dof_IDs after the other stuff, for easier debugging
    if(m_dof_IDs.size() == 0) {
      m_dof_IDs = after_IDs;
      m_dof_subbases.resize(after_IDs.size());
    }
    else {
      Index m;
      for(Index i = 0; i < m_dof_IDs.size(); i++) {
        m = before_IDs.find(m_dof_IDs[i]);
        if(m == before_IDs.size()) {
          std::cerr << "CRITICAL ERROR: In BasisSet::update_dof_IDs(), new IDs are incompatible with current IDs.\n"
                    << "                Exiting...\n";
          assert(0);
          exit(1);
        }
        m_dof_IDs[i] = after_IDs[m];
      }
    }

    for(Index i = 0; i < m_argument.size(); i++)
      m_argument[i]->_update_dof_IDs(before_IDs, after_IDs);

    for(Index i = 0; i < size(); i++) {
      if(at(i))
        at(i)->update_dof_IDs(before_IDs, after_IDs);
    }


    return;
  }

  //*******************************************************************************************

  std::vector< std::set<Index> > BasisSet::independent_sub_bases() const {
    std::vector< std::set<Index> > result;
    std::vector<bool> unclaimed(size(), true);
    for(Index i = 0; i < dof_sub_bases().size(); i++) {
      if(dof_sub_basis(i).size() == 0)
        continue;
      Index j = 0;
      for(j = 0; j < result.size(); j++) {
        if(result[j].find(dof_sub_basis(i)[0]) != result[j].end()) {
          break;
        }
      }
      if(j == result.size())
        result.emplace_back();
      for(Index k = 0; k < dof_sub_basis(i).size(); k++) {
        unclaimed[dof_sub_basis(i)[k]] = false;
        result[j].insert(dof_sub_basis(i)[k]);
      }
    }
    bool has_unclaimed = false;
    for(Index i = 0; i < unclaimed.size(); i++) {
      if(unclaimed[i]) {
        if(!has_unclaimed) {
          result.emplace_back();
          has_unclaimed = true;
        }
        result.back().insert(i);
      }
    }
    return result;
  }

  //*******************************************************************************************

  int BasisSet::register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles) {
    int sum(0);
    for(Index i = 0; i < size(); i++)
      sum += at(i)->register_remotes(dof_name, remote_handles);

    for(Index i = 0; i < m_argument.size(); i++)
      sum += m_argument[i]->register_remotes(dof_name, remote_handles);

    return sum;
  }

  //*******************************************************************************************

  void BasisSet::construct_polynomials_by_order(const Array<BasisSet const * > &tsubs, Index order) {


    _set_arguments(tsubs);
    //Array<BasisSet const *> abset(m_argument);
    PolynomialFunction tpoly(m_argument);
    Array<Index> curr_expon(tpoly.poly_coeffs().depth(), 0);

    IsoCounter<Array<Index> > exp_count(Array<Index>(tsubs[0]->size(), 0), Array<Index>(tsubs[0]->size(), order), 1, order);

    do {
      for(int i = 0; i < exp_count.size(); i++) {
        //std::cout << exp_count()[i] << " ";
        curr_expon[i] = exp_count()[i];
      }
      PolyTrie<double> new_trie(curr_expon.size());
      new_trie(curr_expon) = 1.0;
      push_back(tpoly.copy(new_trie));
      //std::cout << std::endl;
    }
    while(++exp_count);


  }

  //*******************************************************************************************

  void BasisSet::construct_invariant_polynomials(const Array<BasisSet const *> &tsubs, const SymGroup &head_group, Index order, Index min_dof_order) {

    _set_arguments(tsubs);
    //std::cout << "Constructing invariant polynomials from DoFs:\n";
    for(Index i = 0; i < tsubs.size(); ++i) {
      if((tsubs[i]->basis_symrep_ID()).empty()) {
        tsubs[i]->get_symmetry_representation(head_group);
      }
      //std::cout << tsubs[i]->name() << "  ";
    }

    if(tsubs.size() == 0 && min_poly_order() < 1) {
      PolyTrie<double> ttrie(0);
      ttrie(Array<Index>(0)) = 1.0;
      push_back(new PolynomialFunction(m_argument, ttrie));
      return;
    }
    //std::cout << "\n\n";
    //std::cout << "dof_IDs are: " << dof_IDs() << "\n\n";
    //std::cout << "dof_sub_bases are: " << dof_sub_bases() << "\n\n";
    PolynomialFunction *tpoly;
    Array<Index> curr_exp;
    typedef BasisSet::SubBasis SubBasis;
    typedef IsoCounter<Array<Index> > OrderCount;
    typedef MultiCounter<OrderCount> ExpCount;

    // Add constraints that set a minimum combined order for each local dofset
    for(Index i = 0; i < dof_IDs().size() && min_dof_order >= 0; i++) {
      SubBasis big_sub_basis;
      Index offset = 0;
      for(Index j = 0; j < tsubs.size(); j++) {
        //std::cout << "dof_IDs of tsub " << j << " are " << tsubs[j]->dof_IDs() << "\n";
        Index ID_ind = (tsubs[j]->dof_IDs()).find(dof_IDs()[i]);
        if(ID_ind == (tsubs[j]->dof_IDs()).size())
          continue;
        const SubBasis &sub_basis(tsubs[j]->dof_sub_basis(ID_ind));
        for(Index b = 0; b < sub_basis.size(); b++) {
          big_sub_basis.push_back(sub_basis[b] + offset);
        }
        offset += tsubs[j]->size();
      }
      //std::cout << "Adding min constraint: " << big_sub_basis << ":  " << min_dof_order << "\n";
      add_min_poly_constraint(big_sub_basis, min_dof_order);
    }
    //\ end DoFset constraints

    OrderCount::Container initial_order(tsubs.size(), 0), final_order(tsubs.size(), order);
    for(Index i = 0; i < tsubs.size(); i++) {
      if(valid_index(tsubs[i]->min_poly_order()))
        initial_order[i] = tsubs[i]->min_poly_order();
      if(valid_index(tsubs[i]->max_poly_order()))
        final_order[i] = tsubs[i]->max_poly_order();
    }
    //std::cout << "order_count from: " << initial_order << " to final order: " << final_order << "\n\n";
    OrderCount order_count(initial_order, final_order, 1, order);

    ExpCount exp_count;
    for(Index i = 0; i < tsubs.size(); i++) {
      exp_count.push_back(OrderCount(Array<Index>(tsubs[i]->size(), 0), Array<Index>(tsubs[i]->size(), order), 1, order_count[i]));
      curr_exp.append(exp_count.back());
    }

    for(; order_count.valid(); ++order_count) {
      for(Index i = 0; i < exp_count.size(); i++)
        exp_count[i].set_sum_constraint(order_count[i]);
      for(; exp_count.valid(); ++exp_count) {
        bool valid_expon = true;
        for(Index i = 0; i < exp_count.size() && valid_expon; i++)
          valid_expon = tsubs[i]->satisfies_exponent_constraints(exp_count[i]());

        if(!valid_expon)
          continue;

        Index ne = 0;
        ExpCount::const_value_iterator it(exp_count.value_begin()), it_end(exp_count.value_end());
        for(; it != it_end; ++it) {
          curr_exp[ne++] = *it;
        }
        if(!satisfies_exponent_constraints(curr_exp))
          continue;
        //else:
        //std::cout << "Adding exponent " << curr_exp << "\n\n";

        tpoly = new PolynomialFunction(m_argument);
        for(Index i = 0; i < head_group.size(); i++) {
          tpoly->transform_monomial_and_add(1, curr_exp, head_group[i]);
        }
        push_back(tpoly);
      }
    }
    Gram_Schmidt();

    //std::cout << "Result, with " << size() << " invariant functions:\n";
    //for(Index i = 0; i < size(); i++)
    //std::cout << "F_" << i << " = " << at(i)->tex_formula() << "\n\n";
    return;
  }

  //*******************************************************************************************
  // Construct generalized occupant functions that is orthonormal with respect to
  // the inner product defined by gram_mat.
  //   - gram_mat should look like a covariance matrix of random vectors, of dimension
  //              allowed_occs.size().  This means that gram_mat should be positive-definite
  //              and symmetric. Additionally, the variance of the vector of all ones
  //              (i.e., v={1,1,...,1}) should have variance = 1, which indicates the constraint
  //              of conservation of number of occupants.
  //
  // The goal is to find a matrix 'B' of column vectors such that
  //       B.transpose()*gram_mat*B = Identity
  // and such that the first column of 'B' is the vector of all ones (i.e., v={1,1,...,1})
  // The first property is general solved by the matrix
  //        B = gram_mat^(-1/2)*W
  // where gram_mat^(-1/2) is the inverse matrix square root of gram_mat, and 'W' is an arbitrary orthogonal matrix
  // The second property places constraints on 'W'.  We attempt to find a 'B' that is similar to the Chebychev
  // basis in certain limiting cases.

  void BasisSet::construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs,
                                                          const Eigen::MatrixXd &gram_mat,
                                                          Index basis_ind,
                                                          const SymGroup &symgroup) {

    m_argument.clear();
    m_max_poly_order = 1;
    Index N = allowed_occs.size();
    if(N <= 1) {
      m_basis_symrep_ID = SymGroupRepID::identity(0);
      return;
    }

    if(!allowed_occs.is_locked()) {
      set_dof_IDs(Array<Index>(1, allowed_occs.ID()));
      m_dof_subbases[0] = Array<Index>::sequence(0, N - 2);
    }
    //std::cout << "INSIDE construct_orthonormal_discrete_functions and gram_mat is \n";
    //std::cout << gram_mat << "\n\n";
    if(!almost_zero(Eigen::MatrixXd(gram_mat - gram_mat.transpose()))) {
      // we could 'fix' gram_mat, but that could cause mysterious behavior - leave it to the user
      std::cerr << "CRITICAL ERROR: Passed a Gram Matrix to BasisSet::construct_orthonormal_discrete_functions that is not symmetric.\n"
                << "                Gram Matrix is:\n" << gram_mat << "\nExiting...\n";
      exit(1);
    }

    Eigen::VectorXd conc_vec(gram_mat * Eigen::MatrixXd::Ones(N, 1));

    if(!almost_equal(1.0, conc_vec.sum())) {
      // we could 'fix' gram_mat, but that could cause mysterious behavior - leave it to the user
      std::cerr << "CRITICAL ERROR: Passed ill-conditioned Gram Matrix to BasisSet::construct_orthonormal_discrete_functions.\n"
                << "                The sum of the elements of the Gram matrix must be equal to 1.\n"
                << "                Gram Matrix is:\n" << gram_mat << "\nExiting...\n";
      exit(1);
    }


    // ** step 1: find a generic 'B' matrix

    //Use SVD instead of eigendecomposition so that 'U' and 'V' matrices are orthogonal
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> teig(gram_mat, Eigen::ComputeEigenvectors);
    if(teig.eigenvalues().minCoeff() < TOL) {
      // we could 'fix' gram_mat, but that could cause mysterious behavior - leave it to the user
      std::cerr << "CRITICAL ERROR: Passed a Gram Matrix to BasisSet::construct_orthonormal_discrete_functions that is not positive-definite.\n"
                << "                Gram Matrix is:\n" << gram_mat << "\nSmallest Eigenvalue = " << teig.eigenvalues().minCoeff() << "; Exiting...\n";
      exit(1);
    }

    // B matrix is matrix square root of gram_mat.inverse(). Its columns form an orthonormal basis wrt gram_mat
    // In other words,
    //         B = V*(1.0/sqrt(D))*V.inverse()
    // where V is matrix of eigenvectors (as columns) of gram_mat and D is diagonal matrix of eigenvalues of gram_mat
    // B is not ideally oriented for our purposes, so the rest of the algorithm will be focused on fixing the orientation
    Eigen::MatrixXd B = teig.eigenvectors() * (teig.eigenvalues().array().cwiseInverse().cwiseSqrt().matrix().asDiagonal()) * teig.eigenvectors().inverse();


    // ** step 2: Make seed basis. This will be used to seed optimized orientation of 'B'
    Eigen::MatrixXd tseed(Eigen::MatrixXd::Zero(N, N));
    Eigen::MatrixXd::Index max_ind(0);
    if(conc_vec.maxCoeff(&max_ind) < 0.75) {
      // no "outlier" probabilities-> Use Chebychev polynomials as seed
      //Fill cosine table -- columns contain powers of x from 0 to N-1
      Eigen::MatrixXd tcos_table(N, N);
      for(Index i = 0; i < N; i++) {
        tcos_table(i, 0) = 1.0;
        double x = -cos(M_PI * (i + 0.5) / double(N));
        for(Index j = 1; j < N; j++) {
          tcos_table(i, j) = x * tcos_table(i, j - 1);
        }
      }

      // QR decomposition of tcos_table yields Q matrix that holds chebychev basis
      tseed = tcos_table.householderQr().householderQ();
    }
    else {
      // there is an outlier probability --> set seed matrix to occupation basis, with specis 'i==max_ind' as solvent
      Eigen::MatrixXd::Index curr_i(0);
      for(Eigen::MatrixXd::Index i = 0; i < B.rows(); i++) {
        tseed(i, 0) = 1;
        if(i == max_ind)
          continue;
        for(Eigen::MatrixXd::Index j = 1; j < B.cols(); j++) {
          if(curr_i + 1 == j)
            tseed(i, j) = 1;
        }
        curr_i++;
      }
    }

    // ** step 3: use seed matric to find a unitary matrix that rotates 'B' a more optimal form
    // Assume: tseed = B * Q, with unitary Q
    // approximate Q by finding QR decomposition of (B.inverse() * tseed)
    Eigen::MatrixXd Q = (B.inverse() * tseed).householderQr().householderQ();

    // Rotate 'B' by multiplication with 'W'
    // eigen matrix multiplication doesn't alias

    B = B * Q;

    // Columns of B are our basis functions, orthonormal wrt gram_mat
    for(Index i = 1; i < N; i++) {
      int sign_change = 1;
      double max_abs(0.0);
      // The sign of each OccupantFunction is ambiguous, so we use a convention
      // Force sign convention max(function(occupation))=max(abs(function(occupation)))
      // If multiple occupations evaluate to the same abs(phi) and it is the maximal abs(phi),
      // then use convention that the last occurence is positive
      // It IS confusing, but here's a simple example:
      //
      //    phi(occ) = {-1, 0, 1}  is always preferred over phi_alt(occ) = {1, 0, -1}
      //
      // even though they are otherwise both equally valid basis functions
      for(Index j = 0; j < B.rows(); j++) {
        if(std::abs(B(j, i)) > (max_abs - TOL)) {
          max_abs = std::abs(B(j, i));
          sign_change = float_sgn(B(j, i));
        }
      }
      OccupantFunction tOF(allowed_occs, double(sign_change)*B.col(i), size(), basis_ind, allowed_occs.sym_rep_ID());

      push_back(tOF.copy());
    }

    // ** step 4: Calculate BasisSet symmetry representation, based on allowed_occs.sym_rep_ID() && B matrix
    // Q*B.T=B.T*S, where we know S (how to transform a column vector), and we want Q (how to transform row vector)
    // so Q=B.T*S*inv(B.T)
    if(allowed_occs.sym_rep_ID().is_identity())
      m_basis_symrep_ID = SymGroupRepID::identity(N - 1);
    else
      m_basis_symrep_ID = symgroup.master_group().add_transformed_rep(allowed_occs.sym_rep_ID(), Eigen::MatrixXd(B.transpose()));

  }

  //*******************************************************************************************
  // Construct orthonormal basis set of OccupantFunctions for degree of freedom 'allowed_occs'
  //    - 'occ_probs' is an array of probabilities; occ_prob[i] is the probability of allowed_occs
  //                  taking on value allowed_occs[i].
  //                  These are equivalent ot composition, but we make a distinction to avoid confusion
  //                  with global "parameterized compositions"
  //
  // This method finds a Gram matrix that is consistent with the probabilities and then passes it on to
  // the overloaded version that takes a Gram matrix.
  // The convention we use for the gram matrix is kind of arbitrary. Its limiting cases should yield orthonormality
  // of the Chebychev polynomials when the probabilities are equal, and orthonormality of the occupation basis when
  // only one probability is non-zero.

  void BasisSet::construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs,
                                                          const Array<double> &occ_probs,
                                                          Index basis_ind,
                                                          const SymGroup &symgroup) {

    Index N = allowed_occs.size();
    if(allowed_occs.size() != occ_probs.size()) {
      std::cerr << "CRITICAL ERROR: In BasiSet::construct_orthonormal_discrete_functions(), occ_probs and allowed_occs are incompatible!\nExiting...\n";
      assert(0);
      exit(1);
    }

    if(!almost_equal(1.0, occ_probs.sum())) {
      std::cerr << "CRITICAL ERROR: In BasiSet::construct_orthonormal_discrete_functions(), occ_probs must sum to 1 (they specify a probability distributation).\n"
                << "                occ_probs is: " << occ_probs << "\nExiting...\n";
      assert(0);
      exit(1);
    }

    // Build a matrix with N-1 non-zero eigenvalues equal to 1/N.
    // Remaining eigenvalue is zero, and corresponds to vector of all ones
    Eigen::MatrixXd gram_mat(Eigen::MatrixXd::Zero(N, N));

    for(Index i = 0; i < N; i++) {
      gram_mat(i, i) += occ_probs[i] - occ_probs[i] * occ_probs[i];
      for(Index j = 0; j < N; j++) {
        if(i == j) continue;

        gram_mat(i, i) += (occ_probs[i] - occ_probs[j]) * (occ_probs[i] - occ_probs[j]);

        gram_mat(i, j) -= occ_probs[i] * occ_probs[j] + (occ_probs[i] - occ_probs[j]) * (occ_probs[i] - occ_probs[j]);
      }
    }


    // Add in the component corresponding to vector of all ones
    // this is the uncorrelated part of the covariance
    for(Index i = 0; i < N; i++) {
      for(Index j = 0; j < N; j++) {
        gram_mat(i, j) += occ_probs[i] * occ_probs[j];
      }
    }
    construct_orthonormal_discrete_functions(allowed_occs, gram_mat, basis_ind, symgroup);
  }


  //*******************************************************************************************
  void BasisSet::calc_invariant_functions(const SymGroup &head_group) {
    Function *tfunc, *trans_func;
    for(Index nf = 0; nf < size(); nf++) {
      //std::cout << "trying function " << nf << " of " << size() << ":  ";
      //at(nf)->print(std::cout);
      //std::cout << "\n";
      tfunc = at(nf)->copy();
      at(nf)->scale(0.0);
      for(Index ng = 0; ng < head_group.size(); ng++) {
        trans_func = tfunc->sym_copy_coeffs(head_group[ng]);
        at(nf)->plus_in_place(trans_func);
        delete trans_func;
      }
      //std::cout << "Result of Reynold's operator is ";
      //at(nf)->print(std::cout);
      //std::cout << "\n";
      delete tfunc;
    }
    Gram_Schmidt();
    return;
  }

  //*******************************************************************************************
  // Checks block_shape_matrix of BasisSet symmetry representation to see that there are
  // as many blocks as there are irreducible representations

  bool BasisSet::is_normal_basis_for(const SymGroup &head_group) {
    //First do some basic checks to ensure problem is well-defined
    if(m_basis_symrep_ID.empty()) {
      get_symmetry_representation(head_group);
    }
    if(m_basis_symrep_ID.empty()) {
      std::cerr << "CRITICAL ERROR: Inside BasisSet::is_normal_basis_for() and cannot calculate a valid SymGroup representation. Exiting...\n";
      exit(1);
    }
    if(!head_group.size())
      return true;

    SymGroupRep const &t_rep(head_group[0].master_group().representation(m_basis_symrep_ID));

    //Check that block-diagonalization matches number of irreps
    return t_rep.num_blocks(head_group) == (t_rep.num_each_real_irrep(head_group)).sum();

  }

  //*******************************************************************************************

  BasisSet BasisSet::calc_normal_basis(const SymGroup &head_group, Eigen::MatrixXd &trans_mat) const {
    if(m_basis_symrep_ID.empty()) {
      get_symmetry_representation(head_group);
    }
    if(m_basis_symrep_ID.empty()) {
      std::cerr << "CRITICAL ERROR: Inside BasisSet::calc_normal_basis() and cannot calculate a valid SymGroup representation. Exiting...\n";
      exit(1);
    }
    if(!head_group.size()) {
      std::cerr << "CRITICAL ERROR: Inside BasisSet::calc_normal_basis() and and passed empty SymGroup. Exiting...\n";
    }


    SymGroupRep const &t_rep(head_group[0].master_group().representation(m_basis_symrep_ID));
    trans_mat = t_rep.get_irrep_trans_mat(head_group);
    BasisSet normal_basis(transform_copy(trans_mat));
    normal_basis.m_basis_symrep_ID = (t_rep.coord_transformed_copy(trans_mat)).add_copy_to_master();
    return normal_basis;
  }

  //*******************************************************************************************

  void BasisSet::push_back(Function *new_func) {
    Array<Function *>::push_back(new_func);
    m_eval_cache.push_back(0.0);
    m_deval_cache.push_back(0.0);
    if(m_argument.size() == 0 && new_func != nullptr) {
      _set_arguments(new_func->argument_bases());
    }
  }

  //*******************************************************************************************

  Function *BasisSet::_linear_combination(const Eigen::VectorXd &coeffs) const {
    if(!size()) return nullptr;
    if(size() != coeffs.size()) {
      std::cerr << "FATAL ERROR: In BasisSet::_linear_combination, the number of basis functions \n"
                << "does not match the size of the coefficient vector. Exiting...\n";
      exit(1);
    }

    Function *combfunc(nullptr), *tfunc(nullptr);

    for(EigenIndex i = 0; i < coeffs.size(); i++) {
      if(almost_zero(coeffs[i])) continue;
      if(!combfunc) {
        combfunc = at(i)->copy();
        combfunc->scale(coeffs[i]);
        continue;
      }
      tfunc = at(i)->copy();
      tfunc->scale(coeffs[i]);
      combfunc->plus_in_place(tfunc);
      delete tfunc;
    }
    if(!combfunc) {
      combfunc = at(0)->copy();
      combfunc->scale(0.0);
    }
    return combfunc;
  }

  //*******************************************************************************************

  /// Essentially, perform a change of basis on BasisSet as defined by trans_mat.
  /// Returns a BasisSet whos elements are linear combinations of the original BasisSet.
  /// The linear combinations are specified by the ROWS of trans_matx

  BasisSet BasisSet::transform_copy(const Eigen::MatrixXd &trans_mat) const {
    BasisSet copy_basis;
    if(trans_mat.cols() != size()) {
      std::cerr << "In BasisSet::transform_copy(), attempting to transform basis with a transformation\n"
                << "matrix that has incompatible number of columns (has " << trans_mat.cols() << " and needs " << size() << ").\n"
                << "Exiting...\n";
      exit(1);
    }
    for(EigenIndex nc = 0; nc < trans_mat.rows(); nc++) {
      copy_basis.push_back(_linear_combination(trans_mat.row(nc)));
    }

    // We should also transform the SymGroupRep, but we only have the ID, not the MasterSymGroup

    return copy_basis;
  }

  //*******************************************************************************************
  //TODO: Transform functions with mixed dependency layers
  BasisSet &BasisSet::apply_sym(const SymOp &op, int _dependency_layer /* = 1 */) {
    int this_dep_layer = dependency_layer();
    if(this_dep_layer == _dependency_layer) {
      for(Index i = 0; i < size(); i++)
        at(i)->apply_sym_coeffs(op);
    }
    else {
      for(Index i = 0; i < m_argument.size(); i++) {
        m_argument[i]->apply_sym(op, _dependency_layer);
      }
    }
    return *this;
  }

  //*******************************************************************************************

  void BasisSet::_set_arguments(const Array<BasisSet const *> &new_args) {
    if(m_argument.size()) {
      std::cerr << "CRITICAL ERROR: In BasisSet::_set_arguments(), cannot reset arguments of already-initialized BasisSet.\n"
                << "                Exiting...\n";
      exit(1);
    }
    m_argument.reserve(new_args.size());
    for(Index i = 0; i < new_args.size(); i++)
      m_argument.push_back(new_args[i]->shared_copy());
  }

  //*******************************************************************************************
  // Does modified Gram-Schmidt procedure, using tolerance checking for increased speed and stability
  // In worst case, this requires N*(N-1)/2 binary operations

  // In future, may wish to use alternative approach: 1) find Gram matrix G, where G(i,j)=at(i)->dot(at(j));
  //                                                  2) find orthonormal eigenvectors of G, forming columns of the matrix V
  //                                                  3) take the linear combination V.transpose()*(*this)
  // This alternate approach may be more numerically stable, but probably results in less sparse representations, unless there
  // is a way to compute an optimally sparse V matrix

  bool BasisSet::Gram_Schmidt() {
    bool is_unchanged(true);
    Index i, j;
    double tcoeff;
    Function *tfunc(nullptr);

    // loop over functions
    for(i = 0; i < size(); i++) {
      at(i)->small_to_zero(2 * TOL);

      tcoeff = sqrt(at(i)->dot(at(i)));

      if(tcoeff < TOL) {
        is_unchanged = false;
        delete at(i);
        remove(i);
        i--;
        continue;
      }
      else if(!almost_zero(tcoeff - 1.0)) {
        is_unchanged = false;
        at(i)->scale(1.0 / tcoeff);
      }

      // loop from i+1 to end and subtract projection of Function i onto Function j
      for(j = i + 1; j < size(); j++) {
        tcoeff = (at(i)->dot(at(j)));
        if(almost_zero(tcoeff)) {
          continue;
        }

        is_unchanged = false;

        tfunc = at(i)->copy();


        if(!almost_zero(tcoeff - 1)) {
          tfunc->scale(tcoeff);
        }

        at(j)->minus_in_place(tfunc);

        delete tfunc;
      }

    }
    if(!is_unchanged) m_basis_symrep_ID = SymGroupRepID();
    return is_unchanged;

  }

  //*******************************************************************************************

  bool BasisSet::Gaussian_Elim() {
    bool is_unchanged(true);
    Index i(0), j(0);
    Index j_min, i_min, i_temp;
    double tcoeff, min_coeff;
    Function *tfunc;
    while(i < size()) {
      j_min = Index(-1);
      for(i_temp = i; i_temp < size(); i_temp++) {
        tcoeff = at(i_temp)->leading_coefficient(j);
        if(almost_zero(tcoeff)) {
          delete at(i_temp);
          remove(i_temp);
          i_temp--;
          is_unchanged = false;
          continue;
        }
        if(!valid_index(j_min) || j < j_min || (j == j_min && std::abs(tcoeff) > std::abs(min_coeff))) {
          j_min = j;
          i_min = i_temp;
          min_coeff = tcoeff;
        }
      }
      if(i >= size()) break;

      j = j_min;
      if(i != i_min) {
        is_unchanged = false;
        swap_elem(i, i_min);
      }
      if(!almost_zero(min_coeff - 1)) {
        at(i)->scale(1.0 / min_coeff);
      }

      for(i_temp = 0; i_temp < size(); i_temp++) {
        if(i_temp == i) continue;
        tcoeff = at(i_temp)->get_coefficient(j);
        if(almost_zero(tcoeff)) continue;
        is_unchanged = false;
        if(almost_zero(tcoeff - 1)) {
          at(i_temp)->minus_in_place(at(i));
          continue;
        }
        tfunc = at(i)->copy();
        tfunc->scale(tcoeff);
        at(i_temp)->minus_in_place(tfunc);
        delete tfunc;
      }
      i++;
    }
    for(i = 0; i < size(); i++) {
      at(i)->small_to_zero(2 * TOL);
    }

    if(!is_unchanged) m_basis_symrep_ID = SymGroupRepID();

    return is_unchanged;

  }
  //*******************************************************************************************
  void BasisSet::get_symmetry_representation(const SymGroup &head_group) const {
    if(!head_group.size() || !head_group[0].has_valid_master()) return;

    m_basis_symrep_ID = head_group.add_empty_representation();
    Function *tfunct(nullptr);
    Eigen::MatrixXd tRep(size(), size());

    for(Index ng = 0; ng < head_group.size(); ng++) {
      //Get representation for operation head_group[ng]
      //store it in matrix tRep
      for(Index nb1 = 0; nb1 < size(); nb1++) {

        tfunct = at(nb1)->sym_copy_coeffs(head_group[ng]);
        for(Index nb2 = 0; nb2 < size(); nb2++) {
          tRep(nb2, nb1) = tfunct->dot(at(nb2));
        }
        delete tfunct;
      }

      //We have the representation of operation head_group[ng]
      //Make a new SymOpRepresentation out of it and push it back
      head_group[ng].set_rep(m_basis_symrep_ID, SymMatrixXd(tRep));
    }

  }

  //*******************************************************************************************

  bool BasisSet::make_orthogonal_to(const BasisSet &ortho_basis) {
    bool ortho_flag(true);
    for(Index i = 0; i < ortho_basis.size(); i++) {
      ortho_flag = make_orthogonal_to(ortho_basis[i]) && ortho_flag;

    }
    return ortho_flag;
  }

  //*******************************************************************************************

  bool BasisSet::make_orthogonal_to(Function const *ortho_func) {

    bool ortho_flag(true);
    if(!size()) {
      return ortho_flag;
    }

    double tcoeff;
    Function *tfunc(ortho_func->copy());

    //std::cout << "Ortho_func: " << tfunc->tex_formula() << "\n\n";
    tcoeff = tfunc->dot(tfunc);
    //std::cout << "Squared magnitude is " << tcoeff << "\n";

    if(almost_zero(tcoeff)) {
      delete tfunc;
      return ortho_flag;
    }

    if(!almost_zero(tcoeff - 1.0)) {
      tfunc->scale(1.0 / sqrt(tcoeff));
    }
    for(Index i = 0; i < size(); i++) {
      //std::cout << "Func " << i << ":  " << at(i)->tex_formula() << "\n\n";
      tcoeff = (at(i)->dot(tfunc));
      //std::cout << "Projection: " << tcoeff << "\n\n";
      if(almost_zero(tcoeff)) {
        continue;
      }
      ortho_flag = false;
      //You're changing the BasisSet, so the representation is no longer useable!
      m_basis_symrep_ID = SymGroupRepID();

      tfunc->scale(tcoeff);
      at(i)->minus_in_place(tfunc);
      tfunc->scale(1.0 / tcoeff);
      at(i)->normalize();

      //std::cout << "Result is " << at(i)->tex_formula() << "\n\n";
    }

    delete tfunc;


    return ortho_flag;
  }


  //*******************************************************************************************
  //** jsonParser stuff - BasisSet
  //*******************************************************************************************

  jsonParser &BasisSet::to_json(jsonParser &json) const {

    json.put_obj();

    // class BasisSet: public Array<Function *>
    json["basis_functions"].put_array();
    for(Index i = 0; i < size(); i++) {
      json["basis_functions"].push_back(at(i));
    }

    // mutable int m_basis_symrep_ID;
    json["m_basis_symrep_ID"] = m_basis_symrep_ID;

    // Array<BasisSet> subspaces;
    //json["subspaces"] = m_subspaces;

    return json;
  }

  //*******************************************************************************************

  /*
  void BasisSet::from_json(const jsonParser &json) {

    // no reading BasisSet for now

  }
  */


  //*******************************************************************************************

  jsonParser &to_json(const BasisSet &bset, jsonParser &json) {
    return bset.to_json(json);
  }

  //*******************************************************************************************
  /*
  // No reading functions for now
  void from_json(BasisSet &bset, const jsonParser &json) {
    return bset.from_json(json);
  }
  */

  //*******************************************************************************************
  BasisSet operator*(const SymOp &LHS, const BasisSet &RHS) {
    return BasisSet(RHS).apply_sym(LHS);
  }

}

