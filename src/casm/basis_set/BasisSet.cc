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
      m_dof_subbases[ID_ind].reserve(m_dof_subbases[ID_ind].size() + RHS.m_dof_subbases[i].size());
      for(Index j = 0; j < RHS.m_dof_subbases[i].size(); j++) {
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
    if(RHS.min_poly_order() > 0)
      _min_poly_constraints().push_back(PolyConstraint(Array<Index>::sequence(size(), size() + RHS.size() - 1), RHS.min_poly_order()));
    if(valid_index(RHS.max_poly_order()))
      _max_poly_constraints().push_back(PolyConstraint(Array<Index>::sequence(size(), size() + RHS.size() - 1), RHS.max_poly_order()));
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

  int BasisSet::register_remotes(const std::string &dof_name, const Array<DoF::RemoteHandle> &remote_handles) {
    int sum(0);
    for(Index i = 0; i < size(); i++)
      sum += at(i)->register_remotes(dof_name, remote_handles);

    for(Index i = 0; i < m_argument.size(); i++)
      sum += m_argument[i]->register_remotes(dof_name, remote_handles);

    return sum;
  }

  //*******************************************************************************************
  //void BasisSet::construct_polynomials_by_order(const Array<BasisSet const * > &tsubs, Index order) {
  //  std::cout << "beginning of construct_polynomials_by_order" << std::endl;
  //  _set_arguments(tsubs);
  //  PolynomialFunction tpoly(m_argument);
  //  Array<Index> curr_expon(tpoly.num_args(), 0);

  //  typedef IsoCounter<Array<Index> > OrderCount;
  //  typedef MultiCounter<OrderCount> ExpCount;

  //  OrderCount order_count(0, order, 1, order);

  //  ExpCount exp_count;
  //  for(Index i = 0; i < tsubs.size(); i++) {
  //      std::cout << " tsubs.size() : " << tsubs.size() << std::endl;
  //      std::cout << "  i  ::  " << i << std::endl;
  //      std::cout << "in first loop" << std::endl;
  //    exp_count.push_back(OrderCount(Array<Index>(tsubs[i]->size(), 0), Array<Index>(tsubs[i]->size(), order), 1, order_count[i]));
  //      std::cout << "after push_back" << std::endl;
  //    curr_expon.append(exp_count.back());
  //      std::cout << "after " << std::endl;
  //  }

  //  for(; order_count.valid(); ++order_count) {
  //      std::cout << "in second loop" << std::endl;
  //    for(Index i = 0; i < exp_count.size(); i++)
  //      exp_count[i].set_sum_constraint(order_count[i]);
  //    for(; exp_count.valid(); ++exp_count) {
  //      //bool valid_expon = true;
  //      //for(Index i = 0; i < exp_count.size() && valid_expon; i++)
  //      // valid_expon = tsubs[i]->satisfies_exponent_constraints(exp_count[i]());

  //      //if(!valid_expon)
  //      //  continue;

  //      Index ne = 0;
  //      ExpCount::const_value_iterator it(exp_count.value_begin()), it_end(exp_count.value_end());
  //      for(; it != it_end; ++it) {
  //        curr_expon[ne++] = *it;
  //      }
  //      //if(!satisfies_exponent_constraints(curr_expon))
  //      //continue;
  //      //else:
  //      //std::cout << "Adding exponent " << curr_expon << "\n\n";

  //      PolyTrie<double> new_trie(curr_expon.size());
  //      new_trie(curr_expon) = 1.0;
  //      push_back(tpoly.copy(new_trie));

  //    }
  //  }
  //  std::cout << "End of construct_polynomials by order\n";

  //}

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
    std::cout << "Constructing invariant polynomials from DoFs:\n";
    for(Index i = 0; i < tsubs.size(); ++i) {
      if((tsubs[i]->basis_symrep_ID()).empty()) {
        tsubs[i]->get_symmetry_representation(head_group);
      }

      std::cout << tsubs[i]->name() << "  ";

    }
    std::cout << "\n\n";
    std::cout << "dof_IDs are: " << dof_IDs() << "\n\n";
    std::cout << "dof_sub_bases are: " << dof_sub_bases() << "\n\n";
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
        std::cout << "dof_IDs of tsub " << j << " are " << tsubs[j]->dof_IDs() << "\n";
        Index ID_ind = (tsubs[j]->dof_IDs()).find(dof_IDs()[i]);
        if(ID_ind == (tsubs[j]->dof_IDs()).size())
          continue;
        const SubBasis &sub_basis(tsubs[j]->dof_sub_basis(ID_ind));
        for(Index b = 0; b < sub_basis.size(); b++) {
          big_sub_basis.push_back(sub_basis[b] + offset);
        }
        offset += tsubs[j]->size();
      }
      std::cout << "Adding min constraint: " << big_sub_basis << ":  " << min_dof_order << "\n";
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
    std::cout << "order_count from: " << initial_order << " to final order: " << final_order << "\n\n";
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
    accept(SubExpressionLabeler("CCD", "Q_%n[%f]"));
    std::cout << "Result, with " << size() << " invariant functions:\n";
    for(Index i = 0; i < size(); i++)
      std::cout << "F_" << i << " = " << at(i)->tex_formula() << "\n\n";
    return;
  }

  //*******************************************************************************************

  void BasisSet::construct_green_lagrange_dot_prods(const Array<BasisSet const *> &site_disp_dofs, BasisSet const *LG_strain_dofs, const Eigen::MatrixXd &ref_clust) {
    Index num_sites(site_disp_dofs.size());

    Array<BasisSet const *> all_bset(site_disp_dofs);
    all_bset.push_back(LG_strain_dofs);
    _set_arguments(all_bset);
    PolynomialFunction tpoly(m_argument);

    // Initialize strain_inds to number of non-strain DoFs, and then increment to correct values
    Array<Array<Index> > strain_ind(3, Array<Index>(3, num_sites * 3));
    strain_ind[0][0] += 0;
    strain_ind[1][1] += 1;
    strain_ind[2][2] += 2;
    strain_ind[1][2] += 3;
    strain_ind[2][1] += 3;
    strain_ind[0][2] += 4;
    strain_ind[2][0] += 4;
    strain_ind[0][1] += 5;
    strain_ind[1][0] += 5;


    Array<Index> exp_count(tpoly.poly_coeffs().depth(), 0);
    for(Index i = 0; i < num_sites; i++) {
      for(Index j = i; j < num_sites; j++) {
        PolyTrie<double> proto_trie(tpoly.poly_coeffs());
        for(Index a = 0; a < 3; a++) {

          exp_count[a + j * 3]++;                   //disp(a,j)

          //ref_clust(a,i)*disp(a,j)
          proto_trie(exp_count) += ref_clust(a, i);

          exp_count[a + i * 3]++;                   //disp(a,i)

          //disp(a,i)*disp(a,j)
          proto_trie(exp_count) += 1.0;

          exp_count[a + j * 3]--;                   //~disp(a,j)

          //ref_clust(a,j)*disp(a,i)
          proto_trie(exp_count) += ref_clust(a, j);

          exp_count[a + i * 3]--;                   //~disp(a,i)

          for(Index b = 0; b < 3; b++) {
            double norm_factor(a == b ? 1.0 : 1.0 / sqrt(2));
            exp_count[strain_ind[a][b]]++;          //E(a,b)

            exp_count[b + j * 3]++;                 //disp(b,j)

            //2*ref_clust(a,i)*E(a,b)*disp(b,j)
            proto_trie(exp_count) += 2.0 * norm_factor * ref_clust(a, i);

            exp_count[a + i * 3]++;                 //disp(a,i)

            //2*disp(a,i)*E(a,b)*disp(b,j)
            proto_trie(exp_count) += 2.0;

            exp_count[b + j * 3]--;                 //~disp(b,j)

            //2*disp(a,i)*E(a,b)*ref_clust(b,j)
            proto_trie(exp_count) += 2.0 * norm_factor * ref_clust(b, j);

            exp_count[a + i * 3]--;                 //~disp(a,i)

            //2*ref_clust(a,i)*E(a,b)*ref_clust(b,j)
            proto_trie(exp_count) += 2.0 * norm_factor * ref_clust(a, i) * ref_clust(b, j);

            exp_count[strain_ind[a][b]]--;          //~E(a,b)

          }//\End loop over 'b'
        }//\End loop over 'a'

        push_back(new PolynomialFunction(tpoly, proto_trie));
      }
    }
  }


  //*******************************************************************************************

  void BasisSet::construct_disp_grad_dot_prods(const Array<BasisSet const *> &site_disp_dofs, BasisSet const *F_strain_dofs, const Eigen::MatrixXd &ref_clust) {
    Array<BasisSet const *> all_bset(site_disp_dofs);
    all_bset.push_back(F_strain_dofs);

    _set_arguments(all_bset);

    Index num_sites(site_disp_dofs.size());
    PolynomialFunction tpoly(m_argument);

    // ref_covar is a matrix of dot products -> ref_covar(i,j)=site[i].coord.dot(site[j].coord)
    Eigen::MatrixXd ref_covar(ref_clust.transpose()*ref_clust);

    // delta_covar = covar-ref_covar =
    //             = (ref_clust.transpose()*F.transpose()+displacements.transpose())*(F*ref_clust+displacements)-ref_clust.transpose()*ref_clust.transpose()
    //             = RC'*F'*F*RC + RC'*F'*disp + disp'*F*RC + disp'*disp - RC'*RC  // RC == ref_clust
    Array<Index> exp_count(tpoly.poly_coeffs().depth(), 0);
    for(Index i = 0; i < num_sites; i++) {
      for(Index j = i; j < num_sites; j++) {
        PolyTrie<double> proto_trie(tpoly.poly_coeffs());

        proto_trie(exp_count) += ref_covar(i, j);
        for(Index a = 0; a < 3; a++) {

          exp_count[a + j * 3]++;             //disp(a,j)
          exp_count[a + i * 3]++;             //disp(a,i)

          //disp(a,i)*disp(a,j)
          proto_trie(exp_count) += 1.0;

          exp_count[a + j * 3]--;             //~disp(a,j)
          exp_count[a + i * 3]--;             //~disp(a,i)

          for(Index b = 0; b < 3; b++) {

            exp_count[3 * num_sites + 3 * b + a]++; //F(b,a)
            exp_count[b + j * 3]++;           //disp(b,j)

            // ref_clust(a,i)*F(b,a)*disp(b,j)
            proto_trie(exp_count) += ref_clust(a, i);

            exp_count[3 * num_sites + 3 * b + a]--; //~F(b,a)
            exp_count[b + j * 3]--;	        //~disp(b,j)

            exp_count[3 * num_sites + 3 * a + b]++; //F(a,b)
            exp_count[a + i * 3]++;           //disp(a,i)

            // disp(a,i)*F(a,b)*ref_clust(b,j)
            proto_trie(exp_count) += ref_clust(b, j);

            exp_count[3 * num_sites + 3 * a + b]--; //~F(a,b)
            exp_count[a + i * 3]--;           //~disp(a,i)

            for(Index c = 0; c < 3; c++) {
              exp_count[3 * num_sites + 3 * b + a]++; //F(b,a)
              exp_count[3 * num_sites + 3 * b + c]++; //F(b,c)

              // ref_clust(a,i)*F(b,a)*F(b,c)*ref_clust(c,j)
              proto_trie(exp_count) += ref_clust(a, i) * ref_clust(c, j);

              exp_count[3 * num_sites + 3 * b + a]--; //~F(b,a)
              exp_count[3 * num_sites + 3 * b + c]--; //~F(b,c)
            }

          }//\End loop over 'b'
        }//\End loop over 'a'

        push_back(new PolynomialFunction(tpoly, proto_trie));
      }
    }
  }


  //*******************************************************************************************
  void BasisSet::construct_quadratic_ccds(const Array<BasisSet const *> &site_disp_dofs, BasisSet const *strain_vars, const Eigen::MatrixXd &ref_clust,
                                          const Array<SymOp>::X2 &equiv_map, const SymGroupRep &permute_group, double sigma) {

    //Array<BasisSet const *> all_bset(site_disp_dofs);
    //all_bset.push_back(strain_vars);
    //_set_arguments(all_bset);

    if(!equiv_map.size() || !equiv_map[0].size()) {
      std::cerr << "CRITICAL ERROR: In BasisSet::construct_quadratic_ccds(), passed a malformed equivalence_map Array.\n"
                << "                Exiting...\n";
      exit(1);
    }
    std::cout << "Inside construct_quadratic_ccds! Cluster coordinates are:\n" << ref_clust << "\n\n";

    m_name = "CCD";

    // Construct basis from dot products: coord_i.dot(coord_j);
    // ultimate goal is to go from the 'covar' basis to the symmetrized ccd basis, which is simply a
    // (rank-deficient) linear transformation such that
    //    *this = covar2symccd * delta_covar_basis, where covar2symccd is the transformation matrix
    BasisSet delta_covar_basis;

    if((strain_vars->name()) == "STRAIN_GL") {
      delta_covar_basis.construct_green_lagrange_dot_prods(site_disp_dofs, strain_vars, ref_clust);
    }
    else if((strain_vars->name()) == "STRAIN_F") {
      delta_covar_basis.construct_disp_grad_dot_prods(site_disp_dofs, strain_vars, ref_clust);
    }

    // The first row of equiv_map is the clust_group
    SymGroup head_group(equiv_map[0]);

    //Inspect the site_args and global_args
    Index num_sites(site_disp_dofs.size());

    // upper triangle of a NxN matrix, including diagonal
    Index num_prods((num_sites * (num_sites + 1)) / 2);

    // upper triangle of a NxN matrix, excluding diagonal <-- number of distinct (i,j) pairs of cluster sites
    Index num_dists((num_sites * (num_sites - 1)) / 2);

    // ref_covar is a matrix of dot products -> ref_covar(i,j)=site[i].coord.dot(site[j].coord)
    Eigen::MatrixXd ref_covar(ref_clust.transpose()*ref_clust);

    // The covariance matrix of the *CURRENT* deformation state is specified by variables in site_disp_dofs.
    // To ease explanation, we define the following concepts:
    //  - The covariance matrix of the cluster in its 'current' state:
    //          Eigen::Matrix covar(curr_clust.transpose()*curr_clust)
    //
    //  - The covariance matrix of the cluster in its 'current' state:
    //          Eigen::Matrix ref_covar(ref_clust.transpose()*ref_clust);
    //
    //  - The difference of the covariance matrices
    //          Eigen::Matrix delta_covar(covar-ref_covar);

    // dist2covar converts vector
    //   v = {dist_sqrd(0,1), dist_sqrd(0,2), ..., dist_sqrd(0,N), dist_sqrd(1,2), ..., dist_sqrd(N-1,N)}
    // to vector
    //   w = {covar(0,0), covar(0,1), ..., covar(0,N), covar(1,1), ..., covar(N,N)}
    // via w = dist2covar*v
    //   where dist_sqrd(i,j) is the square of the distance from site 'i' to site 'j'
    // dist2covar is equally able to convert the vector
    //   delta_v = {delta_dist_sqrd(0,1), delta_dist_sqrd(0,2), ..., delta_dist_sqrd(0,N), delta_dist_sqrd(1,2), ..., delta_dist_sqrd(N-1,N)}
    // to vector
    //   delta_w = {delta_covar(0,0), delta_covar(0,1), ..., delta_covar(0,N), delta_covar(1,1), ..., delta_covar(N,N)}
    Eigen::MatrixXd dist2covar(Eigen::MatrixXd::Zero(num_prods, num_dists)); //<-- Any xxxx2covar transformation is rank-deficient & can't be inverted directly

    //Array<Array< Index> > prod_inds(num_prods, Array<Index>(2, 0));
    std::cout << "Delta covar functions are:\n";
    delta_covar_basis.accept(VariableLabeler("%s[%n]"));
    Index c(0);
    for(Index i = 0; i < num_sites; i++) {
      for(Index j = i; j < num_sites; j++) {
        std::cout  << "dcovar(" << i << ", " << j << ") = " << delta_covar_basis[c]->formula() << "\n";
        c++;
      }
    }
    std::cout << "\n";

    for(Index i = 0; i < num_sites; i++) {
      for(Index j = i + 1; j < num_sites; j++) {
        //prod_inds[i * num_sites - (i * (i + 1)) / 2 + j][0] = i;
        //prod_inds[i * num_sites - (i * (i + 1)) / 2 + j][1] = j;

        // Fill dist2covar for the {i,j} pair using
        // d_{ij}^2 = covar(i,i)+covar(j,j)-2*covar(i,j)
        dist2covar(i * num_sites - (i * (i + 1)) / 2 + i, i * num_sites - (i * (i + 3)) / 2 + j - 1) = 1.0;  //*{covar(i,i), dist(i,j)}
        dist2covar(j * num_sites - (j * (j + 1)) / 2 + j, i * num_sites - (i * (i + 3)) / 2 + j - 1) = 1.0;  //*{covar(j,j), dist(i,j)}
        dist2covar(i * num_sites - (i * (i + 1)) / 2 + j, i * num_sites - (i * (i + 3)) / 2 + j - 1) = -2.0; //*{covar(i,j), dist(i,j)}
      }
    }

    std::cout << "dist2covar is:\n" << dist2covar << "\n\n";

    // matrix of scalar products
    //    gram_mat(i,j) = scalar_prod(delta_w[i],delta_w[j]);
    // The scalar product is defined as a weighted integral over all the displacement degrees of freedom
    Eigen::MatrixXd gram_mat(Eigen::MatrixXd::Zero(num_prods, num_prods));

    // construct gram matrix for delta_covar_basis.
    // gram_mat(i,j)=delta_covar_basis[i]->scalar_prod(delta_covar_basis[j]); <- We have to construct this by hand, since
    //   Function::scalar_prod() assumes Frobenius inner product in most cases
    // The following loop expands delta_covar(i,j)*delta_covar(i2,j2) component-wise and
    // calculates the gaussian overlap integral:
    //
    //  gram_mat({i,j},{i2,j2})=  \int_{-\infty}^{\infty} gaussian(sigma, ref_clust) * { ref_covar(i,j)*ref_covar(i2,j2) - covar(i,j)*ref_covar(i2,j2)
    //                                                                                 - ref_covar(i,j)*covar(i2,j2) + covar(i,j)*covar(i2,j2) }*
    //                                                                                 dx^{(1)}...dx^{(N)}dy^{(1)}...dy^{(N)}dz^{(1)}...dz^{(N)}
    //
    // where gaussian(sigma, ref_clust) is a 3*N-dimensional Gaussian distribution centered at the coordinates of the reference cluster
    // and with std deviation 'sigma'. Its integral over all of phase space is one.

    // Since delta_covar_basis is the unrolled version of delta_covar, we loop over i<=j and i2<=j2
    Index k(0); // {i,j} matrix index maps to 'k' linear index
    for(Index i = 0; i < num_sites; i++) {
      for(Index j = i; j < num_sites; j++) {

        Index l(0); // {i2,j2} matrix index maps to 'l' linear index
        double gint1, gint2, g2int;

        Index n;
        Array<Array<Index> > var_inds(4, Array<Index>(2, 0));
        for(Index i2 = 0; i2 < num_sites; i2++) {
          for(Index j2 = i2; j2 < num_sites; j2++) {

            // This is the first term in the integral
            gram_mat(k, l) += ref_covar(i, j) * ref_covar(i2, j2);

            // Begin constructing second term in the integral arising from integrand:
            //   -covar(i,j)*ref_covar(i2,j2) - ref_covar(i,j)*covar(i2,j2)

            //initialize gint1 and gint2 with contributions from constant terms
            gint1 = ref_covar(i, j);
            gint2 = ref_covar(i2, j2);

            // add on 3 quadratic terms, if i==j or i2==j2, arising from integrating
            //        coord_i[x]*coord_i[x] + coord_i[y]*coord_i[y] + coord_i[z]*coord_i[z]
            if(i == j)
              gint1 += 3.0 * gaussian_moment(2, sigma);
            if(i2 == j2)
              gint2 += 3.0 * gaussian_moment(2, sigma);

            gram_mat(k, l) += -(gint1 * ref_covar(i2, j2) + ref_covar(i, j) * gint2);
            //\End construction of second term

            //\Begin construction of third term, arising from integrand:
            //        covar(i,j)*covar(i2,j2)

            // expand the product (a fourth-order monomial) into its components by
            // looping over cartesian coordinates 'a' of {i,j} and 'b' of {i2,j2}
            for(Index a = 0; a < 3; a++) {
              for(Index b = 0; b < 3; b++) {
                //var_inds[i][0] is the site index of variable 'i' in the monomial term
                //var_inds[i][1] is the cartesian index of variable 'i' in the monomial term
                var_inds[0][0] = i;
                var_inds[1][0] = j;
                var_inds[2][0] = i2;
                var_inds[3][0] = j2;
                var_inds[0][1] = a;
                var_inds[1][1] = a;
                var_inds[2][1] = b;
                var_inds[3][1] = b;

                // count up how many times each variable appears in the monomial term:
                //    e.g. x^{(1)}*x^{(1)}*y{(3)}*y^{(4)} vs x^{(1)}*x^{(1)}*x{(1)}*x^{(1)}
                //         vs x^{(1)}*x^{(1)}*z{(2)}*z^{(2)} all have different gaussian integrals
                g2int = 1.0;
                for(Index v = 0; v < var_inds.size(); v++) {
                  if(var_inds.find(var_inds[v]) != v) continue; // true if we've already taken care of this term
                  n = var_inds.incidences(var_inds[v]);
                  g2int *= gaussian_moment(n, sigma, ref_clust(var_inds[v][1], var_inds[v][0]));
                }
                gram_mat(k, l) += g2int;
              }
            }
            //\End third (final) term in the integral
            l++;
          }
        }
        k++;
      }
    }//\End constructing gram_mat for delta_covar_basis

    std::cout << "gram_mat for delta_covar_basis:\n" << gram_mat << "\n\n";

    //   using dist2covar we can convert gram_mat^{(covar)} to gram_mat^{(dist)} via
    //
    //      gram_mat^{(dist)} = dist2covar.transpose()*gram_mat^{(covar)}*dist2covar
    //
    //   which restricts it to a translationally-invariant supspace (of intersite distances)
    //   we will now find an orthogonal basis that diagonalizes gram_mat^{(dist)} and then
    //   scale the basis vectors to normalize them.  This will yield the matrix ccd2dist that satisfies
    //
    //        ccd2dist.transpose()*gram_mat^{(dist)}*ccd2dist = Identity <-- the definition of orthonormality w.r.t. our gaussian inner product
    //
    //   It's easy to see that ccd2dist = inverse(matrix_sqrt(gram_mat^{(dist)})), which we can calculate via SVD:
    //        gram_mat^{(dist)} = U*S*V.transpose() <--  NOTE: gram_mat^{(dist)} is symmetric, so U==V
    Eigen::JacobiSVD<Eigen::MatrixXd> tsvd(dist2covar.transpose()*gram_mat * dist2covar, Eigen::ComputeFullU);

    std::cout << "SVD of transformed gram_mat:\n ***U:\n" << tsvd.matrixU() << "\n\n ***S:\n" << tsvd.singularValues() << "\n\n";

    //   dist2ccd converts from delta_dist_sqrd to normalized ccds -- use this to find symmetry representation
    //   i.e., for vector 'delta_v' of delta_dist_sqrd, decompose it into vector 'Q' via
    //        Q = dist2ccd*delta_v
    //   dist2ccd is matrixSqrt(gram_mat^{(dist)})
    Eigen::MatrixXd sqrtS(Eigen::MatrixXd::Zero(num_dists, num_dists));
    sqrtS.diagonal() = tsvd.singularValues().cwiseSqrt();
    //Eigen::MatrixXd dist2ccd(tsvd.matrixU() * sqrtS * tsvd.matrixU().transpose()); <-- THIS IS WRONG
    Eigen::MatrixXd dist2ccd(sqrtS * tsvd.matrixU().transpose());

    //std::cout << "dist2ccd:\n" << dist2ccd  << "\n\n";

    //  ccd2dist.transpose()*gram_Mat^{(dist)}*ccd2dist = Identity <-- ccd basis diagonalizes the gram matrix
    Eigen::MatrixXd ccd2dist(dist2ccd.inverse());

    std::cout << "ccd2dist.transpose()*gram_Mat^{(dist)}*ccd2dist:\n"
              << ccd2dist.transpose()*dist2covar.transpose()*gram_mat *dist2covar *ccd2dist << "\n\n";

    //std::cout << "ccd2dist:\n" << ccd2dist  << "\n\n";

    //ccd2covar converts from covar to normalized ccds -- use this to construct final basis
    Eigen::MatrixXd ccd2covar(dist2covar * dist2ccd.inverse());

    //std::cout << "ccd2covar:\n" << ccd2covar  << "\n\n";

    //  Now we do any additional symmetrization that wasn't encoded into the inner product
    //  Diagonalizing the inner product partitioned the CCD vector space in a way that accounts for the topology
    //  of the cluster. This partitioning has to respect the symmetry of the cluster, but symmetry may allow further partitioning
    //  of the subspaces, or at least allow us to align the basis vectors along high-symmetry directions. It is important
    //  to do the partitioning in this order, since the topological partitioning resolves some ambiguities that arise
    //  when partitioning is done based on symmetry alone. To do this, we first need the ccd symgrouprep
    Array<Index> tperm(num_dists);
    SymGroupRep tccd_rep(head_group[0].master_group());
    //tccd_rep.resize(head_group[0].master_group().size(), nullptr);

    //  The individual site displacements don't know anything about permutation, so we can't use BasisSet::get_symmetry_representation()
    //  We will construct it directly from the permute_group of the cluster instead by inspecting how each {i,j} pair transforms, thus
    //  constructing the symgrouprep for delta_dist_sqrd.  We will then use transformation matrices to convert it to the symccd basis
    //  we could have started by making the covar_basis symrep, but then we would have had to invert a rank-deficient matrix (dist2covar)
    //  this seems easier, since Eigen doesn't currently have built-in pseudoinverse routines
    Index i2, j2;
    for(Index ng = 0; ng < equiv_map[0].size(); ng++) {
      for(Index i = 0; i < num_sites; i++) {
        for(Index j = i + 1; j < num_sites; j++) {

          i2 = ((permute_group.get_permutation(equiv_map[0][ng].index()))->perm_array()).find(i);
          j2 = ((permute_group.get_permutation(equiv_map[0][ng].index()))->perm_array()).find(j);
          if(i2 > j2)
            CASM::swap(i2, j2);

          // (i,j) pair goes to (i2,j2) pair
          //tperm[i2 * num_sites - (i2 * (i2 + 1)) / 2 + j2] = i * num_sites - (i * (i + 1)) / 2 + j;
          tperm[i2 * num_sites - (i2 * (i2 + 3)) / 2 + j2 - 1] = i * num_sites - (i * (i + 3)) / 2 + j - 1;

        }
      }
      //std::cout << "tperm " << ng << " is:  " << tperm << "\n\n";
      //std::cout << "which corresponds to matrix:\n" << (*(SymPermutation(tperm).get_MatrixXd())) << "\n\n";
      //std::cout << "ccd representation " << ng << ":\n" << dist2ccd * (*(SymPermutation(tperm).get_MatrixXd()))*dist2ccd.inverse() << "\n\n";
      tccd_rep.set_rep(equiv_map[0][ng], SymMatrixXd(dist2ccd * (*(SymPermutation(tperm).get_MatrixXd()))*dist2ccd.inverse()));

    }

    // get the representations for the rest of the factor_group using equiv_map
    // transforming by equiv_map[i][0] maps the cluster onto an equivalent but doesn't permute the sites
    // so tccd_rep can be constructed using the multiplication table
    for(Index ne = 1; ne < equiv_map.size(); ne++) { // loop over equivalents (excluding prototype)
      for(Index ns = 0; ns < equiv_map[ne].size(); ns++) { // loop over symops
        //std::cout << "equiv_map[" << ne << "][0].index() = " << equiv_map[ne][0].index() << ", equiv_map[0][" << ns << "].index() = " << equiv_map[0][ns].index()
        //<< ", and prod.index() = " << equiv_map[ne][0].ind_prod(equiv_map[0][ns]) << "\n";
        tccd_rep.set_rep(equiv_map[ne][0].ind_prod(equiv_map[0][ns]), *(tccd_rep[equiv_map[0][ns].index()]));
      }
    }

    // now we have the full representation

    // Use the singular values to partition CCDs into degenerate subspaces.  singly degenerate subspaces are *definitely* irreducible
    // 2-fold degeneracies and higher require additional consideration
    Array<Index> degeneracies(1, 1);
    for(Index i = 1; i < num_dists; i++) {
      if(!almost_equal(tsvd.singularValues()[i], tsvd.singularValues()[i - 1], 0.001 * CASM::min(tsvd.singularValues()[i], tsvd.singularValues()[i - 1]))) {
        degeneracies.push_back(1);
      }
      else {
        degeneracies.back()++;
      }
    }

    // we will symmetrize each subspace independently, and store the transformation matrix in block_trans_mat
    // Then, fill in appropriate block of ccd2symccd matrix using block_trans_mat
    Eigen::MatrixXd block_trans_mat, ccd2symccd(num_dists, num_dists);
    k = 0; // keep track of how many dimensions we have already symmetrized
    for(Index i = 0; i < degeneracies.size(); i++) {
      // Start out with block_trans_mat as a simple rank-reducing matrix that selects out the subspace associated with degeneracies[i]
      block_trans_mat = Eigen::MatrixXd::Zero(degeneracies[i], num_dists);
      block_trans_mat.block(0, k, degeneracies[i], degeneracies[i]) = Eigen::MatrixXd::Identity(degeneracies[i], degeneracies[i]);

      if(degeneracies[i] > 1) { // deal with degeneracies greater than 1
        SymGroupRep block_rep(tccd_rep.coord_transformed_copy(block_trans_mat));
        // find (degeneracies[i] x num_dists) matrix that symmetrizes the subspace
        Eigen::MatrixXd tmat(block_rep.get_irrep_trans_mat(head_group));
        block_trans_mat = tmat * block_trans_mat; //Eigen assumes aliasing for multiplication
      }
      //place block in matrix
      ccd2symccd.block(k, 0, degeneracies[i], num_dists) = block_trans_mat;
      k += degeneracies[i];
    }

    // ccd2covar*symccd2ccd = ccd2covar*ccd2symccd.inverse() = ccd2covar*ccd2symccd.transpose()
    Eigen::MatrixXd symccd2covar(ccd2covar * ccd2symccd.transpose());

    for(EigenIndex nc = 0; nc < symccd2covar.cols(); nc++) {
      push_back(delta_covar_basis._linear_combination(symccd2covar.col(nc)));
      _back()->accept(VariableLabeler("%s[%n]"));
      std::cout << "Q" << nc << " = " << back()->tex_formula() << "\n";
    }

    std::cout << "  symccd2covar:\n" << symccd2covar << "\n\n";

    std::cout << "  ccd2symccd:\n" << ccd2symccd << "\n\n";

    std::cout << "  dist2symccd:\n" << ccd2symccd *dist2ccd << "\n\n";

    //Finally add the symrep for symccds
    m_basis_symrep_ID = head_group[0].master_group().add_representation(tccd_rep.coord_transformed_copy(ccd2symccd));

    // and fill the subbases
    for(Index i = 0; i < m_dof_subbases.size(); ++i) {
      m_dof_subbases[i] = SubBasis::sequence(0, size() - 1);
    }
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
                                                          Index basis_ind) {

    m_argument.clear();
    if(!allowed_occs.is_locked()) {
      set_dof_IDs(Array<Index>(1, allowed_occs.ID()));
    }
    m_max_poly_order = 1;
    Index N = allowed_occs.size();
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
    Eigen::MatrixXd inv_sqrtE(Eigen::MatrixXd::Zero(N, N));
    for(Index i = 0; i < N; i++) {
      if(almost_zero(teig.eigenvalues()[i]) || teig.eigenvalues()[i] < 0) {
        // we could 'fix' gram_mat, but that could cause mysterious behavior - leave it to the user
        std::cerr << "CRITICAL ERROR: Passed a Gram Matrix to BasisSet::construct_orthonormal_discrete_functions that is not positive-definite.\n"
                  << "                Gram Matrix is:\n" << gram_mat << "\nExiting...\n";
        exit(1);
      }
      inv_sqrtE(i, i) = 1.0 / sqrt(teig.eigenvalues()[i]);
    }

    // B matrix is matrix square root of gram_mat.inverse(). Its columns form an orthonormal basis
    // wrt gram_mat but it is ill-suited for our purposes.
    Eigen::MatrixXd B(teig.eigenvectors()*inv_sqrtE * teig.eigenvectors().inverse());

    // ** step 2: Make seed basis. This will be used to seed optimized transformation of 'B'
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
      OccupantFunction tOF(allowed_occs, B.col(i), size(), basis_ind, allowed_occs.sym_rep_ID());

      push_back(tOF.copy());
    }
    // ********* TODO!!!!!!!!! ********
    // Calculate BasisSet symmetry representation here, based on allowed_occs.sym_rep_ID() && B matrix
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
                                                          Index basis_ind) {

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
    construct_orthonormal_discrete_functions(allowed_occs, gram_mat, basis_ind);
  }


  //*******************************************************************************************


  //void BasisSet::construct_harmonic_polynomials_old(const Array<BasisSet const *> &rotation_dofs, Index order, Index min_order){
  //  //rotation_dofs contains n-continuous variables
  //  //goal is to make an orthogonal basis on an nd-sphere
  //  // i.e. must satisfy the constraint that a^2+b^2+c^2+ ... = 1
  //  //BasisSet working_dofs(rotation_dofs[0]);
  //
  //  Array<BasisSet const *> abset(rotation_dofs);
  //  PolynomialFunction tpoly(abset);
  //  PolyTrie<double> mytrie(tpoly.poly_coeffs());
  //  Array<Index> expon(tpoly.poly_coeffs().depth(),0);
  //  Index ndofs(tpoly.poly_coeffs().depth());
  //  Array<BasisSet> myarray;
  //  BasisSet trivial;

  //  //std::cout << "working_dofs.size()" << working_dofs.size() << std::endl;

  //  // Always add first and second order
  //  BasisSet zeroth;
  //  // The zeroth order will always be 1.
  //  //for(Index i = 0; i < expon.size(); i++){
  //  //    expon[i]=0;
  //  //    mytrie(expon)=1;
  //  //}
  //  //PolynomialFunction zeroth_shell(abset);
  //  mytrie(expon) = 1.0;
  //  PolynomialFunction zero(abset,mytrie);
  //  //zero.prune_zeros();
  //  //std::cout << "Initialized zeroth trie formula " <<  zeroth.tex_formula() <<  std::endl;
  //  zeroth.push_back(zero.copy());
  //  myarray.push_back(zeroth);
  //  trivial.push_back(zero.copy());

  //  BasisSet firsth;
  //  // Always want the first order of each rotation
  //  for(Index k = 0; k < expon.size(); k++){
  //      PolyTrie<double> tmptrie(expon.size());
  //          for(Index i = 0; i < expon.size(); i++){
  //              if(i==k){
  //              expon[i]=1;
  //              tmptrie(expon)=1;
  //              expon[i]=0;
  //          }
  //      }
  //     //tmptrie.print_sparse(std::cout);
  //     //PolynomialFunction first_shell(rotation_dofs);
  //     PolynomialFunction first_order(abset,tmptrie);
  //     //first_order.prune_zeros();
  //      firsth.push_back(first_order.copy());
  //     trivial.push_back(first_order.copy());
  //  }
  //  myarray.push_back(firsth);


  //  //std::cout << "TRIVIAL" << std::endl;
  //  //for(Index i = 0; i < trivial.size(); i++) {

  //  //  std::cout << "   V" << i + 1 << " = " << trivial[i]->tex_formula() << '\n';
  //  //}

  //  //probably always want the radial function
  //  // need to figure out how to raise it to power
  //  // initial thought is to use IsoCounter to build up a polytrie
  //  //
  //  PolyTrie<double> rtrie(expon.size());
  //  for(Index i = 0; i < expon.size(); i++){
  //      expon[i]=2;
  //      rtrie(expon)=1;
  //      expon[i]=0;
  //  }
  //  PolynomialFunction radial_func(rotation_dofs,rtrie);


  //  // Here we need to start an iterative process to make as many polynomials
  //  // as desried from the passed 'Index order'
  //
  //  // co = CurrentOrder
  //  for(Index co = 2; co < order+1; co++){

  //      //First want to get the basis that is already orthogonal
  //      // Add all functions in myarray as well as all P_{i}^{m} * r^{n-m}/2 for n-m even
  //      BasisSet rhs_basis;
  //      BasisSet ntriv;
  //      BasisSet xbasis;

  //      // bs = BasisSet
  //      for(Index bs = 0; bs<myarray.size(); bs++){
  //          // bf = BasisFunction
  //          for(Index bf = 0; bf<myarray[bs].size(); bf++){
  //              rhs_basis.push_back(myarray[bs][bf]->copy());
  //              ntriv.push_back(myarray[bs][bf]->copy());


  //              //if((bs == 0 || (co-bs)%2==0) && (co-bs)!=0 ){
  //              if((co-bs)%2==0 ){
  //
  //                   // i.e. if n-m is even
  //                   // now we need to make radial function
  //                   PolyTrie<double> rtrie(ndofs);
  //                   IsoCounter<Array<Index> > ex_vecs(Array<Index>(ndofs,0), Array<Index>(ndofs,co),2,co-bs);
  //                   do{
  //                      rtrie(ex_vecs())=1;
  //                   } while(++ex_vecs);

  //                   PolynomialFunction rfunc(abset,rtrie);
  //                   if(bs==0){
  //                       xbasis.push_back(rfunc.copy());
  //                   }
  //                   else {
  //                       Function* p_times_r;
  //                       //rfunc.prune_zeros();
  //                       p_times_r = rfunc.multiply(myarray[bs][bf]);
  //                       //myarray[bs][bf]->prune_zeros();
  //                       xbasis.push_back(p_times_r->copy());
  //                       }
  //              } // end if
  //          }
  //      }

  //      for(Index xx = 0; xx < xbasis.size();xx++){
  //          rhs_basis.push_back(xbasis[xx]->copy());
  //      }
  //      //std::cout << " what is my rhs_basis??????? See below: \n";
  //      //for(Index ww = 0; ww<rhs_basis.size();ww++){
  //      //    std::cout << "RHS" << ww << ": " << rhs_basis[ww]->tex_formula() << std::endl;
  //      //}
  //
  //      BasisSet lhs_basis;
  //      lhs_basis.construct_polynomials_by_order(abset,co);
  //      //std::cout << " what is my lhs_basis??????? See below: \n";
  //      //for(Index ww = 0; ww<lhs_basis.size();ww++){
  //      //    std::cout << "LHS" << ww << ": " << lhs_basis[ww]->tex_formula() << std::endl;
  //      //}

  //      //lhs_basis.make_orthogonal_to(ntriv);
  //      lhs_basis.make_orthogonal_to(rhs_basis);
  //      //std::cout << " lhs_basis after make_orthogonal_to See below: \n";
  //      //for(Index ww = 0; ww<lhs_basis.size();ww++){
  //      //    std::cout << "LHS" << ww << ": " << lhs_basis[ww]->tex_formula() << std::endl;
  //      //}
  //      lhs_basis.Gram_Schmidt();
  //      //std::cout << " lhs_basis after Gram_Schmidt See below: \n";
  //      //for(Index ww = 0; ww<lhs_basis.size();ww++){
  //      //    std::cout << "LHS" << ww << ": " << lhs_basis[ww]->tex_formula() << std::endl;
  //      //}
  //      myarray.push_back(lhs_basis);
  //  }

  //
  //  for(Index bs = 0; bs < myarray.size(); bs++){
  //      for(Index bf = 0; bf < myarray[bs].size(); bf++){
  //          push_back(myarray[bs][bf]->copy());
  //      }
  //  }

  //  Gram_Schmidt();
  //
  //}

  void BasisSet::construct_harmonic_polynomials(const Array<BasisSet const *> &rotation_dofs, Index order, Index min_order, bool even) {

    std::cout << "beginning of construct_harmonic_polynomials" << std::endl;
    _set_arguments(rotation_dofs);
    PolynomialFunction tpoly(m_argument);
    Array<Index> expon(tpoly.num_args(), 0);
    BasisSet even_basis, odd_basis;
    PolyTrie<double> rtrie(expon.size());
    for(Index i = 0; i < expon.size(); i++) {
      expon[i] = 2;
      rtrie(expon) = 1;
      expon[i] = 0;
    }
    PolynomialFunction Rsqrd(m_argument, rtrie);

    //std::cout << " RSQRD 22222222 ===========" << Rsqrd.tex_formula() << std::endl;
    //PolynomialFunction Rsqrd(/*initialize to R^2 = (a^2+b^2+c^2+...) */)
    //initialize even basis with s-orbital
    even_basis.construct_polynomials_by_order(rotation_dofs, 0);

    for(Index i = 0; i < even_basis.size(); i++) {
      //std::cout << "EVEN BASIS" << i << ": " << even_basis[i]->tex_formula() << std::endl;
    }


    //initialize odd basis with p-orbitals
    odd_basis.construct_polynomials_by_order(rotation_dofs, 1);
    //std::cout << " construct_polys by order went OK" << std::endl;
    for(Index i = 0; i < odd_basis.size(); i++) {
      //std::cout << "ODD BASIS" << i << ": " << odd_basis[i]->tex_formula() << std::endl;
    }

    // add even_basis to (*this)
    if(even) append(even_basis);

    // add odd_basis to (*this)
    if(!even) append(odd_basis);

    for(Index i = 2; i <= order; i++) {
      //Get a reference to either the odd or even basis
      if(even && i % 2 != 0) continue;
      if(!even && i % 2 == 0) continue;

      BasisSet &orthog_basis(i % 2 == 0 ? even_basis : odd_basis);

      // Multiply orthog_basis by R^2
      for(Index j = 0; j < orthog_basis.size(); j++) {
        // In lieu of PolynomialFunction::multiply_by() :
        //std::cout << "about to multiply_by() shenanigans" << std::endl;
        Function *tfunc = orthog_basis._at(j);
        //std::cout << "TFUNC NOW :::::::: " << tfunc->tex_formula() << std::endl;
        //Function const *temp;
        //temp = tfunc->multiply(&Rsqrd);
        //std::cout << " TEMP TEST " << temp->tex_formula() << std::endl;
        orthog_basis._at(j) = (tfunc->multiply(&Rsqrd));
        //*(orthog_basis._at(j)) = *temp;
        orthog_basis._at(j)->set_arguments(orthog_basis.m_argument);
        //std::cout << " NEW ORTHOG MULTI " << orthog_basis[j]->tex_formula() << std::endl;
        // Function const *tfunc=orthog_basis[j];
        // orthog_basis._at(j)= tfunc->multiply(&Rsqrd);
        // deleting tfunc doesn't work. is this a problem?
        delete tfunc;
      }
      //std::cout << " after inner for loop" << std::endl;
      BasisSet curr_order_basis;
      curr_order_basis.construct_polynomials_by_order(rotation_dofs, i);
      for(Index i = 0; i < curr_order_basis.size(); i++) {
        //std::cout << "CURR ORDER BASIS" << i << ": " << curr_order_basis[i]->tex_formula() << std::endl;
      }
      for(Index i = 0; i < orthog_basis.size(); i++) {
        //std::cout << "Orthog BASIS" << i << ": " << orthog_basis[i]->tex_formula() << std::endl;
      }

      //std::cout << " after construct_polynomials" << std::endl;
      curr_order_basis.make_orthogonal_to(orthog_basis);
      for(Index i = 0; i < curr_order_basis.size(); i++) {
        //std::cout << "CURR ORDER BASIS" << i << ": " << curr_order_basis[i]->tex_formula() << std::endl;
      }
      //std::cout << " after make_orthogonal_to" << std::endl;
      curr_order_basis.Gram_Schmidt();
      for(Index i = 0; i < curr_order_basis.size(); i++) {
        //std::cout << "CURR ORDER BASIS" << i << ": " << curr_order_basis[i]->tex_formula() << std::endl;
      }
      //std::cout << " after Gram_Schmidt" << std::endl;
      orthog_basis.append(curr_order_basis);
      //std::cout << " after append" << std::endl;
      // add basis for current order to (*this)
      append(curr_order_basis);
      //std::cout << " after one instance of outer for loop " << std::endl;
    }
    /*you could get the symmetry representation here if you took a SymGroup*/
    //get_symmetry_representation(the_sym_group);
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

