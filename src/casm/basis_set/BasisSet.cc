#include "casm/basis_set/BasisSet.hh"

#include <algorithm>
#include <functional>

#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/OccupantFunction.hh"
#include "casm/basis_set/PolynomialFunction.hh"
#include "casm/basis_set/Variable.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/container/Counter.hh"
#include "casm/container/IsoCounter.hh"
#include "casm/container/MultiCounter.hh"
#include "casm/container/Permutation.hh"
#include "casm/container/multivector.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"
#include "casm/symmetry/SymPermutation.hh"
#include "casm/symmetry/SymRepTools.hh"

namespace CASM {
namespace BasisSet_impl {

void ArgList::add(std::vector<BasisSet> const &BB) {
  for (auto const &B : BB) push_back(&B);
}

void ArgList::add(std::vector<BasisSet const *> const &BB) {
  for (auto B : BB) push_back(B);
}

void ArgList::add(Array<BasisSet> const &BB) {
  for (auto const &B : BB) push_back(&B);
}

void ArgList::add(Array<BasisSet const *> const &BB) {
  for (auto B : BB) push_back(B);
}

}  // namespace BasisSet_impl
//*******************************************************************************************

BasisSet::BasisSet(const BasisSet &init_basis)
    : Array<Function *>(0),
      m_basis_symrep_ID(init_basis.basis_symrep_ID()),
      m_name(init_basis.name()),
      m_basis_ID(init_basis.m_basis_ID),
      m_min_poly_order(init_basis.m_min_poly_order),
      m_max_poly_order(init_basis.m_max_poly_order),
      // m_eval_cache & m_deval_cache are taken care of by BasisSet::push_back
      // m_eval_cache(init_basis.m_eval_cache),
      // m_deval_cache(init_basis.m_deval_cache),
      m_dof_IDs(init_basis.m_dof_IDs),
      m_dof_subbases(init_basis.m_dof_subbases),
      m_min_poly_constraints(init_basis.min_poly_constraints()),
      m_max_poly_constraints(init_basis.max_poly_constraints()) {
  // std::cout << "Constructing B BSet " << this << ", copy of " << &init_basis
  // << std::endl;

  for (Index i = 0; i < init_basis.m_argument.size(); i++) {
    m_argument.push_back(init_basis.m_argument[i]->shared_copy());
  }

  for (Index i = 0; i < init_basis.size(); i++) {
    if (!init_basis[i])
      push_back(nullptr);
    else {
      push_back(init_basis[i]->copy());
    }
  }
}

//*******************************************************************************************

const BasisSet &BasisSet::operator=(const BasisSet &RHS) {
  // std::cout << "Copying to C BSet " << this << ", copy of " << &RHS <<
  // std::endl;
  if (this == &RHS) {
    return *this;
  }
  clear();
  m_basis_symrep_ID = RHS.basis_symrep_ID();
  m_min_poly_order = RHS.min_poly_order();
  m_max_poly_order = RHS.max_poly_order();
  // m_eval_cache & m_deval_cache are taken care of by BasisSet::push_back
  // m_eval_cache = RHS.m_eval_cache;
  // m_deval_cache = RHS.m_deval_cache;
  m_eval_cache.clear();
  m_deval_cache.clear();
  m_name = RHS.name();
  m_dof_IDs = RHS.m_dof_IDs;
  m_dof_subbases = RHS.m_dof_subbases;
  _min_poly_constraints() = RHS.min_poly_constraints();
  _max_poly_constraints() = RHS.max_poly_constraints();

  m_argument.clear();
  for (Index i = 0; i < RHS.m_argument.size(); i++) {
    m_argument.push_back(RHS.m_argument[i]->shared_copy());
  }

  for (Index i = 0; i < RHS.size(); i++) {
    if (!RHS[i])
      push_back(nullptr);
    else {
      push_back(RHS[i]->copy());
    }
  }
  return *this;
}

//*******************************************************************************************

BasisSet::~BasisSet() {
  // std::cout << "Deleting BSet " << this << std::endl;
  clear();
}

//*******************************************************************************************

void BasisSet::clear() {
  for (Index i = 0; i < size(); i++) {
    if (at(i)) delete at(i);
  }
  Array<Function *>::clear();
  _refresh_ID();
}

//*******************************************************************************************

void BasisSet::append(const BasisSet &RHS,
                      std::function<Function *(Function *)> const &transform) {
  // Before appending functions, copy over  DoF IDs and subbasis info
  if (dof_IDs().size()) {
    for (Index i = 0; i < RHS.m_dof_IDs.size(); i++) {
      Index ID_ind = find_index(m_dof_IDs, RHS.m_dof_IDs[i]);
      if (ID_ind == m_dof_IDs.size()) {
        throw std::runtime_error(
            "In BasisSet::append(), it is unsafe to append a BasisSet whose "
            "dependencies differ from (this).");
        // m_dof_IDs.push_back(RHS.m_dof_IDs[i]);
        // m_dof_subbases.push_back(SubBasis());
      }

      for (Index j = 0; j < RHS.m_dof_subbases[i].size(); j++) {
        // std::cout << "*** Push back " << j << " to subbases " << ID_ind << "
        // value " << RHS.m_dof_subbases[i][j] + size() << "\n\n";
        m_dof_subbases[ID_ind].push_back(RHS.m_dof_subbases[i][j] + size());
      }
    }
  } else
    set_dof_IDs(RHS.dof_IDs());

  // Convert polynomial constraints
  _min_poly_constraints().reserve(min_poly_constraints().size() +
                                  RHS.min_poly_constraints().size());
  for (Index i = 0; i < RHS.min_poly_constraints().size(); i++) {
    _min_poly_constraints().push_back(RHS.min_poly_constraints()[i]);
    for (Index j = 0; j < min_poly_constraints().back().first.size(); j++)
      _min_poly_constraints().back().first[j] += size();
  }

  _max_poly_constraints().reserve(max_poly_constraints().size() +
                                  RHS.max_poly_constraints().size());
  for (Index i = 0; i < RHS.max_poly_constraints().size(); i++) {
    _max_poly_constraints().push_back(RHS.max_poly_constraints()[i]);
    for (Index j = 0; j < max_poly_constraints().back().first.size(); j++)
      _max_poly_constraints().back().first[j] += size();
  }

  // Convert RHS.min_poly_order() and RHS.max_poly_order into polynomial
  // constraints
  if (RHS.min_poly_order() > 0) {
    _min_poly_constraints().push_back(
        PolyConstraint(Array<Index>::sequence(size(), size() + RHS.size() - 1),
                       RHS.min_poly_order()));
  }
  if (valid_index(RHS.max_poly_order())) {
    _max_poly_constraints().push_back(
        PolyConstraint(Array<Index>::sequence(size(), size() + RHS.size() - 1),
                       RHS.max_poly_order()));
  }
  for (Index i = 0; i < RHS.size(); i++) {
    if (!RHS[i])
      push_back(nullptr);
    else {
      push_back(transform(RHS[i]->copy()));
    }
  }
  _refresh_ID();
}

//*******************************************************************************************

BasisSet BasisSet::poly_quotient_set(const Function *divisor) const {
  BasisSet new_set;
  for (Index i = 0; i < size(); i++) {
    if (!at(i))
      new_set.push_back(nullptr);
    else
      new_set.push_back(at(i)->poly_quotient(divisor));
  }
  return new_set;
}

//*******************************************************************************************

bool BasisSet::satisfies_exponent_constraints(
    const Array<Index> &expons) const {
  Index exp_sum;
  for (Index i = 0; i < max_poly_constraints().size(); i++) {
    exp_sum = 0;
    for (Index j = 0; j < max_poly_constraints()[i].first.size(); j++) {
      exp_sum += expons[max_poly_constraints()[i].first[j]];
    }
    if (max_poly_constraints()[i].second < exp_sum) return false;
  }
  for (Index i = 0; i < min_poly_constraints().size(); i++) {
    exp_sum = 0;
    for (Index j = 0; j < min_poly_constraints()[i].first.size(); j++) {
      exp_sum += expons[min_poly_constraints()[i].first[j]];
    }
    if (exp_sum < min_poly_constraints()[i].second) return false;
  }
  exp_sum = expons.sum();
  if (valid_index(min_poly_order()) && exp_sum < min_poly_order()) return false;
  if (valid_index(max_poly_order()) && max_poly_order() < exp_sum) return false;
  return true;
}
//*******************************************************************************************
bool BasisSet::accept(const FunctionVisitor &visitor) {
  bool is_changed(false), tchanged;
  // std::cout << "Accepting visitor of type " << visitor.type_name() << "...
  // \n";
  for (Index i = 0; i < m_argument.size(); i++) {
    tchanged = m_argument[i]->accept(visitor);
    is_changed = tchanged || is_changed;
  }
  // std::cout << "BasisSet " << name() << " (" << this << ") and is_changed is
  // " << is_changed << "\n";

  for (Index i = 0; i < size(); i++) {
    if (at(i)) {
      if (is_changed) at(i)->clear_formula();
      tchanged = at(i)->accept(visitor, this);
      is_changed = tchanged || is_changed;
    }
  }
  // if(is_changed)_refresh_ID(); //should we do this here?
  return is_changed;
}

//*******************************************************************************************
/// Remotely evaluate each basis function and add it to the respective value in
/// cumulant
void BasisSet::remote_eval_and_add_to(Array<double> &cumulant) const {
  assert(size() == cumulant.size());
  for (Index i = 0; i < m_argument.size(); i++) m_argument[i]->_eval_to_cache();

  for (Index i = 0; i < size(); i++) {
    if (at(i)) cumulant[i] += at(i)->cache_eval();
  }
}

//*******************************************************************************************
/// Remotely evaluate derivative of each basis function (w.r.t. dvar) and add it
/// to the respective value in cumulant
void BasisSet::remote_deval_and_add_to(Array<double> &cumulant,
                                       const DoF::RemoteHandle &dvar) const {
  assert(size() == cumulant.size());
  for (Index i = 0; i < m_argument.size(); i++)
    m_argument[i]->_deval_to_cache(dvar);

  for (Index i = 0; i < size(); i++) {
    if (at(i)) cumulant[i] += at(i)->cache_deval(dvar);
  }
}

//*******************************************************************************************

bool BasisSet::compare(const BasisSet &RHS) const {
  if (m_argument.size() != RHS.m_argument.size() || size() != RHS.size())
    return false;

  for (Index i = 0; i < m_argument.size(); i++) {
    if (!m_argument[i]->compare(*(RHS.m_argument[i]))) return false;
  }

  for (Index i = 0; i < size(); i++) {
    if (!(at(i)->shallow_compare(RHS[i]))) return false;
  }
  return true;
}

//*******************************************************************************************

int BasisSet::dependency_layer() const {
  int result = -1;
  for (Index i = 0; i < m_argument.size(); i++)
    result = CASM::max(result, m_argument[i]->dependency_layer());
  return ++result;
}
//*******************************************************************************************
/// Remotely evaluate each basis function and add it to the respective value in
/// cumulant
void BasisSet::_eval_to_cache() const {
  remote_eval_to(m_eval_cache.begin(), m_eval_cache.end());
}
//*******************************************************************************************
/// Remotely evaluate each basis function and add it to the respective value in
/// cumulant
void BasisSet::_deval_to_cache(const DoF::RemoteHandle &dvar) const {
  remote_eval_to(m_eval_cache.begin(), m_eval_cache.end());
  remote_deval_to(m_deval_cache.begin(), m_deval_cache.end(), dvar);
  // std::cout << "eval_cache: " << m_eval_cache << "\ndeval_cache: " <<
  // m_deval_cache << "\n\n";
}
//*******************************************************************************************
void BasisSet::set_variable_basis(DoFSet const &_dof_set) {
  set_name(_dof_set.type_name());
  m_argument.clear();
  m_basis_symrep_ID = _dof_set.symrep_ID();
  // std::vector<Index> tdof_IDs;
  for (Index i = 0; i < _dof_set.size(); i++) {
    if (!_dof_set[i].is_locked()) {
      Index t_ind = CASM::find_index(m_dof_IDs, _dof_set[i].ID());
      if (t_ind == m_dof_IDs.size()) {
        m_dof_IDs.push_back(_dof_set[i].ID());
        m_dof_subbases.push_back(Array<Index>());
      }
      m_dof_subbases[t_ind].push_back(i);
    }
  }
  // set_dof_IDs(tdof_IDs);
  for (Index i = 0; i < _dof_set.size(); i++) {
    push_back(new Variable(_dof_set, i));
  }
  _refresh_ID();
}
//*******************************************************************************************
void BasisSet::set_dof_IDs(const std::vector<Index> &new_IDs) {
  _update_dof_IDs(m_dof_IDs, new_IDs);
}
//*******************************************************************************************
// Pass before_IDs by value to avoid aliasing issues when called from
// BasisSet::set_dof_IDs()
bool BasisSet::_update_dof_IDs(std::vector<Index> before_IDs,
                               const std::vector<Index> &after_IDs) {
  bool changed = false;

  // std::cout << "BasisSet " << this << " --> m_dof_IDs is " << m_dof_IDs << ";
  // before_IDs: " << before_IDs << "; after_ID: " << after_IDs << "\n" <<
  // std::endl;
  if (before_IDs.size() != after_IDs.size() && size() > 0) {
    throw std::runtime_error(
        "In BasisSet::update_dof_IDs(), new IDs are incompatible with current "
        "IDs.");
  }
  // std::cout << "***_update_dof_IDs on basisset: " << name()<< "\n";
  // update m_dof_IDs after the other stuff, for easier debugging
  if (m_dof_IDs.size() == 0) {
    m_dof_IDs = after_IDs;
    m_dof_subbases.resize(after_IDs.size());
    changed = true;
  } else {
    Index m;
    for (Index i = 0; i < m_dof_IDs.size(); i++) {
      m = find_index(before_IDs, m_dof_IDs[i]);
      if (m == before_IDs.size())
        throw std::runtime_error(
            "In BasisSet::update_dof_IDs(), new IDs are incompatible with "
            "current IDs.");

      m_dof_IDs[i] = after_IDs[m];
    }
  }

  for (Index i = 0; i < m_argument.size(); i++)
    changed = m_argument[i]->_update_dof_IDs(before_IDs, after_IDs) || changed;

  // std::cout << "Working on basisset: " << name() << " changed is " << changed
  // << "\n"
  //        << "before_IDs: " << before_IDs << "; after_IDs: " << after_IDs <<
  //        "\n";

  for (Index i = 0; i < size(); i++) {
    // std::cout << "Looping " << i+1 << " of " << size() << "; value: " <<
    // at(i) << "; changed: " << changed <<"\n";
    if (changed && at(i)) {
      // std::cout << "clearing formula " << i << "\n";
      at(i)->clear_formula();
    }
  }

  for (Index i = 0; i < size(); i++) {
    if (at(i))
      changed = at(i)->update_dof_IDs(before_IDs, after_IDs) || changed;
  }

  return changed;
}

//*******************************************************************************************

std::vector<std::set<Index>> BasisSet::independent_sub_bases() const {
  std::vector<std::set<Index>> result;
  std::vector<bool> unclaimed(size(), true);
  for (Index i = 0; i < dof_sub_bases().size(); i++) {
    if (dof_sub_basis(i).size() == 0) continue;
    Index j = 0;
    for (j = 0; j < result.size(); j++) {
      if (result[j].find(dof_sub_basis(i)[0]) != result[j].end()) {
        break;
      }
    }
    if (j == result.size()) result.emplace_back();
    for (Index k = 0; k < dof_sub_basis(i).size(); k++) {
      unclaimed[dof_sub_basis(i)[k]] = false;
      result[j].insert(dof_sub_basis(i)[k]);
    }
  }
  bool has_unclaimed = false;
  for (Index i = 0; i < unclaimed.size(); i++) {
    if (unclaimed[i]) {
      if (!has_unclaimed) {
        result.emplace_back();
        has_unclaimed = true;
      }
      result.back().insert(i);
    }
  }
  return result;
}

//*******************************************************************************************

int BasisSet::register_remotes(
    const std::vector<DoF::RemoteHandle> &remote_handles) {
  int sum(0);
  for (Index i = 0; i < size(); i++)
    sum += at(i)->register_remotes(remote_handles);

  for (Index i = 0; i < m_argument.size(); i++)
    sum += m_argument[i]->register_remotes(remote_handles);

  return sum;
}

//*******************************************************************************************

void BasisSet::construct_polynomials_by_order(BasisSet::ArgList const &tsubs,
                                              Index order) {
  _set_arguments(tsubs);
  // Array<BasisSet const *> abset(m_argument);
  PolynomialFunction tpoly(m_argument);
  Array<Index> curr_expon(tpoly.poly_coeffs().depth(), 0);

  IsoCounter<Array<Index>> exp_count(Array<Index>(tsubs[0]->size(), 0),
                                     Array<Index>(tsubs[0]->size(), order), 1,
                                     order);

  do {
    for (int i = 0; i < exp_count.size(); i++) {
      // std::cout << exp_count()[i] << " ";
      curr_expon[i] = exp_count()[i];
    }
    PolyTrie<double> new_trie(curr_expon.size());
    new_trie(curr_expon) = 1.0;
    push_back(tpoly.copy(new_trie));
    // std::cout << std::endl;
  } while (++exp_count);
}

//*******************************************************************************************

void BasisSet::construct_invariant_polynomials(BasisSet::ArgList const &tsubs,
                                               const SymGroup &head_group,
                                               Index order,
                                               Index min_dof_order) {
  // std::cout << "Constructing invariant polynomials of order " << order <<
  // "\n";
  _set_arguments(tsubs);
  // std::cout << "Constructing invariant polynomials from DoFs:\n";
  for (Index i = 0; i < tsubs.size(); ++i) {
    if ((tsubs[i]->basis_symrep_ID()).empty()) {
      tsubs[i]->get_symmetry_representation(head_group);
    }
    // std::cout << i+1 << ") " <<  tsubs[i]->name() << " : ";
    // for(Index j=0; j<tsubs[i]->size(); ++j){
    // std::cout << (*tsubs[i])[j]->tex_formula() << "\n";
    //}

    // std::cout << "DoF_IDs: " << tsubs[i]->dof_IDs() << "\n";
    // std::cout << "subbases: ";
    // for(auto const &subbasis : tsubs[i]->dof_sub_bases())
    // std::cout << " + " << subbasis << "\n";
  }

  // Base case, zero-order function with no arguments -> unit function
  if (tsubs.size() == 0 && order == 0 && min_poly_order() < 1) {
    PolyTrie<double> ttrie(0);
    ttrie(Array<Index>(0)) = 1.0;
    push_back(new PolynomialFunction(m_argument, ttrie));
    return;
  }
  // std::cout << "\n\n";
  // std::cout << "dof_IDs are: " << dof_IDs() << "\n\n";
  // std::cout << "dof_sub_bases are: " << dof_sub_bases() << "\n\n";
  PolynomialFunction *tpoly;
  Array<Index> curr_exp;
  typedef BasisSet::SubBasis SubBasis;
  typedef IsoCounter<Array<Index>> OrderCount;
  typedef MultiCounter<OrderCount> ExpCount;

  // Add constraints that set a minimum combined order for each local dofset
  for (Index i = 0; i < dof_IDs().size() && min_dof_order >= 0; i++) {
    SubBasis big_sub_basis;
    Index offset = 0;
    for (Index j = 0; j < tsubs.size(); offset += tsubs[j]->size(), j++) {
      // std::cout << "tsub " << j << ": size " << tsubs[j]->size() << ",
      // dof_IDs " << tsubs[j]->dof_IDs() << "\n";
      Index ID_ind = find_index(tsubs[j]->dof_IDs(), dof_IDs()[i]);
      if (ID_ind == (tsubs[j]->dof_IDs()).size()) continue;
      const SubBasis &sub_basis(tsubs[j]->dof_sub_basis(ID_ind));
      // std::cout << "sub_basis[" << j << "](" << ID_ind << "): " << sub_basis
      // << "\n";
      for (Index b = 0; b < sub_basis.size(); b++) {
        big_sub_basis.push_back(sub_basis[b] + offset);
      }
    }
    // std::cout << "Adding min constraint: " << big_sub_basis << ":  " <<
    // min_dof_order << "\n";
    add_min_poly_constraint(big_sub_basis, min_dof_order);
  }
  //\ end DoFset constraints

  OrderCount::Container initial_order(tsubs.size(), 0),
      final_order(tsubs.size(), order);
  for (Index i = 0; i < tsubs.size(); i++) {
    if (valid_index(tsubs[i]->min_poly_order()))
      initial_order[i] = tsubs[i]->min_poly_order();
    if (valid_index(tsubs[i]->max_poly_order()))
      final_order[i] = tsubs[i]->max_poly_order();
  }
  // std::cout << "order_count from: " << initial_order << " to final order: "
  // << final_order << "\n\n";
  OrderCount order_count(initial_order, final_order, 1, order);

  ExpCount exp_count;
  for (Index i = 0; i < tsubs.size(); i++) {
    exp_count.push_back(OrderCount(Array<Index>(tsubs[i]->size(), 0),
                                   Array<Index>(tsubs[i]->size(), order), 1,
                                   order_count[i]));
    curr_exp.append(exp_count.back());
  }

  for (; order_count.valid(); ++order_count) {
    // std::cout << "New order_count: ";
    for (Index i = 0; i < exp_count.size(); i++) {
      exp_count[i].set_sum_constraint(order_count[i]);
      // std::cout << order_count[i] << "  ";
    }
    exp_count.reset();
    // std::cout << "\n";

    for (; exp_count.valid(); ++exp_count) {
      bool valid_expon = true;
      for (Index i = 0; i < exp_count.size() && valid_expon; i++) {
        // std::cout << "exp_count[" << i << "] : " <<  exp_count[i]() << " --
        // valid? ";
        valid_expon = tsubs[i]->satisfies_exponent_constraints(exp_count[i]());
        // std::cout << valid_expon << "\n";
      }

      if (!valid_expon) continue;

      Index ne = 0;
      ExpCount::const_value_iterator it(exp_count.value_begin()),
          it_end(exp_count.value_end());
      // std::cout << "exp: ";
      for (; it != it_end; ++it) {
        curr_exp[ne++] = *it;
        // std::cout << *it << " ";
      }
      if (!satisfies_exponent_constraints(curr_exp)) {
        // std::cout << "failure to add " << curr_exp << "\n\n";
        continue;
      }
      // else:
      // std::cout << "Adding exponent " << curr_exp << "\n\n";

      tpoly = new PolynomialFunction(m_argument);
      for (Index i = 0; i < head_group.size(); i++) {
        tpoly->transform_monomial_and_add(1, curr_exp, head_group[i]);
      }
      push_back(tpoly);
    }
  }
  Gram_Schmidt();

  // std::cout << "Result, with " << size() << " invariant functions:\n";
  // for(Index i = 0; i < size(); i++)
  // std::cout << "F_" << i << " = " << at(i)->tex_formula() << "\n\n";
  return;
}

/// \brief Directly specify site basis functions
///
/// \param site_basis_functions The non-constant site basis functions. Of size
///     (allowed_occs.size()-1 x allowed_occs.size()). Use
///     site_basis_functions[i] as phi_i.
///
void BasisSet::construct_discrete_functions(
    const DiscreteDoF &allowed_occs,
    std::vector<std::vector<double>> site_basis_functions, Index basis_ind,
    const SymGroup &symgroup) {
  m_argument.clear();
  set_name(allowed_occs.type_name());
  m_max_poly_order = 1;
  Index N = allowed_occs.size();
  if (N <= 1) {
    m_basis_symrep_ID = SymGroupRepID::identity(0);
    return;
  }

  if (!allowed_occs.is_locked()) {
    set_dof_IDs(std::vector<Index>(1, allowed_occs.ID()));
    m_dof_subbases[0] = Array<Index>::sequence(0, N - 2);
  }

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(N, N);
  B.col(0) = Eigen::VectorXd::Constant(N, 1.0);
  if (site_basis_functions.size() != N - 1) {
    throw std::runtime_error(
        "Error in BasisSet::construct_discrete_functions: number of site basis "
        "functions is incorrect");
  }
  for (Index j = 0; j < site_basis_functions.size(); ++j) {
    if (site_basis_functions[j].size() != N) {
      throw std::runtime_error(
          "Error in BasisSet::construct_discrete_functions: site basis "
          "function size is incorrect");
    }
    for (Index i = 0; i < site_basis_functions[j].size(); ++i) {
      B(i, j + 1) = site_basis_functions[j][i];
    }
  }

  // Columns of B are our basis functions
  for (Index i = 1; i < N; i++) {
    OccupantFunction tOF(allowed_occs, B.col(i), size(), basis_ind,
                         allowed_occs.symrep_ID());

    push_back(tOF.copy());
  }

  // Calculate BasisSet symmetry representation, based on
  // allowed_occs.symrep_ID() && B matrix Q*B.T=B.T*S, where we know S (how to
  // transform a column vector), and we want Q (how to transform row vector) so
  // Q=B.T*S*inv(B.T)
  if (allowed_occs.symrep_ID().is_identity())
    m_basis_symrep_ID = SymGroupRepID::identity(N - 1);
  else
    m_basis_symrep_ID = symgroup.master_group().add_transformed_rep(
        allowed_occs.symrep_ID(),
        Eigen::MatrixXd(B.rightCols(N - 1).transpose()));
}

//*******************************************************************************************
// Construct generalized occupant functions that is orthonormal with respect to
// the inner product defined by gram_mat.
//   - gram_mat should look like a covariance matrix of random vectors, of
//   dimension
//              allowed_occs.size().  This means that gram_mat should be
//              positive-definite and symmetric. Additionally, the variance of
//              the vector of all ones (i.e., v={1,1,...,1}) should have
//              variance = 1, which indicates the constraint of conservation of
//              number of occupants.
//
// The goal is to find a matrix 'B' of column vectors such that
//       B.transpose()*gram_mat*B = Identity
// and such that the first column of 'B' is the vector of all ones (i.e.,
// v={1,1,...,1}) The first property is general solved by the matrix
//        B = gram_mat^(-1/2)*W
// where gram_mat^(-1/2) is the inverse matrix square root of gram_mat, and 'W'
// is an arbitrary orthogonal matrix The second property places constraints on
// 'W'.  We attempt to find a 'B' that is similar to the Chebychev basis in
// certain limiting cases.

void BasisSet::construct_orthonormal_discrete_functions(
    const DiscreteDoF &allowed_occs, const Eigen::MatrixXd &gram_mat,
    Index basis_ind, const SymGroup &symgroup) {
  m_argument.clear();
  set_name(allowed_occs.type_name());
  m_max_poly_order = 1;
  Index N = allowed_occs.size();
  if (N <= 1) {
    m_basis_symrep_ID = SymGroupRepID::identity(0);
    return;
  }

  if (!allowed_occs.is_locked()) {
    set_dof_IDs(std::vector<Index>(1, allowed_occs.ID()));
    m_dof_subbases[0] = Array<Index>::sequence(0, N - 2);
  }

  // std::cout << "INSIDE construct_orthonormal_discrete_functions and gram_mat
  // is \n"; std::cout << gram_mat << "\n\n";
  if (!almost_zero(Eigen::MatrixXd(gram_mat - gram_mat.transpose()))) {
    // we could 'fix' gram_mat, but that could cause mysterious behavior - leave
    // it to the user
    throw std::runtime_error(
        "Passed a Gram Matrix to "
        "BasisSet::construct_orthonormal_discrete_functions that is not "
        "symmetric.");
  }

  Eigen::VectorXd conc_vec(gram_mat * Eigen::MatrixXd::Ones(N, 1));

  if (!almost_equal(1.0, conc_vec.sum())) {
    std::stringstream ss;
    // we could 'fix' gram_mat, but that could cause mysterious behavior - leave
    // it to the user
    ss << "Passed ill-conditioned Gram Matrix to "
          "BasisSet::construct_orthonormal_discrete_functions."
       << "The sum of the elements of the Gram matrix must be equal to 1."
       << "Gram Matrix is:\n"
       << gram_mat;
    throw std::runtime_error(ss.str());
  }

  // ** step 1: find a generic 'B' matrix

  // Use SVD instead of eigendecomposition so that 'U' and 'V' matrices are
  // orthogonal
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> teig(
      gram_mat, Eigen::ComputeEigenvectors);
  if (teig.eigenvalues().minCoeff() < TOL) {
    std::stringstream ss;
    // we could 'fix' gram_mat, but that could cause mysterious behavior - leave
    // it to the user
    ss << "Gram Matrix in BasisSet::construct_orthonormal_discrete_functions "
          "is not positive-definite."
       << "Gram Matrix is:\n"
       << gram_mat
       << "\nSmallest Eigenvalue = " << teig.eigenvalues().minCoeff();
    throw std::runtime_error(ss.str());
  }

  // B matrix is matrix square root of gram_mat.inverse(). Its columns form an
  // orthonormal basis wrt gram_mat In other words,
  //         B = V*(1.0/sqrt(D))*V.inverse()
  // where V is matrix of eigenvectors (as columns) of gram_mat and D is
  // diagonal matrix of eigenvalues of gram_mat B is not ideally oriented for
  // our purposes, so the rest of the algorithm will be focused on fixing the
  // orientation
  Eigen::MatrixXd B = teig.eigenvectors() *
                      (teig.eigenvalues()
                           .array()
                           .cwiseInverse()
                           .cwiseSqrt()
                           .matrix()
                           .asDiagonal()) *
                      teig.eigenvectors().inverse();

  // ** step 2: Make seed basis. This will be used to seed optimized orientation
  // of 'B'
  Eigen::MatrixXd tseed(Eigen::MatrixXd::Zero(N, N));
  Eigen::MatrixXd::Index max_ind(0);
  if (conc_vec.maxCoeff(&max_ind) < 0.75) {
    // no "outlier" probabilities-> Use Chebychev polynomials as seed
    // Fill cosine table -- columns contain powers of x from 0 to N-1
    Eigen::MatrixXd tcos_table(N, N);
    for (Index i = 0; i < N; i++) {
      tcos_table(i, 0) = 1.0;
      double x = -cos(M_PI * (i + 0.5) / double(N));
      for (Index j = 1; j < N; j++) {
        tcos_table(i, j) = x * tcos_table(i, j - 1);
      }
    }

    // QR decomposition of tcos_table yields Q matrix that holds chebychev basis
    tseed = tcos_table.householderQr().householderQ();
  } else {
    // there is an outlier probability --> set seed matrix to occupation basis,
    // with specis 'i==max_ind' as solvent
    Eigen::MatrixXd::Index curr_i(0);
    for (Eigen::MatrixXd::Index i = 0; i < B.rows(); i++) {
      tseed(i, 0) = 1;
      if (i == max_ind) continue;
      for (Eigen::MatrixXd::Index j = 1; j < B.cols(); j++) {
        if (curr_i + 1 == j) tseed(i, j) = 1;
      }
      curr_i++;
    }
  }

  // ** step 3: use seed matric to find a unitary matrix that rotates 'B' a more
  // optimal form Assume: tseed = B * Q, with unitary Q approximate Q by finding
  // QR decomposition of (B.inverse() * tseed)
  Eigen::MatrixXd Q = (B.inverse() * tseed).householderQr().householderQ();

  // Rotate 'B' by multiplication with 'W'
  // eigen matrix multiplication doesn't alias

  B = B * Q;

  // Columns of B are our basis functions, orthonormal wrt gram_mat
  for (Index i = 1; i < N; i++) {
    int sign_change = 1;
    double max_abs(0.0);
    // The sign of each OccupantFunction is ambiguous, so we use a convention
    // Force sign convention
    // max(function(occupation))=max(abs(function(occupation))) If multiple
    // occupations evaluate to the same abs(phi) and it is the maximal abs(phi),
    // then use convention that the last occurence is positive
    // It IS confusing, but here's a simple example:
    //
    //    phi(occ) = {-1, 0, 1}  is always preferred over phi_alt(occ) = {1, 0,
    //    -1}
    //
    // even though they are otherwise both equally valid basis functions
    for (Index j = 0; j < B.rows(); j++) {
      if (std::abs(B(j, i)) > (max_abs - TOL)) {
        max_abs = std::abs(B(j, i));
        sign_change = float_sgn(B(j, i));
      }
    }

    B.col(i) *= double(sign_change);
    OccupantFunction tOF(allowed_occs, B.col(i), size(), basis_ind,
                         allowed_occs.symrep_ID());

    push_back(tOF.copy());
  }

  // ** step 4: Calculate BasisSet symmetry representation, based on
  // allowed_occs.symrep_ID() && B matrix Q*B.T=B.T*S, where we know S (how to
  // transform a column vector), and we want Q (how to transform row vector) so
  // Q=B.T*S*inv(B.T)
  if (allowed_occs.symrep_ID().is_identity())
    m_basis_symrep_ID = SymGroupRepID::identity(N - 1);
  else
    m_basis_symrep_ID = symgroup.master_group().add_transformed_rep(
        allowed_occs.symrep_ID(),
        Eigen::MatrixXd(B.rightCols(N - 1).transpose()));
}

//*******************************************************************************************
// Construct orthonormal basis set of OccupantFunctions for degree of freedom
// 'allowed_occs'
//    - 'occ_probs' is an array of probabilities; occ_prob[i] is the probability
//    of allowed_occs
//                  taking on value allowed_occs[i].
//                  These are equivalent ot composition, but we make a
//                  distinction to avoid confusion with global "parameterized
//                  compositions"
//
// This method finds a Gram matrix that is consistent with the probabilities and
// then passes it on to the overloaded version that takes a Gram matrix. The
// convention we use for the gram matrix is kind of arbitrary. Its limiting
// cases should yield orthonormality of the Chebychev polynomials when the
// probabilities are equal, and orthonormality of the occupation basis when only
// one probability is non-zero.

void BasisSet::construct_orthonormal_discrete_functions(
    const DiscreteDoF &allowed_occs, const std::vector<double> &occ_probs,
    Index basis_ind, const SymGroup &symgroup) {
  Index N = allowed_occs.size();
  if (allowed_occs.size() != occ_probs.size()) {
    std::stringstream ss;
    ss << "Error in BasiSet::construct_orthonormal_discrete_functions(): "
          "occ_probs ([";
    for (Index i = 0; i < occ_probs.size(); ++i) {
      if (i != 0) ss << ",";
      ss << occ_probs[i];
    }
    ss << "]) and allowed_occs (size=" + std::to_string(allowed_occs.size()) +
              ") are incompatible"
              " on sublattice " +
              std::to_string(basis_ind) + ".";
    throw std::runtime_error(ss.str());
  }

  double prob_sum(0.);
  for (double dub : occ_probs) prob_sum += dub;

  if (!almost_equal(1.0, prob_sum)) {
    throw std::runtime_error(
        "In BasiSet::construct_orthonormal_discrete_functions(), occ_probs "
        "must sum to 1 (they specify a probability distributation).\n");
    //<< "                occ_probs is: " << occ_probs << "\nExiting...\n";
  }

  // Build a matrix with N-1 non-zero eigenvalues equal to 1/N.
  // Remaining eigenvalue is zero, and corresponds to vector of all ones
  Eigen::MatrixXd gram_mat(Eigen::MatrixXd::Zero(N, N));

  for (Index i = 0; i < N; i++) {
    gram_mat(i, i) += occ_probs[i] - occ_probs[i] * occ_probs[i];
    for (Index j = 0; j < N; j++) {
      if (i == j) continue;

      gram_mat(i, i) +=
          (occ_probs[i] - occ_probs[j]) * (occ_probs[i] - occ_probs[j]);

      gram_mat(i, j) -=
          occ_probs[i] * occ_probs[j] +
          (occ_probs[i] - occ_probs[j]) * (occ_probs[i] - occ_probs[j]);
    }
  }

  // Add in the component corresponding to vector of all ones
  // this is the uncorrelated part of the covariance
  for (Index i = 0; i < N; i++) {
    for (Index j = 0; j < N; j++) {
      gram_mat(i, j) += occ_probs[i] * occ_probs[j];
    }
  }
  construct_orthonormal_discrete_functions(allowed_occs, gram_mat, basis_ind,
                                           symgroup);
}

//*******************************************************************************************

void BasisSet::construct_harmonic_polynomials(BasisSet::ArgList const &tsubs,
                                              Index order, Index min_order,
                                              bool even_only) {
  _set_arguments(tsubs);
  // std::cout <<"IN HARMONIC POLYNOMIALS:\n";

  // std::cout << "dof_IDs(): " << dof_IDs() << "\n";
  // std::cout << "tsubs->dof_IDs(): " << tsubs[0]->dof_IDs() << "\n";
  // std::cout << "subbases:\n";
  // for(auto const& sub : m_dof_subbases){
  // std::cout << " + "<< sub << "\n";
  //}
  PolynomialFunction tpoly(m_argument);
  Array<Index> expon(tpoly.num_args(), 0);
  BasisSet even_basis, odd_basis;
  PolyTrie<double> rtrie(expon.size());
  for (Index i = 0; i < expon.size(); i++) {
    expon[i] = 2;
    rtrie(expon) = 1;
    expon[i] = 0;
  }
  PolynomialFunction Rsqrd(m_argument, rtrie);

  // PolynomialFunction Rsqrd(/*initialize to R^2 = (a^2+b^2+c^2+...) */)
  // initialize even basis with s-orbital
  even_basis.construct_polynomials_by_order(tsubs, 0);

  // initialize odd basis with p-orbitals
  odd_basis.construct_polynomials_by_order(tsubs, 1);
  // std::cout << " construct_polys by order went OK" << std::endl;

  // add even_basis to (*this)
  if (min_order <= 0) append(even_basis);

  // add odd_basis to (*this)
  if (!even_only && min_order <= 1) append(odd_basis);

  for (Index i = 2; i <= order; i++) {
    // Get a reference to either the odd or even basis
    if (even_only && i % 2 != 0) continue;
    //      if(!even && i%2==0) continue;

    BasisSet &orthog_basis(i % 2 == 0 ? even_basis : odd_basis);

    // Multiply orthog_basis by R^2
    for (Index j = 0; j < orthog_basis.size(); j++) {
      // In lieu of PolynomialFunction::multiply_by() :
      Function *tfunc = orthog_basis._at(j);
      // Function const *temp;
      // temp = tfunc->multiply(&Rsqrd);
      orthog_basis._at(j) = (tfunc->multiply(&Rsqrd));
      //*(orthog_basis._at(j)) = *temp;
      orthog_basis._at(j)->set_arguments(orthog_basis.m_argument);
      // Function const *tfunc=orthog_basis[j];
      // orthog_basis._at(j)= tfunc->multiply(&Rsqrd);
      // deleting tfunc doesn't work. is this a problem?
      delete tfunc;
    }
    BasisSet curr_order_basis;
    curr_order_basis.construct_polynomials_by_order(tsubs, i);
    curr_order_basis.make_orthogonal_to(orthog_basis);
    curr_order_basis.Gram_Schmidt();

    orthog_basis.append(curr_order_basis);
    // add basis for current order to (*this)
    if (min_order <= i) append(curr_order_basis);
  }
  // and fill the sub_bases
  for (Index i = 0; i < m_dof_subbases.size(); ++i) {
    m_dof_subbases[i] = SubBasis::sequence(0, size() - 1);
  }

  // add polynomial constraint to prevent multiplication of harmonic functions
  // with each other.
  if (size() > 0)
    _max_poly_constraints().push_back(
        PolyConstraint(Array<Index>::sequence(0, size() - 1), 1));
}

//*******************************************************************************************
void BasisSet::calc_invariant_functions(const SymGroup &head_group) {
  Function *tfunc, *trans_func;
  for (Index nf = 0; nf < size(); nf++) {
    // std::cout << "trying function " << nf << " of " << size() << ":  ";
    // at(nf)->print(std::cout);
    // std::cout << "\n";
    tfunc = at(nf)->copy();
    at(nf)->scale(0.0);
    for (Index ng = 0; ng < head_group.size(); ng++) {
      trans_func = tfunc->sym_copy_coeffs(head_group[ng]);
      at(nf)->plus_in_place(trans_func);
      delete trans_func;
    }
    // std::cout << "Result of Reynold's operator is ";
    // at(nf)->print(std::cout);
    // std::cout << "\n";
    delete tfunc;
  }
  Gram_Schmidt();
  return;
}

//*******************************************************************************************

BasisSet BasisSet::calc_normal_basis(const SymGroup &head_group,
                                     Eigen::MatrixXd &trans_mat) const {
  if (m_basis_symrep_ID.empty()) {
    get_symmetry_representation(head_group);
  }
  if (m_basis_symrep_ID.empty()) {
    std::cerr
        << "CRITICAL ERROR: Inside BasisSet::calc_normal_basis() and cannot "
           "calculate a valid SymGroup representation. Exiting...\n";
    exit(1);
  }
  if (!head_group.size()) {
    std::cerr << "CRITICAL ERROR: Inside BasisSet::calc_normal_basis() and and "
                 "passed empty SymGroup. Exiting...\n";
  }

  SymGroupRep const &t_rep(
      head_group[0].master_group().representation(m_basis_symrep_ID));
  trans_mat = irrep_trans_mat(t_rep, head_group);
  BasisSet normal_basis(transform_copy(trans_mat));
  normal_basis.m_basis_symrep_ID = t_rep.master_group().add_representation(
      coord_transformed_copy(t_rep, trans_mat));
  return normal_basis;
}

//*******************************************************************************************

void BasisSet::push_back(Function *new_func) {
  Array<Function *>::push_back(new_func);
  m_eval_cache.push_back(0.0);
  m_deval_cache.push_back(0.0);

  // Nothing having to do with DoF_IDs should be touched here.
  // That should happen before pushing back onto the BasisSet
  if (new_func != nullptr) {
    if (m_argument.size() == 0) _set_arguments(new_func->argument_bases());

    _back()->set_arguments(m_argument);
  }
}

//*******************************************************************************************

Function *BasisSet::_linear_combination(const Eigen::VectorXd &coeffs) const {
  if (!size()) return nullptr;
  if (size() != coeffs.size()) {
    std::cerr
        << "FATAL ERROR: In BasisSet::_linear_combination, the number of basis "
           "functions \n"
        << "does not match the size of the coefficient vector. Exiting...\n";
    exit(1);
  }

  Function *combfunc(nullptr), *tfunc(nullptr);

  for (EigenIndex i = 0; i < coeffs.size(); i++) {
    if (almost_zero(coeffs[i])) continue;
    if (!combfunc) {
      combfunc = at(i)->copy();
      combfunc->scale(coeffs[i]);
      continue;
    }
    tfunc = at(i)->copy();
    tfunc->scale(coeffs[i]);
    combfunc->plus_in_place(tfunc);
    delete tfunc;
  }
  if (!combfunc) {
    combfunc = at(0)->copy();
    combfunc->scale(0.0);
  }
  return combfunc;
}

//*******************************************************************************************

/// Essentially, perform a change of basis on BasisSet as defined by trans_mat.
/// Returns a BasisSet whos elements are linear combinations of the original
/// BasisSet. The linear combinations are specified by the ROWS of trans_matx

BasisSet BasisSet::transform_copy(const Eigen::MatrixXd &trans_mat) const {
  BasisSet copy_basis;
  if (trans_mat.cols() != size()) {
    std::cerr << "In BasisSet::transform_copy(), attempting to transform basis "
                 "with a transformation\n"
              << "matrix that has incompatible number of columns (has "
              << trans_mat.cols() << " and needs " << size() << ").\n"
              << "Exiting...\n";
    exit(1);
  }
  for (EigenIndex nc = 0; nc < trans_mat.rows(); nc++) {
    copy_basis.push_back(_linear_combination(trans_mat.row(nc)));
  }

  // We should also transform the SymGroupRep, but we only have the ID, not the
  // MasterSymGroup

  return copy_basis;
}

//*******************************************************************************************
// TODO: Transform functions with mixed dependency layers
BasisSet &BasisSet::apply_sym(const SymOp &op,
                              int _dependency_layer /* = 1 */) {
  int this_dep_layer = dependency_layer();
  if (this_dep_layer == _dependency_layer) {
    if (op.has_valid_master() && !m_basis_symrep_ID.empty() &&
        !m_basis_symrep_ID.is_identity()) {
      auto ptr = op.get_matrix_rep(m_basis_symrep_ID);
      if (ptr == nullptr)
        m_basis_symrep_ID = SymGroupRepID();
      else if (!ptr->isIdentity()) {
        m_basis_symrep_ID =
            op.master_group().add_transformed_rep(m_basis_symrep_ID, *ptr);
      }
    }
    for (Index i = 0; i < size(); i++) at(i)->apply_sym_coeffs(op);
  } else {
    for (Index i = 0; i < m_argument.size(); i++) {
      m_argument[i]->apply_sym(op, _dependency_layer);
    }
  }
  return *this;
}

//*******************************************************************************************

void BasisSet::_set_arguments(BasisSet::ArgList const &new_args) {
  if (m_argument.size()) {
    throw std::runtime_error(
        "In BasisSet::_set_arguments(), cannot reset arguments of "
        "already-initialized BasisSet.");
  }
  m_argument.reserve(new_args.size());
  Index t_ind;
  for (Index i = 0; i < new_args.size(); i++) {
    for (Index ID : new_args[i]->dof_IDs()) {
      t_ind = CASM::find_index(dof_IDs(), ID);
      if (t_ind == dof_IDs().size()) {
        m_dof_IDs.push_back(ID);
        m_dof_subbases.push_back(Array<Index>());
      }
    }
    m_argument.push_back(new_args[i]->shared_copy());
  }
}

//*******************************************************************************************
// Does modified Gram-Schmidt procedure, using tolerance checking for increased
// speed and stability In worst case, this requires N*(N-1)/2 binary operations

// In future, may wish to use alternative approach: 1) find Gram matrix G, where
// G(i,j)=at(i)->dot(at(j));
//                                                  2) find orthonormal
//                                                  eigenvectors of G, forming
//                                                  columns of the matrix V 3)
//                                                  take the linear combination
//                                                  V.transpose()*(*this)
// This alternate approach may be more numerically stable, but probably results
// in less sparse representations, unless there is a way to compute an optimally
// sparse V matrix

bool BasisSet::Gram_Schmidt() {
  bool is_unchanged(true);
  Index i, j;
  double tcoeff;
  Function *tfunc(nullptr);

  // loop over functions
  for (i = 0; i < size(); i++) {
    at(i)->small_to_zero(2 * TOL);

    tcoeff = sqrt(at(i)->dot(at(i)));

    if (tcoeff < TOL) {
      is_unchanged = false;
      delete at(i);
      remove(i);
      i--;
      continue;
    } else if (!almost_zero(tcoeff - 1.0)) {
      is_unchanged = false;
      at(i)->scale(1.0 / tcoeff);
    }

    // loop from i+1 to end and subtract projection of Function i onto Function
    // j
    for (j = i + 1; j < size(); j++) {
      tcoeff = (at(i)->dot(at(j)));
      if (almost_zero(tcoeff)) {
        continue;
      }

      is_unchanged = false;

      tfunc = at(i)->copy();

      if (!almost_zero(tcoeff - 1)) {
        tfunc->scale(tcoeff);
      }

      at(j)->minus_in_place(tfunc);

      delete tfunc;
    }
  }
  if (!is_unchanged) m_basis_symrep_ID = SymGroupRepID();
  return is_unchanged;
}

//*******************************************************************************************

bool BasisSet::Gaussian_Elim() {
  bool is_unchanged(true);
  Index i(0), j(0);
  Index j_min(0), i_min(0), i_temp;
  double tcoeff, min_coeff(1e10);
  Function *tfunc;
  while (i < size()) {
    j_min = Index(-1);
    for (i_temp = i; i_temp < size(); i_temp++) {
      tcoeff = at(i_temp)->leading_coefficient(j);
      if (almost_zero(tcoeff)) {
        delete at(i_temp);
        remove(i_temp);
        i_temp--;
        is_unchanged = false;
        continue;
      }
      if (!valid_index(j_min) || j < j_min ||
          (j == j_min && std::abs(tcoeff) > std::abs(min_coeff))) {
        j_min = j;
        i_min = i_temp;
        min_coeff = tcoeff;
      }
    }
    if (i >= size()) break;

    j = j_min;
    if (i != i_min) {
      is_unchanged = false;
      swap_elem(i, i_min);
    }
    if (!almost_zero(min_coeff - 1)) {
      at(i)->scale(1.0 / min_coeff);
    }

    for (i_temp = 0; i_temp < size(); i_temp++) {
      if (i_temp == i) continue;
      tcoeff = at(i_temp)->get_coefficient(j);
      if (almost_zero(tcoeff)) continue;
      is_unchanged = false;
      if (almost_zero(tcoeff - 1)) {
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
  for (i = 0; i < size(); i++) {
    at(i)->small_to_zero(2 * TOL);
  }

  if (!is_unchanged) m_basis_symrep_ID = SymGroupRepID();

  return is_unchanged;
}
//*******************************************************************************************
void BasisSet::get_symmetry_representation(const SymGroup &head_group) const {
  if (!head_group.size() || !head_group[0].has_valid_master()) return;

  m_basis_symrep_ID = head_group.allocate_representation();
  Function *tfunct(nullptr);
  Eigen::MatrixXd tRep(size(), size());

  for (Index ng = 0; ng < head_group.size(); ng++) {
    // Get representation for operation head_group[ng]
    // store it in matrix tRep
    for (Index nb1 = 0; nb1 < size(); nb1++) {
      tfunct = at(nb1)->sym_copy_coeffs(head_group[ng]);
      for (Index nb2 = 0; nb2 < size(); nb2++) {
        tRep(nb2, nb1) = tfunct->dot(at(nb2));
      }
      delete tfunct;
    }

    // We have the representation of operation head_group[ng]
    // Make a new SymOpRepresentation out of it and push it back
    head_group[ng].set_rep(m_basis_symrep_ID, SymMatrixXd(tRep));
  }
}

//*******************************************************************************************

bool BasisSet::make_orthogonal_to(const BasisSet &ortho_basis) {
  bool ortho_flag(true);
  for (Index i = 0; i < ortho_basis.size(); i++) {
    ortho_flag = make_orthogonal_to(ortho_basis[i]) && ortho_flag;
  }
  return ortho_flag;
}

//*******************************************************************************************

bool BasisSet::make_orthogonal_to(Function const *ortho_func) {
  bool ortho_flag(true);
  if (!size()) {
    return ortho_flag;
  }

  double tcoeff;
  Function *tfunc(ortho_func->copy());

  // std::cout << "Ortho_func: " << tfunc->tex_formula() << "\n\n";
  tcoeff = tfunc->dot(tfunc);
  // std::cout << "Squared magnitude is " << tcoeff << "\n";

  if (almost_zero(tcoeff)) {
    delete tfunc;
    return ortho_flag;
  }

  if (!almost_zero(tcoeff - 1.0)) {
    tfunc->scale(1.0 / sqrt(tcoeff));
  }
  for (Index i = 0; i < size(); i++) {
    // std::cout << "Func " << i << ":  " << at(i)->tex_formula() << "\n\n";
    tcoeff = (at(i)->dot(tfunc));
    // std::cout << "Projection: " << tcoeff << "\n\n";
    if (almost_zero(tcoeff)) {
      continue;
    }
    ortho_flag = false;
    // You're changing the BasisSet, so the representation is no longer useable!
    m_basis_symrep_ID = SymGroupRepID();

    tfunc->scale(tcoeff);
    at(i)->minus_in_place(tfunc);
    tfunc->scale(1.0 / tcoeff);
    at(i)->normalize();

    // std::cout << "Result is " << at(i)->tex_formula() << "\n\n";
  }

  delete tfunc;

  return ortho_flag;
}

BasisSet direct_sum(BasisSet::ArgList const &tsubs) {
  if (tsubs.empty()) return BasisSet();
  multivector<Index>::X<2> compat_maps;
  compat_maps.reserve(tsubs.size());
  BasisSet::ArgList tot_args;
  std::set<Index> dof_ids;
  std::string new_name = tsubs.size() ? tsubs[0]->name() : "";

  for (auto const &sub : tsubs) {
    if ((sub->name()).find(new_name) == std::string::npos)
      new_name += ("+" + sub->name());
    dof_ids.insert(sub->dof_IDs().begin(), sub->dof_IDs().end());
    compat_maps.push_back({});
    for (auto const &arg : sub->arguments()) {
      Index i = 0;
      for (; i < tot_args.size(); ++i) {
        if (tot_args[i]->compare(*arg)) break;
      }
      if (i == tot_args.size()) tot_args.add(*arg);
      compat_maps.back().push_back(i);
    }
  }

  std::vector<Index> dof_ids_vec(dof_ids.begin(), dof_ids.end());
  BasisSet result(new_name, tot_args);
  result.set_dof_IDs(dof_ids_vec);
  for (Index i = 0; i < tsubs.size(); ++i) {
    Index l = 0;
    result.append(*(tsubs[i]), [&](Function *f) -> Function * {
      f->set_arguments(result.arguments(), compat_maps[i]);
      if (f->identifier('l') == "?")
        f->set_identifier('l', std::to_string(l++));
      return f;
    });
  }
  return result;
}
//*******************************************************************************************
//** jsonParser stuff - BasisSet
//*******************************************************************************************

jsonParser &BasisSet::to_json(jsonParser &json) const {
  json.put_obj();

  // class BasisSet: public Array<Function *>
  json["basis_functions"].put_array();
  for (Index i = 0; i < size(); i++) {
    json["basis_functions"].push_back(at(i));
  }

  // mutable int m_basis_symrep_ID;
  json["m_basis_symrep_ID"] = m_basis_symrep_ID;

  // Array<BasisSet> subspaces;
  // json["subspaces"] = m_subspaces;

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

}  // namespace CASM
