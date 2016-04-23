#include "casm/basis_set/PolynomialFunction.hh"

#include "casm/container/Counter.hh"
#include "casm/container/IsoCounter.hh"
#include "casm/container/MultiCounter.hh"

#include "casm/symmetry/SymOp.hh"

#include "casm/basis_set/OccupantFunction.hh"
#include "casm/basis_set/Variable.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/basis_set/FunctionVisitor.hh"
namespace CASM {


  PolynomialFunction::PolynomialFunction(const std::vector<std::shared_ptr<BasisSet> > &_args) : Function(_args), m_coeffs(0) {
    Index i;
    //i counts over subspaces (i.e., BasisSets in init_args)
    for(i = 0; i < _args.size(); i++) {
      //m_subspaces.push_back(Array<Index>(init_args[i].size()));
      //m_sub_sym_reps.push_back(init_args[i].basis_symrep_ID());
      if((_args[i]->basis_symrep_ID()).empty()) {
        std::cerr << "WARNING: Initializing a PolynomialFunction without knowing how to apply symmetry to it. \n"
                  << "         Something bad will probably happen; consider this your warning!\n";
      }
      if(!_args[i]->size()) {
        std::cerr << "WARNING: In PolynomialFunction constructor, initial arguments are ill-defined. Continuing..\n";
        return;
      }
    }

    m_coeffs.redefine(num_args());
  }

  //*******************************************************************************************

  PolynomialFunction::PolynomialFunction(const std::vector<std::shared_ptr<BasisSet> > &_args, const PolyTrie<double> &_coeffs)
    : Function(_args), m_coeffs(_coeffs) {
    Index i;
    //i counts over subspaces (i.e., BasisSets in init_args)
    for(i = 0; i < _args.size(); i++) {
      //m_subspaces.push_back(Array<Index>(init_args[i].size()));
      //m_sub_sym_reps.push_back(init_args[i].basis_symrep_ID());
      if((_args[i]->basis_symrep_ID()).empty()) {
        std::cerr << "WARNING: Initializing a PolynomialFunction without knowing how to apply symmetry to it. \n"
                  << "         Something bad will probably happen; consider this your warning!\n";
      }
      if(!_args[i]->size()) {
        std::cerr << "WARNING: In PolynomialFunction constructor, initial arguments are ill-defined. Continuing..\n";
        return;
      }
    }

    if(m_coeffs.depth() != num_args()) {
      std::cerr << "CRITICAL ERROR: PolynomialFunction initialized with coefficients that are incompatible with its argument list.\n"
                << "                Exiting...\n";
      exit(1);
    }
  }

  //*******************************************************************************************

  PolynomialFunction::PolynomialFunction(const PolynomialFunction &LHS, const PolynomialFunction &RHS) : m_coeffs(0) {

    Index i, j;
    //RHS_ind tracks how exponents of RHS map into resulting product
    Array<Index> RHS_ind;
    for(i = 0; i < LHS.m_argument.size(); i++) {
      m_argument.push_back(LHS.m_argument[i]);
    }
    m_arg2sub = LHS.m_arg2sub;
    m_arg2fun = LHS.m_arg2fun;
    Index linear_ind;

    // check each argument set from RHS to see if it already exists in m_argument
    for(j = 0; j < RHS.m_argument.size(); j++) {
      bool add_subspace(true);
      linear_ind = 0;
      for(i = 0; i < LHS.m_argument.size() && add_subspace; i++) {
        if(LHS.m_argument[i]->compare(*RHS.m_argument[j])) {
          add_subspace = false;
          break;
        }
        linear_ind += LHS.m_argument[i]->size();
      }

      // RHS.m_argument[j] does not already exist in the produt, so we add it
      if(add_subspace) {
        m_argument.push_back(RHS.m_argument[j]);
        for(i = 0; i < RHS.m_argument[j]->size(); i++) {
          m_arg2sub.push_back(m_argument.size() - 1);
          m_arg2fun.push_back(i);
        }
      }

      // linear_ind is the index in (*this) that corresponds to argument (*RHS.m_argument[j])[0]
      RHS_ind.append(Array<Index>::sequence(linear_ind, linear_ind + RHS.m_argument[j]->size() - 1));
    }

    m_coeffs.redefine(num_args());
    Array<Index> LHS_ind(num_args()), tot_ind(num_args());

    PolyTrie<double>::const_iterator LHS_it(LHS.m_coeffs.begin()), LHS_end(LHS.m_coeffs.end());
    double t_coeff;
    for(; LHS_it != LHS_end; ++LHS_it) {
      for(i = 0; i < LHS_it.key().size(); i++) {
        LHS_ind[i] = LHS_it.key()[i];
      }
      PolyTrie<double>::const_iterator RHS_it(RHS.m_coeffs.begin()), RHS_end(RHS.m_coeffs.end());
      for(; RHS_it != RHS_end; ++RHS_it) {
        t_coeff = (*RHS_it) * (*LHS_it);
        if(almost_zero(t_coeff))
          continue;

        tot_ind = LHS_ind;
        for(i = 0; i < RHS_it.key().size(); i++) {
          tot_ind[RHS_ind[i]] += RHS_it.key()[i];
        }
        m_coeffs.at(tot_ind) += t_coeff;
      }
    }
    return;
  }

  //*******************************************************************************************

  PolynomialFunction::PolynomialFunction(const PolynomialFunction &RHS) :
    Function(RHS), m_coeffs(RHS.m_coeffs) {
  }

  //*******************************************************************************************

  PolynomialFunction::PolynomialFunction(const PolynomialFunction &RHS, const PolyTrie<double> &_coeffs) :
    Function(RHS), m_coeffs(_coeffs) {
    if(m_coeffs.depth() != RHS.m_coeffs.depth()) {
      std::cerr << "WARNING: In PolynomialFunction::PolynomialFunction(const PolynomialFunction&, const PolyTrie<double>&),\n"
                << "         the new PolyTrie is incompatible with with the number of arguments. Initializing to zero instead.\n";
      m_coeffs.redefine(RHS.m_coeffs.depth());
    }
  }

  //*******************************************************************************************

  Function *PolynomialFunction::copy() const {
    return new PolynomialFunction(*this);
  }

  //*******************************************************************************************

  Function *PolynomialFunction::copy(const PolyTrie<double> &new_coeffs) const {
    return new PolynomialFunction(*this, new_coeffs);
  }

  //*******************************************************************************************

  bool PolynomialFunction::_accept(const FunctionVisitor &visitor, BasisSet const *home_basis_ptr/*=nullptr*/) {
    return visitor.visit(*this, home_basis_ptr);
  }

  //*******************************************************************************************
  bool PolynomialFunction::depends_on(const Function *test_func) const {
    Index sub_ind(0), arg_ind(0), linear_ind(0);
    for(sub_ind = 0; sub_ind < m_argument.size(); sub_ind++) {
      arg_ind = m_argument[sub_ind]->find(const_cast<Function *const>(test_func));
      linear_ind += arg_ind;
      if(arg_ind < m_argument[sub_ind]->size()) {
        break;
      }
    }
    if(sub_ind == m_argument.size())
      return false;

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(almost_zero(*it))
        continue;

      if(it.key()[linear_ind])
        return true;
    }
    return false;

  }
  //*******************************************************************************************
  bool PolynomialFunction::is_zero() const {
    //Could check to see if arguments are zero.  For now, assume this is done at time of construction
    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(!almost_zero(*it))
        return false;
    }
    return true;
  }

  //*******************************************************************************************
  void PolynomialFunction::small_to_zero(double tol) {
    m_coeffs.prune_zeros(tol);
  }

  //*******************************************************************************************
  Index PolynomialFunction::num_terms()const {
    Index np(0);
    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(!almost_zero(*it))
        np++;
    }
    return np;
  }

  //*******************************************************************************************

  double PolynomialFunction::leading_coefficient()const {
    Index t_ind;
    return leading_coefficient(t_ind);
  }

  //*******************************************************************************************

  double PolynomialFunction::leading_coefficient(Index &index)const {
    index = 0;
    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(!almost_zero(*it))
        return *it;
      index++;
    }
    return 0.0;
  }

  //*******************************************************************************************

  double PolynomialFunction::get_coefficient(Index i)const {
    Index index = 0;
    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(!almost_zero(*it)) {
        if(index == i)
          return *it;

        index++;
      }
    }
    return 0.0;
  }

  //*******************************************************************************************

  void PolynomialFunction::make_formula()const {
    PolyTrie<double> tcoeffs(m_coeffs);
    tcoeffs.sort_leaves(ComparePTLeaf::CustomOrder());
    m_formula.clear();
    m_tex_formula.clear();
    //std::cout << "Making PolynomialFunction Formula:\n";
    //tcoeffs.print_sparse(std::cout);
    //std::cout << std::endl;
    std::stringstream tformula, ttex;
    Index np;
    bool is_zero(true);
    Array<Array<Index> > unique_product;
    Array<double> prefactor;

    PolyTrie<double>::const_iterator it(tcoeffs.begin()), it_end(tcoeffs.end());
    for(; it != it_end; ++it) {
      if(almost_zero(*it))
        continue;

      unique_product.push_back(it.key());
      prefactor.push_back(*it);
      is_zero = false;
    }

    // Comment out following block to turn off monomial sorting
    //Array<Index> iperm;
    //unique_product.sort(iperm);
    //prefactor.permute(iperm);
    // \end monomial sort

    if(is_zero) {
      m_formula = "0";
      m_tex_formula = "0";
      return;
    }


    double func_scale(prefactor[0]);
    if(!almost_zero(func_scale - 1)) {
      if(almost_zero(func_scale + 1)) {
        ttex << '-';
      }
      else {
        ttex << irrational_to_tex_string(func_scale, num_args()*num_args());
      }
    }

    if(unique_product.size() > 1) {
      tformula << '(';
      ttex << '(';
    }

    for(np = 0; np < unique_product.size(); np++) {

      if(np > 0 && prefactor[np] > 0.0) {
        tformula << '+';
      }
      if(almost_zero(prefactor[np] + 1)) {
        tformula << '-';
      }
      else if(!almost_zero(prefactor[np] - 1)) {
        tformula << prefactor[np] << '*';
      }

      if(np > 0 && prefactor[np] / func_scale > 0.0) {
        ttex << '+';
      }
      if(almost_zero(prefactor[np] / func_scale + 1)) {
        ttex << '-';
      }
      else if(!almost_zero(prefactor[np] / func_scale - 1)) {
        ttex << irrational_to_tex_string(prefactor[np] / func_scale, num_args()*num_args());
      }


      //Loop over arguments of unique_product[np]
      int tot_pow(0);
      for(Index linear_ind = 0; linear_ind < unique_product[np].size(); linear_ind++) {
        if(!unique_product[np][linear_ind]) continue;
        if(tot_pow > 0) {
          tformula << '*';
        }
        tot_pow += unique_product[np][linear_ind];

        if(unique_product[np][linear_ind] > 1) {
          tformula << "pow(" << _argument(linear_ind)->formula() << ", " << unique_product[np][linear_ind] << ")";
          ttex << _argument(linear_ind)->tex_formula() << "^{" << unique_product[np][linear_ind] << "} ";
        }
        else {
          tformula << _argument(linear_ind)->formula();
          ttex << _argument(linear_ind)->tex_formula();
        }
      }
      if(tot_pow == 0 && almost_zero(std::abs(prefactor[np] / func_scale) - 1)) {
        tformula << '1';
        ttex << '1';
      }

    }

    if(unique_product.size() > 1) {
      tformula << ')';
      ttex << ')';
    }
    m_tex_formula = ttex.str();
    m_formula = tformula.str();
    //std::cout << "Formula is " << m_formula << '\n';
    //std::cout << "TeX Formula is " << m_tex_formula << '\n';
  }

  //*******************************************************************************************

  void PolynomialFunction::fill_dispatch_table() {
    Function::inner_prod_table[sclass_ID()][sclass_ID()] = new BasicPolyPolyScalarProd();
    Function::operation_table[sclass_ID()][sclass_ID()] = new PolyPolyOperation();

    // IMPORTANT: Do
    //       Function::operation_table[OTHER::sclass_ID()][sclass_ID()] = nullptr;
    // Before
    //       Function::operation_table[sclass_ID()][OTHER::sclass_ID()] = Whatever;

    Function::operation_table[Variable::sclass_ID()][sclass_ID()] = nullptr;
    Function::operation_table[OccupantFunction::sclass_ID()][sclass_ID()] = nullptr;
    Function::operation_table[sclass_ID()][Variable::sclass_ID()] = new PolyVarOperation();
    Function::operation_table[sclass_ID()][OccupantFunction::sclass_ID()] = new PolyOccOperation();

    //BasicPolyVarProduct does not yet exist
    //Function::inner_prod_table[sclass_ID()][Variable::sclass_ID()]=new BasicPolyVarProduct();
  }

  //*******************************************************************************************

  SparseTensor<double> const *PolynomialFunction::get_coeffs() const {
    return nullptr;
  }

  //*******************************************************************************************

  void PolynomialFunction::scale(double scale_factor) {
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    m_coeffs *= scale_factor;
  }

  //*******************************************************************************************
  Function *PolynomialFunction::_apply_sym(const SymOp &op) {
    //std::cout << "Applying symmetry operation to coefficients:\n";
    //m_coeffs.print_sparse(std::cout);
    //m_coeffs.sort_leaves(ComparePTLeaf::CustomOrder());
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    PolyTrie<double> t_trie(m_coeffs.depth());
    m_coeffs.swap(t_trie);

    PolyTrie<double>::const_iterator it(t_trie.begin()), it_end(t_trie.end());
    for(; it != it_end; ++it)
      transform_monomial_and_add(*it, it.key(), op);

    //std::cout << "\n\n And result is:\n";
    //m_coeffs.print_sparse(std::cout);
    return this;
  }

  //*******************************************************************************************
  Function *PolynomialFunction::_apply_sym(const SymOp &op, const std::vector<bool> &transform_flags) {
    //std::cout << "Applying symmetry operation to coefficients:\n";
    //m_coeffs.print_sparse(std::cout);
    //m_coeffs.sort_leaves(ComparePTLeaf::CustomOrder());
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    PolyTrie<double> t_trie(m_coeffs.depth());
    m_coeffs.swap(t_trie);

    PolyTrie<double>::const_iterator it(t_trie.begin()), it_end(t_trie.end());
    for(; it != it_end; ++it)
      transform_monomial_and_add_new(*it, it.key(), op, transform_flags);

    //std::cout << "\n\n And result is:\n";
    //m_coeffs.print_sparse(std::cout);
    return this;
  }

  //*******************************************************************************************
  Function *PolynomialFunction::transform_monomial_and_add(double prefactor, const Array<Index> &ind, const SymOp &op) {
    assert(ind.size() == m_coeffs.depth() && "\'ind\' Array is not compatible with PolynomialFunction in PolynomialFunction::transform_monomial_and_add");

    Array<Eigen::MatrixXd const *> rep_mats(op.get_matrix_reps(_sub_sym_reps()));
    Array<Array<Array<Index> > > exp_counter(m_argument.size(), Array<Array<Index> >());
    Array<Array<Array<Index> > > nz_terms(m_argument.size(), Array<Array<Index> >());
    Array<Array<Array<double> > > nz_coeffs(m_argument.size(), Array<Array<double> >());

    // nz_terms and nz_coeffs keep track of how each individual variable transforms into
    // a linear combination of vectors
    // exp_counter will be used to expand the resulting multinomial expression
    // by counting over the multinomial terms
    // for example, if x^2 transforms to (x-y)^2/2, the end result will be a PolynomialFunction that encodes
    // x^2/2-x*y+y^2/2
    Index linear_offset(0);
    for(Index ns = 0; ns < m_argument.size(); ns++) {
      for(Index na1 = 0; na1 < m_argument[ns]->size(); na1++) {
        nz_terms[ns].push_back(Array<Index>());
        exp_counter[ns].push_back(Array<Index>());
        nz_coeffs[ns].push_back(Array<double>());

        if(!ind[linear_offset + na1])
          continue;

        for(Index na2 = 0; na2 < m_argument[ns]->size(); na2++) {
          if(rep_mats[ns] && !almost_zero((*rep_mats[ns])(na2, na1))) {
            nz_terms[ns][na1].push_back(linear_offset + na2);
            exp_counter[ns][na1].push_back(0);
            nz_coeffs[ns][na1].push_back((*rep_mats[ns])(na2, na1));
          }
          else if(!rep_mats[ns] && na1 == na2) { // assume identity if no symrep exists
            nz_terms[ns][na1].push_back(linear_offset + na2);
            exp_counter[ns][na1].push_back(0);
            nz_coeffs[ns][na1].push_back(1.0);
          }
        }
        exp_counter[ns][na1][0] = ind[linear_offset + na1];
      }
      linear_offset += m_argument[ns]->size();
    }
    /*
    for(int ns=0; ns<exp_counter.size(); ns++){
      for(int na1=0; na1<exp_counter[ns].size(); na1++){
    for(int na2=0; na2<exp_counter[ns][na1].size(); na2++){
    std::cout << "exp_counter[" << ns << "][" << na1 << "][" << na2 << "] = " << exp_counter[ns][na1][na2] << '\n';
    std::cout << "nz_terms[" << ns << "][" << na1 << "][" << na2 << "] = " << nz_terms[ns][na1][na2] << '\n';
    std::cout << "nz_coeffs[" << ns << "][" << na1 << "][" << na2 << "] = " << nz_coeffs[ns][na1][na2] << '\n';
    }
      }
    }
    */
    Array<Index> out_ind(ind.size(), 0);
    double out_coeff;
    Index ns, na1, na2;
    bool cflag(true);
    while(cflag) {
      //std::cout << "Still inside while loop and exp_counter is " << exp_counter << '\n';
      out_coeff = prefactor;
      out_ind.resize(ind.size(), 0);
      for(ns = 0; ns < exp_counter.size(); ns++) {
        for(na1 = 0; na1 < exp_counter[ns].size(); na1++) {
          for(na2 = 0; na2 < exp_counter[ns][na1].size(); na2++) {
            out_ind[nz_terms[ns][na1][na2]] += exp_counter[ns][na1][na2];
            out_coeff *= pow(nz_coeffs[ns][na1][na2], exp_counter[ns][na1][na2]);
          }
          out_coeff *= multinomial_coeff(exp_counter[ns][na1]);
        }
      }
      if(out_ind.sum() != ind.sum()) {
        std::cerr << "WARNING: Starting from " << ind << " a portion of the result is at " << out_ind << '\n';
      }
      /*
        if(ind != out_ind || !almost_equal(prefactor, out_coeff)) {
      std::cout << "Monomial term: " << ind << " : " << prefactor << "\n";
      std::cout << "Transforms to: " << out_ind << " : " << out_coeff << " <--";
      if(ind != out_ind || !almost_zero(prefactor + out_coeff))
      std::cout << "******";
      std::cout << "\n";
        }*/

      //this is where the coefficient gets transforemd
      m_coeffs.at(out_ind) += out_coeff;

      //Increment exponent counters
      for(ns = 0; ns < exp_counter.size(); ns++) {
        for(na1 = 0; na1 < exp_counter[ns].size(); na1++) {
          if(exp_counter[ns][na1].size() <= 1)
            continue;

          for(na2 = 0; na2 < exp_counter[ns][na1].size() - 1; na2++) {
            if(exp_counter[ns][na1][na2] == 0)
              continue;

            exp_counter[ns][na1][na2]--;
            exp_counter[ns][na1][na2 + 1]++;
            if(na2 > 0) {
              exp_counter[ns][na1][0] = exp_counter[ns][na1][na2];
              exp_counter[ns][na1][na2] = 0;
            }
            break;
          }
          if(na2 == exp_counter[ns][na1].size() - 1) {
            exp_counter[ns][na1][0] = exp_counter[ns][na1][na2];
            exp_counter[ns][na1][na2] = 0;
          }
          else {
            break;
          }
        }
        if(na1 < exp_counter[ns].size())
          break;
      }
      //\done incrementing exponent counters

      cflag = ns < exp_counter.size();
    }

    //std::cout << "Transformed PolyTrie: \n";
    //m_coeffs.print_sparse(std::cout);
    //std::cout << '\n';
    return this;
  }

  //*******************************************************************************************
  // Improved by incorporating MultiCounter and IsoCounter-- not yet tested

  Function *PolynomialFunction::transform_monomial_and_add_new(double prefactor,
                                                               const Array<Index> &ind,
                                                               const SymOp &op,
                                                               const std::vector<bool> &transform_flags) {
    assert(ind.size() == m_coeffs.depth() && "\'ind\' Array is not compatible with PolynomialFunction in PolynomialFunction::transform_monomial_and_add");
    typedef IsoCounter<Array<Index> > TermCounter;
    typedef MultiCounter<TermCounter> SubspaceCounter;
    typedef MultiCounter<SubspaceCounter> ExpCounter;

    assert(transform_flags.size() == m_argument.size());
    Array<Eigen::MatrixXd const *> rep_mats(op.get_matrix_reps(_sub_sym_reps()));
    for(Index i = 0; i < transform_flags.size(); i++) {
      if(!transform_flags[i])
        rep_mats[i] = nullptr;
    }
    ExpCounter exp_counter;
    Array<Array<Array<Index> > > nz_terms(m_argument.size(), Array<Array<Index> >());
    Array<Array<Array<double> > > nz_coeffs(m_argument.size(), Array<Array<double> >());

    // nz_terms and nz_coeffs keep track of how each individual variable transforms into
    // a linear combination of vectors
    // exp_counter will be used to expand the resulting multinomial expression
    // by counting over the multinomial terms
    // for example, if x^2 transforms to (x-y)^2/2, the end result will be a genericpolynomialfunction that encodes
    // x^2/2-x*y+y^2/2

    // nz_terms and nz_coeffs is basically a sparse depiction of the symrep matrices.  Instead of constructing them here, they
    // should probably be stored in the symreps themselves and accessed via lazy evaluation.
    Index linear_offset(0);
    for(Index ns = 0; ns < m_argument.size(); ns++) {
      exp_counter.push_back(SubspaceCounter());
      for(Index na1 = 0; na1 < m_argument[ns]->size(); na1++) {
        nz_terms[ns].push_back(Array<Index>());
        nz_coeffs[ns].push_back(Array<double>());

        int pow_to_distribute(ind[linear_offset + na1]);

        if(pow_to_distribute == 0) {
          exp_counter[ns].push_back(TermCounter());
          continue;
        }


        for(Index na2 = 0; na2 < m_argument[ns]->size(); na2++) {
          if(rep_mats[ns] && !almost_zero((*rep_mats[ns])(na2, na1))) {
            //nz_flag = true; // indicate that this variable can transform into something else
            nz_terms[ns][na1].push_back(linear_offset + na2);
            nz_coeffs[ns][na1].push_back((*rep_mats[ns])(na2, na1));
          }
          else if(!rep_mats[ns] && na1 == na2) { // assume identity if no symrep exists
            nz_terms[ns][na1].push_back(linear_offset + na2);
            nz_coeffs[ns][na1].push_back(1.0);
          }
        }
        exp_counter[ns].push_back(TermCounter(Array<Index>(nz_terms[ns][na1].size(), 0),
                                              Array<Index>(nz_terms[ns][na1].size(), pow_to_distribute),
                                              1, pow_to_distribute));
      }
      linear_offset += m_argument[ns]->size();
    }
    /*
    for(Index ns=0; ns<exp_counter.size(); ns++){
      for(Index na1=0; na1<exp_counter[ns].size(); na1++){
    for(Index na2=0; na2<exp_counter[ns][na1].size(); na2++){
    std::cout << "exp_counter[" << ns << "][" << na1 << "][" << na2 << "] = " << exp_counter[ns][na1][na2] << '\n';
    std::cout << "nz_terms[" << ns << "][" << na1 << "][" << na2 << "] = " << nz_terms[ns][na1][na2] << '\n';
    std::cout << "nz_coeffs[" << ns << "][" << na1 << "][" << na2 << "] = " << nz_coeffs[ns][na1][na2] << '\n';
    }
      }
    }
    */
    Array<Index> out_ind(ind.size(), 0);
    double out_coeff;
    Index ns, na1, na2;

    while(exp_counter.valid()) {
      //std::cout << "Still inside while loop and exp_counter is " << exp_counter << '\n';
      out_coeff = prefactor;
      out_ind.resize(ind.size(), 0);
      for(ns = 0; ns < exp_counter.size(); ns++) {
        for(na1 = 0; na1 < exp_counter[ns].size(); na1++) {
          for(na2 = 0; na2 < exp_counter[ns][na1].size(); na2++) {
            out_ind[nz_terms[ns][na1][na2]] += exp_counter[ns][na1][na2];
            out_coeff *= pow(nz_coeffs[ns][na1][na2], exp_counter[ns][na1][na2]);
          }
          out_coeff *= multinomial_coeff(exp_counter[ns][na1]());
        }
      }
      m_coeffs.at(out_ind) += out_coeff;

      exp_counter++;
    }
    return this;
  }

  //*******************************************************************************************

  double PolynomialFunction::remote_eval() const {
    double t_sum(0.0);
    double t_prod;
    Array<double> arg_states(num_args());

    for(Index i = 0; i < num_args(); i++)
      arg_states[i] = _argument(i)->remote_eval();

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      t_prod = *it;

      for(Index i = 0; i < it.key().size(); i++) {
        t_prod *= pow(arg_states[i], it.key()[i]);
      }

      t_sum += t_prod;
    }

    return t_sum;
  }

  //*******************************************************************************************

  double PolynomialFunction::remote_deval(const DoF::RemoteHandle &dvar) const {
    double t_sum(0.0);
    double t_prod;

    //Array<bool> arg_calc(m_argument.size(), false);
    Array<double> arg_states(num_args());
    Array<double> arg_derivs(num_args());

    for(Index i = 0; i < num_args(); i++) {
      arg_states[i] = _argument(i)->remote_eval();
      arg_derivs[i] = _argument(i)->remote_deval(dvar);
    }

    Index const *k_begin, *k_end;

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(almost_zero(*it))
        continue;
      k_begin = it.key().begin();
      k_end = it.key().end();

      //calculate monomial derivative
      for(Index const *i = k_begin; i < k_end; i++) {
        if(*i == 0) continue;
        t_prod = (*it) * arg_derivs[i - k_begin];

        for(Index const *j = k_begin; j < k_end; j++) {
          for(Index p = 0; p < (*j) - int(i == j); p++)
            t_prod *= arg_states[j - k_begin];
        }

        t_sum += t_prod;
      }
    }

    return t_sum;
  }

  //*******************************************************************************************

  double PolynomialFunction::cache_eval() const {
    double t_sum(0.0);
    double t_prod;

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      t_prod = *it;

      for(Index i = 0; i < it.key().size(); i++) {
        t_prod *= pow(_arg_eval_cache(i), it.key()[i]);
      }

      t_sum += t_prod;
    }

    return t_sum;
  }

  //*******************************************************************************************

  double PolynomialFunction::cache_deval(const DoF::RemoteHandle &dvar) const {
    double t_sum(0.0);
    double t_prod;

    //Array<bool> arg_calc(m_argument.size(), false);
    //Array<double> arg_states(num_args());
    //Array<double> arg_derivs(num_args());

    //for(Index i = 0; i < num_args(); i++) {
    //arg_states[i] = _argument(i)->remote_eval();
    //arg_derivs[i] = _argument(i)->remote_deval(dvar);
    //}

    Index const *k_begin, *k_end;

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(almost_zero(*it))
        continue;
      k_begin = it.key().begin();
      k_end = it.key().end();

      //calculate monomial derivative
      for(Index const *i = k_begin; i < k_end; i++) {
        if(*i == 0) continue;
        t_prod = (*it) * _arg_deval_cache(i - k_begin);

        for(Index const *j = k_begin; j < k_end; j++) {
          for(Index p = 0; p < (*j) - int(i == j); p++)
            t_prod *= _arg_eval_cache(j - k_begin);
        }

        t_sum += t_prod;
      }
    }

    return t_sum;
  }

  //*******************************************************************************************

  double PolynomialFunction::eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const {
    if(var_states.size() != num_args()) {
      std::cerr << "WARNING: Evalaution error in PolynomialFunction::eval(). Argument list has\n"
                << "         incompatible number of variables for evaluation.\n";
      return NAN;
    }
    Array<double> eval_table;
    eval_table.reserve(num_args());
    for(Index i = 0; i < num_args(); i++) {
      eval_table.push_back(_argument(i)->eval(dof_IDs, var_states));
    }
    return poly_eval(eval_table);
  }

  //*******************************************************************************************

  double PolynomialFunction::eval(const Array<Index> &dof_IDs, const Array<double> &var_states) const {
    if(var_states.size() != num_args()) {
      std::cerr << "WARNING: Evalaution error in PolynomialFunction::eval(). Argument list has\n"
                << "         incompatible number of variables for evaluation.\n";
      return NAN;
    }
    Array<double> eval_table;
    eval_table.reserve(num_args());
    for(Index i = 0; i < num_args(); i++) {
      eval_table.push_back(_argument(i)->eval(dof_IDs, var_states));
    }

    return poly_eval(eval_table);

  }

  //*******************************************************************************************

  double PolynomialFunction::poly_eval(const Array<double> &arg_states) const {
    if(arg_states.size() != num_args()) {
      std::cerr << "WARNING: Evalaution error in PolynomialFunction::eval(). Argument list has\n"
                << "         incompatible number of variables for evaluation.\n";
      return NAN;
    }

    double t_sum(0.0);
    double t_prod;

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      t_prod = *it;

      for(Index i = 0; i < it.key().size(); i++) {
        t_prod *= pow(arg_states[i], it.key()[i]);
      }

      t_sum += t_prod;
    }

    return t_sum;
  }

  //*******************************************************************************************

  double PolynomialFunction::frobenius_scalar_prod(const PolynomialFunction &RHS) const {
    // This assumes that all the arguments of *this and RHS are the same and mutually orthogonal
    // This is not generally the case, but it is unclear how to resolve this in the case of the
    // Frobenius product (unlike an inner product based on integration over some domain).
    double prod_result(0.0), tprod(0.0);
    Array<Index> texp;
    //std::cout << "Scalar product of depth " << m_coeffs.depth() << " and depth " << RHS.m_coeffs.depth() << "\n";
    //std::cout << "               arg_size " << num_args() << " and arg_size " << RHS.num_args() << "\n";
    PolyTrie<double>::const_iterator it(RHS.m_coeffs.begin()), it_end(RHS.m_coeffs.end());
    std::vector<std::vector<std::set<Index> > > indep;
    indep.reserve(m_argument.size());
    for(Index i = 0; i < m_argument.size(); i++)
      indep.push_back(m_argument[i]->independent_sub_bases());

    for(; it != it_end; ++it) {
      tprod = (*it) * m_coeffs.get(it.key());
      if(almost_zero(tprod, TOL * TOL))
        continue;
      Index l = 0;
      for(Index i = 0; i < indep.size(); i++) {
        for(auto const &indep_set : indep[i]) {
          texp.clear();
          for(Index const &func_ind : indep_set) {
            texp.push_back(it.key()[l + func_ind]);
          }
          tprod /= multinomial_coeff(texp);
        }
        l += m_argument[i]->size();
      }
      prod_result += tprod;
    }
    return prod_result;

  }


  //*******************************************************************************************

  double PolynomialFunction::box_integral_scalar_prod(const PolynomialFunction &RHS, double edge_length)const {

    // This assumes that all the arguments of *this and RHS are the same and mutually orthogonal
    // This is not generally the case, and perhaps we should add a check

    // returns integral of (*this) x RHS using limits -edge_length to +edge_length (a box of specified edge_length centered at origin)
    // divided by volume of box (so that '1' is always normalized)

    double tprod(0), tintegral(0);
    int texp;
    //volume of the box
    double tvol(pow(edge_length, RHS.num_args()));

    if(m_coeffs.depth() != RHS.m_coeffs.depth()) {
      std::cerr << "WARNING!!! Attempting to get scalar_product between incompatible PolynomialFunctions. Assuming that they are orthogonal...\n";
      return 0;
    }

    PolyTrie<double>::const_iterator RHS_it(RHS.m_coeffs.begin()), RHS_end(RHS.m_coeffs.end());
    PolyTrie<double>::const_iterator LHS_it(m_coeffs.begin()), LHS_end(m_coeffs.end());
    if(RHS_it == RHS_end || LHS_it == LHS_end) return 0;

    for(; RHS_it != RHS_end; ++RHS_it) {
      if(almost_zero(*RHS_it))
        continue;

      for(LHS_it = m_coeffs.begin(); LHS_it != LHS_end; ++LHS_it) {
        tintegral = (*RHS_it) * (*LHS_it);
        for(Index i = 0; i < LHS_it.key().size(); i++) {
          texp = LHS_it.key()[i] + RHS_it.key()[i];

          //Check to see if texp is odd, in which case integral is zero
          if(texp != 2 * (texp / 2)) {
            tintegral = 0;
            break;
          }

          tintegral *= 2 * pow(edge_length / 2, texp + 1) / double(texp + 1);
        }
        tprod += tintegral;
      }
    }
    return tprod / tvol;

  }

  //*******************************************************************************************

  Function *PolynomialFunction::minus_equals(const PolynomialFunction *RHS) {
    m_coeffs -= RHS->m_coeffs;
    return this;
  }

  //*******************************************************************************************
  bool PolynomialFunction::compare(const PolynomialFunction *RHS)const {
    // Function::compare() checks that  arguments are equivalent.


    //std::cout << " In compare ... " << std::endl;
    //std::cout << " LHS " <<   std::endl;
    //m_coeffs.print_sparse(std::cout);
    //std::cout << " RHS " <<  std::endl;
    //(RHS->m_coeffs).print_sparse(std::cout);
    return m_coeffs.compare(RHS->m_coeffs, TOL);
  }
  //*******************************************************************************************
  bool PolynomialFunction::prune_zeros() {
    return m_coeffs.prune_zeros();
  }
  //*******************************************************************************************
  Function *PolynomialFunction::plus_equals(const PolynomialFunction *RHS) {
    m_coeffs += RHS->m_coeffs;
    return this;
  }

  //*******************************************************************************************
  Function *PolynomialFunction::apply_sym_coeffs(const SymOp &op, int dependency_layer) {
    std::vector<bool> transform_flags(m_argument.size());
    int tdep;
    for(Index i = 0; i < m_argument.size(); i++) {
      tdep = (m_argument[i]->dependency_layer()) + 1;
      transform_flags[i] = (tdep == dependency_layer);
      if(dependency_layer < tdep)
        m_argument[i]->apply_sym(op, dependency_layer);
    }
    return _apply_sym(op, transform_flags);
  }

  //*******************************************************************************************

  Function *PolynomialFunction::poly_quotient(const Variable *RHS) const {
    Index i = 0;
    for(i = 0; i < num_args(); i++) {
      if(_argument(i)->compare(RHS))
        break;
    }
    if(i == num_args()) {
      return nullptr;
    }

    PolynomialFunction *tpoly(new PolynomialFunction(*this, PolyTrie<double>(m_coeffs.depth())));

    Array<Index> new_key;

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(!it.key()[i])
        continue;

      new_key = it.key();
      new_key[i]--;

      (tpoly->m_coeffs).set(new_key, *it);
    }

    return tpoly;
  }

  //*******************************************************************************************

  Function *PolynomialFunction::poly_remainder(const Variable *RHS) const {
    Index i = 0;
    for(i = 0; i < num_args(); i++) {
      if(_argument(i)->compare(RHS))
        break;
    }
    if(i == num_args()) {
      return copy();
    }


    PolynomialFunction *tpoly(new PolynomialFunction(*this, PolyTrie<double>(m_coeffs.depth())));

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(it.key()[i])
        continue;

      (tpoly->m_coeffs).set(it.key(), *it);
    }

    return tpoly;

  }

  //*******************************************************************************************

  Function *PolynomialFunction::poly_quotient(const OccupantFunction *RHS) const {
    Index i = 0;
    for(i = 0; i < num_args(); i++) {
      if(_argument(i)->compare(RHS))
        break;
    }
    if(i == num_args()) {
      return nullptr;
    }

    PolynomialFunction *tpoly(new PolynomialFunction(*this, PolyTrie<double>(m_coeffs.depth())));

    Array<Index> new_key;
    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(!it.key()[i])
        continue;

      new_key = it.key();
      new_key[i]--;

      (tpoly->m_coeffs).set(new_key, *it);
    }

    return tpoly;
  }

  //*******************************************************************************************

  Function *PolynomialFunction::poly_remainder(const OccupantFunction *RHS) const {
    Index i = 0;
    for(i = 0; i < num_args(); i++) {
      if(_argument(i)->compare(RHS))
        break;
    }
    if(i == num_args()) {
      return copy();
    }

    PolynomialFunction *tpoly(new PolynomialFunction(*this, PolyTrie<double>(m_coeffs.depth())));

    PolyTrie<double>::const_iterator it(m_coeffs.begin()), it_end(m_coeffs.end());
    for(; it != it_end; ++it) {
      if(it.key()[i])
        continue;

      (tpoly->m_coeffs).set(it.key(), *it);
    }

    return tpoly;

  }

  //*******************************************************************************************

  double BasicPolyPolyScalarProd::dot(Function const *LHS, Function const *RHS) const {
    PolynomialFunction const *PLHS(static_cast<PolynomialFunction const *>(LHS)),
                       *PRHS(static_cast<PolynomialFunction const *>(RHS));

    return PLHS->frobenius_scalar_prod(*PRHS);

  }

  //*******************************************************************************************
  bool PolyPolyOperation::compare(Function const *LHS, Function const *RHS) const {
    PolynomialFunction const *PLHS(static_cast<PolynomialFunction const *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    return PLHS->compare(PRHS);
  }

  //*******************************************************************************************
  Function *PolyPolyOperation::multiply(Function const *LHS, Function const *RHS) const {

    PolynomialFunction const *PLHS(static_cast<PolynomialFunction const *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));

    return new PolynomialFunction(*PLHS, *PRHS);

  }

  //*******************************************************************************************

  Function *PolyPolyOperation::multiply_by(Function *LHS, Function const *RHS) const {
    std::cerr << "WARNING: In-place multiplication of one PolynomialFunction with another is not yet implemented!!\n"
              << "         Exiting...\n";
    exit(1);
    return nullptr;
  }

  //*******************************************************************************************
  Function *PolyPolyOperation::subtract(Function const *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS_copy(static_cast<PolynomialFunction *>(LHS->copy()));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS_copy->minus_equals(PRHS);

    PLHS_copy->refresh_ID();

    return PLHS_copy;
  }

  //*******************************************************************************************
  Function *PolyPolyOperation::subtract_from(Function *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS(static_cast<PolynomialFunction *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS->minus_equals(PRHS);

    PLHS->refresh_ID();

    return PLHS;
  }

  //*******************************************************************************************
  Function *PolyPolyOperation::add(Function const *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS_copy(static_cast<PolynomialFunction *>(LHS->copy()));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS_copy->plus_equals(PRHS);

    PLHS_copy->refresh_ID();

    return PLHS_copy;
  }

  //*******************************************************************************************
  Function *PolyPolyOperation::add_to(Function *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS(static_cast<PolynomialFunction *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS->plus_equals(PRHS);

    PLHS->refresh_ID();

    return PLHS;
  }

  //*******************************************************************************************

  Function *PolyVarOperation::poly_quotient(Function const *LHS, Function const *RHS) const {
    return static_cast<PolynomialFunction const *>(LHS)->poly_quotient(static_cast<Variable const *>(RHS));
  }

  //*******************************************************************************************

  Function *PolyVarOperation::poly_remainder(Function const *LHS, Function const *RHS) const {
    return static_cast<PolynomialFunction const *>(LHS)->poly_remainder(static_cast<Variable const *>(RHS));
  }

  //*******************************************************************************************

  Function *PolyOccOperation::poly_quotient(Function const *LHS, Function const *RHS) const {
    return static_cast<PolynomialFunction const *>(LHS)->poly_quotient(static_cast<OccupantFunction const *>(RHS));
  }

  //*******************************************************************************************

  Function *PolyOccOperation::poly_remainder(Function const *LHS, Function const *RHS) const {
    return static_cast<PolynomialFunction const *>(LHS)->poly_remainder(static_cast<OccupantFunction const *>(RHS));
  }

}
