#include "casm/basis_set/PolynomialFunction.hh"

#include "casm/container/IsoCounter.hh"
#include "casm/container/MultiCounter.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/basis_set/OccupantFunction.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {


  PolynomialFunction::PolynomialFunction(const Array<BasisSet> &init_args) : m_coeffs(0) {
    Index i, j;
    Array<Index> dim_array;
    //i counts over subspaces (i.e., BasisSets in init_args)
    for(i = 0; i < init_args.size(); i++) {
      m_subspaces.push_back(Array<Index>(init_args[i].size()));
      m_sub_sym_reps.push_back(init_args[i].basis_symrep_ID());
      if(m_sub_sym_reps.back() == Index(-1)) {
        std::cerr << "WARNING: Initializing a PolynomialFunction without knowing how to apply symmetry to it. \n"
                  << "         Something bad will probably happen; consider this your warning!\n";
      }
      if(!init_args[i].size()) {
        std::cerr << "WARNING: In PolynomialFunction constructor, initial arguments are ill-defined. Continuing..\n";
        return;
      }

      for(j = 0; j < init_args[i].size(); j++) {
        m_subspaces[i][j] = m_argument.size();
        m_argument.push_back(init_args[i][j]->copy());
        m_arg2sub.push_back(i);
      }

    }

    m_coeffs.redefine(m_argument.size());
  }

  //********************************************************

  PolynomialFunction::PolynomialFunction(const Array<BasisSet const *> &init_args) : m_coeffs(0) {
    Index i, j;
    Array<Index> dim_array;
    //i counts over subspaces (i.e., BasisSets in init_args)
    for(i = 0; i < init_args.size(); i++) {
      m_subspaces.push_back(Array<Index>(init_args[i]->size()));
      m_sub_sym_reps.push_back(init_args[i]->basis_symrep_ID());
      if(m_sub_sym_reps.back() == Index(-1)) {
        std::cerr << "WARNING: Initializing a PolynomialFunction without knowing how to apply symmetry to it. \n"
                  << "         Something bad will probably happen; consider this your warning!\n";
      }
      if(!init_args[i]->size()) {
        std::cerr << "WARNING: In PolynomialFunction constructor, initial arguments are ill-defined. Continuing..\n";
        return;
      }

      for(j = 0; j < init_args[i]->size(); j++) {
        m_subspaces[i][j] = m_argument.size();
        m_argument.push_back((init_args[i]->at(j))->copy());
        m_arg2sub.push_back(i);
      }

    }

    m_coeffs.redefine(m_argument.size());
  }

  //********************************************************

  PolynomialFunction::PolynomialFunction(const PolynomialFunction &LHS, const PolynomialFunction &RHS) : m_coeffs(0) {

    Index i, j, ii, k;
    Index LID, RID;
    Array<Index> dim_array;
    Array<Index> RHS_dim;

    for(i = 0; i < LHS.m_argument.size(); i++) {
      m_argument.push_back(LHS.m_argument[i]->copy());
    }
    m_sub_sym_reps = LHS.m_sub_sym_reps;
    m_subspaces = LHS.m_subspaces;
    m_arg2sub = LHS.m_arg2sub;

    for(j = 0; j < RHS.m_subspaces.size(); j++) {
      bool add_subspace(true);
      for(i = 0; i < LHS.m_subspaces.size() && add_subspace; i++) {
        if(LHS.m_subspaces[i].size() != RHS.m_subspaces[j].size()) {
          continue;
        }
        for(ii = 0; ii < LHS.m_subspaces[i].size(); ii++) {
          LID = (LHS.m_argument[LHS.m_subspaces[i][ii]]->ID());
          RID = (RHS.m_argument[RHS.m_subspaces[j][ii]]->ID());
          if(LID != RID) {
            break;
          }
        }
        if(ii == LHS.m_subspaces[i].size()) { //subspaces are equivalent
          RHS_dim.append(LHS.m_subspaces[i]);
          add_subspace = false;
        }
      }
      if(add_subspace) {
        m_subspaces.push_back(Array<Index>());
        m_sub_sym_reps.push_back(RHS.m_sub_sym_reps[j]);
        for(k = 0; k < RHS.m_subspaces[j].size(); k++) {
          m_subspaces.back().push_back(m_argument.size());
          m_argument.push_back(RHS.m_argument[RHS.m_subspaces[j][k]]->copy());
          m_arg2sub.push_back(m_subspaces.size() - 1);
        }
        RHS_dim.append(m_subspaces.back());
      }
    }

    m_coeffs.redefine(m_argument.size());
    Array<Index> LHS_ind(m_argument.size()), tot_ind(m_argument.size());

    PTLeaf<double> const *LHS_leaf(LHS.m_coeffs.begin());
    double t_coeff;
    while(LHS_leaf) {
      for(i = 0; i < (LHS_leaf->key()).size(); i++) {
        LHS_ind[i] = (LHS_leaf->key())[i];
      }
      PTLeaf<double> const *RHS_leaf(RHS.m_coeffs.begin());
      while(RHS_leaf) {
        t_coeff = (RHS_leaf->val()) * (LHS_leaf->val());
        if(almost_zero(t_coeff)) {
          RHS_leaf = RHS_leaf->next();
          continue;
        }
        tot_ind = LHS_ind;
        for(i = 0; i < (RHS_leaf->key()).size(); i++) {
          tot_ind[RHS_dim[i]] += (RHS_leaf->key())[i];
        }
        m_coeffs.at(tot_ind) += t_coeff;

        RHS_leaf = RHS_leaf->next();
      }
      LHS_leaf = LHS_leaf->next();
    }
    return;
  }

  //********************************************************

  PolynomialFunction::PolynomialFunction(const PolynomialFunction &RHS) :
    Function(RHS), m_sub_sym_reps(RHS.m_sub_sym_reps), m_subspaces(RHS.m_subspaces),
    m_arg2sub(RHS.m_arg2sub), m_coeffs(RHS.m_coeffs) {
    for(Index i = 0; i < RHS.m_argument.size(); i++) {
      m_argument.push_back(RHS.m_argument[i]->copy());
    }
  }

  //********************************************************

  PolynomialFunction::PolynomialFunction(const PolynomialFunction &RHS, const PolyTrie<double> &_coeffs) :
    Function(RHS), m_sub_sym_reps(RHS.m_sub_sym_reps), m_subspaces(RHS.m_subspaces),
    m_arg2sub(RHS.m_arg2sub), m_coeffs(_coeffs) {
    for(Index i = 0; i < RHS.m_argument.size(); i++) {
      m_argument.push_back(RHS.m_argument[i]->copy());
    }
    if(m_coeffs.depth() != RHS.m_coeffs.depth()) {
      std::cerr << "WARNING: In PolynomialFunction::PolynomialFunction(const PolynomialFunction&, const PolyTrie<double>&),\n"
                << "         the new PolyTrie is incompatible with with the number of arguments. Initializing to zero instead.\n";
      m_coeffs.redefine(RHS.m_coeffs.depth());
    }
  }

  //********************************************************

  Function *PolynomialFunction::copy() const {
    return new PolynomialFunction(*this);
  }

  //********************************************************

  bool PolynomialFunction::accept(const FunctionVisitor &visitor) {
    bool is_updated(Function::accept(visitor));

    return visitor.visit(*this) || is_updated;
  }

  //********************************************************
  bool PolynomialFunction::depends_on(const Function *test_func) const {
    PTLeaf<double> const *current(m_coeffs.begin());
    Index arg_ind = m_argument.find(const_cast<Function *const>(test_func));
    if(arg_ind == m_argument.size())
      return false;

    while(current) {
      if(almost_zero(current->val()))
        continue;

      if((current->key())[arg_ind])
        return true;

      current = current->next();
    }
    return false;

  }
  //********************************************************
  bool PolynomialFunction::is_zero() const {
    //Could check to see if arguments are zero.  For now, assume this is done at time of construction
    PTLeaf<double> const *current(m_coeffs.begin());
    while(current) {
      if(!almost_zero(current->val())) {
        return false;
      }
      current = current->next();
    }
    return true;
  }

  //********************************************************
  void PolynomialFunction::small_to_zero(double tol) {
    m_coeffs.prune_zeros(tol);
  }

  //********************************************************
  Index PolynomialFunction::num_terms()const {
    Index np(0);
    PTLeaf<double> const *current(m_coeffs.begin());
    while(current) {
      if(!almost_zero(current->val())) {
        np++;
      }
      current = current->next();
    }
    return np;
  }

  //********************************************************

  double PolynomialFunction::leading_coefficient()const {
    PTLeaf<double> const *current(m_coeffs.begin());
    while(current) {
      if(!almost_zero(current->val())) {
        return current->val();
      }
      current = current->next();
    }
    return 0.0;
  }

  //********************************************************

  double PolynomialFunction::leading_coefficient(Index &index)const {
    index = 0;
    PTLeaf<double> const *current(m_coeffs.begin());
    while(current) {
      if(!almost_zero(current->val())) {
        return current->val();
      }
      index++;
      current = current->next();
    }
    return 0.0;
  }

  //********************************************************

  double PolynomialFunction::get_coefficient(Index i)const {

    Index index = 0;
    PTLeaf<double> const *current(m_coeffs.begin());
    while(current) {
      if(!almost_zero(current->val())) {
        if(index == i)
          return current->val();
        else
          index++;
      }
      current = current->next();
    }
    return 0.0;
  }

  //********************************************************

  void PolynomialFunction::make_formula()const {
    m_formula.clear();
    m_tex_formula.clear();
    //std::cout << "Making PolynomialFunction Formula:\n";
    //m_coeffs.print_sparse(std::cout);
    //std::cout << std::endl;
    std::stringstream tformula, ttex;
    tformula.precision(10);
    Index np;
    PTLeaf<double> const *current(m_coeffs.begin());
    bool is_zero(true);
    Array<Array<Index> > unique_product;
    Array<double> prefactor;

    while(current) {
      if(almost_zero(current->val())) {
        current = current->next();
        continue;
      }
      unique_product.push_back(current->key());
      prefactor.push_back(current->val());
      is_zero = false;

      current = current->next();

    }

    // Comment out following block to turn off monomial sorting
    Array<Index> iperm;
    unique_product.sort(iperm);
    prefactor.permute(iperm);
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
        ttex << irrational_to_tex_string(func_scale, m_argument.size()*m_argument.size());
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
        ttex << irrational_to_tex_string(prefactor[np] / func_scale, m_argument.size()*m_argument.size());
      }


      //Loop over arguments of unique_product[np]
      int tot_pow(0);
      for(Index na = 0; na < unique_product[np].size(); na++) {
        if(!unique_product[np][na]) continue;
        if(tot_pow > 0) {
          tformula << '*';
        }
        tot_pow += unique_product[np][na];

        if(unique_product[np][na] > 1) {
          tformula << "pow(" << m_argument[na]->formula() << ", " << unique_product[np][na] << ")";
          ttex << m_argument[na]->tex_formula() << "^{" << unique_product[np][na] << "} ";
        }
        else {
          tformula << m_argument[na]->formula();
          ttex << m_argument[na]->tex_formula();
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

  //********************************************************

  void PolynomialFunction::fill_dispatch_table() {
    Function::inner_prod_table[sclass_ID()][sclass_ID()] = new BasicPolyPolyScalarProd();
    Function::operation_table[sclass_ID()][sclass_ID()] = new PolyPolyOperation();

    // IMPORTANT: Do
    //       Function::operation_table[OTHER::sclass_ID()][sclass_ID()] = NULL;
    // Before
    //       Function::operation_table[sclass_ID()][OTHER::sclass_ID()] = Whatever;

    Function::operation_table[OccupantFunction::sclass_ID()][sclass_ID()] = NULL;

    Function::operation_table[sclass_ID()][OccupantFunction::sclass_ID()] = new PolyOccOperation();

    //BasicPolyVarProduct does not yet exist
    //Function::inner_prod_table[sclass_ID()][Variable::sclass_ID()]=new BasicPolyVarProduct();
  }

  //********************************************************

  void PolynomialFunction::scale(double scale_factor) {
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    m_coeffs *= scale_factor;
  }

  //********************************************************
  Function *PolynomialFunction::apply_sym(const SymOp &op) {
    m_formula.clear();
    m_tex_formula.clear();
    refresh_ID();
    PolyTrie<double> t_trie(m_coeffs.depth());
    m_coeffs.swap(t_trie);
    PTLeaf<double> *current(t_trie.begin());
    while(current) {
      transform_monomial_and_add(current->val(), current->key(), op);
      current = current->next();
    }
    return this;
  }

  //********************************************************
  Function *PolynomialFunction::transform_monomial_and_add(double prefactor, const Array<Index> &ind, const SymOp &op) {
    assert(ind.size() == m_coeffs.depth() && "\'ind\' Array is not compatible with PolynomialFunction in PolynomialFunction::transform_monomial_and_add");

    Array<Eigen::MatrixXd const *> rep_mats(op.get_matrix_reps(m_sub_sym_reps));
    Array<Array<Array<Index> > > exp_counter(m_subspaces.size(), Array<Array<Index> >());
    Array<Array<Array<Index> > > nz_terms(m_subspaces.size(), Array<Array<Index> >());
    Array<Array<Array<double> > > nz_coeffs(m_subspaces.size(), Array<Array<double> >());

    // nz_terms and nz_coeffs keep track of how each individual variable transforms into
    // a linear combination of vectors
    // exp_counter will be used to expand the resulting multinomial expression
    // by counting over the multinomial terms
    // for example, if x^2 transforms to (x-y)^2/2, the end result will be a genericpolynomialfunction that encodes
    // x^2/2-x*y+y^2/2
    for(Index ns = 0; ns < m_subspaces.size(); ns++) {
      for(Index na1 = 0; na1 < m_subspaces[ns].size(); na1++) {
        nz_terms[ns].push_back(Array<Index>());
        exp_counter[ns].push_back(Array<Index>());
        nz_coeffs[ns].push_back(Array<double>());

        if(!ind[m_subspaces[ns][na1]])
          continue;

        for(Index na2 = 0; na2 < m_subspaces[ns].size(); na2++) {
          if(rep_mats[ns] && !almost_zero((*rep_mats[ns])(na2, na1))) {
            //nz_flag = true; // indicate that this variable can transform into something else
            nz_terms[ns][na1].push_back(na2);
            exp_counter[ns][na1].push_back(0);
            nz_coeffs[ns][na1].push_back((*rep_mats[ns])(na2, na1));
          }
          else if(!rep_mats[ns] && na1 == na2) { // assume identity if no symrep exists
            nz_terms[ns][na1].push_back(na2);
            exp_counter[ns][na1].push_back(0);
            nz_coeffs[ns][na1].push_back(1.0);
          }
        }
        exp_counter[ns][na1][0] = ind[m_subspaces[ns][na1]];
      }
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
            out_ind[m_subspaces[ns][nz_terms[ns][na1][na2]]] += exp_counter[ns][na1][na2];
            out_coeff *= pow(nz_coeffs[ns][na1][na2], exp_counter[ns][na1][na2]);
          }
          out_coeff *= multinomial_coeff(exp_counter[ns][na1]);
        }
      }
      if(out_ind.sum() != ind.sum()) {
        std::cerr << "WARNING: Starting from " << ind << " a portion of the result is at " << out_ind << '\n';
      }
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
    // std::cout << "Transformed PolyTrie: \n";
    // m_coeffs.print_sparse(std::cout);
    // std::cout << '\n';
    return this;
  }

  //********************************************************
  // Improved by incorporating MultiCounter and IsoCounter-- not yet tested

  Function *PolynomialFunction::transform_monomial_and_add_new(double prefactor, const Array<Index> &ind, const SymOp &op) {
    assert(ind.size() == m_coeffs.depth() && "\'ind\' Array is not compatible with PolynomialFunction in PolynomialFunction::transform_monomial_and_add");
    typedef IsoCounter<Array<Index> > TermCounter;
    typedef MultiCounter<TermCounter> SubspaceCounter;
    typedef MultiCounter<MultiCounter<TermCounter> > ExpCounter;

    Array<Eigen::MatrixXd const *> rep_mats(op.get_matrix_reps(m_sub_sym_reps));

    ExpCounter exp_counter;
    Array<Array<Array<Index> > > nz_terms(m_subspaces.size(), Array<Array<Index> >());
    Array<Array<Array<double> > > nz_coeffs(m_subspaces.size(), Array<Array<double> >());

    // nz_terms and nz_coeffs keep track of how each individual variable transforms into
    // a linear combination of vectors
    // exp_counter will be used to expand the resulting multinomial expression
    // by counting over the multinomial terms
    // for example, if x^2 transforms to (x-y)^2/2, the end result will be a genericpolynomialfunction that encodes
    // x^2/2-x*y+y^2/2

    // nz_terms and nz_coeffs is basically a sparse depiction of the symrep matrices.  Instead of constructing them here, they
    // should probably be stored in the symreps themselves and accessed via lazy evaluation.
    for(Index ns = 0; ns < m_subspaces.size(); ns++) {
      exp_counter.push_back(SubspaceCounter());
      for(Index na1 = 0; na1 < m_subspaces[ns].size(); na1++) {
        nz_terms[ns].push_back(Array<Index>());
        nz_coeffs[ns].push_back(Array<double>());

        int pow_to_distribute(ind[m_subspaces[ns][na1]]);

        if(pow_to_distribute == 0)
          continue;

        for(Index na2 = 0; na2 < m_subspaces[ns].size(); na2++) {
          if(rep_mats[ns] && !almost_zero((*rep_mats[ns])(na2, na1))) {
            //nz_flag = true; // indicate that this variable can transform into something else
            nz_terms[ns][na1].push_back(na2);
            nz_coeffs[ns][na1].push_back((*rep_mats[ns])(na2, na1));
          }
          else if(!rep_mats[ns] && na1 == na2) { // assume identity if no symrep exists
            nz_terms[ns][na1].push_back(na2);
            nz_coeffs[ns][na1].push_back(1.0);
          }
        }
        exp_counter[ns].push_back(TermCounter(Array<Index>(nz_terms[ns][na1].size(), 0),
                                              Array<Index>(nz_terms[ns][na1].size(), pow_to_distribute),
                                              1, pow_to_distribute));
      }
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
            out_ind[m_subspaces[ns][nz_terms[ns][na1][na2]]] += exp_counter[ns][na1][na2];
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

  //********************************************************

  double PolynomialFunction::remote_eval() const {
    double t_sum(0.0);
    double t_prod;
    PTLeaf<double> const *current(m_coeffs.begin());
    Array<double> arg_states(m_argument.size());

    for(Index i = 0; i < m_argument.size(); i++)
      arg_states[i] = m_argument[i]->remote_eval();

    while(current) {
      t_prod = current->val();

      for(Index i = 0; i < (current->key()).size(); i++) {
        t_prod *= pow(arg_states[i], (current->key())[i]);
      }

      t_sum += t_prod;

      current = current->next();
    }

    return t_sum;
  }

  //********************************************************

  double PolynomialFunction::eval(const Array<Index> &dof_IDs, const Array<Index> &var_states) const {
    if(var_states.size() != m_subspaces.size()) {
      std::cerr << "WARNING: Evalaution error in PolynomialFunction::eval(). Argument list has\n"
                << "         incompatible number of variables for evaluation.\n";
      return NAN;
    }
    Array<double> eval_table(m_argument.size(), 0);
    for(Index i = 0; i < m_subspaces.size(); i++) {
      for(Index j = 0; j < m_subspaces[i].size(); j++) {
        eval_table[m_subspaces[i][j]] = m_argument[m_subspaces[i][j]]->eval(dof_IDs, var_states);
      }
    }
    return poly_eval(eval_table);
  }

  //********************************************************

  double PolynomialFunction::eval(const Array<Index> &dof_IDs, const Array<double> &var_states) const {
    if(var_states.size() != m_subspaces.size()) {
      std::cerr << "WARNING: Evalaution error in PolynomialFunction::eval(). Argument list has\n"
                << "         incompatible number of variables for evaluation.\n";
      return NAN;
    }
    Array<double> eval_table(m_argument.size(), 0);
    for(Index i = 0; i < m_subspaces.size(); i++) {
      for(Index j = 0; j < m_subspaces[i].size(); j++) {
        eval_table[m_subspaces[i][j]] = m_argument[m_subspaces[i][j]]->eval(dof_IDs, var_states);
      }
    }
    return poly_eval(eval_table);

  }

  //********************************************************

  double PolynomialFunction::poly_eval(const Array<double> &arg_states) const {
    if(arg_states.size() != m_argument.size()) {
      std::cerr << "WARNING: Evalaution error in PolynomialFunction::eval(). Argument list has\n"
                << "         incompatible number of variables for evaluation.\n";
      return NAN;
    }

    double t_sum(0.0);
    double t_prod;
    PTLeaf<double> const *current(m_coeffs.begin());

    while(current) {
      t_prod = current->val();

      for(Index i = 0; i < (current->key()).size(); i++) {
        t_prod *= pow(arg_states[i], (current->key())[i]);
      }

      t_sum += t_prod;

      current = current->next();
    }

    return t_sum;
  }

  //********************************************************

  double PolynomialFunction::frobenius_scalar_prod(const PolynomialFunction &RHS) const {
    // This assumes that all the arguments of *this and RHS are the same and mutually orthogonal
    // This is not generally the case, but it is unclear how to resolve this in the case of the
    // Frobenius product (unlike an inner product based on integration over some domain).
    double tprod(0);
    PTLeaf<double>const *current(RHS.m_coeffs.begin());
    while(current) {
      tprod += (current->val()) * m_coeffs.get(current->key());
      current = current->next();
    }
    return tprod;

  }


  //********************************************************

  double PolynomialFunction::box_integral_scalar_prod(const PolynomialFunction &RHS, double edge_length)const {

    // This assumes that all the arguments of *this and RHS are the same and mutually orthogonal
    // This is not generally the case, and perhaps we should add a check

    // returns integral of (*this) x RHS using limits -edge_length to +edge_length (a box of specified edge_length centered at origin)
    // divided by volume of box (so that '1' is always normalized)

    double tprod(0), tintegral(0);
    int texp;
    PTLeaf<double>const *current_RHS(RHS.m_coeffs.begin());
    PTLeaf<double>const *current_LHS(m_coeffs.begin());
    if(!current_RHS || !current_LHS) return 0;

    //volume of the box
    double tvol(pow(edge_length, RHS.m_argument.size()));

    while(current_RHS) {
      current_LHS = m_coeffs.begin();
      if((current_LHS->key()).size() != (current_RHS->key()).size()) {
        std::cerr << "WARNING!!! Attempting to get scalar_product between incompatible PolynomialFunctions. Assuming that they are orthogonal...\n";
        return 0;
      }
      if(almost_zero(current_RHS->val())) {
        current_RHS = current_RHS->next();
        continue;
      }
      while(current_LHS) {
        tintegral = (current_LHS->val()) * (current_RHS->val());
        for(Index i = 0; i < (current_LHS->key()).size(); i++) {
          texp = (current_LHS->key())[i] + (current_RHS->key())[i];

          //Check to see if texp is odd, in which case integral is zero
          if(texp != 2 * (texp / 2)) {
            tintegral = 0;
            break;
          }

          tintegral *= 2 * pow(edge_length / 2, texp + 1) / double(texp + 1);
        }
        tprod += tintegral;

        current_LHS = current_LHS->next();
      }
      current_RHS = current_RHS->next();
    }
    return tprod / tvol;

  }

  //********************************************************

  Function *PolynomialFunction::minus_equals(const PolynomialFunction *RHS) {
    m_coeffs -= RHS->m_coeffs;
    return this;
  }

  //********************************************************

  Function *PolynomialFunction::plus_equals(const PolynomialFunction *RHS) {
    m_coeffs += RHS->m_coeffs;
    return this;
  }

  //********************************************************

  Function *PolynomialFunction::poly_quotient(const OccupantFunction *RHS) const {
    Index i = 0;
    for(i = 0; i < m_argument.size(); i++) {
      if(m_argument[i]->compare(RHS))
        break;
    }
    if(i == m_argument.size()) {
      return NULL;
    }

    PolynomialFunction *tpoly(new PolynomialFunction(*this, PolyTrie<double>(m_coeffs.depth())));

    PTLeaf<double>const *current(m_coeffs.begin());
    Array<Index> new_key;
    while(current) {
      if(!(current->key())[i]) {
        current = current->next();
        continue;
      }

      new_key = current->key();
      new_key[i]--;

      (tpoly->m_coeffs).set(new_key, current->val());
      current = current->next();
    }

    return tpoly;
  }

  //********************************************************


  Function *PolynomialFunction::poly_remainder(const OccupantFunction *RHS) const {
    Index i = 0;
    for(i = 0; i < m_argument.size(); i++) {
      if(m_argument[i]->compare(RHS))
        break;
    }
    if(i == m_argument.size()) {
      return copy();
    }


    PolynomialFunction *tpoly(new PolynomialFunction(*this, PolyTrie<double>(m_coeffs.depth())));

    PTLeaf<double>const *current(m_coeffs.begin());
    while(current) {
      if((current->key())[i]) {
        current = current->next();
        continue;
      }

      (tpoly->m_coeffs).set(current->key(), current->val());
      current = current->next();
    }

    return tpoly;

  }

  //********************************************************

  double BasicPolyPolyScalarProd::dot(Function const *LHS, Function const *RHS) const {
    PolynomialFunction const *PLHS(static_cast<PolynomialFunction const *>(LHS)),
                       *PRHS(static_cast<PolynomialFunction const *>(RHS));

    return PLHS->frobenius_scalar_prod(*PRHS);

  }

  //********************************************************

  Function *PolyPolyOperation::multiply(Function const *LHS, Function const *RHS) const {

    PolynomialFunction const *PLHS(static_cast<PolynomialFunction const *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));

    return new PolynomialFunction(*PLHS, *PRHS);

  }

  //********************************************************

  Function *PolyPolyOperation::multiply_by(Function *LHS, Function const *RHS) const {
    std::cerr << "WARNING: In-place multiplication of one PolynomialFunction with another is not yet implemented!!\n"
              << "         Exiting...\n";
    exit(1);
    return NULL;
  }

  //********************************************************
  Function *PolyPolyOperation::subtract(Function const *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS_copy(static_cast<PolynomialFunction *>(LHS->copy()));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS_copy->minus_equals(PRHS);

    PLHS_copy->refresh_ID();

    return PLHS_copy;
  }

  //********************************************************
  Function *PolyPolyOperation::subtract_from(Function *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS(static_cast<PolynomialFunction *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS->minus_equals(PRHS);

    PLHS->refresh_ID();

    return PLHS;
  }

  //********************************************************
  Function *PolyPolyOperation::add(Function const *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS_copy(static_cast<PolynomialFunction *>(LHS->copy()));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS_copy->plus_equals(PRHS);

    PLHS_copy->refresh_ID();

    return PLHS_copy;
  }

  //********************************************************
  Function *PolyPolyOperation::add_to(Function *LHS, Function const *RHS) const {
    PolynomialFunction *PLHS(static_cast<PolynomialFunction *>(LHS));
    PolynomialFunction const *PRHS(static_cast<PolynomialFunction const *>(RHS));
    PLHS->plus_equals(PRHS);

    PLHS->refresh_ID();

    return PLHS;
  }

  //********************************************************

  Function *PolyOccOperation::poly_quotient(Function const *LHS, Function const *RHS) const {
    return static_cast<PolynomialFunction const *>(LHS)->poly_quotient(static_cast<OccupantFunction const *>(RHS));
  }

  //********************************************************

  Function *PolyOccOperation::poly_remainder(Function const *LHS, Function const *RHS) const {
    return static_cast<PolynomialFunction const *>(LHS)->poly_remainder(static_cast<OccupantFunction const *>(RHS));
  }

}

