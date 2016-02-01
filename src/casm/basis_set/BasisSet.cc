#include "casm/basis_set/BasisSet.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/container/Permutation.hh"
#include "casm/container/IsoCounter.hh"
#include "casm/container/Counter.hh"
#include "casm/container/MultiCounter.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/PolynomialFunction.hh"
#include "casm/basis_set/OccupantFunction.hh"
#include "casm/symmetry/SymGroupRep.hh"
#include "casm/symmetry/SymMatrixXd.hh"

namespace CASM {

  //******************************************************************************

  BasisSet::BasisSet(const BasisSet &init_basis) :
    Array<Function * >(0), m_basis_symrep_ID(init_basis.basis_symrep_ID()),
    m_min_poly_order(init_basis.m_min_poly_order),
    m_max_poly_order(init_basis.m_max_poly_order),
    subspaces(init_basis.subspaces) {
    for(Index i = 0; i < init_basis.size(); i++) {
      if(!init_basis[i])
        push_back(nullptr);
      else
        push_back(init_basis[i]->copy());
    }

  }

  //******************************************************************************

  const BasisSet &BasisSet::operator=(const BasisSet &RHS) {
    if(this == &RHS) {
      return *this;
    }
    clear();
    m_basis_symrep_ID = RHS.basis_symrep_ID();
    subspaces = RHS.subspaces;
    m_min_poly_order = RHS.min_poly_order();
    m_max_poly_order = RHS.max_poly_order();

    for(Index i = 0; i < RHS.size(); i++) {
      if(!RHS[i])
        push_back(nullptr);
      else
        push_back(RHS[i]->copy());
    }
    return *this;
  }


  //******************************************************************************

  void BasisSet::append(const BasisSet &RHS) {
    for(Index i = 0; i < RHS.size(); i++) {
      if(!RHS[i])
        push_back(nullptr);
      else
        push_back(RHS[i]->copy());
    }
  }

  //******************************************************************************

  void BasisSet::clear() {
    for(Index i = 0; i < size(); i++) {
      if(at(i))
        delete at(i);
    }
    Array<Function *>::clear();
  }

  //******************************************************************************

  BasisSet::~BasisSet() {
    clear();
  }

  //******************************************************************************
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

  //******************************************************************************
  void BasisSet::accept(const FunctionVisitor &visitor) {
    //std::cout << "Accepting visitor of type " << visitor.type_name() << "... \n";
    for(Index i = 0; i < size(); i++) {
      if(at(i))
        at(i)->accept(visitor);
    }
  }

  //******************************************************************************
  /// Remotely evaluate each basis function and add it to the respective value in cumulant
  void BasisSet::remote_eval_and_add_to(Array<double> &cumulant)const {
    assert(size() == cumulant.size());
    for(Index i = 0; i < size(); i++) {
      if(at(i))
        cumulant[i] += at(i)->remote_eval();
    }
  }

  //******************************************************************************

  void BasisSet::update_dof_IDs(const Array<Index> &before_IDs, const Array<Index> &after_IDs) {
    for(Index i = 0; i < size(); i++) {
      if(at(i))
        at(i)->update_dof_IDs(before_IDs, after_IDs);
    }

    return;
  }

  //******************************************************************************

  void BasisSet::construct_invariant_cluster_polynomials(const Array<Array<BasisSet const *> > &site_args,
                                                         const Array<BasisSet const *> &global_args,
                                                         const SymGroup &head_group,
                                                         const SymGroupRep &permute_group,
                                                         Index max_poly_order) {

    construct_invariant_cluster_polynomials(site_args, global_args, head_group, permute_group, Array<Index>(site_args.size(), 1), max_poly_order);

  }


  //******************************************************************************
  // construct_invariant_cluster_polynomials is specialized for clusters with multiple degrees of freedom per site
  //    - site_args[s][f] is a pointer to the set 'f' of degrees of freedom on site 's' in the cluster
  //    - global_args[f] is a pointer to the set 'f' of degrees of freedom not attached to a site
  //                     which may include cluster DoFs (e.g., CDAs, which aren't attached to specific sites) or crystal DoFs (e.g., strains)
  //    - head_group is the cluster point group
  //    - permute_group is the permutation representation of the cluster group (which does not get stored as Master representation, but perhaps should)
  //    - min_site_order specifies the minimum order allowed for all DoFs of a particular site. Overloaded (see above) so that default value is (1,...,1)
  //                     This default is correct for most cases, but not when cluster DoFs exist. More complicated situations may require an overloaded version
  //                     that accepts some sort of filter functor.
  //    - max_poly_order specifies the overall maximum polynomial order allowed for this basis set
  //
  // We assume that each polynomial must include at least one DoF from each

  void BasisSet::construct_invariant_cluster_polynomials(const Array<Array<BasisSet const *> > &site_args,
                                                         const Array<BasisSet const *> &global_args,
                                                         const SymGroup &head_group,
                                                         const SymGroupRep &permute_group,
                                                         const Array<Index> &min_site_order,
                                                         Index max_poly_order) {
    typedef IsoCounter<CASM::Array<Index> > OrderCount;
    typedef MultiCounter<OrderCount> MultiOrderCount;
    typedef MultiCounter<MultiOrderCount> ExponCount;

    //Inspect the site_args and global_args
    Index N_args(0);
    Array<BasisSet const *> all_bset;
    Array<Index> main_min(site_args.size() + 1, 0), main_max(site_args.size() + 1, 0), block_dims(site_args.size() + 1, 0);
    Array<Array<Index> > sub_mins(site_args.size() + 1), sub_maxs(site_args.size() + 1);
    for(Index i = 0; i < site_args.size(); i++) {
      sub_mins[i].resize(site_args[i].size());
      sub_maxs[i].resize(site_args[i].size());
      //std::cout << "site_args[" << i << "].size() is " << site_args[i].size() << '\n';
      for(Index j = 0; j < site_args[i].size(); j++) {
        block_dims[i] += site_args[i][j]->size();
        N_args += site_args[i][j]->size();

        main_min[i] += site_args[i][j]->min_poly_order();
        main_max[i] += site_args[i][j]->max_poly_order();

        sub_mins[i][j] = site_args[i][j]->min_poly_order();
        sub_maxs[i][j] = site_args[i][j]->max_poly_order();

        all_bset.push_back(site_args[i][j]);
      }
      main_min[i] = CASM::max(main_min[i], min_site_order[i]);
    }

    for(Index j = 0; j < global_args.size(); j++) {
      block_dims.back() += global_args[j]->size();
      N_args += global_args[j]->size();

      main_min.back() += global_args[j]->min_poly_order();
      main_max.back() += global_args[j]->max_poly_order();

      sub_mins.back()[j] = global_args[j]->min_poly_order();
      sub_maxs.back()[j] = global_args[j]->max_poly_order();

      all_bset.push_back(global_args[j]);
    }
    //std::cout << "BLOCK_DIMS: " << block_dims << '\n';
    //std::cout << "MAIN_MIN: " << main_min << '\n';
    //std::cout << "MAIN_MAX: " << main_max << '\n';
    //std::cout << "SUB_MINS: " << sub_mins << '\n';
    //std::cout << "SUB_MAXS: " << sub_maxs << '\n';

    // Make some IsoCounters
    // main_partition ensures that each
    OrderCount main_partition(main_min, main_max, 1, main_min.sum());

    MultiOrderCount sub_partitions;
    ExponCount exp_counter;
    for(Index i = 0; i < sub_mins.size(); i++) {
      sub_partitions.push_back(OrderCount(sub_mins[i], sub_maxs[i], 1, main_partition[i]));
    }

    for(Index i = 0; i < site_args.size(); i++) {
      exp_counter.push_back(MultiOrderCount());
      for(Index j = 0; j < site_args[i].size(); j++) {
        exp_counter[i].push_back(OrderCount(Array<Index>(site_args[i][j]->size(), 0),
                                            Array<Index>(site_args[i][j]->size(), sub_maxs[i][j]),
                                            1, sub_partitions[i][j]));
      }
    }

    exp_counter.push_back(MultiOrderCount());
    for(Index j = 0; j < global_args.size(); j++) {
      exp_counter.back().push_back(OrderCount(Array<Index>(global_args[j]->size(), 0),
                                              Array<Index>(global_args[j]->size(), sub_maxs.back()[j]),
                                              1, sub_partitions.back()[j]));
    }

    // Make exponent permutations out of perm_group
    Array<Permutation> exp_perm;
    Permutation tperm(0);
    exp_perm.reserve(permute_group.size());
    for(Index i = 0; i < permute_group.size(); i++) {
      tperm = *(permute_group.get_permutation(i));
      tperm.append_fixed_points(1);
      exp_perm.push_back(tperm.make_block_permutation(block_dims));
    }

    Index poly_order = main_min.sum();
    max_poly_order = CASM::min(main_max.sum(), max_poly_order);
    PolynomialFunction *tpoly;
    Array<Index> curr_exp(N_args, 0);
    Index ne;
    //std::cout << "poly_order is " << poly_order << " and max_poly_order is " << max_poly_order << '\n';
    for(; poly_order <= max_poly_order; poly_order++) {
      main_partition.set_sum_constraint(poly_order);
      ////std::cout << "INIT main_partition: " << main_partition() << '\n';
      while(main_partition.valid()) {
        //std::cout << "main_partition: " << main_partition() << '\n';
        for(Index i = 0; i < sub_partitions.size(); i++) {
          sub_partitions[i].set_sum_constraint(main_partition[i]);
        }
        //std::cout << "INIT sub_partitions: " << sub_partitions() << '\n';
        while(sub_partitions.valid()) {
          //std::cout << "sub_partitions: " << sub_partitions() << '\n';
          for(Index i = 0; i < exp_counter.size(); i++) {
            for(Index j = 0; j < exp_counter[i].size(); j++) {
              exp_counter[i][j].set_sum_constraint(sub_partitions[i][j]);
            }
          }
          //std::cout << "INIT exp_counter: " << exp_counter() << '\n';
          while(exp_counter.valid()) {
            //std::cout << "exp_counter: " << exp_counter() << '\n';
            ne = 0;
            for(Index i = 0; i < exp_counter.size(); i++) {
              for(Index j = 0; j < exp_counter[i].size(); j++) {
                for(Index k = 0; k < exp_counter[i][j].size(); k++) {
                  curr_exp[ne++] = exp_counter[i][j][k];
                }
              }
            }
            //std::cout << "curr_exp is" << curr_exp << "\n";
            tpoly = new PolynomialFunction(all_bset);
            for(Index i = 0; i < head_group.size(); i++) {
              tpoly->transform_monomial_and_add(1, exp_perm[i].permute(curr_exp), head_group[i]);
            }
            push_back(tpoly);


            exp_counter++;
          }

          sub_partitions++;
        }

        main_partition++;
      }
    }
    // If Gram_Schmidt() is too slow for large basis sets, we can change it to work on each order separately
    Gram_Schmidt();
    return;
  }



  //John G 011013

  //********************************************************
  /**	Fills up your basis with occupation basis polynomials.
   *    For a site of N components (0,1,2,3...N-1) there are
   *    N-1 basis functions. Polynomial i evaluates to 1 if
   *    if component i+1 is occupying the site. If that's not
   *    the case the polynomial evaluates to 0.
   */
  //********************************************************

  void BasisSet::construct_discrete_occupations(const DiscreteDoF &allowed_occs, Index basis_ind, Index sym_rep_ind) {
    m_max_poly_order = 1;
    //Start by making the strings for the individual formula pieces (e.g. "1", "p_Ni", "p_Si"...)
    //Array<std::string> tformula_bits;
    Index N = allowed_occs.size();


    //We create a table such that multiplying it with an array of coefficients gives
    //the right evaluation for the basis functions. Identity works.
    //Eigen::MatrixXd teval_table = Eigen::MatrixXd::Identity(N, N);

    //Start making occupation polynomials excluding the 0th one
    //We just need to make the coefficient vectors now that give the right values
    //when multiplied by the eval_table (e.g. [0100],[0010],[0001] for quaternary)

    Eigen::VectorXd tcoeffs(N);
    for(Index i = 1; i < N; i++) {
      tcoeffs.setZero();

      tcoeffs[i] = 1;
      OccupantFunction tOF(allowed_occs, tcoeffs, size(), basis_ind, sym_rep_ind);

      //tOF.formula_bits = tformula_bits;

      push_back(tOF.copy());
    }
    if(sym_rep_ind == Index(-2))
      m_basis_symrep_ID = -2;

    return;
  }

  //*******************************************************************************
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

  void BasisSet::construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs, const Eigen::MatrixXd &gram_mat, Index basis_ind, Index sym_rep_ind) {
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
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> t_es(gram_mat);
    if(t_es.eigenvalues().minCoeff() < TOL) {
      // we could 'fix' gram_mat, but that could cause mysterious behavior - leave it to the user
      std::cerr << "CRITICAL ERROR: Passed a Gram Matrix to BasisSet::construct_orthonormal_discrete_functions that is not positive-definite.\n"
                << "                Gram Matrix is:\n" << gram_mat << "\nSmallest Eigenvalue = " << t_es.eigenvalues().minCoeff() << "; Exiting...\n";
      exit(1);
    }

    // B matrix is matrix square root of gram_mat.inverse(). Its columns form an orthonormal basis wrt gram_mat
    // In other words,
    //         B = V*(1.0/sqrt(D))*V.inverse()
    // where V is matrix of eigenvectors (as columns) of gram_mat and D is diagonal matrix of eigenvalues of gram_mat
    // B is not ideally oriented for our purposes, so the rest of the algorithm will be focused on fixing the orientation
    Eigen::MatrixXd B = t_es.eigenvectors() * (t_es.eigenvalues().array().cwiseInverse().cwiseSqrt().matrix().asDiagonal()) * t_es.eigenvectors().inverse();

    /*
    std::cout << "Eigenvalues: " << t_es.eigenvalues().transpose() << std::endl
              << "Operated: " << t_es.eigenvalues().array().cwiseInverse().cwiseSqrt().matrix().transpose() << std::endl
              << "Eigenvectors: \n" << t_es.eigenvectors() << std::endl
              << "B: \n" << B << std::endl;
    */


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

    // ** step 3: use seed matrix to find a unitary matrix that rotates 'B' a more 'intuitive' orientation
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
          sign_change = sgn(B(j, i));
        }
      }
      OccupantFunction tOF(allowed_occs, double(sign_change)*B.col(i), size(), basis_ind, sym_rep_ind);

      push_back(tOF.copy());
    }
    if(sym_rep_ind == Index(-2))
      m_basis_symrep_ID = -2;
  }

  //*******************************************************************************
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

  void BasisSet::construct_orthonormal_discrete_functions(const DiscreteDoF &allowed_occs, const Array<double> &occ_probs, Index basis_ind, Index sym_rep_ind) {
    Index N = allowed_occs.size();
    if(allowed_occs.size() != occ_probs.size()) {
      std::cerr << "CRITICAL ERROR: In BasiSet::construct_orthonormal_discrete_functions(), occ_probs and allowed_occs are incompatible!\nExiting...\n";
      exit(1);
    }

    if(!almost_equal(1.0, occ_probs.sum())) {
      std::cerr << "CRITICAL ERROR: In BasiSet::construct_orthonormal_discrete_functions(), occ_probs must sum to 1 (they definite a probability distributation).\n"
                << "                occ_probs is: " << occ_probs << "\nExiting...\n";
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
    construct_orthonormal_discrete_functions(allowed_occs, gram_mat, basis_ind, sym_rep_ind);
  }


  //******************************************************************************

  void BasisSet::calc_invariant_functions(const SymGroup &head_group) {
    Function *tfunc, *trans_func;
    for(Index nf = 0; nf < size(); nf++) {
      // std::cout << "trying function " << nf << " of " << size() << ":  ";
      // at(nf)->print(std::cout);
      // std::cout << "\n";
      tfunc = at(nf)->copy();
      at(nf)->scale(0.0);
      for(Index ng = 0; ng < head_group.size(); ng++) {
        trans_func = tfunc->sym_copy(head_group[ng]);
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

  //******************************************************************************

  Function *BasisSet::linear_combination(const Eigen::VectorXd &coeffs) const {
    if(!size()) return nullptr;
    if(size() != coeffs.size()) {
      std::cerr << "FATAL ERROR: In BasisSet::linear_combination, the number of basis functions \n"
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

  //******************************************************************************

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
      copy_basis.push_back(linear_combination(trans_mat.row(nc)));
    }

    // We should also transform the SymGroupRep, but we only have the ID, not the MasterSymGroup

    return copy_basis;
  }

  //******************************************************************************

  BasisSet &BasisSet::apply_sym(const SymOp &op) {
    for(Index i = 0; i < size(); i++) {
      at(i)->apply_sym(op);
    }

    return *this;
  }

  //******************************************************************************
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
    if(!is_unchanged) m_basis_symrep_ID = -1;
    return is_unchanged;

  }

  //******************************************************************************

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

    if(!is_unchanged) m_basis_symrep_ID = -1;

    return is_unchanged;

  }
  //******************************************************************************
  void BasisSet::get_symmetry_representation(const SymGroup &head_group) const {
    if(!head_group.size() || !head_group[0].has_valid_master()) return;

    m_basis_symrep_ID = head_group.make_empty_representation();
    Function *tfunct(nullptr);
    Eigen::MatrixXd tRep(size(), size());

    for(Index ng = 0; ng < head_group.size(); ng++) {
      //Get representation for operation head_group[ng]
      //store it in matrix tRep
      for(Index nb1 = 0; nb1 < size(); nb1++) {

        tfunct = at(nb1)->sym_copy(head_group[ng]);
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

  //******************************************************************************

  bool BasisSet::make_orthogonal_to(const BasisSet &ortho_basis) {
    bool ortho_flag(true);
    for(Index i = 0; i < ortho_basis.size(); i++) {
      ortho_flag = make_orthogonal_to(ortho_basis[i]) && ortho_flag;
    }
    return ortho_flag;
  }

  //******************************************************************************

  bool BasisSet::make_orthogonal_to(Function const *ortho_func) {

    Index i;
    double tcoeff;
    Function *tfunc(ortho_func->copy());

    bool ortho_flag(true);
    if(!size()) {
      return ortho_flag;
    }

    tcoeff = tfunc->dot(tfunc);

    if(almost_zero(tcoeff)) {
      delete tfunc;
      return ortho_flag;
    }

    if(!almost_zero(tcoeff - 1.0)) {
      tfunc->scale(1.0 / sqrt(tcoeff));
    }
    for(i = 0; i < size(); i++) {
      tcoeff = (at(i)->dot(tfunc));

      if(almost_zero(tcoeff)) {
        continue;
      }
      ortho_flag = false;
      //You're changing the BasisSet, so the representation is no longer useable!
      m_basis_symrep_ID = -1;

      tfunc->scale(tcoeff);
      at(i)->minus_in_place(tfunc);
      tfunc->scale(1.0 / tcoeff);
      at(i)->normalize();
    }

    delete tfunc;


    return ortho_flag;
  }


  //********************************************************
  //** jsonParser stuff - BasisSet
  //********************************************************

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
    json["subspaces"] = subspaces;

    return json;
  }

  //********************************************************

  /*
  void BasisSet::from_json(const jsonParser &json) {

    // no reading BasisSet for now

  }
  */

  //********************************************************

  jsonParser &to_json(const BasisSet &bset, jsonParser &json) {
    return bset.to_json(json);
  }

  /*
  // No reading functions for now
  void from_json(BasisSet &bset, const jsonParser &json) {
    return bset.from_json(json);
  }
  */

}

