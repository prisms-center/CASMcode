#include "casm/crystallography/SupercellEnumerator.hh"

#include "casm/external/boost.hh"
#include "casm/external/Eigen/Dense"

#include "casm/crystallography/Structure.hh"

namespace CASM {

  /**
   * If you really need a n x m (vs n x n) matrix, or to allow
   * zeros along the longest diagonal, this class won't work.
   */

  HermiteCounter::HermiteCounter(int init_start_determinant, int init_dim):
    m_pos(0) {
    m_diagonal = Eigen::VectorXi::Ones(init_dim);
    m_diagonal(0) = init_start_determinant;

    m_upper_tri = HermiteCounter_impl::_upper_tri_counter(m_diagonal);
  }

  HermiteCounter::Index HermiteCounter::position() const {
    return m_pos;
  }

  Eigen::VectorXi HermiteCounter::diagonal() const {
    return m_diagonal;
  }

  Eigen::MatrixXi HermiteCounter::current() const {
    return HermiteCounter_impl::_zip_matrix(m_diagonal, m_upper_tri.current());
  }

  HermiteCounter::Index HermiteCounter::dim() const {
    return m_diagonal.size();
  }

  Eigen::MatrixXi HermiteCounter::operator()() const {
    return current();
  }

  HermiteCounter::value_type HermiteCounter::determinant() const {
    return m_diagonal.prod();
  }

  void HermiteCounter::reset_current() {
    jump_to_determinant(determinant());
    return;
  }

  void HermiteCounter::next_determinant() {
    jump_to_determinant(determinant() + 1);
    return;
  }

  Index HermiteCounter::_increment_diagonal() {
    m_pos = HermiteCounter_impl::next_spill_position(m_diagonal, m_pos);
    return m_pos;
  }

  HermiteCounter &HermiteCounter::operator++() {
    //Vary elements above diagonal
    m_upper_tri++;

    //Nothing to permute for identity
    if(determinant() == 1) {
      jump_to_determinant(2);
    }

    //Still other matrices available for current diagonal
    else if(m_upper_tri.valid()) {
      return *this;
    }

    //Find next diagonal, reset elements above diagonal
    else if(_increment_diagonal() < dim()) {
      m_upper_tri = HermiteCounter_impl::_upper_tri_counter(m_diagonal);
    }

    //Reset diagonal with next determinant value, reset elements above diagonal
    else {
      jump_to_determinant(determinant() + 1);
    }

    return *this;
  }

  void HermiteCounter::jump_to_determinant(value_type new_det) {
    if(new_det < 1) {
      throw std::runtime_error("Determinants of Hermite normal form matrices must be greater than 0");
    }

    m_diagonal = Eigen::VectorXi::Ones(dim());
    m_diagonal(0) = new_det;

    m_upper_tri = HermiteCounter_impl::_upper_tri_counter(m_diagonal);
    m_pos = 0;

    return;
  }

  namespace HermiteCounter_impl {
    HermiteCounter::Index _spill_factor(Eigen::VectorXi &diag, HermiteCounter::Index position, HermiteCounter::value_type attempt) {
      //If you fall into these traps you're using this wrong
      assert(attempt <= diag(position));
      assert(position < diag.size() - 1);
      assert(diag(position + 1) == 1);
      assert(diag(position) > 1);
      assert(attempt >= 2);


      //Use the given attempt as a starting guess for factorizing the element at position,
      //but if that doesn't work, keep increasing the value of the attempt until it does.
      while(diag(position) % attempt != 0) {
        attempt++;
      }

      diag(position) = diag(position) / attempt;
      diag(position + 1) = attempt;

      position++;

      return position;
    }

    HermiteCounter::Index next_spill_position(Eigen::VectorXi &diag, HermiteCounter::Index position) {
      HermiteCounter::value_type attempt = 2;

      //If you reached the end of the diagonal, backtrack to nearest non-1 value
      //and perform the next spill
      if(position == diag.size() - 1) {
        do {
          position--;
        }
        while(position >= 0 && diag(position) == 1);

        //You're at the last spill already. Your diagonal is 1 1 ... 1 1 det. No more increments are possible.
        if(position < 0) {
          return diag.size();
        }

        //Flush everything to the right of position with ones, and reset the value at position
        //with the next attempt for factorization
        attempt = diag(diag.size() - 1) + 1;
        diag(position) = diag(position) * diag(diag.size() - 1);
        diag(diag.size() - 1) = 1;
      }

      return _spill_factor(diag, position, attempt);
    }

    HermiteCounter::Index upper_size(HermiteCounter::Index init_dim) {
      assert(init_dim > 0);

      HermiteCounter::Index tritotal = 0;

      for(int i = 0; i < init_dim; i++) {
        tritotal += i;
      }

      return tritotal;
    }

    EigenVectorXiCounter _upper_tri_counter(const Eigen::VectorXi &current_diag) {
      //Find out how many slots you need for all the elements above the diagonal
      HermiteCounter::Index uppersize = upper_size(current_diag.size());

      //Start the counter with zero everywhere
      Eigen::VectorXi begincount(Eigen::VectorXi::Zero(uppersize));

      //Increments always go one at a time
      Eigen::VectorXi stepcount(Eigen::VectorXi::Ones(uppersize));

      //The m_upper_tri is unrolled left to right, top to bottom
      Eigen::VectorXi endcount(Eigen::VectorXi::Zero(uppersize));
      HermiteCounter::Index slot = 0;
      for(HermiteCounter::Index i = 0; i < current_diag.size(); i++) {
        for(HermiteCounter::Index j = i + 1; j < current_diag.size(); j++) {
          //The counter value should always be smaller than the diagonal
          endcount(slot) = current_diag(i) - 1;
          slot++;
        }
      }

      return EigenVectorXiCounter(begincount, endcount, stepcount);
    }

    Eigen::MatrixXi _zip_matrix(const Eigen::VectorXi &current_diag, const Eigen::VectorXi &current_upper_tri) {
      assert(current_upper_tri.size() == upper_size(current_diag.size()));

      Eigen::MatrixXi hmat(current_diag.asDiagonal());

      Index slot = 0;
      for(Index i = 0; i < current_diag.size() - 1; i++) {
        for(Index j = i + 1; j < current_diag.size(); j++) {
          hmat(i, j) = current_upper_tri(slot);
          slot++;
        }
      }

      return hmat;
    }

    /**
     * A particular Hermit normal form matrix can have its dimensions expanded by inserting
     * extra 0s and a 1 in the diagonal. Specify the "active" indexes (n values of 1) and
     * "inactive" indexes (m-n values of 0) in a vector.
     *
     * For example, you can expand a 4 x 4 matrix
     *  2 1 0 0
     *  0 5 2 2
     *  0 0 6 1
     *  0 0 0 3
     *
     * into a 6 x 6 matrix by specifying the vector [101110]
     *  2 0 1 0 0 0
     *  0 1 0 0 0 0
     *  0 0 5 2 2 0
     *  0 0 0 6 1 0
     *  0 0 0 0 3 0
     *  0 0 0 0 0 1
     *
     * Rows/columns where a 0 was specified are all zero except for the diagonal element, which is 1.
     */

    Eigen::MatrixXi _expand_dims_old(const Eigen::MatrixXi &hermit_mat, const Eigen::VectorXi &active_dims) {
      assert(hermit_mat.rows() == active_dims.sum() && hermit_mat.cols() == active_dims.sum());
      assert(active_dims.maxCoeff() == 1 && active_dims.minCoeff() == 0);

      Eigen::MatrixXi expanded(Eigen::MatrixXi::Identity(active_dims.size(), active_dims.size()));

      HermiteCounter::Index i, j, si, sj;
      si = -1;
      sj = -1;

      for(i = 0; i < expanded.rows(); i++) {
        if(active_dims(i) == 0) {
          continue;
        }

        else {
          si++;
        }


        for(j = 0; j < expanded.cols(); j++) {
          if(active_dims(j) == 0) {
            continue;
          }

          else {
            sj++;
          }

          expanded(i, j) = hermit_mat(si, sj);
        }

        j = 0;
        sj = -1;
      }

      return expanded;
    }

    /*
     * This is a generalized way to insert new dimensions into a HNF matrix.
     * If there is a n x n HNF matrix H, it will be expanded to a new non-HNF
     * m x m matrix E through a m x m transformation matrix T.
     *
     * The first n columns of T will be multiplied with the values of H, while
     * the last m-n columns of T will remain unaffected. For example, if you
     * are counting over a 2x2 H matrix and currently have
     *
     * 2 1
     * 0 3
     *
     * Then you can specify that you want a 3x3 matrix E so that the values of
     * H work on your lattice vectors a and c, but not b with a T matrix
     *
     * 1 0 0
     * 0 0 1
     * 0 1 0
     *
     * The resulting expanded matrix E is now
     *
     * 2 1 0
     * 0 0 1
     * 0 3 0
     *
     * This is achieved by converting H to a block matrix of dimensions m x m,
     * which is an identity matrix with the upper left block equal to H. For
     * the example above B would be
     *
     * 2 1 0
     * 0 3 0
     * 0 0 1
     *
     * This way E=T*B
     *
     * You may specify arbitrary combinations of vectors in the columns of T that
     * H will work on.
     *
     * Note that the resulting matrix will probably *NOT* retain it's Hermite normal
     * form. For use in the SupercellEnumerator class, the order within the first n
     * vectors and the order within the last m-n vectors will not affect your enumerations.
     */

    Eigen::MatrixXi _expand_dims(const Eigen::MatrixXi &H, const Eigen::MatrixXi &T) {
      assert(H.rows() == H.cols());
      assert(T.rows() == T.cols());
      assert(T.rows() >= H.rows());

      Index n = H.rows();
      Index m = T.rows();

      //First convert H into a block matrix with dimensions m x m, the H block is on the upper left
      Eigen::MatrixXi B(m, m);
      Eigen::MatrixXi I_block(Eigen::MatrixXi::Identity(m - n, m - n));
      Eigen::MatrixXi Z_block(Eigen::MatrixXi::Zero(n, m - n));
      B << H, Z_block, Z_block.transpose(), I_block;

      return T * B;
    }


    /**
     * When comparing two Hermite normal form matrices, the canonical one is determined on which one
     * has larger elements along a particular order. This routine unrolls a matrix in the following
     * manner:
     *
     * 1 6 5
     * 0 2 4
     * 0 0 3
     *
     * But also works for higher dimensional matrices, such as
     *
     * 1  12 11 10 9
     * 0  2  13 15 8
     * 0  0  3  14 7
     * 0  0  0  4  6
     * 0  0  0  0  5
     */

    Eigen::VectorXi _canonical_unroll(const Eigen::MatrixXi &hermit_mat) {
      assert(hermit_mat.rows() == hermit_mat.cols());

      int dims = hermit_mat.rows();
      int unrolldim = 0;

      for(int i = dims; i > 0; i--) {
        unrolldim += i;
      }

      Eigen::VectorXi unrolled_herm(unrolldim);

      Index lasti, lastj;
      lasti = lastj = -1;

      enum STEP {DOWN, UP, LEFT};
      STEP curr_step = DOWN;

      Index unroll_ind = 0;

      for(Index i = 0; i < dims; i++) {
        for(Index j = 0; j < dims - i; j++) {
          switch(curr_step) {
          case DOWN:
            lasti = lasti + 1;
            lastj = lastj + 1;
            break;
          case UP:
            lasti = lasti - 1;
            break;
          case LEFT:
            lastj = lastj - 1;
            break;
          }

          unrolled_herm(unroll_ind) = hermit_mat(lasti, lastj);
          unroll_ind++;
        }

        switch(curr_step) {
        case DOWN:
          curr_step = UP;
          break;
        case UP:
          curr_step = LEFT;
          break;
        case LEFT:
          curr_step = DOWN;
          break;
        }
      }

      return unrolled_herm;
    }

  }


  //******************************************************************************************************************//


  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    double tol,
                                                    size_type begin_volume,
                                                    size_type end_volume,
                                                    int init_dims,
                                                    Eigen::Matrix3i init_trans_mat) :
    m_unit(unit),
    m_lat(unit),
    m_begin_volume(begin_volume),
    m_end_volume(end_volume),
    m_dims(init_dims),
    m_trans_mat(init_trans_mat) {

    if(m_trans_mat.determinant() < 1) {
      throw std::runtime_error("The transformation matrix to expand into a 3x3 matrix must have a positive determinant!");
    }

    m_lat.generate_point_group(m_point_group, tol);

  }

  template<>
  SupercellEnumerator<Lattice>::SupercellEnumerator(Lattice unit,
                                                    const SymGroup &point_grp,
                                                    size_type begin_volume,
                                                    size_type end_volume,
                                                    int init_dims,
                                                    Eigen::Matrix3i init_trans_mat) :
    m_unit(unit),
    m_lat(unit),
    m_point_group(point_grp),
    m_begin_volume(begin_volume),
    m_end_volume(end_volume),
    m_dims(init_dims),
    m_trans_mat(init_trans_mat)

  {
    if(m_trans_mat.determinant() < 1) {
      throw std::runtime_error("The transformation matrix to expand into a 3x3 matrix must have a positive determinant!");
    }

  }



  template<>
  Eigen::Matrix3i enforce_min_volume<Lattice>(
    const Lattice &unit,
    const Eigen::Matrix3i &T,
    const SymGroup &point_grp,
    Index volume,
    bool fix_shape) {

    if(fix_shape) {
      auto init_vol = T.determinant();
      Index m = 1;
      while(m * m * m * init_vol < volume) {
        ++m;
      }

      return m * Eigen::Matrix3i::Identity();
    }
    else {
      auto init_vol = T.determinant();
      Index m = 1;
      while(m * init_vol < volume) {
        ++m;
      }

      auto compactness = [](const Lattice & lat) {
        Eigen::Matrix3d L = lat.lat_column_mat();
        return (L.transpose() * L).trace();
      };

      auto compare = [&](const Lattice & A, const Lattice & B) {
        return compactness(A) < compactness(B);
      };

      SupercellEnumerator<Lattice> scel(unit, point_grp, m, m + 1);
      auto best_it = std::min_element(scel.begin(), scel.end(), compare);
      return best_it.matrix();
    }

  }


  /// \brief Return canonical hermite normal form of the supercell matrix, and op used to find it
  ///
  /// \returns std::pair<Eigen::Matrix3i, Eigen::Matrix3d> of H in canonical form, and op used to find it from T
  ///
  /// \param T a supercell matrix (Eigen::Matrix3i), such that S = U*T,
  ///          where S is the superlattice and U the unit lattice, as column vector matrices
  /// \param unitcell the unit BasicStructure<Site>
  ///
  /// Canonical form is such that T is in hermite normal form (as from casm::hermite_normal_form),
  ///  and the unrolled coefficients
  /// \code
  /// [a f e]
  /// [0 b d] -> abcdef
  /// [0 0 c]
  /// \endcode
  /// form the highest lexicographic order when considering equivalent superlattices by point group operations.
  ///
  /// - Equivalent superlattices can be obtained using point group operations: S' = op*S = U*H*V,
  ///   where V is integer and has determinant +/- 1
  /// - Substituting S = U*T, we have op*U*T = U*H*V.
  /// - Or H*V = U.inverse*op*U*T, which is the hermite normal form of U.inverse*op*U*T
  /// - So T is canonical if it is in hermite normal form and for all point group operations
  ///   it has a higher lexicographic order than the resulting H
  ///
  ///
  /// \relatesalso Lattice
  ///
  std::pair<Eigen::MatrixXi, Eigen::MatrixXd> canonical_hnf(const Eigen::MatrixXi &T, const SymGroup &effective_pg, const Lattice &ref_lattice) {

    Eigen::Matrix3d lat = ref_lattice.lat_column_mat();

    Index i_canon = -1;

    //get T in hermite normal form
    //H is the canonical form of the initial T matrix
    const Eigen::Matrix3i H = hermite_normal_form(T).first;

    //H_best will be the most canonical version and is returned
    Eigen::Matrix3i H_best;
    H_best = H;

    Eigen::VectorXi unrolled_H_best = HermiteCounter_impl::_canonical_unroll(H_best);

    for(Index i = 0; i < effective_pg.size(); i++) {
      Eigen::Matrix3i transformed = iround(lat.inverse() * effective_pg[i].matrix() * lat) * H;
      Eigen::Matrix3i H_transformed = hermite_normal_form(transformed).first;
      Eigen::VectorXi unrolled_H_transformed = HermiteCounter_impl::_canonical_unroll(H_transformed);

      // canonical only if H_best is '>=' H, for all H, so if H '>' m_current, make H the H_best
      for(Index j = 0; j < unrolled_H_best.size(); j++) {
        if(unrolled_H_transformed(j) > unrolled_H_best(j)) {
          i_canon = i;
          H_best = H_transformed;
          unrolled_H_best = HermiteCounter_impl::_canonical_unroll(H_best);
          break;
        }
      }
    }

    return std::make_pair<Eigen::MatrixXi, Eigen::MatrixXd>(H_best, effective_pg[i_canon].matrix());
  }

}
