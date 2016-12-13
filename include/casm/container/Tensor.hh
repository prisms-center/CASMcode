#ifndef TENSOR_HH
#define TENSOR_HH

#include <iostream>
#include <cmath>

#include "casm/container/Permutation.hh"
#include "casm/container/Template_Algorithms.hh"
#include "casm/container/Counter.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  template<class T>
  class DynamicMatrix;

  template<class T>
  class ReturnTensor;

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  /** \defgroup Tensor
   *
   *  \ingroup Container
   *  \brief Data structure representing tensors
   *  @{
  */

  template<class T>
  class Tensor : public Array<T> {
    Index Nrank;
    Array<Index> Ndim, idim, max_ind;

  public:
    using Array<T> :: resize;
    using Array<T> :: size;
    using Array<T> :: at;

    //These constructors take two arguments even though they really just need one, since rank==dim.size() if
    //you're doing it right.
    Tensor(Index tNrank, Array<Index> tNdim) : Nrank(tNrank), Ndim(tNdim) {
      Index N = 1;
      idim.resize(Nrank);
      max_ind.resize(Nrank);
      for(Index i = 0; i < Nrank; i++) {
        max_ind[i] = Ndim[i] - 1;
        idim[i] = N;
        N *= Ndim[i];
      }

      resize(N);
    };

    Tensor(Index tNrank, Array<Index> tNdim, T init_val) : Nrank(tNrank), Ndim(tNdim) {
      Index N = 1;
      idim.resize(Nrank);
      max_ind.resize(Nrank);
      for(Index i = 0; i < Nrank; i++) {
        max_ind[i] = Ndim[i] - 1;
        idim[i] = N;
        N *= Ndim[i];
      }

      resize(N, init_val);
    };

    Tensor(Index tNrank = 0) : Nrank(tNrank), Ndim(tNrank, 3) {
      Index N = 1;
      idim.resize(Nrank);
      max_ind.resize(Nrank);
      for(Index i = 0; i < Nrank; i++) {
        max_ind[i] = Ndim[i] - 1;
        idim[i] = N;
        N *= Ndim[i];
      }
      resize(N);
    };

    Tensor(const Eigen::Matrix<T, 3, 3> &mat_init) : Nrank(2), Ndim(2, 3) {
      Index N = 1;
      idim.resize(Nrank);
      max_ind.resize(Nrank);
      for(Index i = 0; i < Nrank; i++) {
        max_ind[i] = Ndim[i] - 1;
        idim[i] = N;
        N *= Ndim[i];
      }
      resize(N);
      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
          (*this)(i, j) = mat_init(i, j);
    }

    Tensor(ReturnTensor<T> &init_cont) {
      swap(init_cont);
    };

    Tensor &operator=(ReturnTensor<T> &RHS) {
      swap(RHS);
      return *this;
    };

    Tensor &operator=(const T &RHS) {
      for(Index i = 0; i < size(); i++)
        at(i) = RHS;
      return *this;
    }

    template<int dim1, int dim2, int flag1>
    Tensor &operator=(const Eigen::Matrix<T, dim1, dim2, flag1> &LHS);

    static ReturnTensor<T> identity(Index trank, Index tdim);

    void swap(Tensor &RHS);
    void redefine(Array<Index> _Ndim);
    void redefine(Array<Index> _Ndim, T fill_val);

    Tensor &dim_permute(const Array<Index> &iperm);
    Tensor &dim_permute(const Permutation &perm);

    Tensor &dim_ipermute(const Array<Index> &perm);
    Tensor &dim_ipermute(const Permutation &perm);

    Tensor &dim_unpermute();
    bool next_permute();
    /// Symmetrize tensor with respect to index permutation
    void permute_symmetrize();
    /// Symmetrize tensor with respect to index permutation if not all dimensions are quivalent
    /// unique_dim lists each subset of indices which are mutually equivalent
    void permute_symmetrize(const Array< Array<Index> > &unique_dim);
    void reset_idim();

    Counter<Array<Index> > element_counter() const;

    Index rank() const;
    const Array<Index> &ind_max() const;
    const Array<Index> &dim()const;
    Index dim(Index i) const;

    const Array<Index> &mult_array() const;

    const T &get(const Array<Index> &inds) const;
    T &at(const Array<Index> &inds);
    const T &at(const Array<Index> &inds) const;
    T &operator()(const Array<Index> &inds);
    const T &operator()(const Array<Index> &inds) const;
    T &operator()(Index i, Index j);
    const T &operator()(Index i, Index j) const;

    Tensor &operator+=(const Tensor &RHS);
    Tensor &operator-=(const Tensor &RHS);
    Tensor &operator*=(T RHS);
    Tensor &operator/=(T RHS);

    ReturnTensor<T> operator-();
    ReturnTensor<T> operator+(const Tensor &RHS);
    ReturnTensor<T> operator-(const Tensor &RHS);

    //We should find a way to do tensor contraction/multiplication

    T scalar_prod(const Tensor &RHS) const;
    ReturnTensor<T> tensor_prod(const Tensor &RHS);
    Tensor &normalize();

    Tensor &transform(const Eigen::MatrixXd &op);

    //rep_IDs specify the representations that transform each dimension of tensor
    Tensor &transform(const SymOp &op, Array<SymGroupRepID> rep_IDs);

    Tensor &apply_sym(const SymOp &op);

    Tensor &symmetrize_index(const SymOp &op, Array<SymGroupRepID> rep_IDs, Array<Index> inner_ind);

    Tensor slice(Array<Index> slice_Ndim, Array<Index> slice_ind);

    void read(std::istream &stream); //Added by Ivy -- TODO: Make more general for not just rank 2 tensors

    //Tensor<T> matrix_multiply(Vector3<T> &tvec);

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> convert_to_Eigen() const;

  };


  template<class T>
  Tensor<T> operator*(const T &LHS, const Tensor<T> &RHS);

  template<class T>
  Tensor<T> operator*(const SymOp &LHS, const Tensor<T> &RHS);

  template<class T>
  std::ostream &operator << (std::ostream &stream, const Tensor<T> &RHS);

  template<class T>
  std::istream &operator >> (std::istream &stream, Tensor<T> &RHS);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<class T>
  class ReturnTensor : public Tensor<T> {

  public:
    using Tensor<T>::swap;

    ReturnTensor(Tensor<T> &init_tens) {
      swap(init_tens);
    }

    ReturnTensor &operator=(Tensor<T> &RHS) {
      swap(RHS);
      return *this;
    };

  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  template<class T>
  class TensorBasis : public Array< Tensor<T> > {
  public:
    using Array< Tensor<T> > :: size;
    using Array< Tensor<T> > :: at;
    using Array< Tensor<T> > :: clear;

    Array<double> coeffs;

    double eci(Index i) const {
      return coeffs[i];
    };
    double &eci(Index i) {
      return coeffs[i];
    }; //Added by Ivy 09/25/2012

    void generate_basis(Index Nrank, const SymGroup &sym_group);
    void generate_basis(Index Nrank, const SymGroup &sym_group, Index Rep_ID);
    void generate_basis(Index Nrank, const SymGroup &sym_group, const SymGroupRep &perm_group);

    void make_orthogonal_to(const TensorBasis &ortho_basis);
    TensorBasis &apply_sym(const SymOp &op);
    void normalize();
    void idealize();
    bool read(std::istream &stream); //Added by Ivy
  };

  //************************************************************

  template<class T>
  void Tensor<T>::swap(Tensor<T> &RHS) {
    Index trank(Nrank);
    Nrank = RHS.Nrank;
    RHS.Nrank = trank;
    Ndim.swap(RHS.Ndim);
    idim.swap(RHS.idim);
    max_ind.swap(RHS.max_ind);
    Array<T>::swap(RHS);
    return;
  }

  //************************************************************
  template< class T >
  template<int dim1, int dim2, int flag1>
  Tensor<T> &Tensor<T>::operator=(const Eigen::Matrix<T, dim1, dim2, flag1> &LHS) {
    Array<Index> tdim(2);
    if(rank() != 2 || dim(0) != LHS.rows() || dim(1) != LHS.cols()) {
      tdim[0] = LHS.rows();
      tdim[1] = LHS.cols();
      redefine(tdim);
    }
    for(int i = 0; i < LHS.rows(); i++) {
      tdim[0] = i;
      for(int j = 0; j < LHS.rows(); j++) {
        tdim[1] = j;
        at(tdim) = LHS(i, j);
      }
    }
    return (*this);
  }


  //************************************************************

  template<class T>
  void Tensor<T>::redefine(Array<Index> _Ndim) {
    Nrank = _Ndim.size();
    Ndim = _Ndim;
    Index N = 1;
    idim.resize(Nrank);
    max_ind.resize(Nrank);
    for(Index i = 0; i < Nrank; i++) {
      max_ind[i] = Ndim[i] - 1;
      idim[i] = N;
      N *= Ndim[i];
    }

    resize(N);
  }

  //************************************************************

  template<class T>
  void Tensor<T>::redefine(Array<Index> _Ndim, T fill_val) {
    Nrank = _Ndim.size();
    Ndim = _Ndim;
    Index N = 1;
    idim.resize(Nrank);
    max_ind.resize(Nrank);
    for(Index i = 0; i < Nrank; i++) {
      max_ind[i] = Ndim[i] - 1;
      idim[i] = N;
      N *= Ndim[i];
    }

    resize(N, fill_val);
  }

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::dim_permute(const Array<Index> &iperm) {
    if(iperm.size() != Nrank) {
      // temporarily changed this to cout
      std::cout << "WARNING: In Tensor<T>::dim_permute, permutation array is incompatible with tensor rank! Continuing without permuting indeces...\n";
      return *this;
    }

    idim.permute(iperm);
    Ndim.permute(iperm);
    max_ind.permute(iperm);

    return *this;
  }

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::dim_permute(const Permutation &perm) {
    return dim_permute(perm.perm_array());
  }
  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::dim_ipermute(const Array<Index> &iperm) {
    if(iperm.size() != Nrank) {
      std::cerr << "WARNING: In Tensor<T>::dim_permute, permutation array is incompatible with tensor rank! Continuing without permuting indeces...\n";
      return *this;
    }

    idim.ipermute(iperm);
    Ndim.ipermute(iperm);
    max_ind.ipermute(iperm);

    return *this;
  }

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::dim_ipermute(const Permutation &perm) {
    return dim_ipermute(perm.perm_array());
  }
  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::dim_unpermute() {
    Array<Index> tperm;

    idim.sort(tperm);
    Ndim.permute(tperm);
    max_ind.permute(tperm);

    return *this;
  }

  //************************************************************

  template<class T>
  bool Tensor<T>::next_permute() {
    Index i = Nrank - 2;
    bool is_valid = true;
    Array<Index> iperm;
    for(Index k = 0; k < Nrank; k++)
      iperm[k] = k;

    while(valid_index(i) && idim[i + 1] <= idim[i])
      i--;


    Index j = Nrank - 1;

    if(valid_index(i)) {
      while(idim[j] <= idim[i])
        j--;

      idim.swap_elem(i, j);
      Ndim.swap_elem(i, j);
      max_ind.swap_elem(i, j);

      i++;
      j = Nrank - 1;
    }

    else {
      is_valid = false;
      i = 0;
    }

    while(i < j && valid_index(j)) {
      idim.swap_elem(i, j);
      Ndim.swap_elem(i, j);
      max_ind.swap_elem(i, j);

      i++;
      j--;
    }

    return is_valid;
  }

  //************************************************************

  template<class T>
  void Tensor<T>::reset_idim() {
    if(idim.is_ascending())
      return;

    Tensor ttens(Nrank, Ndim);
    Counter<Array<Index> > icount(ttens.element_counter());

    do {
      ttens(icount()) = at(icount());
    }
    while((++icount).valid());

    swap(ttens);

    return;
  }

  //************************************************************

  template<class T>
  void Tensor<T>::permute_symmetrize() {
    reset_idim();
    Tensor<T> ttens(*this);
    Index nperm = 1;
    while(ttens.next_permute()) {
      (*this) += ttens;
      nperm++;
    }
    (*this) /= nperm;
    return;
  }

  //************************************************************

  template<class T>
  void Tensor<T>::permute_symmetrize(const Array< Array<Index> > &unique_dim) {
    reset_idim();
    Array<Array<Index> > tunique(unique_dim);
    for(Index i = 0; i < tunique.size(); i++)
      tunique[i].sort();
    Array<Index> tperm = Ndim;
    Tensor<T> ttens(*this);
    Index iset(0), i(0), j(0), nperm(1);
    (*this) = 0;
    while(iset < tunique.size()) {
      for(i = 0; i < tunique.size(); i++)
        for(j = 0; j < tunique[i].size(); j++)
          tperm[unique_dim[i][j]] = tunique[i][j];
      (*this) += ttens.dim_permute(tperm);
      ttens.dim_unpermute();
      nperm++;
      iset = 0;
      while(iset < tunique.size() && !tunique[iset].next_permute())
        iset++;
    }
    (*this) /= double(nperm);

  }

  //************************************************************

  template<class T>
  Counter<Array<Index> > Tensor<T>::element_counter() const {
    return Counter<Array<Index> >(Array<Index>(Nrank, 0), max_ind, Array<Index>(Nrank, 1));
  }

  //************************************************************

  template<class T>
  const Array<Index> &Tensor<T>::ind_max()const {
    return max_ind;
  }

  //************************************************************

  template<class T>
  const Array<Index> &Tensor<T>::dim()const {
    return Ndim;
  }

  //************************************************************

  template<class T>
  Index Tensor<T>::dim(Index i) const {
    return Ndim[i];
  }

  //************************************************************

  template<class T>
  const Array<Index> &Tensor<T>::mult_array() const {
    return idim;
  }

  //************************************************************

  template<class T>
  Index Tensor<T>::rank() const {
    return Nrank;
  }

  //************************************************************

  template<class T>
  const T &Tensor<T>::get(const Array<Index> &inds) const {
    Index lin_index = 0;
    for(Index i = 0; i < Nrank; i++)
      lin_index += inds[i] * idim[i];

    return at(lin_index);
  };

  //************************************************************

  template<class T>
  T &Tensor<T>::at(const Array<Index> &inds) {
    Index lin_index = 0;
    for(Index i = 0; i < Nrank; i++)
      lin_index += inds[i] * idim[i];

    return at(lin_index);
  };

  //************************************************************

  template<class T>
  const T &Tensor<T>::at(const Array<Index> &inds) const {
    Index lin_index = 0;
    for(Index i = 0; i < Nrank; i++)
      lin_index += inds[i] * idim[i];

    return at(lin_index);
  };

  //************************************************************

  template<class T>
  T &Tensor<T>::operator()(const Array<Index> &inds) {
    Index lin_index = 0;
    for(Index i = 0; i < Nrank; i++)
      lin_index += inds[i] * idim[i];

    return at(lin_index);
  };

  //************************************************************

  template<class T>
  const T &Tensor<T>::operator()(const Array<Index> &inds) const {
    Index lin_index = 0;
    for(Index i = 0; i < Nrank; i++)
      lin_index += inds[i] * idim[i];

    return at(lin_index);
  };
  //************************************************************

  template<class T>
  T &Tensor<T>::operator()(Index i, Index j) {
    return at(i * idim[0] + j * idim[1]);
  };

  //************************************************************

  template<class T>
  const T &Tensor<T>::operator()(Index i, Index j) const {
    return at(i * idim[0] + j * idim[1]);
  };

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::operator+=(const Tensor &RHS) {
    if(Ndim != RHS.Ndim) {
      std::cerr << "Attempting to add two incompatible tensors!\n";
      exit(1);
    }

    if(idim == RHS.idim) {
      for(Index i = 0; i < size(); i++)
        at(i) += RHS[i];

      return *this;
    }


    Counter<Array<Index> > icount(Array<Index>(Nrank, 0), max_ind, Array<Index>(Nrank, 1));
    do {
      at(icount()) += RHS(icount());

    }
    while((++icount).valid());

    return *this;
  }

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::operator-=(const Tensor &RHS) {
    if(Ndim != RHS.Ndim) {
      std::cerr << "Attempting to add two incompatible tensors!\n";
      exit(1);
    }

    if(idim == RHS.idim) {
      for(Index i = 0; i < size(); i++)
        at(i) -= RHS[i];

      return *this;
    }

    Counter<Array<Index> > icount(Array<Index>(Nrank, 0), max_ind, Array<Index>(Nrank, 1));
    do {
      at(icount()) -= RHS(icount());
    }
    while((++icount).valid());


    return *this;
  }

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::operator*=(T RHS) {
    for(Index i = 0; i < size(); i++)
      at(i) *= RHS;

    return *this;

  }

  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::operator/=(T RHS) {
    for(Index i = 0; i < size(); i++)
      at(i) /= RHS;

    return *this;

  }

  //************************************************************

  template<class T>
  ReturnTensor<T> Tensor<T>::operator-() {
    Tensor<T> ttens(*this);
    for(Index i = 0; i < ttens.size(); i++)
      ttens[i] = -at(i);
    return ReturnTensor<T>(ttens);
  }

  //************************************************************

  template<class T>
  ReturnTensor<T> Tensor<T>::operator+(const Tensor &RHS) {
    return ReturnTensor<T>(Tensor<T>(*this) += RHS);
  }

  //************************************************************

  template<class T>
  ReturnTensor<T> Tensor<T>::operator-(const Tensor &RHS) {
    return ReturnTensor<T>(Tensor<T>(*this) -= RHS);
  }

  //************************************************************

  template<class T>
  T Tensor<T>::scalar_prod(const Tensor &RHS) const {
    if(Ndim != RHS.Ndim) {
      std::cerr << "Attempting to take scalar product of two incompatible tensors!\n";
      assert(0);
      //exit(1);
    }

    T tprod = 0;

    if(idim == RHS.idim) {
      for(Index i = 0; i < size(); i++)
        tprod += at(i) * RHS[i];
    }

    else {
      Counter<Array<Index> > icount(Array<Index>(Nrank, 0), max_ind, Array<Index>(Nrank, 1));
      do {
        tprod += at(icount()) * RHS(icount());
      }
      while((++icount).valid());
    }

    return tprod;

  }

  //************************************************************

  template<class T>
  ReturnTensor<T> Tensor<T>::tensor_prod(const Tensor &RHS) {
    Tensor<T> ttens(rank() + RHS.rank(), array_cat(dim(), RHS.dim()), 0);

    Counter<Array<Index> > icount(element_counter()), jcount(RHS.element_counter());

    Index i, j, tot_ind;
    do {
      jcount.reset();
      do {
        tot_ind = 0;
        for(i = 0; i < icount().size(); i++)
          tot_ind += icount[i] * ttens.idim[i];

        for(j = 0; j < jcount().size(); j++)
          tot_ind += jcount[j] * ttens.idim[i + j];

        ttens[tot_ind] = at(icount()) * RHS(jcount());

      }
      while((++jcount).valid());
    }
    while((++icount).valid());

    return ReturnTensor<T>(ttens);
  }


  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::normalize() {

    return (*this) /= sqrt(scalar_prod(*this));
  }

  //************************************************************

  template<class T>
  ReturnTensor<T> Tensor<T>::identity(Index trank, Index tdim) {
    Tensor<T> ttens(trank, Array<Index>(trank, tdim), 0);
    for(Index i = 0; i < tdim; i++)
      ttens(Array<Index>(trank, i)) = 1;

    return ReturnTensor<T>(ttens);
  }

  //************************************************************
  template<class T>
  Tensor<T> &Tensor<T>::transform(const SymOp &op, Array<SymGroupRepID> rep_IDs) {
    if(!rank()) return *this;
    Array<Eigen::MatrixXd const *> rep_mats(op.get_matrix_reps(rep_IDs));

    Tensor<T> ttens(Nrank, Ndim);
    Index nd;

    T tval;
    Counter<Array<Index> > icount(element_counter());
    Counter<Array<Index> > jcount(icount);
    do {
      ttens(icount()) = 0;
      jcount.reset();
      do {
        tval = at(jcount());
        for(nd = 0; nd < Nrank; nd++)
          tval *= (*rep_mats[nd])(icount[nd], jcount[nd]);
        ttens(icount()) += tval;

      }
      while((++jcount).valid());

    }
    while((++icount).valid());
    swap(ttens);
    return *this;
  }

  //************************************************************
  template<class T>
  Tensor<T> &Tensor<T>::symmetrize_index(const SymOp &op, Array<SymGroupRepID> rep_IDs, Array<Index> inner_ind) {
    if(!rank()) return *this;
    Array<Eigen::MatrixXd const *> rep_mats(op.get_matrix_reps(rep_IDs));

    Index nd;

    T tval;
    Counter<Array<Index> > icount(element_counter());
    do {
      tval = 1.0;
      for(nd = 0; nd < Nrank; nd++)
        tval *= (*rep_mats[nd])(icount[nd], inner_ind[nd]);
      at(icount()) += tval;

    }
    while((++icount).valid());
    return *this;
  }

  //************************************************************
  template<class T>
  Tensor<T> &Tensor<T>::apply_sym(const SymOp &op) {
    if(!rank()) return *this;

    Tensor<T> ttens(Nrank, Ndim);
    Index nd;

    T tval;
    Counter<Array<Index> > icount(element_counter());
    Counter<Array<Index> > jcount(icount);
    do {
      ttens(icount()) = 0;
      jcount.reset();
      do {
        tval = at(jcount());
        for(nd = 0; nd < Nrank; nd++)
          tval *= op.matrix()(icount[nd], jcount[nd]);
        ttens(icount()) += tval;

      }
      while((++jcount).valid());

    }
    while((++icount).valid());
    swap(ttens);
    return *this;
  }

  //************************************************************
  //************************************************************
  template<class T>
  void Tensor<T>::read(std::istream &stream) {

    double val;
    Array<Index> min(rank(), 0), max(ind_max()), inc(rank(), 1);
    Counter<Array<Index> > icount(min, max, inc);

    do {
      stream >> val;
      at(icount()) = val;
    }
    while((++icount).valid());

    return;
  };

  //************************************************************
  /**
   * A temporary hack for multiplying tensors -- This only
   * works for a rank 2 tensor of dimension 3x3 multiplying
   * a Vector3.  It returns a rank 2 tensor of dimension 3x1
   */
  //************************************************************
  /*template<class T>
  Tensor<T> Tensor<T>::matrix_multiply(Vector3<T> &tvec) {

    Index j = 0;
    Array<Index> tensordim(2, 0);
    tensordim[0] = 3;
    tensordim[1] = 1;
    Tensor<T> result(2, tensordim, 0);

    for(Index i = 0; i < size(); i++) {
      if(i % 3 == 0) {
        j++;
      }
      result.at(j - 1) += at(i) * tvec.at(i % 3);
    }

    return result;

  };
  */
  //************************************************************

  template<class T>
  Tensor<T> &Tensor<T>::transform(const Eigen::MatrixXd &mat) {
    if(!rank()) return *this;

    Tensor<T> ttens(Nrank, Ndim);
    Index nd;

    T tval;
    Counter<Array<Index> > icount(element_counter());
    Counter<Array<Index> > jcount(icount);
    do {
      ttens(icount()) = 0;
      jcount.reset();
      do {
        tval = at(jcount());
        for(nd = 0; nd < Nrank; nd++)
          tval *= mat(icount[nd], jcount[nd]);
        ttens(icount()) += tval;

      }
      while((++jcount).valid());

    }
    while((++icount).valid());
    swap(ttens);
    return *this;
  }

  //************************************************************

  template<class T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Tensor<T>::convert_to_Eigen() const {
    if(rank() != 2) {
      std::cerr << "WARNING: Attempting to assign Eigen::Matrix object using a Tensor of rank " << rank()
                << "!\nTensor must have rank=2 for assignment to be valid. Exiting...\n";
      exit(1);
    }
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmat(dim(0), dim(1));
    tmat.setZero();
    for(Index i = 0; i < tmat.rows(); i++) {
      for(Index j = 0; j < tmat.rows(); j++) {
        tmat(i, j) = (*this)(i, j);
      }
    }
    return tmat;
  }

  //************************************************************

  template<class T>
  Tensor<T> operator*(const T &LHS, const Tensor<T> &RHS) {
    return Tensor<T>(RHS) *= LHS;
  }

  //************************************************************

  template<class T>
  Tensor<T> operator*(const SymOp &LHS, const Tensor<T> &RHS) {
    return Tensor<T>(RHS).apply_sym(LHS);
  }

  //************************************************************

  template<class T>
  T dot(const Tensor<T> &LHS, const Tensor<T> &RHS) {
    return LHS.scalar_prod(RHS);
  }

  //************************************************************

  template<class T>
  T norm(const Tensor<T> &ttens) {
    return sqrt(ttens.scalar_prod(ttens));
  }


  //************************************************************

  template<class T>
  std::ostream &operator<<(std::ostream &stream, const Tensor<T> &RHS) {
    std::ios::fmtflags old_flags = stream.flags();
    stream.flags(std::ios::showpoint | std::ios::fixed | std::ios::right);
    int twide = stream.width();

    if(!RHS.rank()) {
      stream << RHS[0] << '\n';
      stream.flags(old_flags);
      return stream;
    }

    Array<Index> min(RHS.rank(), 0), max(RHS.ind_max()), inc(RHS.rank(), 1);
    Counter<Array<Index> > icount(min, max, inc);
    do {
      if(RHS.rank() > 2 && icount[0] == 0 && icount[1] == 0) {
        stream << "Block (:,:," << icount[2];
        for(Index j = 3; j < RHS.rank(); j++)
          stream  << "," << icount[j];
        stream << ") :\n";
      }

      stream << std::setw(0) << "  " << std::setw(twide) << RHS(icount());
      if(icount[0] == RHS.ind_max()[0]) {
        stream << '\n';

        // if(RHS.rank() > 1 && icount[1] == RHS.ind_max()[1])
        //   stream << "\n"; //commented out by Ivy 07/01/13 so that an extra new line isn't printed.
      }

    }
    while((++icount).valid());
    stream.flags(old_flags);
    return stream;
  }



  //************************************************************

  template<class T>
  std::istream &operator >> (std::istream &stream, Tensor<T> &RHS) {
    std::cerr << "WARNING: Input stream operator has not been implemented for Tensor class!  Exiting...\n";
    exit(1);
  }



  //************************************************************

  template<class T>
  void TensorBasis<T>::generate_basis(Index Nrank, const SymGroup &sym_group) {
    Index i, ns;
    clear();

    Tensor<T> ttens(Nrank, Array<Index>(Nrank, 3), 0);

    for(i = 0; i < ttens.size(); i++) {
      Tensor<T> avtens(Nrank, Array<Index>(Nrank, 3), 0);
      ttens[i] = 1;
      for(ns = 0; ns < sym_group.size(); ns++)
        avtens += sym_group[ns] * ttens;
      if(!this->contains(avtens))
        this->push_back(avtens);
      ttens[i] = 0;
    }
    idealize();

    coeffs.resize(size(), NAN); //Changed by Ivy 09/25/2012

  }

  //************************************************************

  template<class T>
  void TensorBasis<T>::generate_basis(Index Nrank, const SymGroup &sym_group, Index Rep_ID) {
    if(!sym_group.size()) {
      std::cerr << "WARNING: Attempting to generate tensor basis from an invalid symmetry group. Exiting...";
      exit(0);
    }
    Index i, ns;
    clear();

    Tensor<T> ttens(Nrank, Array<Index>(Nrank, 3), 0);

    for(i = 0; i < ttens.size(); i++) {
      Tensor<T> avtens(Nrank, Array<Index>(Nrank, 3), 0);
      ttens[i] = 1;
      for(ns = 0; ns < sym_group.size(); ns++)
        avtens += sym_group[ns] * ttens;
      if(!contains(avtens))
        push_back(avtens);
      ttens[i] = 0;
    }
    idealize();

    coeffs.resize(size(), NAN); //Changed by Ivy 09/25/2012

  }


  //************************************************************

  template<class T>
  void TensorBasis<T>::generate_basis(Index Nrank, const SymGroup &sym_group, const SymGroupRep &perm_group) {
    Index i, ns;
    clear();

    if(sym_group.size() != perm_group.size()) {
      std::cerr << "WARNING: Attempting to generate tensor basis, but symmetry group and permutation group are not compatible!\n";
    }

    Tensor<T> ttens(Nrank, Array<Index>(Nrank, 3), 0);

    for(i = 0; i < ttens.size(); i++) {
      Tensor<T> avtens(Nrank, Array<Index>(Nrank, 3), 0);
      ttens[i] = 1;
      for(ns = 0; ns < sym_group.size(); ns++)
        avtens += (sym_group[ns] * ttens).dim_permute(*(perm_group[ns]->get_permutation()));
      if(!this->contains(avtens))
        this->push_back(avtens);
      ttens[i] = 0;
    }
    idealize();

    coeffs.resize(size(), NAN); //Changed by Ivy 09/25/2012
    return;
  }

  //************************************************************

  template<class T>
  void TensorBasis<T>::make_orthogonal_to(const TensorBasis<T> &ortho_basis) {
    Index i, j;
    for(i = 0; i < size(); i++) {
      for(j = 0; j < ortho_basis.size(); j++) {
        at(i) -= (at(i).scalar_prod(ortho_basis[j]) / norm(ortho_basis[j])) * ortho_basis[j];
      }
    }
    return;
  }

  //************************************************************

  template<class T>
  TensorBasis<T> &TensorBasis<T>::apply_sym(const SymOp &op) {
    for(Index i = 0; i < size(); i++) {
      at(i).apply_sym(op);
    }
    return *this;
  }

  //************************************************************

  template<class T>
  void TensorBasis<T>::normalize() {
    for(Index i = 0; i < size(); i++)
      at(i).normalize();
    return;
  }

  //************************************************************

  template<class T>
  void TensorBasis<T>::idealize() {
    //GaussJordan<Tensor<T>, T>::eliminate(*this);

    //    for(Index i=0; i<size(); i++)
    //std::cout << "Basis Tensor " << i << ":  \n" << at(i);

    GramSchmidt<Tensor<T>, T>::orthogonalize(*this);
    return;
  }

  //************************************************************
  template<class T>
  bool TensorBasis<T>::read(std::istream &stream) {

    char ch;
    Index num_elem, rank;
    Array<Index> dim;
    double tcoeff;

    ch = stream.peek();
    while(ch < '0' || ch > '9') {
      stream.ignore(256, '\n');
      ch = stream.peek();
    }
    stream >> num_elem;
    if(!num_elem) {
      stream.ignore(256, '\n');
      return false;
    }
    stream >> rank;

    dim.resize(rank);

    for(Index i = 0; i < rank; i++) {
      stream >> dim[i];
    }

#ifdef DEBUG
    std::cout << "dim is " << dim << "\n";
#endif //DEBUG

    stream.ignore(1000, '\n');

    coeffs.resize(num_elem, NAN);

    Tensor<T> ttensor(rank, dim);

    for(Index j = 0; j < num_elem; j++) {
      ttensor.redefine(dim);

      //Check for "<" to find "<ECI>"
      ch = stream.peek();
      while(ch != char(60)) {
        if(((ch < '0') || (ch > '9')) && (ch != '-')) {
          stream.ignore(1000, '\n');
        }
        else {
          break;
        }
        ch = stream.peek();
      }

      //Reading in coefficient
      if(((ch != char(60)) && ((ch >= '0') || (ch <= '9'))) || (ch == '-')) {
        stream >> tcoeff;
#ifdef DEBUG
        std::cout << "Coefficient is " << tcoeff << "\n";
#endif //DEBUG
        eci(j) = tcoeff;
      }

      stream.ignore(1000, '\n');

      ttensor.read(stream);
      this->push_back(ttensor);
    }


    return true;
  }

  /** @} */
};

#endif
