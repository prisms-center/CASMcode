#ifndef SYMPERMUTATION_HH
#define SYMPERMUTATION_HH

#include <iostream>
#include <cmath>
#include "casm/symmetry/SymOpRepresentation.hh"
#include "casm/container/Permutation.hh"

namespace CASM {
  class MasterSymGroup;
  class Lattice;

  /** \ingroup Symmetry
   *  @{
   */

  ///\brief  SymPermutation describes how a symmetry operation permutes a list of 'things'
  /// For example, Coordinates in a Cluster, basis atoms in a Structure, Clusters in an Orbit, etc.
  class SymPermutation: public SymOpRepresentation {
  public:
    /// Initialize a SymPermutation with the permutation array.
    /// The corresponding matrix is generated automatically
    SymPermutation(const Array<Index> &init_permute) : m_permute(init_permute), m_has_mat(false) {

    }

    /// Initialize a SymPermutation with the permutation array.
    /// The corresponding matrix is generated automatically
    SymPermutation(const Permutation &init_permute) : m_permute(init_permute), m_has_mat(false) {

    }

    double character() const override;

    /// Return pointer to a copy of this SymPermutation
    SymOpRepresentation *copy() const override {
      return new SymPermutation(*this);
    }

    /// Access the permutation array 'm_permute'
    Permutation const *get_permutation() const override {
      return &m_permute;
    }

    /// Access the permutation matrix
    Eigen::MatrixXd const *get_MatrixXd() const override {
      if(!m_has_mat)
        _calc_mat();
      return &m_mat;
    }

    jsonParser &to_json(jsonParser &json) const override;

    void from_json(const jsonParser &json) override;

  private:
    /// Array of indices, of length 'n'. An index 'm_permute[j]' before application of symmetry
    /// resides at index 'j' after application of symmetry
    /// example: For an 'Array<THINGS> my_array', transforms as my_array.permute(m_permute);
    Permutation m_permute;

    // true if m_mat is initialized
    bool m_has_mat;

    /// Matrix of ones and zeroes that reorders elements of a vector
    /// Matrix is nxn, where 'n' is the number of things that are permuted
    mutable Eigen::MatrixXd m_mat;

    /// Generate the matrix of permutation, when m_permute is known
    void _calc_mat() const;

  };

  jsonParser &to_json(const SymPermutation &sym, jsonParser &json);
  void from_json(SymPermutation &sym, const jsonParser &json);

  /** @{ */
}
#endif
