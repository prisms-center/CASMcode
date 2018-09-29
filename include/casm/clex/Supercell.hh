#ifndef CASM_Supercell
#define CASM_Supercell

#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/Comparisons.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/clex/SupercellTraits.hh"
#include "casm/clex/HasPrimClex.hh"
#include "casm/clex/HasCanonicalForm.hh"
#include "casm/database/Named.hh"

namespace CASM {

  template<typename T, typename U> class ConfigIterator;
  class PermuteIterator;
  class PrimClex;
  class Clexulator;
  class Configuration;
  template<typename T> class BasicStructure;
  class Site;
  class SuperNeighborList;
  class Structure;

  namespace DB {
    template<typename T> class DatabaseIterator;
  }

  struct ConfigMapCompare {
    bool operator()(const Configuration *A, const Configuration *B) const;
  };

  /** \defgroup Supercell
   *  \ingroup Clex
   *  \brief Represents a supercell of the primitive parent crystal structure
   *  @{
   */

  /// \brief Represents a supercell of the primitive parent crystal structure
  ///
  class Supercell : public
    HasPrimClex <
    DB::Named <
    Comparisons <
    SupercellCanonicalForm <
    CRTPBase<Supercell >>> >> {

  public:

    using permute_const_iterator = SupercellSymInfo::permute_const_iterator;

    // **** Constructors ****

    Supercell(const Supercell &RHS);
    Supercell(const PrimClex *_prim, const Lattice &superlattice);
    Supercell(const PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &superlattice_matrix);

    ~Supercell();

    // **** Coordinates ****

    /// \brief Return the sublattice index for a linear index
    Index sublat(Index linear_index) const;

    /// \brief Given a Coordinate and tolerance, return linear index into Configuration
    Index linear_index(const Coordinate &coord, double tol = TOL) const;

    /// \brief Return the linear index corresponding to integral coordinates
    Index linear_index(const UnitCellCoord &bijk) const;

    /// \brief Return the linear index corresponding to integral coordinates
    Coordinate coord(Index linear_index) const;

    /// \brief Return the integral coordinates corresponding to a linear index
    UnitCellCoord uccoord(Index linear_index) const;


    // returns maximum allowed occupation bitstring -- used for initializing enumeration counters
    std::vector<int> max_allowed_occupation() const;

    Configuration configuration(const BasicStructure<Site> &structure_to_config, double tol = TOL) const;

    /// return Structure corresponding to this supercell
    ///    w/ basis site occupation as per primclex.prim
    Structure superstructure() const;

    /// return Structure corresponding to this supercell
    ///    w/ basis site occupation as per config
    Structure superstructure(const Configuration &config) const;

    // **** Accessors ****

    const PrimClex &primclex() const;

    const PrimGrid &prim_grid() const;

    ///Return number of primitive cells that fit inside of *this
    Index volume() const;

    Index basis_size() const;

    Index num_sites() const;

    SymGroupRep const &permutation_symrep() const;

    Eigen::Matrix3i transf_mat() const;

    /// \brief The super lattice
    const Lattice &lattice() const;

    /// \brief Returns the SuperNeighborList
    const SuperNeighborList &nlist() const;

    // Factor group of this supercell
    const SymGroup &factor_group() const;

    // SymInfo object of this supercell
    const SupercellSymInfo &sym_info() const;

    // Returns the permutation representation of the i'th element of m_factor_group
    //const Permutation &factor_group_permute(Index i) const;

    // Returns the i'th element of m_trans_permute
    //const Permutation &translation_permute(Index i) const;

    // Const access of m_trans_permute
    //const std::vector<Permutation> &translation_permute() const;

    /// \brief Begin iterator over pure translational permutations
    //permute_const_iterator translate_begin() const;

    /// \brief End iterator over pure translational permutations
    //permute_const_iterator translate_end() const;

    // begin and end iterators for iterating over translation and factor group permutations
    //permute_const_iterator permute_begin() const;
    //permute_const_iterator permute_end() const;
    //permute_const_iterator permute_it(Index fg_index, Index trans_index) const;
    //permute_const_iterator permute_it(Index fg_index, UnitCell trans) const;


    bool operator<(const Supercell &B) const;

    /// \brief Insert the canonical form of this into the database
    std::pair<DB::DatabaseIterator<Supercell>, bool> insert() const;

    // **** Other ****

    bool is_supercell_of(const Structure &structure) const;
    bool is_supercell_of(const Structure &structure, Eigen::Matrix3d &multimat) const;
    std::vector<int> vacant()const;

    void write_pos() const;
    std::ostream &write_pos(std::ostream &sout) const;

  private:

    friend Comparisons<SupercellCanonicalForm<CRTPBase<Supercell>>>;
    friend DB::Named<Comparisons<SupercellCanonicalForm<CRTPBase<Supercell>>>>;

    bool eq_impl(const Supercell &B) const;

    // **** Generating functions ****

    // Populate m_factor_group -- probably should be private
    void _generate_factor_group() const;

    // Populate m_trans_permute -- probably should be private
    void _generate_permutations() const;

    std::string generate_name_impl() const;


    const PrimClex *m_primclex;

    // lattice of supercell in real space
    Lattice m_lattice;

    SupercellSymInfo m_sym_info;

    // superlattice arithmetic
    //PrimGrid m_prim_grid;

    // m_perm_symrep_ID is the ID of the SymGroupRep of prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       the the Supercel::factor_group_permute(int) method, so that you don't encounter the gaps
    //       OR, see note for Supercell::permutation_symrep() below.
    //mutable SymGroupRepID m_perm_symrep_ID;

    // m_factor_group is factor group of the super cell, found by identifying the subgroup of
    // (*this).prim().factor_group() that leaves the supercell lattice vectors unchanged
    // if (*this).prim() is actually primitive, then m_factor_group.size() <= 48
    // NOTE: This is different from the SymGroup found by doing (*this).superstruc().factor_group()
    //       if Tprim is the translation group formed by the primitive cell lattice vectors, then
    //       m_factor_group is the group formed by the cosets of Tprim in the supercell space group
    //       if Tsuper is the translation group formed by the supercell lattice vectors, then,
    //       m_occupation(init_config.occupation()),
    //       m_displacement(init_config.displacement()),
    //       m_strain(init_config.supercell().strain
    //       (*this).superstruc().factor_group() is the group formed by the cosets of Tsuper in the supercell space group
    //mutable SymGroup m_factor_group;

    /// SuperNeighborList, mutable for lazy construction
    mutable notstd::cloneable_ptr<SuperNeighborList> m_nlist;

    /// Store size of PrimNeighborList at time of construction of SuperNeighborList
    /// to enable checking if SuperNeighborList should be re-constructed
    mutable Index m_nlist_size_at_construction;

  };

  /// \brief Get canonical supercell from name. If not yet in database, construct and insert.
  const Supercell &make_supercell(const PrimClex &primclex, std::string name);

  SupercellSymInfo make_supercell_sym_info(Structure const &_prim, Lattice const &_slat);

  /// \brief Construct non-canonical supercell from name. Uses equivalent niggli lattice.
  std::shared_ptr<Supercell> make_shared_supercell(const PrimClex &primclex, std::string name);

  Supercell &apply(const SymOp &op, Supercell &scel);

  Supercell copy_apply(const SymOp &op, const Supercell &scel);

  Eigen::Matrix3i transf_mat(const PrimClex &primclex, const Lattice &super_lat);

  std::string generate_name(const Eigen::Matrix3i &transf_mat);

  std::string scelname(const PrimClex &primclex, const Lattice &lat);

  std::string canonical_scelname(const PrimClex &primclex, const Lattice &lat);

  /// \brief Construct the subgroup of permutations that leaves an element unchanged
  template<typename Element>
  std::vector<PermuteIterator> make_invariant_subgroup(const Element &element, const Supercell &scel);

  /** @} */
}
#endif
