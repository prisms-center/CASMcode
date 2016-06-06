#ifndef SUPERCELL_HH
#define SUPERCELL_HH

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/PrimGrid.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {


  template<typename T, typename U> class ConfigIterator;
  class PermuteIterator;
  class PrimClex;
  class Clexulator;

  /// \brief Represents a supercell of a PrimClex
  ///
  /// \ingroup Clex
  ///
  class Supercell {

  public:

    typedef boost::container::stable_vector<Configuration> ConfigList;

    typedef ConfigIterator<Configuration, PrimClex> config_iterator;
    typedef ConfigIterator<const Configuration, const PrimClex> config_const_iterator;
    typedef PermuteIterator permute_const_iterator;

  private:
    // pointer to Primcell containing all the cluster expansion data
    PrimClex *m_primclex;

    // lattice of supercell in real space
    Lattice m_real_super_lattice;

    // superlattice arithmetic
    PrimGrid m_prim_grid;

    // m_perm_symrep_ID is the ID of the SymGroupRep of prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       the the Supercel::factor_group_permute(int) method, so that you don't encounter the gaps
    //       OR, see note for Supercell::permutation_symrep() below.
    mutable SymGroupRepID m_perm_symrep_ID;

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
    mutable SymGroup m_factor_group;

    /// unique name of the supercell based on hermite normal form (see generate_name() )
    std::string m_name;

    /// SuperNeighborList, mutable for lazy construction
    mutable notstd::cloneable_ptr<SuperNeighborList> m_nlist;

    /// Store size of PrimNeighborList at time of construction of SuperNeighborList
    /// to enable checking if SuperNeighborList should be re-constructed
    mutable Index m_nlist_size_at_construction;


    // Could hold either enumerated configurations or any 'saved' configurations
    ConfigList m_config_list;

    Eigen::Matrix3i m_transf_mat;

    /// index into PrimClex::supercell_list, used for configuration iterators
    Index m_id;

  public:

    // **** Constructors ****

    Supercell(const Supercell &RHS);
    Supercell(PrimClex *_prim, const Lattice &superlattice);
    Supercell(PrimClex *_prim, const Eigen::Ref<const Eigen::Matrix3i> &superlattice_matrix);


    // **** Coordinates ****

    /// \brief Return the sublattice index for a linear index
    ///
    /// Linear indices are grouped by sublattice, then ordered as determined by
    /// PrimGrid. This function is equivalent to:
    /// \code
    /// linear_index / volume();
    /// \endcode
    Index sublat(Index linear_index) const {
      return linear_index / volume();
    }

    /// \brief Given a Coordinate and tolerance, return linear index into Configuration
    ///
    ///   This may be slow, first converts Coordinate -> UnitCellCoord,
    ///   then gets linear_index from UnitCellCoord
    ///
    /// Implementation:
    /// \code
    /// Coordinate tcoord(coord);
    /// tcoord.within();
    /// return linear_index(UnitCellCoord(prim(), coord, tol));
    /// \endcode
    Index linear_index(const Coordinate &coord, double tol = TOL) const {
      Coordinate tcoord(coord);
      tcoord.within();
      return linear_index(UnitCellCoord(prim(), coord, tol));
    };

    /// \brief Return the linear index corresponding to integral coordinates
    ///
    /// Linear indices are grouped by sublattice, then ordered as determined by
    /// PrimGrid. This function is equivalent to:
    /// \code
    /// bijk[0] * volume() + m_prim_grid.find(bijk.unitcell());
    /// \endcode
    Index linear_index(const UnitCellCoord &bijk) const {
      return bijk[0] * volume() + m_prim_grid.find(bijk.unitcell());
    }

    /// \brief Return the linear index corresponding to integral coordinates
    ///
    /// Equivalent to:
    /// \code
    /// uccoord(linear_index).coordinate()
    /// \endcode
    Coordinate coord(Index linear_index) const {
      return uccoord(linear_index).coordinate();
    }

    /// \brief Return the integral coordinates corresponding to a linear index
    ///
    /// Linear indices are grouped by sublattice, then ordered as determined by
    /// PrimGrid. This function is equivalent to:
    /// \code
    /// UnitCellCoord(prim(), sublat(linear_index), m_prim_grid.unitcell(linear_index % volume()))
    /// \endcode
    UnitCellCoord uccoord(Index linear_index) const {
      return UnitCellCoord(prim(), sublat(linear_index), m_prim_grid.unitcell(linear_index % volume()));
    };


    // returns maximum allowed occupation bitstring -- used for initializing enumeration counters
    ReturnArray<int> max_allowed_occupation() const;

    Configuration configuration(const BasicStructure<Site> &structure_to_config, double tol = TOL);

    // return Structure corresponding to this supercell
    //    w/ basis site occupation as per primclex.prim
    Structure superstructure() const;
    //    w/ basis site occupation as per config
    Structure superstructure(const Configuration &config) const;        //This should be private
    ///Returns a structure corresponding to the specified configuration.
    Structure superstructure(Index config_index) const;
    //    w/ basis site occupation as per config
    //    and itself as prim, not primclex->prim


    // **** Accessors ****

    const PrimClex &primclex() const {
      return *m_primclex;
    }

    const PrimGrid &prim_grid() const {
      return m_prim_grid;
    }

    const Structure &prim() const;

    ///Return number of primitive cells that fit inside of *this
    Index volume()const {
      return m_prim_grid.size();
    };

    Index basis_size() const {
      return prim().basis.size();
    }

    Index num_sites()const {
      return volume() * basis_size();
    };

    // the permutation_symrep is the SymGroupRep of prim().factor_group() that describes how
    // operations of m_factor_group permute sites of the Supercell.
    // NOTE: The permutation representation is for (*this).prim().factor_group(), which may contain
    //       more operations than m_factor_group, so the Permutation SymGroupRep may have 'gaps' at the
    //       operations that aren't in m_factor_group. You should access elements of the SymGroupRep using
    //       SymGroupRep::get_representation(m_factor_group[i]) or SymGroupRep::get_permutation(m_factor_group[i]),
    //       so that you don't encounter the gaps (i.e., the representation can be indexed using the
    //       SymOps of m_factor_group
    SymGroupRepID permutation_symrep_ID()const {
      if(m_perm_symrep_ID.empty())
        generate_permutations();
      return m_perm_symrep_ID;
    }

    SymGroupRep const &permutation_symrep() const {
      return prim().factor_group().representation(permutation_symrep_ID());
    }

    const Eigen::Matrix3i &transf_mat() const {
      return m_transf_mat;
    };

    const Lattice &real_super_lattice() const {
      return m_real_super_lattice;
    };

    /// \brief Returns the SuperNeighborList
    const SuperNeighborList &nlist() const;


    ConfigList &config_list() {
      return m_config_list;
    };

    const ConfigList &config_list() const {
      return m_config_list;
    };

    const Configuration &config(Index i) const {
      return m_config_list[i];
    };

    Configuration &config(Index i) {
      return m_config_list[i];
    }

    // begin and end iterators for iterating over configurations
    config_iterator config_begin();
    config_iterator config_end();

    // begin and end const_iterators for iterating over configurations
    config_const_iterator config_cbegin() const;
    config_const_iterator config_cend() const;

    Index id() const {
      return m_id;
    }

    void set_id(Index id) {
      m_id = id;
    }

    std::string name() const {
      return m_name;
    };

    // Populates m_factor_group (if necessary) and returns it.
    const SymGroup &factor_group() const;

    // Returns the permutation representation of the i'th element of m_factor_group
    const Permutation &factor_group_permute(Index i) const;

    // Returns the i'th element of m_trans_permute
    // Populates m_trans_permute if needed
    const Permutation &translation_permute(Index i) const;

    // Const access of m_trans_permute
    // Populates m_trans_permute if needed
    const Array<Permutation> &translation_permute() const;

    // begin and end iterators for iterating over translation and factor group permutations
    permute_const_iterator permute_begin() const;
    permute_const_iterator permute_end() const;

    ///Return path to supercell directory
    fs::path path() const;

    ///Count how many configs are selected in *this
    Index amount_selected() const;


    // **** Generating functions ****

    // Populate m_factor_group -- probably should be private
    void generate_factor_group() const;

    // Populate m_trans_permute -- probably should be private
    void generate_permutations() const;

    //\John G 070713
    void generate_name();


    // **** Enumerating functions ****

    bool contains_config(const Configuration &config) const;
    bool contains_config(const Configuration &config, Index &index) const;
    bool add_config(const Configuration &config);
    bool add_config(const Configuration &config, Index &index, Supercell::permute_const_iterator &permute_it);
    bool add_canon_config(const Configuration &config, Index &index);
    void read_config_list(const jsonParser &json);

    template<typename ConfigIterType>
    void add_unique_canon_configs(ConfigIterType it_begin, ConfigIterType it_end);

    template<typename ConfigIterType>
    void add_configs(ConfigIterType it_begin, ConfigIterType it_end);

    // **** Other ****
    // Reads a relaxed structure and calculates the strains and stretches using the reference structure
    void read_relaxed_structure(Index configNum, const Lattice &home_lattice);
    void read_relaxed_structure(Index configNum);
    void read_clex_relaxations(const Lattice &home_lattice);

    bool is_supercell_of(const Structure &structure) const;
    bool is_supercell_of(const Structure &structure, Eigen::Matrix3d &multimat) const;
    ReturnArray<int> vacant()const;

    // **** Printing ****

    void print_bijk(std::ostream &stream);
    //   void print_clex_correlations(std::ostream &corrFile);
    ///Old CASM style corr.in output for all the configurations in *this supercell
    //   void print_global_correlations_simple(std::ostream &corrstream) const;
    void print_sublat_to_comp(std::ostream &stream);

    ///Call Configuration::write out every configuration in supercell
    jsonParser &write_config_list(jsonParser &json);

  };


  //*******************************************************************************
  // Warning: Assumes configurations are in canonical form
  template<typename ConfigIterType>
  void Supercell::add_unique_canon_configs(ConfigIterType it_begin, ConfigIterType it_end) {
    // Remember existing configs, to avoid duplicates
    //   Enumerated configurations are added after existing configurations
    Index N_existing = m_config_list.size();
    Index N_existing_enumerated = 0;
    //std::cout << "ADDING CONFIGS TO SUPERCELL; N_exiting: " << N_existing << " N_enumerated: " << N_existing_enumerated << "\n";
    //std::cout << "beginning iterator: " << it_begin->occupation() << "\n";
    // Loops through all possible configurations
    for(; it_begin != it_end; ++it_begin) {
      //std::cout << "Attempting to add configuration: " << it_begin->occupation() << "\n";
      // Adds the configuration to the list, if not among previously existing configurations

      bool add = true;
      if(N_existing_enumerated != N_existing) {
        for(Index i = 0; i < N_existing; i++) {
          if(m_config_list[i].configdof() == it_begin->configdof()) {
            m_config_list[i].push_back_source(it_begin->source());
            add = false;
            N_existing_enumerated++;
            break;
          }
        }
      }
      if(add) {
        m_config_list.push_back(*it_begin);
        // get source info from enumerator
        //m_config_list.back().set_source(it_begin.source());
        m_config_list.back().set_id(m_config_list.size() - 1);
        m_config_list.back().set_selected(false);
      }
    }

  }

  //*******************************************************************************

  template<typename ConfigIterType>
  void Supercell::add_configs(ConfigIterType it_begin, ConfigIterType it_end) {
    if(ConfigIterType::is_canonical_iter()) {
      for(; it_begin != it_end; ++it_begin) {
        add_canon_config(*it_begin);
      }
    }
    else {
      for(; it_begin != it_end; ++it_begin) {
        add_config(*it_begin);
      }
    }
  }

}
#endif
