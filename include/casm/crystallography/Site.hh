#ifndef SITE_HH
#define SITE_HH

#include <iostream>

#include "casm/crystallography/Molecule.hh"

namespace CASM {

  class SymOp;

  /** \ingroup Coordinate
   *  @{
   */

  class Site : public Coordinate {
  public:
    explicit Site(const Lattice &init_home);
    Site(const Coordinate &init_pos, const std::string &occ_name);

    /// \brief Construct site with initial position and the allowed Molecule
    Site(const Coordinate &init_pos, std::initializer_list<Molecule> site_occ);

    const MoleculeOccupant &site_occupant() const {
      return m_site_occupant;
    }

    const Array<ContinuousDoF> &displacement() const {
      return m_displacement;
    }

    void update_data_members(const Site &_ref_site);

    /// Checks if current occupant is a vacancy
    bool is_vacant() const;

    ///access m_nlist_ind
    Index nlist_ind() const;

    /// Name of current occupant (name of molecule, but for single atom, molecule name is species name)

    std::string occ_name()const;

    /// Const reference to occupying molecule. ***WARNING*** only use if you are certain the occupant has been set.
    /// If you only need to know occupant name or whether site is vacant, use Site::is_vacant() or Site::occ_name() instead
    const Molecule &occ() const;

    bool compare(const Coordinate &test_coord, double compare_tol = TOL) const;
    bool compare(const Site &test_site, double compare_tol = TOL) const; //Ivy
    bool compare(const Site &test_site, const Coordinate &shift, double compare_tol = TOL) const;
    bool compare_type(const Site &test_site) const; //Ivy
    bool operator==(const Site &test_site) const;
    bool almost_equal(const Site &test_site, double tol) const;

    //checks to see if species with name 'name' is allowed at site.
    bool contains(const std::string &name) const;
    bool contains(const std::string &name, int &index) const;

    void set_lattice(const Lattice &new_lat, COORD_TYPE mode);//John G

    void set_site_occupant(const MoleculeOccupant &new_dof) {
      m_site_occupant = new_dof;
      m_type_ID = -1;
    }

    void set_occ_value(int new_val) {
      m_site_occupant.set_value(new_val);
    }

    void set_occ(const Molecule &new_occ) {
      m_site_occupant.set_current_state(new_occ);
    }


    Array<std::string> allowed_occupants() const;

    /// set basis_ind of site and its occupant functions
    void set_basis_ind(Index);

    /// set m_nlist_ind of Site and its DoFs
    void set_nlist_ind(Index);

    Site &apply_sym(const SymOp &op);
    Site &apply_sym_no_trans(const SymOp &op);

    void read(std::istream &stream, bool SD_is_on = false);
    void read(std::istream &stream, std::string &elem, bool SD_is_on);

    void print(std::ostream &stream) const;
    void print_occ(std::ostream &stream) const;
    void print_mol(std::ostream &stream, int spaces, char delim, bool SD_is_on = false)const;


    Site &operator+=(const Coordinate &translation);
    Site &operator-=(const Coordinate &translation);

    jsonParser &to_json(jsonParser &json) const;
    void from_json(const jsonParser &json);

  private:
    static Array<Site> &_type_prototypes() {
      static Array<Site> m_type_prototypes;
      return m_type_prototypes;
    }

    //Index into PrimClex neighbor list
    Index m_nlist_ind;     //John G 230913
    mutable Index m_type_ID;

    // displacement degrees of freedom of the molecule.
    // These may be x,y,z, or a subspace (e.g., displacement only in the x--y plane).
    Array<ContinuousDoF> m_displacement;

    // Configuration state is fundamentally different from most other degrees of freedom,
    // so we'll treat it separately. 'occupant' is the discrete degree of freedom associated
    // with the molecule that occupies the site
    MoleculeOccupant m_site_occupant;

    //============

    bool _compare_type_no_ID(const Site &test_site) const;
    Index _type_ID() const;


  };

  jsonParser &to_json(const Site &value, jsonParser &json);
  void from_json(Site &value, const jsonParser &json);

  std::ostream &operator<< (std::ostream &stream, const Site &site);

  Site operator*(const SymOp &LHS, const Site &RHS);
  Site operator+(const Site &LHS, const Coordinate &RHS);
  Site operator+(const Coordinate &LHS, const Site &RHS);

  /** @} */
}

#endif
