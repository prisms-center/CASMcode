#ifndef SITE_HH
#define SITE_HH

#include "casm/misc/cloneable_ptr.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/basis_set/DoFSet.hh"

#include <iostream>
#include <string>

namespace CASM {
  class jsonParser;
  class SymOp;

  template<typename OccType>
  class OccupantDoF;
  class DoFSet;

  namespace xtal {


    class Molecule;


    /** \ingroup Coordinate
     *  @{
     */

    class Site : public Coordinate {
    public:
      explicit Site(const Lattice &init_home);

      Site(const Coordinate &init_pos, const std::string &occ_name);

      /// \brief Construct site with initial position and the allowed Molecule
      Site(const Coordinate &init_pos, std::initializer_list<Molecule> site_occ);

      ~Site();

      const OccupantDoF<Molecule> &occupant_dof() const;

      DoFSet const &dof(std::string const &_dof_type) const;

      Index dof_size() const;

      bool has_dof(std::string const &_dof_type) const;

      std::vector<std::string> dof_types() const;

      bool time_reversal_active() const;

      /// Checks if current occupant is a vacancy
      bool is_vacant() const;

      ///access m_label;
      Index label() const;

      /// Name of current occupant (name of molecule, but for single atom, molecule name is species name)

      std::string occ_name()const;

      /// Const reference to occupying molecule. ***WARNING*** only use if you are certain the occupant has been set.
      /// If you only need to know occupant name or whether site is vacant, use Site::is_vacant() or Site::occ_name() instead
      const Molecule &occ() const;

      bool compare(const Coordinate &test_coord) const;
      bool compare(const Site &test_site) const; //Ivy
      bool compare(const Site &test_site, const Coordinate &shift) const;
      bool compare_type(const Site &test_site) const; //Ivy
      bool operator==(const Site &test_site) const;
      bool almost_equal(const Site &test_site) const;

      //checks to see if species with name 'name' is allowed at site.
      bool contains(const std::string &name) const;
      bool contains(const std::string &name, int &index) const;

      void set_allowed_occupants(std::vector<Molecule> const &_occ_domain);

      void set_occ_value(int new_val);

      void set_occ(const Molecule &new_occ);

      void set_dofs(std::map<std::string, DoFSet> _dofs);

      std::map<std::string, DoFSet> const &dofs() const {
        return m_dof_map;
      }

      std::vector<std::string> allowed_occupants() const;

      /// set basis_ind of site and its occupant functions
      void set_basis_ind(Index _basis_ind);

      /// set m_label of Site
      void set_label(Index _new_label);

      Site &apply_sym(const SymOp &op);
      Site &apply_sym_no_trans(const SymOp &op);

      void read(std::istream &stream, bool SD_is_on = false);
      void read(std::istream &stream, std::string &elem, bool SD_is_on);

      void print(std::ostream &stream, Eigen::IOFormat format = Eigen::IOFormat(7, 12)) const;
      void print_occ(std::ostream &stream, Eigen::IOFormat format = Eigen::IOFormat(7, 12)) const;
      void print_mol(std::ostream &stream, int spaces, char delim, bool SD_is_on = false)const;


      Site &operator+=(const Coordinate &translation);
      Site &operator-=(const Coordinate &translation);

    private:
      static std::vector<Site> &_type_prototypes() {
        static std::vector<Site> m_type_prototypes;
        return m_type_prototypes;
      }

      Site &_apply_sym_attributes(const SymOp &op);

      /// Integer label used to differentiate sites of otherwise identical type
      Index m_label;

      mutable Index m_type_ID;

      // Configuration state is fundamentally different from most other degrees of freedom,
      // so we'll treat it separately. 'occupant' is the discrete degree of freedom associated
      // with the molecule that occupies the site
      notstd::cloneable_ptr<OccupantDoF<Molecule>> m_occupant_dof;

      /// additional continuous degrees of freedom
      std::map <std::string, DoFSet > m_dof_map;

      //============

      bool _compare_type_no_ID(const Site &test_site) const;
      Index _type_ID() const;


    };

    std::ostream &operator<< (std::ostream &stream, const Site &site);

    Site operator*(const SymOp &LHS, const Site &RHS);
    Site operator+(const Site &LHS, const Coordinate &RHS);
    Site operator+(const Coordinate &LHS, const Site &RHS);

    /** @} */
  }
}

#endif
