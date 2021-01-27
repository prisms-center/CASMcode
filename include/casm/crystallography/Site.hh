#ifndef SITE_HH
#define SITE_HH

#include <iostream>
#include <string>
#include <vector>

#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/DoFSet.hh"
#include "casm/crystallography/Molecule.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {
namespace xtal {

struct SymOp;
class Molecule;

/** \ingroup Coordinate
 *  @{
 */

class Site : public Coordinate {
 public:
  explicit Site(const Lattice &init_home);

  Site(const Coordinate &init_pos, const std::string &occ_name);

  /// \brief Construct site with initial position, allowed molecules (occupants)
  Site(const Coordinate &init_pos, const std::vector<Molecule> &site_occ);

  /// \brief Construct site with initial position, allowed molecules
  /// (occupants), and local degrees of freedom
  Site(const Coordinate &init_pos, const std::vector<Molecule> &site_occ,
       const std::map<std::string, SiteDoFSet> &site_dof);

  /// \brief Construct site with initial position, allowed molecules
  /// (occupants), and local degrees of freedom
  Site(const Coordinate &init_pos, const std::vector<Molecule> &site_occ,
       const std::vector<SiteDoFSet> &site_dof);

  ~Site();

  // TODO: Decide where this should actually go
  static void print_occupant_dof(const std::vector<Molecule> &allowed_occupants,
                                 std::ostream &out_stream);

  const std::vector<Molecule> &occupant_dof() const;

  SiteDoFSet const &dof(std::string const &_dof_type) const;

  Index dof_size() const;

  bool has_dof(std::string const &_dof_type) const;

  std::vector<std::string> dof_types() const;

  bool time_reversal_active() const;

  /// access m_label;
  Index label() const;

  // TODO: Make comparators independent classes?
  bool compare(const Coordinate &test_coord) const;
  bool compare(const Site &test_site) const;       // Ivy
  bool compare_type(const Site &test_site) const;  // Ivy
  bool operator==(const Site &test_site) const;
  bool almost_equal(const Site &test_site) const;

  // checks to see if species with name 'name' is allowed at site.
  bool contains(const std::string &name) const;
  bool contains(const std::string &name, int &index) const;

  void set_allowed_occupants(const std::vector<Molecule> &new_occ_domain);

  void set_dofs(std::map<std::string, SiteDoFSet> _dofs);

  std::map<std::string, SiteDoFSet> const &dofs() const { return m_dof_map; }

  // TODO: Change this to allowed_occupant_names or something
  std::vector<std::string> allowed_occupants() const;

  /// Set label of Site. The label is used to distinguish between otherwise
  /// identical sites.
  void set_label(Index _new_label);

  void read(std::istream &stream, bool SD_is_on = false);
  void read(std::istream &stream, std::string &elem, bool SD_is_on);

  void print(std::ostream &stream,
             Eigen::IOFormat format = Eigen::IOFormat(7, 12)) const;

  Site &operator+=(const Coordinate &translation);
  Site &operator-=(const Coordinate &translation);

 private:
  // TODO: What is this?
  static std::vector<Site> &_type_prototypes() {
    static std::vector<Site> m_type_prototypes;
    return m_type_prototypes;
  }

  /// Integer label used to differentiate sites of otherwise identical type
  Index m_label;

  // TODO: What is this?
  mutable Index m_type_ID;

  // Configuration state is fundamentally different from most other degrees of
  // freedom, so we'll treat it separately. 'occupant' is the discrete degree of
  // freedom associated with the molecule that occupies the site
  std::vector<Molecule> m_occupant_dof;

  /// additional continuous degrees of freedom
  std::map<std::string, SiteDoFSet> m_dof_map;

  //============

  bool _compare_type_no_ID(const Site &test_site) const;
  Index _type_ID() const;
};

std::ostream &operator<<(std::ostream &stream, const Site &site);

/// Apply the symmetry operation to the site, and return a new transformed value
Site operator*(const SymOp &LHS, const Site &RHS);
template <typename ExternSymOp>
Site operator*(const ExternSymOp &LHS, const Site &RHS) {
  return adapter::Adapter<SymOp, ExternSymOp>()(LHS) * RHS;
}

Site operator+(const Site &LHS, const Coordinate &RHS);
Site operator+(const Coordinate &LHS, const Site &RHS);

/** @} */
}  // namespace xtal

namespace sym {
/// \brief Apply SymOp to a Site. Transforms the Coordinate as well as the
/// allowed occupants (Molecules) and other degrees of freedom.
xtal::Site &apply(const xtal::SymOp &op, xtal::Site &mutating_site);
/// \brief Copy and apply SymOp to a Site. Transforms the Coordinate as well as
/// the allowed occupants (Molecules) and other degrees of freedom.
xtal::Site copy_apply(const xtal::SymOp &op, xtal::Site site);
}  // namespace sym
}  // namespace CASM

#endif
