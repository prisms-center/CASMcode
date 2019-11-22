#include "casm/crystallography/Molecule.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
  namespace xtal {

    //****************************************************
    //
    //****************************************************
    /*
    void AtomPosition::print(std::ostream &stream,
                             Eigen::Ref<const Eigen::Vector3d> const &translation,
                             Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
                             int spaces,
                             bool print_sd_flags) const {
      for(int i = 0; i < spaces; i++) {
        stream << ' ';
      }

      stream << (cart2frac * (cart() + translation)).transpose();

      if(print_sd_flags) {
        //for(int i = 0; i < 3; i++) {
        //if(m_sd_flag[i]) stream << "  T";
        //else stream << "  F";
        //}
      }

      stream << "   " << name();

      return;
      }*/

    //****************************************************
    //
    //****************************************************

    AtomPosition &AtomPosition::apply_sym(const SymOp &op) {
      m_position = get_matrix(op) * m_position;
      for(auto it = m_attribute_map.begin(); it != m_attribute_map.end(); ++it)
        (it->second).apply_sym(op);

      return *this;
    }

    //****************************************************
    //
    //****************************************************

    bool AtomPosition::identical(AtomPosition const &RHS, double _tol) const {

      return almost_equal(cart(), RHS.cart(), _tol) && compare_type(*this, RHS, _tol);
    }

    //****************************************************

    bool compare_type(AtomPosition const &A, AtomPosition const &B, double tol) {
      // compare number of attributes
      if(A.attributes().size() != B.attributes().size())
        return false;
      if(A.name() != B.name())
        return false;
      // compare attributes
      auto it_A(A.attributes().cbegin()), end_it(A.attributes().cend());
      for(; it_A != end_it; ++it_A) {
        auto it_B = B.attributes().find(it_A->first);
        if(it_B == B.attributes().cend() || !(it_A->second).identical(it_B->second, tol))
          return false;
      }
      return true;
    }

    //****************************************************

    bool Molecule::is_atomic() const {
      if(size() != 1)
        return false;
      if(!attributes().empty())
        return false;
      for(AtomPosition const &atom : atoms()) {
        if(atom.cart().norm() > TOL)
          return false;
        if(!atom.attributes().empty())
          return false;
      }
      return true;
    }
    //****************************************************

    bool Molecule::is_vacancy() const {
      //return m_atoms.empty();
      return ::CASM::xtal::is_vacancy(m_atoms[0].name());
    }

    //****************************************************
    // Applies symmetry operation to a Molecule
    //****************************************************

    Molecule &Molecule::apply_sym(SymOp const &op) {
      for(Index i = 0; i < size(); i++)
        m_atoms[i].apply_sym(op);
      for(auto it = m_attribute_map.begin(); it != m_attribute_map.end(); ++it)
        (it->second).apply_sym(op);
      return *this;
    }

    //****************************************************
    //
    //****************************************************

    bool Molecule::identical(Molecule const &RHS, double _tol) const {
      // compare number of attributes
      if(m_attribute_map.size() != RHS.m_attribute_map.size())
        return false;

      // compare number of atoms
      if(size() != RHS.size())
        return false;

      // compare atoms, irrespective of order
      for(Index i = 0; i < RHS.size(); i++) {
        Index j = 0;
        for(j = 0; j < size(); j++) {
          if(atom(i).identical(RHS.atom(j), _tol))
            break;
        }
        if(j == size())
          return false;
      }

      // compare attributes
      auto it(m_attribute_map.cbegin()), end_it(m_attribute_map.cend());
      for(; it != end_it; ++it) {
        auto it_RHS = RHS.m_attribute_map.find(it->first);
        if(it_RHS == RHS.m_attribute_map.cend() || !(it->second).identical(it_RHS->second, _tol))
          return false;
      }

      return true;
    }

    //****************************************************
    //
    //****************************************************

    bool Molecule::contains(std::string const &_name) const {
      for(Index i = 0; i < size(); i++)
        if(atom(i).name() == _name)
          return true;
      return false;
    }

    //****************************************************
    //
    //****************************************************
    /*
    void Molecule::print(std::ostream &stream,
                         Eigen::Ref<const Eigen::Vector3d> const &translation,
                         Eigen::Ref<const Eigen::Matrix3d> const &cart2frac,
                         int spaces,
                         char delim,
                         bool print_sd_flags ) const {
      for(Index i = 0; i < size(); i++) {
        atom(i).print(stream, translation, cart2frac, spaces, print_sd_flags);
        stream << delim;
      }
      return;
    }*/

    //****************************************************
    /// \brief Return an atomic Molecule with specified name and Lattice
    //****************************************************

    Molecule Molecule::make_atom(std::string const &atom_name) {
      return Molecule(atom_name, {AtomPosition(0., 0., 0., atom_name)});
    }

    //****************************************************
    //
    //****************************************************

    Molecule Molecule::make_vacancy() {
      //return Molecule("Va", {});
      return make_atom("Va");
    }

  }
}
