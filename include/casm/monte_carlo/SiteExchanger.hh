#ifndef CASM_SiteExchanger_HH
#define CASM_SiteExchanger_HH

#include <vector>
#include "casm/CASM_global_definitions.hh"

namespace CASM {

  class Supercell;

  class SiteExchanger {

  public:

    /// \brief Constructor determine possible swaps in the given Supercell
    SiteExchanger(const Supercell &scel);

    /// \brief std::vector of indices into occupation array of ConfigDoF that have more than one allowed occupant
    const std::vector<Index> &variable_sites() const {
      return m_variable_sites;
    }

    /// \brief sublat()[i] is the sublattice index for site variable_sites()[i]
    const std::vector<int> &sublat() const {
      return m_sublat;
    }

    /// \brief For a given sublattice with a particular occupant, show what OTHER occupants it could be
    ///
    /// possible_swap()[sublat][curr_occupant index][i] -> other occupant index (i, for each possible other occupant)
    const std::vector< std::vector< std::vector<int> > > &possible_swap() const {
      return m_possible_swap;
    }

    /// \brief Map the integer values from the possible swaps or variable sites arrays into actual species
    ///
    /// sublat_to_mol()[sublat][occupant_index] -> species_type_index
    const std::vector< std::vector<int> > &sublat_to_mol() const {
      return m_sublat_to_mol;
    }

  private:

    /// \brief std::vector of indices into occupation array of m_confdof that have more than one allowed occupant
    std::vector<Index> m_variable_sites;

    /// \brief m_sublat[i] is the sublattice index for site m_variable_sites[i]
    std::vector<int> m_sublat;

    /// \brief For a given sublattice with a particular occupant, show what OTHER occupants it could be
    /// m_possible_swap[sublat][curr_occupant index][i] -> other occupant index (i, for each possible other occupant)
    std::vector< std::vector< std::vector<int> > > m_possible_swap;

    /// \brief Map the integer values from the possible swaps or variable sites arrays into actual species
    /// m_sublat_to_mol[sublat][occupant_index] -> species_type_index
    std::vector< std::vector<int> > m_sublat_to_mol;

  };


}

#endif


