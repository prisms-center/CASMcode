#include "casm/monte_carlo/SiteExchanger.hh"

#include "casm/clex/Supercell.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {


  // ---- SiteExchanger Definitions ---------------------------------


  /// \brief Constructor determine possible swaps in the given Supercell
  ///
  ///For Monte Carlo simulations such as Grand Canonical that involve changing the occupants of sites, you need to know:
  ///     -Which sites are allowed to change                      -> variable_sites   std::vector<long int>
  ///     -Which sublattices those sites are on                   -> sublat           std::vector<int>
  ///     -What sites that are allowed to change can change to    -> possible_swap    std::vector<std::vector<std::vector<int> > >
  ///     -How to interpret the occupant indices                  -> sublat_to_mol    std::vector<std::vector<int>
  ///
  /// The Monte Carlo class works on a Configuration, which has an std::vector<int> that describes what's living at each site (m_occupation).
  /// Remember this array is organized in blocks of primitive basis sites. This means that for a PRIM with 3 basis sites in a SUPERCELL
  /// of size four, the first four values in the occupancy array would correspond to the first basis site of the PRIM. Put another way,
  /// the array is organized in bijk order.
  ///
  /// Suppose site 0 can hold type A, D
  ///              1 can hold type A
  ///              2 can hold type A, B, C
  ///
  /// PRIM:   /------------\                              SUPERCELL:  /-------------------------\
  ///         |         2  |                                          |         5  |         11 |
  ///         |     1      |                                          |     4      |     10     |
  ///         | 0          |                                          | 3          | 9          |
  ///         \------------/                                          |------------+------------|
  ///                                                                 |         2  |         8  |
  ///                                                                 |     1      |     7      |
  ///                                                                 | 0          | 6          |
  ///                                                                 \-------------------------/
  ///
  /// Though the order the example of SUPERCELL shows each PRIM is may be off, it is unimportant for this routine. The occupancy array
  /// of the Configurations with parent SUPERCELL will represent sites in the order:
  ///
  /// OCCUPATION:              |   0   3   6   9   |   1   4   7   10  |   2   5   8   11  |       (!)
  ///
  /// (!)      The actual occupancy array holds the values of the current occupant. What's shown above is simply meant to show which
  /// sites correspond to each index.
  ///
  /// Where '|' represents the beginning/end of a primitive basis block (not found in the actual std::vector<int>).
  /// variable_sites is simply a subset of the occupancy array indices, containing only the blocks for sites that
  /// can hold more than one occupant. We don't care about the others because the indices in variable_sites are selected
  /// at random to change the occupation at a site. For this example, this means we only care about sites that came
  /// from primitive sites 0 (A, D) and 2 (A, B, C).
  ///
  /// VARIABLE_SITES:          [   0   1   2   3   8   9   10   11  ]                               (!!)
  ///
  /// (!!)     Unlike the OCCUPATION example, VARIABLE_SITES shown above would hold these values.
  ///
  ///
  /// We also store the sublattice that each variable_site corresponds to
  ///
  /// SUBLAT:          [   0   0   0   0   2   2   2   2  ]                               (!!)
  ///
  /// (!!)     Unlike the OCCUPATION example, SUBLAT shown above would hold these values.
  ///
  /// The result is an array from which we pick a random value, which we use as an index to access a site in the occupancy array
  /// of a Configuration. Whatever site that corresponds to will have more than one allowed occupant. We want to change the current
  /// occupant to something else, but need to know what values that site can go up to. In order to know this we use possible_swap, which
  /// for a given site (outer array) for a given current occupant (middle array) holds all values that are NOT the
  /// current occupant (inner array).
  ///
  /// POSSIBLE_SWAP:           [   1   ][   1   2   ]
  ///                          [   0   ][   0   2   ]
  ///                                   [   0   1   ]
  ///
  /// The example above shows m_possible_swap, with each column of [] representing a sublattice. For a particular
  /// sublattice, we can find all the values the site can change to by looking at the innermost array with the index of
  /// the current occupant. For example, we randomly pick index 6 of VARIABLE_SITES. This corresponds to site 8 in SUPERCELL. This site
  /// can be occupied by A, B, or C, which is represented simply as 0, 1, 2 in the occupancy array of the Configuration. By
  /// accessing said array at index 8, we might find the value to be 2. POSSIBLE_SWAP[6][2] shows the array [  0   1   ], which
  /// are values that SUPERCELL site 8 is allowed to be changed to (currently C, but can change to A, B). When doing Monte Carlo
  /// we can again randomly pick from this innermost array to change the current occupant.
  ///     The last thing needed is a way to know what the integers in the occupancy array mean in terms of actual species. For
  /// example, a value of 1 on SUPERCELL site 3 represents species D, but a value of 1 on SUPERCELL site 8 represents species B.
  /// To know what means what we use site_to_mol. For a particular site (outer array) with a particular current occupant (inner array)
  /// we get the index into the array of species (type std::string) returned by ParamComposition::get_components().
  ///
  /// Suppose the array of allowed components is returned as   [   A   B   C   D   ]
  ///
  /// SUBLAT_TO_MOL:      [   0   3   ]
  ///                     [   0   1   2   ]
  ///
  /// Using the same example as before, we can read that a value of 1 on sublattice 0 corresponds to D, while a value of 1 on
  /// sublattice 1 represents species B:
  ///     -VARIABLE_SITES[1] gives the index into the occupancy array of Configuration to SUPERCELL site 3
  ///     -SUBLAT[1] gives the sublattice index for SUPERCELL site 3, which is 0
  ///     -POSSIBLE_SWAP[SUBLAT[1]][1] (assuming the current occupant is 1) shows only 0 as an alternative allowed occupant
  ///     -SUBLAT_TO_MOL[SUBLAT[1]][1] gives value 3, which corresponds to D in the array of allowed components
  ///
  ///     -VARIABLE_SITES[6] gives the index into the occupancy array of Configuration to SUPERCELL site 8
  ///     -SUBLAT[6] gives the sublattice index for SUPERCELL site 8, which is 2
  ///     -POSSIBLE_SWAP[SUBLAT[6]][1] (assuming the current occupant is 1) shows 0 and 2 as an alternative allowed occupants
  ///     -SUBLAT_TO_MOL[SUBLAT[6]][1] gives value 1, which corresponds to B in the array of allowed components
  ///
  SiteExchanger::SiteExchanger(const Supercell &scel) {

    //The occupants in a Configuration are ordered in blocks of basis sites.
    int scel_volume = scel.volume();

    //int scel_basis = scel->basis_size();
    int prim_basis = scel.get_prim().basis.size();
    std::vector<std::string> allowed_components = scel.get_primclex().composition_axes().components();

    //Count over sites in prim basis.
    for(Index prim_basis_site = 0; prim_basis_site < prim_basis; prim_basis_site++) {

      //If the site we're working allows multiple occupants, we're interested in filling up values for it.
      const auto &allowed = scel.get_prim().basis[prim_basis_site].allowed_occupants();
      std::vector<std::string> site_allowed_occ(allowed.cbegin(), allowed.cend());

      if(site_allowed_occ.size() > 1) {

        //This is the center array of possible_swap and changes with each prim_basis_site with more than one allowed occupant
        std::vector<std::vector<int> > single_possible_swap;

        //This is a site_to_mol inner array and changes with each prim_basis_site with more than one allowed occupant
        std::vector<int> single_site_to_mol;

        //For every current occupant we could have, we need an std::vector<int> for single_possible_swap and an int for single_site_to_mol
        for(int possible_curr_occ = 0; possible_curr_occ < site_allowed_occ.size(); possible_curr_occ++) {
          //This is the inner array of possible_swap and tells you all the OTHER values the current occupant could have been for the site
          std::vector<int> swaps;

          //Loop over the same thing, and only keep different values of the current occupant as swaps that could occur
          for(int other_curr_occ = 0; other_curr_occ < site_allowed_occ.size(); other_curr_occ++) {
            if(possible_curr_occ != other_curr_occ) {
              swaps.push_back(other_curr_occ);
            }
          }

          //At this point we have all the available swaps for one site of a given current occupant, so we save them
          single_possible_swap.push_back(swaps);

          //For site_to_mol we want a way to translate what possible_curr_occ means in terms of allowed components at the PrimClex level
          int mol_ind = std::find(allowed_components.cbegin(), allowed_components.cend(), site_allowed_occ[possible_curr_occ]) - allowed_components.cbegin();

          //This should never ever happen
          if(mol_ind == allowed_components.size()) {
            std::cerr << "ERROR in Monte::populate_occ_exchange_tables" << std::endl;
            std::cerr << "The possible components are " << allowed_components << std::endl;
            std::cerr << "Could not find " << allowed_components[possible_curr_occ] << " in the allowed components of the PrimClex prim." << std::endl;
            exit(9000);
          }

          //Save the index into allowed_components
          single_site_to_mol.push_back(mol_ind);
        }

        //possible_swap has blocks of identical double arrays single_possible_swap. They get pushed back the appropriate amount of times here
        m_possible_swap.push_back(single_possible_swap);

        //site_to_mol has blocks of identical single arrays single_site_to_mol. They get pushed back the appropriate amount of times here
        m_sublat_to_mol.push_back(single_site_to_mol);

        // This loop should happen scel_volume times and takes care of repeating values in m_variable_sites and m_sublat
        // the appropriate amount of times. Both m_variable_sites and m_sublat have the same length on the outside: one
        // slot for each site in the Configuration that can hold more than one occupant.
        for(int variable_site = prim_basis_site * scel_volume; variable_site < prim_basis_site * scel_volume + scel_volume; variable_site++) {

          // variable_sites determined by the counter. Contains indexes of sites in the Configuration that allow
          // more than one occupant
          m_variable_sites.push_back(variable_site);

          // store the sublat index
          m_sublat.push_back(prim_basis_site);

        }
      }

      //Here we populate m_possible_swap and m_sublat_to_mol for sublattices that only have a single allowed occupant.
      else {

        //In this else block, single_site_to_mol can only have a single value in the array
        std::vector<int> single_site_to_mol;

        //The one and only allowed occupant at prim_basis_site must be
        std::string only_site_occ = scel.get_prim().basis[prim_basis_site].allowed_occupants()[0];

        //And the corresponding index for that occupant in terms of allowed_components is
        int mol_ind = std::find(allowed_components.cbegin(), allowed_components.cend(), only_site_occ) - allowed_components.cbegin();
        if(mol_ind == allowed_components.size()) {
          std::cerr << "ERROR in Monte::populate_occ_exchange_tables" << std::endl;
          std::cerr << "The possible components are " << allowed_components << std::endl;
          std::cerr << "Could not find " << only_site_occ << " in the allowed components of the PrimClex prim." << std::endl;
          exit(9000);
        }

        // No swaps are possible on this sublattice, so we push back empty std::vector<std::vector<int> >
        std::vector<std::vector<int> > single_possible_swap;
        m_possible_swap.push_back(single_possible_swap);

        //Put that single value we found into the array
        single_site_to_mol.push_back(mol_ind);

        m_sublat_to_mol.push_back(single_site_to_mol);

      }
    }

    return;
  }

}

