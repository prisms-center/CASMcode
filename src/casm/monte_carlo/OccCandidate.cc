#include "casm/monte_carlo/OccCandidate.hh"
#include "casm/monte_carlo/Conversions.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/clex/Supercell.hh"

namespace CASM {

  namespace Monte {

    jsonParser &to_json(const OccCandidate &cand, const Conversions &convert, jsonParser &json) {
      json.put_obj();
      json["asym"] = cand.asym;
      json["spec"] = convert.species_name(cand.species_index);
      return json;
    }
  }

  Monte::OccCandidate jsonConstructor<Monte::OccCandidate>::from_json(const jsonParser &json, const Monte::Conversions &convert) {
    return Monte::OccCandidate(
             json["asym"].get<Index>(),
             convert.species_index(json["spec"].get<std::string>()));
  }

  namespace Monte {

    std::ostream &operator<<(std::ostream &sout, std::pair<const OccCandidate &, const Conversions &> value) {
      sout << "(" << value.second.species_name(value.first.species_index) << ", "
           << value.first.asym << ")";
      return sout;
    }

    jsonParser &to_json(const OccSwap &swap, const Conversions &convert, jsonParser &json) {
      jsonParser tmp;
      json.put_array();
      json.push_back(to_json(swap.cand_a, convert, tmp));
      json.push_back(to_json(swap.cand_b, convert, tmp));
      return json;
    }
  }

  Monte::OccSwap jsonConstructor<Monte::OccSwap>::from_json(const jsonParser &json, const Monte::Conversions &convert) {
    return Monte::OccSwap(
             jsonConstructor<Monte::OccCandidate>::from_json(json[0], convert),
             jsonConstructor<Monte::OccCandidate>::from_json(json[1], convert));
  }

  namespace Monte {

    std::ostream &operator<<(std::ostream &sout, std::pair<const OccSwap &, const Conversions &> value) {
      sout << std::pair<const OccCandidate &, const Conversions &>(value.first.cand_a, value.second)
           << " <-> " << std::pair<const OccCandidate &, const Conversions &>(value.first.cand_b, value.second);
      return sout;
    }


    OccCandidateList::OccCandidateList(const Conversions &convert) {

      // create set of 'candidate' asym / species pairs
      m_candidate.clear();
      for(Index asym = 0; asym < convert.asym_size(); ++asym) {

        // hard code allowed sublattices: >1 allowed occupant
        if(convert.occ_size(asym) < 2) {
          continue;
        }

        // add candidates
        for(Index i = 0; i < convert.occ_size(asym); ++i) {
          m_candidate.push_back(OccCandidate(asym, convert.species_index(asym, i)));
        }
      }

      // create lookup table of asym, species_index -> candidate index,
      //   will return {Nasym, Nspecies} if {asym, species_index} not allowed
      Index Nspecies = convert.species_size();
      Index Nasym = convert.asym_size();
      m_end = m_candidate.size();
      std::vector<Index> unallowed(Nspecies, m_end);
      m_species_to_cand_index = std::vector<std::vector<Index> >(Nasym, unallowed);

      Index index = 0;
      for(const auto &cand : m_candidate) {
        m_species_to_cand_index[cand.asym][cand.species_index] = index;
        ++index;
      }

      // make canonical and grand canonical swaps
      _make_possible_swaps(convert);
    }

    /// \brief Construct m_canonical_swaps, m_grand_canonical_swaps
    ///
    /// - Currently settings is not used, but we could add restrictions
    void OccCandidateList::_make_possible_swaps(const Conversions &convert) {

      // construct canonical and grand canonical swaps
      m_canonical_swap.clear();
      m_grand_canonical_swap.clear();

      // check that species are different and allowed on both sites
      auto allowed_canonical_swap = [&](OccCandidate cand_a, OccCandidate cand_b) {
        return cand_a.species_index != cand_b.species_index &&
               convert.species_allowed(cand_a.asym, cand_b.species_index) &&
               convert.species_allowed(cand_b.asym, cand_a.species_index);
      };

      // check that asym is the same and species_index is different
      auto allowed_grand_canonical_swap = [&](OccCandidate cand_a, OccCandidate cand_b) {
        return cand_a.asym == cand_b.asym &&
               cand_a.species_index != cand_b.species_index;
      };

      // for each pair of candidates, check if they are allowed to swap
      for(const auto &cand_a : m_candidate) {
        for(const auto &cand_b : m_candidate) {

          // don't repeat a->b, b->a
          // and check that cand_b's species is allowed on cand_a's sublat && vice versa
          if(cand_a < cand_b && allowed_canonical_swap(cand_a, cand_b)) {
            m_canonical_swap.push_back(OccSwap(cand_a, cand_b));
          }

          // allow a->b, b->a
          // check that asym is the same and species_index is different
          if(allowed_grand_canonical_swap(cand_a, cand_b)) {
            m_grand_canonical_swap.push_back(OccSwap(cand_a, cand_b));
          }
        }
      }
    }

    jsonParser &to_json(const OccCandidateList &list, const Conversions &convert, jsonParser &json) {
      jsonParser tmp;

      json.put_obj();

      json["candidate"].put_array();
      for(auto it = list.begin(); it != list.end(); ++it) {
        json["candidate"].push_back(to_json(*it, convert, tmp));
      }

      json["canonical_swap"].put_array();
      for(auto it = list.canonical_swap().begin(); it != list.canonical_swap().end(); ++it) {
        json["candidate_swap"].push_back(to_json(*it, convert, tmp));
      }

      json["grand_canonical_swap"].put_array();
      for(auto it = list.canonical_swap().begin(); it != list.canonical_swap().end(); ++it) {
        json["grand_candidate_swap"].push_back(to_json(*it, convert, tmp));
      }

      return json;
    }

    std::ostream &operator<<(std::ostream &sout, std::pair<const OccCandidateList &, const Conversions &> value) {

      typedef std::pair<const OccCandidate &, const Conversions &> cand_pair;
      typedef std::pair<const OccSwap &, const Conversions &> swap_pair;
      const Conversions &convert = value.second;
      const OccCandidateList &list = value.first;

      sout << "Unit cell for determining equivalent swaps: \n" <<
           convert.unit_scel().get_transf_mat() << "\n\n";

      sout << "Asymmetric Unit: " << std::endl;
      for(Index asym = 0; asym != convert.asym_size(); ++asym) {
        sout << "  " << asym << ": ";
        for(Index i = 0; i != convert.occ_size(asym); ++i) {
          sout << convert.species_name(convert.species_index(asym, i)) << " ";
        }
        sout << "\n";

        const auto &set = convert.asym_to_unitl(asym);
        for(auto it = set.begin(); it != set.end(); ++it) {
          sout << "    " << convert.unitl_to_bijk(*it) << "\n";
        }
      }
      sout << "\n";

      sout << "Candidates: (Species, AsymUnit)" << std::endl;
      for(auto it = list.begin(); it != list.end(); ++it) {
        sout << "  " << cand_pair(*it, convert) << "\n";
      }
      sout << "\n";

      sout << "Canonical swaps: " << std::endl;
      for(auto it = list.canonical_swap().begin(); it != list.canonical_swap().end(); ++it) {
        sout << "  " << swap_pair(*it, convert) << "\n";
      }
      sout << "\n";

      sout << "Grand canonical swaps: " << std::endl;
      for(auto it = list.grand_canonical_swap().begin(); it != list.grand_canonical_swap().end(); ++it) {
        sout << "  " << cand_pair(it->cand_a, convert) << " -> " << cand_pair(it->cand_b, convert) << "\n";
      }
      sout << "\n";
      return sout;
    }

  }

}
