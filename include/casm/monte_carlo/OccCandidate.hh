#ifndef CASM_OccCandidate_HH
#define CASM_OccCandidate_HH

#include <vector>
#include <tuple>
#include <utility>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {

  class jsonParser;
  template<typename T> struct jsonConstructor;

  namespace Monte {

    class Conversions;

    struct OccCandidate : public Comparisons<OccCandidate> {

      OccCandidate(Index _asym, Index _species_index) :
        asym(_asym),
        species_index(_species_index) {}

      Index asym;
      Index species_index;

      bool operator<(OccCandidate B) const {
        if(asym != B.asym) {
          return asym < B.asym;
        }
        return species_index < B.species_index;
      }
    };

    jsonParser &to_json(const OccCandidate &cand, const Conversions &convert, jsonParser &json);

  }

  template<>
  struct jsonConstructor<Monte::OccCandidate> {
    static Monte::OccCandidate from_json(const jsonParser &json, const Monte::Conversions &convert);
  };

  namespace Monte {

    std::ostream &operator<<(std::ostream &sout, std::pair<const OccCandidate &, const Conversions &> value);


    /// \brief Store swap type, mutating sites, and info for keeping OccLocation up-to-date
    class OccSwap : public Comparisons<OccSwap> {

    public:

      OccSwap(const OccCandidate &_cand_a, const OccCandidate &_cand_b) :
        cand_a(_cand_a),
        cand_b(_cand_b) {}

      OccCandidate cand_a;
      OccCandidate cand_b;

      void reverse() {
        using std::swap;
        std::swap(cand_a, cand_b);
      }

      OccSwap &sort() {
        OccSwap B(*this);
        B.reverse();

        if(B._lt(*this)) {
          *this = B;
        }
        return *this;
      }

      OccSwap sorted() const {
        OccSwap res(*this);
        res.sort();
        return res;
      }

      bool operator<(const OccSwap &B) const {
        return this->sorted()._lt(B.sorted());
      }


    private:

      bool _lt(const OccSwap &B) const {
        return this->tuple() < B.tuple();
      }

      typedef std::tuple<OccCandidate, OccCandidate> tuple_type;

      tuple_type tuple() const {
        return std::make_tuple(cand_a, cand_b);
      }

    };

    jsonParser &to_json(const OccSwap &swap, const Conversions &convert, jsonParser &json);
  }

  template<>
  struct jsonConstructor<Monte::OccSwap> {
    static Monte::OccSwap from_json(const jsonParser &json, const Monte::Conversions &convert);
  };

  namespace Monte {

    std::ostream &operator<<(std::ostream &sout, std::pair<const OccSwap &, const Conversions &> value);


    /// List of asym / species_index pairs indicating allowed variable occupation dof
    class OccCandidateList {

    public:

      typedef std::vector<OccCandidate>::const_iterator const_iterator;

      OccCandidateList() {}

      OccCandidateList(const Conversions &convert);

      /// Return index into std::vector<OccCandidate>, or _candidate.size() if not allowed
      Index index(const OccCandidate &cand) const {
        return m_species_to_cand_index[cand.asym][cand.species_index];
      }

      /// Return index into std::vector<OccCandidate>, or _candidate.size() if not allowed
      Index index(Index asym, Index species_index) const {
        return m_species_to_cand_index[asym][species_index];
      }

      const OccCandidate &operator[](Index candidate_index) const {
        return m_candidate[candidate_index];
      }

      const_iterator begin() const {
        return m_candidate.begin();
      }

      const_iterator end() const {
        return m_candidate.end();
      }

      Index size() const {
        return m_end;
      }

      const std::vector<OccSwap> &canonical_swap() const {
        return m_canonical_swap;
      }

      const std::vector<OccSwap> &grand_canonical_swap() const {
        return m_grand_canonical_swap;
      }

    private:

      /// \brief Construct m_canonical_swaps, m_grand_canonical_swaps
      void _make_possible_swaps(const Conversions &convert);

      /// m_converter[asym][species_index] -> candidate_index
      std::vector<std::vector<Index> > m_species_to_cand_index;

      std::vector<OccCandidate> m_candidate;

      /// Number of allowed candidates, what is returned if a candidate is not allowed
      Index m_end;

      /// vector of allowed canonical swaps
      std::vector<OccSwap> m_canonical_swap;

      /// vector of allowed grand canonical swaps
      std::vector<OccSwap> m_grand_canonical_swap;

    };

    jsonParser &to_json(const OccCandidateList &list, const Conversions &convert, jsonParser &json);

    std::ostream &operator<<(std::ostream &sout, std::pair<const OccCandidateList &, const Conversions &> value);
  }
}

#endif
