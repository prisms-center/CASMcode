#ifndef CASM_Orbit_impl
#define CASM_Orbit_impl

#include "casm/symmetry/Orbit.hh"

namespace CASM {

  namespace {

    /// \brief Returns vector containing sorted indices of SymOp in first row of equivalence_map
    ///
    /// This is what the first row of equivalence_map will look like if Element
    /// 'proto' is the prototype.
    std::vector<Index> _sorter(
      const Element &proto,
      const std::vector<Element> &equiv,
      const SymGroup &g) {
      int count = 0;
      std::vector<Index> sorter(-1, equiv.size());
      for(const auto &op : g) {
        Index i = find_index(equiv, prepare(copy_apply(op, proto)), equal));
        if(sorter[i] == -1) {
          sorter[i] = op.index();
          count++;
          if(count == equiv.size()) {
            std::sort(sorter.begin(), sorter.end());
            return sorter;
          }
        }
      }
      throw std::runtime_error("Error generating equivalence map");
    }

    struct _EqMapCompare {
      bool operator()(const std::pair<Element, std::vector<Index> > &A,
                      const std::pair<Element, std::vector<Index> > &B) {
        return lexicographical_compare(A.second.begin(), A.second.end(), B.second.begin(), B.second.end());
      }
    };

  }

  /* -- Orbit Definitions ------------------------------------- */

  /// \brief Construct an Orbit from a generating_element Element, using provided Group
  ///
  /// \param generating_element The element used to generate equivalents
  /// \param generating_group The group used for generating the orbit
  /// \param sym_compare Binary functor that implements symmetry properties
  ///
  template<typename _Element, typename _SymCompareType>
  Orbit<_Element, _SymCompareType>::Orbit(Element generating_element,
                                          const SymGroup &generating_group,
                                          const SymCompare<Element> &sym_compare) :
    m_sym_compare(sym_compare) {

    const SymGroup &g = generating_group;

    auto prepare = [&](const Element & A) {
      return m_sym_compare->prepare(A);
    };
    auto compare = [&](const Element & A, const Element & B) {
      return m_sym_compare->intra_orbit_compare(A, B);
    };
    auto equal = [&](const Element & A, const Element & B) {
      return m_sym_compare->intra_orbit_equal(A, B);
    };

    // generate equivalents
    std::set<Element, decltype(compare)> t_equiv(compare);
    for(const auto &op : g) {
      t_equiv.insert(prepare(copy_apply(op, generating_element)));
    }

    // sort element using each element's first equivalence map row to find prototype
    std::set<std::pair<Element, std::vector<Index> >, _EqMapCompare> _set;
    for(const auto &e : t_equiv) {
      _set.insert(std::make_pair(e, _sorter(e, t_equiv, g)));
    }

    // use _set.begin()->first for prototype, use _set.begin()->second to generate equiv
    for(auto op_index : _set.begin()->second) {
      m_element.push_back(prepare(copy_apply(g[op_index], _set.begin()->first)));
    }

    // generate equivalence map
    m_equivalence_map.resize(m_element.size());
    for(const auto &op : g) {
      Index i = find_index(m_element, prepare(copy_apply(op, m_element[0])), equal));
      m_equivalence_map[i].push_back(op);
    }

  }

  /// \brief Apply symmetry to Orbit
  template<typename _Element, typename _SymCompareType>
  Orbit &Orbit<_Element, _SymCompareType>::apply_sym(const SymOp &op) {

    // transform elements
    for(auto it = m_element.begin(); it != m_element.end(); ++it) {
      it->apply_sym(op);
    }

    // transform equivalence map: std::vector<std::vector<SymOp> >
    for(auto it = m_equivalence_map.begin(); it != m_equivalence_map.end(); ++it) {
      for(auto op_it = it->begin(); op_it != it->end(); ++op_it) {
        op_it->apply_sym(op);
      }
    }

    // transform sym_compare functor
    m_sym_compare->apply_sym(op);

    return *this;
  }


  /// \brief Find orbit containing an element in a range of Orbit<ClusterType>
  ///
  /// \param begin,end Range of Orbit
  /// \param e Element to find
  ///
  /// \returns Iterator to Orbit containing e, or end if not found
  ///
  /// - Expects `typename std::iterator_traits<OrbitIterator>::value_type` to
  ///   be the Orbit<Element, SymCompareType> type
  /// - Expects `e` to be `prepared` via `SymCompareType::prepare`
  /// - Assume range of orbit is sorted according to `SymCompareType::invariants_compare`
  /// - Uses `SymCompareType::intra_orbit_compare` (via `Orbit::contains`) to check for
  ///   element in orbit
  template<typename OrbitIterator, typename Element>
  OrbitIterator find_orbit(OrbitIterator begin, OrbitIterator end, Element e) {

    typedef typename std::iterator_traits<OrbitIterator>::value_type orbit_type;
    const auto &sym_compare = begin->sym_compare();

    // first find range of possible orbit by checking invariants
    auto compare = [&](const Element & A, const Element & B) {
      return sym_compare.invariants_compare(A, B);
    };
    auto _range = std::equal_range(begin, end, e, compare);

    // find if any of the orbits in range [_range.first, _range.second) contain equivalent
    auto contains = [&](const orbit_type & orbit) {
      return orbit.contains(e);
    };
    auto res = std::find_if(_range.first, _range.second, contains);
    if(res == _range.second) {
      return end;
    }
    return res;
  }


}

#endif
