#ifndef CASM_Orbit_impl
#define CASM_Orbit_impl

#include <set>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/misc/algorithm.hh"
#include "casm/kinetics/PrimPeriodicDiffTransOrbitTraits.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/database/Named.hh"
#include "casm/database/Database.hh"
#include "casm/app/DirectoryStructure.hh"

namespace CASM {

  namespace Orbit_impl {

    /// \brief Returns vector containing sorted indices of SymOp in first column of equivalence_map
    ///
    /// - This is what the first column of equivalence_map will look like if Element
    /// 'proto' is the prototype.
    /// - Uses index into generating_group
    template<typename Element, typename SymCompareType, typename EquivContainer>
    std::vector<Index> _sorter(
      const Element &proto,
      const EquivContainer &equiv,
      const SymGroup &g,
      const SymCompareType &sym_compare) {

      auto equal = [&](const Element & A, const Element & B) {
        return sym_compare.equal(A, B);
      };

      int count = 0;
      std::vector<Index> sorter(equiv.size(), -1);

      for(Index op_index = 0; op_index != g.size(); ++op_index) {
        Index i = find_index(equiv, sym_compare.prepare(copy_apply(g[op_index], proto)), equal);
        if(sorter[i] == -1) {
          sorter[i] = op_index;
          count++;
          if(count == equiv.size()) {
            std::sort(sorter.begin(), sorter.end());
            return sorter;
          }
        }
      }
      throw std::runtime_error("Error generating equivalence map");
    }

    template<typename Element>
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
  GenericOrbit<_Element, _SymCompareType>::GenericOrbit(Element generating_element,
                                                        const SymGroup &generating_group,
                                                        const _SymCompareType &sym_compare) :
    m_sym_compare(sym_compare) {

    const SymGroup &g = generating_group;

    auto prepare = [&](const Element & A) {
      return m_sym_compare.prepare(A);
    };
    auto compare = [&](const Element & A, const Element & B) {
      return m_sym_compare.compare(A, B);
    };
    auto equal = [&](const Element & A, const Element & B) {
      return m_sym_compare.equal(A, B);
    };
    // generate equivalents
    std::set<Element, decltype(compare)> t_equiv(compare);
    for(const auto &op : g) {
      t_equiv.insert(prepare(copy_apply(op, generating_element)));
    }
    // sort element using each element's first equivalence map column to find prototype
    std::set<std::pair<Element, std::vector<Index> >, Orbit_impl::_EqMapCompare<Element> > _set;
    for(const auto &e : t_equiv) {
      _set.insert(std::make_pair(e, Orbit_impl::_sorter(e, t_equiv, g, m_sym_compare)));
    }

    // use _set.begin()->first for prototype, use _set.begin()->second to generate equiv
    for(auto op_index : _set.begin()->second) {
      m_element.push_back(prepare(copy_apply(g[op_index], _set.begin()->first)));
    }
    // generate equivalence map
    m_equivalence_map.resize(m_element.size());
    for(const auto &op : g) {
      Index i = find_index(m_element, prepare(copy_apply(op, m_element[0])), equal);
      m_equivalence_map[i].push_back(op);
    }

  }

  /// \brief Apply symmetry to Orbit
  template<typename _Element, typename _SymCompareType>
  GenericOrbit<_Element, _SymCompareType> &GenericOrbit<_Element, _SymCompareType>::apply_sym(const SymOp &op) {

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
    m_sym_compare.apply_sym(op);

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
  /// - Uses `SymCompareType::compare` (via `Orbit::contains`) to check for
  ///   element in orbit
  template<typename OrbitIterator, typename Element>
  OrbitIterator find_orbit(OrbitIterator begin, OrbitIterator end, Element e) {

    typedef typename std::iterator_traits<OrbitIterator>::value_type orbit_type;
    const auto &sym_compare = begin->sym_compare();

    // first find range of possible orbit by checking invariants
    auto compare = [&](const Element & A, const Element & B) {
      return sym_compare.invariants_compare(A.invariants(), B.invariants());
    };
    auto _range = std::equal_range(prototype_iterator(begin), prototype_iterator(end), e, compare);

    // find if any of the orbits in range [_range.first, _range.second) contain equivalent
    auto contains = [&](const orbit_type & orbit) {
      return orbit.contains(e);
    };
    auto res = std::find_if(_range.first.base(), _range.second.base(), contains);
    if(res == _range.second.base()) {
      return end;
    }
    return res;
  }

  /// \brief Construct an Orbit from a generating_element Element, using provided symmetry group
  template<typename _Element, typename _SymCompareType>
  DatabaseTypeOrbit<_Element, _SymCompareType>::DatabaseTypeOrbit(Element generating_element,
                                                                  const SymGroup &generating_group,
                                                                  const SymCompareType &sym_compare,
                                                                  const PrimClex *_primclex) :
    GenericOrbit<_Element, _SymCompareType>(generating_element, generating_group, sym_compare),
    m_primclex(_primclex) {}

  template<typename _Element, typename _SymCompareType>
  void DatabaseTypeOrbit<_Element, _SymCompareType>::write_pos() const {
    const auto &dir = primclex().dir();
    try {
      fs::create_directories(dir.configuration_dir(this->name()));
    }
    catch(const fs::filesystem_error &ex) {
      std::cerr << "Error in DatabaseTypeOrbit::write_pos(): could not create_directories" << std::endl;
      std::cerr << ex.what() << std::endl;
    }

    fs::ofstream file(dir.POS(this->name()));
    write_pos(file);
  }

  template<typename _Element, typename _SymCompareType>
  void DatabaseTypeOrbit<_Element, _SymCompareType>::write_pos(std::ostream &sout) const {
    OrbitTraits<_Element, _SymCompareType>::write_pos(*this, sout);
  }

  template<typename _Element, typename _SymCompareType>
  const PrimClex &DatabaseTypeOrbit<_Element, _SymCompareType>::primclex() const {
    return *m_primclex;
  }

  template<typename _Element, typename _SymCompareType>
  std::string DatabaseTypeOrbit<_Element, _SymCompareType>::generate_name_impl() const {
    return OrbitTraits<_Element, _SymCompareType>::generate_name_impl(*this);
  }

  template<typename _Element, typename _SymCompareType>
  void DatabaseTypeOrbit<_Element, _SymCompareType>::set_primclex(const PrimClex *_primclex) {
    m_primclex = _primclex;
  }

}

#endif
