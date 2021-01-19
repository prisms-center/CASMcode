// #ifndef CASM_ConfigSubOrbits_impl
// #define CASM_ConfigSubOrbits_impl
//
// #include "casm/symmetry/ElementSymApply.hh"
// #include "casm/symmetry/ConfigSubOrbits.hh"
//
// #include "casm/symmetry/InvariantSubgroup_impl.hh"
// #include "casm/symmetry/SubOrbits_impl.hh"
// #include "casm/symmetry/ScelSubOrbits_impl.hh"
//
// namespace CASM {
//
//   template<typename OrbitType, typename ElementOutputIterator>
//   ElementOutputIterator MakeConfigSubOrbitGenerators::operator()(
//     const OrbitType &orbit,
//     ElementOutputIterator result) const {
//
//     typedef typename OrbitType::Element Element;
//
//     // make prim config suborbit generators
//     std::vector<Element> prim_config_suborbit_generators;
//     _slow(orbit, std::back_inserter(prim_config_suborbit_generators));
//
//     // make suborbit generators in potentially non-prim config
//     for(const auto &el : prim_config_suborbit_generators) {
//       result = make_suborbit_generators(
//                  el,
//                  m_prim_config.supercell(),
//                  m_prim_config_fg.begin(), m_prim_config_fg.end(),
//                  m_config_subgroup.begin(), m_config_subgroup.end(),
//                  result);
//     }
//     return result;
//   }
//
//   template<typename OrbitIterator, typename ElementOutputIterator>
//   ElementOutputIterator MakeConfigSubOrbitGenerators::operator()(
//     OrbitIterator begin,
//     OrbitIterator end,
//     ElementOutputIterator result) const {
//
//     while(begin != end) {
//       result = (*this)(*begin, result);
//       ++begin;
//     }
//     return result;
//   }
//
//   template<typename OrbitType, typename ElementOutputIterator>
//   ElementOutputIterator MakeConfigSubOrbitGenerators::_slow(
//     const OrbitType &orbit,
//     ElementOutputIterator result) const {
//
//     typedef typename OrbitType::Element Element;
//
//     // get generating elements for prim->supercell orbit splitting
//     std::vector<Element> scel_suborbit_generators;
//
//     make_suborbit_generators(
//       orbit,
//       m_prim_config.prim().factor_group(),
//       m_prim_config.supercell().factor_group(),
//       typename
//       traits<Element>::copy_apply_f_type(this->m_config.supercell().primclex().shared_prim()),
//       std::back_inserter(scel_suborbit_generators));
//
//     // get generating elements for supercell->config orbit splitting
//     for(const auto &el : scel_suborbit_generators) {
//       result = make_suborbit_generators(
//                  el,
//                  m_prim_config.supercell(),
//                  m_prim_config_fg.begin(),
//                  m_prim_config_fg.end(),
//                  result);
//     }
//
//     return result;
//   }
//
//   /// \brief Output the orbit generators necessary to construct the
//   sub-orbits
//   /// corresponding to Prim Structure -> Configuration symmetry breaking
//   template<typename OrbitType, typename ElementOutputIterator>
//   ElementOutputIterator make_suborbit_generators(
//     OrbitType orbit,
//     const Configuration &config,
//     ElementOutputIterator result) {
//
//     return MakeConfigSubOrbitGenerators{config}(orbit, result);
//   }
//
//   /// \brief Output the orbit generators necessary to construct the
//   sub-orbits
//   /// corresponding to Prim Structure -> Configuration symmetry breaking
//   template<typename OrbitIterator, typename ElementOutputIterator>
//   ElementOutputIterator make_suborbit_generators(
//     OrbitIterator begin,
//     OrbitIterator end,
//     const Configuration &config,
//     ElementOutputIterator result) {
//
//     return MakeConfigSubOrbitGenerators{config}(begin, end, result);
//   }
//
// }
//
// #endif
