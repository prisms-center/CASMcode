// #ifndef CASM_ConfigSubOrbits
// #define CASM_ConfigSubOrbits
//
// #include "casm/clex/Configuration.hh"
// #include "casm/symmetry/PermuteIterator.hh"
//
// namespace CASM {
//
//   /// \brief Output the orbit generators necessary to construct the sub-orbits
//   /// corresponding to Prim Structure -> Configuration symmetry breaking
//   class MakeConfigSubOrbitGenerators {
//   public:
//
//     MakeConfigSubOrbitGenerators(const Configuration &_config);
//
//     template<typename OrbitType, typename ElementOutputIterator>
//     ElementOutputIterator operator()(const OrbitType &orbit, ElementOutputIterator result) const;
//
//     template<typename OrbitIterator, typename ElementOutputIterator>
//     ElementOutputIterator operator()(
//       OrbitIterator begin,
//       OrbitIterator end,
//       ElementOutputIterator result) const;
//
//   private:
//
//     template<typename OrbitType, typename ElementOutputIterator>
//     ElementOutputIterator _slow(const OrbitType &orbit, ElementOutputIterator result) const;
//
//     Configuration m_config;
//     Configuration m_prim_config;
//
//     // operations that map m_prim_config onto itself
//     std::vector<PermuteIterator> m_prim_config_fg;
//
//     // subgroup of m_prim_config_fg that maps m_config supercell onto itself
//     std::vector<PermuteIterator> m_config_subgroup;
//
//   };
//
//   /// \brief Output the orbit generators necessary to construct the sub-orbits
//   /// corresponding to Prim Structure -> Configuration symmetry breaking
//   template<typename OrbitType, typename ElementOutputIterator>
//   ElementOutputIterator make_suborbit_generators(
//     OrbitType orbit,
//     const Configuration &config,
//     ElementOutputIterator result);
//
//   /// \brief Output the orbit generators necessary to construct the sub-orbits
//   /// corresponding to Prim Structure -> Configuration symmetry breaking
//   template<typename OrbitIterator, typename ElementOutputIterator>
//   ElementOutputIterator make_suborbit_generators(
//     OrbitIterator begin,
//     OrbitIterator end,
//     const Configuration &config,
//     ElementOutputIterator result);
//
// }
//
// #endif
