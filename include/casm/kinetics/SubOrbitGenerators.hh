/*
#ifndef CASM_SubOrbitGenerators
#define CASM_SubOrbitGenerators

#include "casm/container/InputEnumerator.hh"
#include "casm/clex/Supercell.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace Kinetics {

    template<typename Element, typename SymCompareType, typename ElementOutputIterator>
    ElementOutputIterator sub_orbit_generators_prim_config(
      const Orbit<Element, SymCompareType> &orbit,
      const Configuration &prim_config,
      ElementOutputIterator result) {
      // Find unique orbit elements in prim_config
      const Supercell scel = prim_config.supercell();
      const std::vector<CASM::PermuteIterator> prim_config_fg = prim_config.factor_group();

      for(auto i = 0; i < orbit.size(); i++) {
        for(auto k = 0; k < scel.prim_grid().size(); k++) {
          SymOp Sik = scel.prim_grid().sym_op(k) * orbit.equivalence_map()[i][0];

          auto lambda = [&](const PermuteIterator & it) {
            SymRepIndexCompare compare;
            return compare(Sik, it.sym_op() * Sik);
          };

          if(std::none_of(prim_config_fg.begin(), prim_config_fg.end(), lambda)) {
            Element Eik = copy_apply(Sik, orbit.prototype());
            //result.push_back(Eik);
            *result++ = Eik;
          }
        }
      }
    }

    template<typename Element, typename SymCompareType, typename ElementOutputIterator>
    ElementOutputIterator sub_orbit_generators(
      const Orbit<Element, SymCompareType> &orbit,
      const Configuration &config,
      ElementOutputIterator result) {
      // Find unique orbit elements in prim_config
      const Configuration prim_config = config.primitive();
      std::vector<Element> unique_prim_config_elements;
      sub_orbit_generators_prim_config <Element, SymCompareType, ElementOutputIterator> (orbit, prim_config, std::back_inserter(unique_prim_config_elements));

      // Map the unique orbit elements found above onto config by using SymOps in primcell that are not present in supercell
      const SymGroup config_scel_fg = config.supercell().factor_group();
      const SymGroup prim_config_scel_fg = prim_config.supercell().factor_group();

      std::vector<SymOp> unique_prim_op;

      for(auto prim_op : prim_config_scel_fg) {

        auto lambda = [&](const SymOp & config_op) {
          SymRepIndexCompare compare;
          return compare(prim_op, config_op * prim_op);
        };

        if(std::none_of(config_scel_fg.begin(), config_scel_fg.end(), lambda)) {
          unique_prim_op.push_back(prim_op);
        }
      }

      //std::vector<Element> unique_config_elements;
      for(auto e : unique_prim_config_elements) {
        for(auto op : unique_prim_op) {
          *result++ = copy_apply(op, e);
        }
      }
    }
  }
}
#endif
*/
