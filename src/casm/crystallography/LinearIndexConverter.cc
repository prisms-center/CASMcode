#include "casm/crystallography/LinearIndexConverter.hh"
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace CASM {
  namespace xtal {

    void LinearIndexConverter::_throw_if_bad_basis_sites_in_prim(int basis_sites_in_prim) {
      if(basis_sites_in_prim < 1) {
        throw std::runtime_error("UnitCellCoords require at least one basis site in the tiling unit, but you specified " +
                                 std::to_string(basis_sites_in_prim));
      }
      return;
    }

    void LinearIndexConverter::_throw_if_incompatible_index(Index ix) const {
      if(ix < 0 || ix >= this->total_sites()) {
        throw std::runtime_error("The specified index is out of range. There are " + std::to_string(this->total_sites()) +
                                 " availabel sites, but you specified index " + std::to_string(ix));
      }
      return;
    }

    void LinearIndexConverter::_throw_if_incompatible_bijk(const UnitCellCoord &bijk) const {
      if(bijk.sublattice() >= m_basis_sites_in_prim) {
        throw std::runtime_error("The given UnitCellCoord has a sublattice index that exceeds the expected value");
      }

      if(!m_automatically_bring_bijk_within && m_bijk_to_linear_index.find(bijk) == m_bijk_to_linear_index.end()) {
        throw std::runtime_error("The given UnitCellCoord lands outside of the superlattice, and you requested to not allow this");
      }

      return;
    }

    std::vector<UnitCellCoord> LinearIndexConverter::_make_all_ordered_bijk_values(const OrderedLatticePointGenerator &make_point,
                                                                                   int basis_sites_in_prim) {
      std::vector<UnitCellCoord> all_bijk_values;
      auto total_lattice_points = make_point.size();

      for(int b = 0; b < basis_sites_in_prim; ++b) {
        for(Index ijk_ix = 0; ijk_ix < total_lattice_points; ++ijk_ix) {
          auto ijk = make_point(ijk_ix);
          all_bijk_values.emplace_back(b, std::move(ijk));
        }
      }

      return all_bijk_values;
    }

    void LinearIndexConverter::never_bring_within() {
      m_automatically_bring_bijk_within = false;
    }

    void LinearIndexConverter::always_bring_within() {
      m_automatically_bring_bijk_within = true;
    }

    UnitCellCoord LinearIndexConverter::bring_within(const UnitCellCoord &bijk) const {
      //equivalent to requesting index of bijk, then using the index to get the
      //bijk within the superlattice
      return UnitCellCoord(bijk.sublattice(), this->m_bring_within_f(bijk.unitcell()));
    }

    UnitCellCoord LinearIndexConverter::operator[](Index ix) const {
      _throw_if_incompatible_index(ix);
      return m_linear_index_to_bijk[ix];
    }

    Index LinearIndexConverter::operator[](const UnitCellCoord &bijk) const {
      //Make sure the UntiCellCoord is allowed
      this->_throw_if_incompatible_bijk(bijk);

      //If only requesting UnitCellCoord conversions that are within, you already have the value
      if(!m_automatically_bring_bijk_within) {
        return m_bijk_to_linear_index.at(bijk);
      }

      //Otherwise you have to bring the UnitCellCoord within the superlattice before checking,
      //but maybe you already hit it
      if(m_bijk_to_linear_index_outside_of_superlattice.find(bijk) != m_bijk_to_linear_index_outside_of_superlattice.end()) {
        return m_bijk_to_linear_index_outside_of_superlattice.at(bijk);
      }

      //You've never seen this UnitCellCoord before. Bring it within the superlattice, figure out its index,
      //and cache the result for the future
      auto bijk_within = m_bring_within_f(bijk);
      auto ix_within = m_bijk_to_linear_index.at(bijk_within);

      m_bijk_to_linear_index_outside_of_superlattice[bijk] = ix_within;
      return ix_within;
    }
  } // namespace xtal
} // namespace CASM
