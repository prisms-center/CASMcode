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

    void LinearIndexConverter::dont_bring_within() {
      m_automatically_bring_within = false;
    }

    void LinearIndexConverter::do_bring_within() {
      m_automatically_bring_within = true;
    }

    UnitCellCoord LinearIndexConverter::operator[](Index ix) const {
      return m_linear_index_to_bijk[ix];
    }
  } // namespace xtal
} // namespace CASM
