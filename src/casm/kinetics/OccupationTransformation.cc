#include "casm/kinetics/OccupationTransformation.hh"
#include "casm/clex/Configuration.hh"

namespace CASM {

  namespace Kinetics {

    // OccupationTransformation

    OccupationTransformation::OccupationTransformation(
      const UnitCellCoord &uccoord,
      Index from_value,
      Index to_value) :
      SiteDoFTransformation(uccoord),
      from_value(from_value),
      to_value(to_value) {}

    OccupationTransformation &OccupationTransformation::apply_sym(const SymOp &op) {
      this->apply_sym_impl(op);
      return *this;
    }

    std::unique_ptr<OccupationTransformation> OccupationTransformation::clone() const {
      return std::unique_ptr<OccupationTransformation>(this->_clone());
    }

    Configuration &OccupationTransformation::apply_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), to_value);
      return config;
    }

    Configuration &OccupationTransformation::apply_reverse_to_impl(Configuration &config) const {
      config.set_occ(config.linear_index(uccoord), from_value);
      return config;
    }

    void OccupationTransformation::apply_sym_impl(const SymOp &op) {
      static_cast<SiteDoFTransformation &>(*this).apply_sym(op);

      //MOLECULE_SUPPORT: apply permutation to from_value & to_value
    }

    void OccupationTransformation::reverse_impl() {
      using std::swap;
      swap(from_value, to_value);
    }

    OccupationTransformation *OccupationTransformation::_clone() const {
      return new OccupationTransformation(*this);
    }
  }
}

