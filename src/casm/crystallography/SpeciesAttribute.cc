#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/crystallography/SpeciesAttribute.hh"
#include "casm/symmetry/SymOp.hh"

namespace CASM {
  namespace xtal {
    namespace SpeciesAttribute_impl {

      SpeciesAttribute BasicTraits::copy_apply(SymOp const &_op, SpeciesAttribute const &_attr) const {
        return _attr;
      }

    }

    //*******************************************************************
    SpeciesAttribute &SpeciesAttribute::apply_sym(SymOp const &_op) {
      return *this = traits().copy_apply(_op, *this);
    }

    //*******************************************************************

    bool SpeciesAttribute::identical(SpeciesAttribute const &other, double _tol) const {
      return name() == other.name() && almost_equal(value(), other.value(), _tol);
    }
  }
  template<>
  ParsingDictionary<xtal::SpeciesAttribute::BasicTraits>  make_parsing_dictionary<xtal::SpeciesAttribute::BasicTraits>() {
    ParsingDictionary<xtal::SpeciesAttribute::BasicTraits> dict;
    return dict;
  }

}