#ifndef CASM_SpeciesProperty
#define CASM_SpeciesProperty

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/ParsingDictionary.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"

namespace CASM {
class jsonParser;
namespace xtal {
struct SymOp;
class SpeciesProperty;

namespace SpeciesProperty_impl {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/// \brief  Parsing dictionary for obtaining the correct SpeciesProperty given
/// a name
using TraitsDictionary = ParsingDictionary<AnisoValTraits>;

}  // namespace SpeciesProperty_impl

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SpeciesProperty {
 public:
  using BasicTraits = AnisoValTraits;
  using KeyType = std::string;

  BasicTraits const &traits(KeyType const &key);

  SpeciesProperty(AnisoValTraits const &_traits,
                  Eigen::Ref<const Eigen::VectorXd> const &_value)
      : m_traits(_traits), m_value(_value) {}

  SpeciesProperty(AnisoValTraits const &_traits)
      : SpeciesProperty(_traits, Eigen::VectorXd::Zero(_traits.dim())) {}

  std::string const &name() const { return traits().name(); }

  Eigen::VectorXd const &value() const { return m_value; }

  void set_value(Eigen::Ref<const Eigen::VectorXd> const &_value) {
    m_value = _value;
  }

  bool identical(SpeciesProperty const &other, double _tol) const;

  BasicTraits const &traits() const { return m_traits; }

 private:
  BasicTraits m_traits;
  Eigen::VectorXd m_value;
};

}  // namespace xtal
}  // namespace CASM

namespace CASM {
namespace sym {
xtal::SpeciesProperty &apply(const xtal::SymOp &op,
                             xtal::SpeciesProperty &mutating_attribute);
xtal::SpeciesProperty copy_apply(const xtal::SymOp &op,
                                 xtal::SpeciesProperty attribute);

template <typename ExternSymOp>
xtal::SpeciesProperty copy_apply(const ExternSymOp &op,
                                 xtal::SpeciesProperty attribute) {
  return sym::copy_apply(adapter::Adapter<xtal::SymOp, ExternSymOp>()(op),
                         attribute);
}
}  // namespace sym
}  // namespace CASM
#endif
