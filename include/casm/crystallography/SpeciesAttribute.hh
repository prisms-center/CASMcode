#ifndef CASM_SpeciesAttribute
#define CASM_SpeciesAttribute

#include "casm/crystallography/Adapter.hh"
#include "casm/global/eigen.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/misc/ParsingDictionary.hh"
#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
  class jsonParser;
  namespace xtal {
    struct SymOp;
    class SpeciesAttribute;

    namespace SpeciesAttribute_impl {


      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      /// \brief  Parsing dictionary for obtaining the correct SpeciesAttribute given a name
      using TraitsDictionary = ParsingDictionary<AnisoValTraits>;

    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    class SpeciesAttribute {
    public:
      using BasicTraits = AnisoValTraits;
      using KeyType = std::string;

      BasicTraits const &traits(KeyType const &key);

      SpeciesAttribute(AnisoValTraits const &_traits,
                       Eigen::Ref<const Eigen::VectorXd> const &_value) :
        m_traits(_traits),
        m_value(_value) {

      }

      SpeciesAttribute(AnisoValTraits const &_traits) :
        SpeciesAttribute(_traits, Eigen::VectorXd::Zero(_traits.dim())) {

      }

      std::string const &name() const {
        return traits().name();
      }

      Eigen::VectorXd const &value() const {
        return m_value;
      }

      void set_value(Eigen::Ref<const Eigen::VectorXd> const &_value) {
        m_value = _value;
      }

      bool identical(SpeciesAttribute const &other, double _tol) const;

      BasicTraits const &traits() const {
        return m_traits;
      }
    private:
      BasicTraits m_traits;
      Eigen::VectorXd m_value;
    };

  }
}

namespace CASM {
  namespace sym {
    xtal::SpeciesAttribute &apply(const xtal::SymOp &op, xtal::SpeciesAttribute &mutating_attribute);
    xtal::SpeciesAttribute copy_apply(const xtal::SymOp &op, xtal::SpeciesAttribute attribute);

    template <typename ExternSymOp>
    xtal::SpeciesAttribute copy_apply(const ExternSymOp &op, xtal::SpeciesAttribute attribute) {
      return sym::copy_apply(adapter::Adapter<xtal::SymOp, ExternSymOp>()(op), attribute);
    }
  }
}
#endif
