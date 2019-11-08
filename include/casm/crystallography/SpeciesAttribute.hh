#ifndef CASM_SpeciesAttribute
#define CASM_SpeciesAttribute

#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/misc/ParsingDictionary.hh"
#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
  class jsonParser;
  namespace xtal {
    class SymOp;
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

      SpeciesAttribute &apply_sym(SymOp const &op);

      BasicTraits const &traits() const {
        return m_traits;
      }
    private:
      BasicTraits m_traits;
      Eigen::VectorXd m_value;
    };

  }
}
#endif
