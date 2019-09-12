#ifndef CASM_SpeciesAttribute
#define CASM_SpeciesAttribute

#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/misc/ParsingDictionary.hh"
#include "casm/crystallography/AnisoValTraits.hh"

namespace CASM {
  class jsonParser;
  class MasterSymGroup;
  class SpeciesAttribute;
  class SymOp;

  namespace SpeciesAttribute_impl {


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief  Parsing dictionary for obtaining the correct SpeciesAttribute given a name
    using TraitsDictionary = ParsingDictionary<AnisoValTraits>;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class SpeciesAttribute {
  public:
    using BasicTraits = SpeciesAttribute_impl::BasicTraits;
    using KeyType = std::string;

    BasicTraits const &traits(KeyType const &key);

    SpeciesAttribute(std::string const &_name) :
      m_name(_name) {
      //_load_traits();
    }

    SpeciesAttribute(std::string const &_name,
                     Eigen::Ref<const Eigen::VectorXd> const &_value):
      m_name(_name),
      m_value(_value) {
      //_load_traits();
    }

    std::string const &name() const {
      return m_name;
    }

    Eigen::VectorXd const &value() const {
      return m_value;
    }

    bool identical(SpeciesAttribute const &other, double _tol) const;

    SpeciesAttribute &apply_sym(SymOp const &op);

    jsonParser &to_json(jsonParser &json) const {
      return traits().to_json(*this, json);
    }

    void from_json(const jsonParser &json) {
      traits().from_json(*this, json);
      return;
    }

    BasicTraits const &traits() const {
      return *m_traits_ptr;
    }
  private:

    std::string m_name;
    Eigen::VectorXd m_value;
    mutable notstd::cloneable_ptr<const BasicTraits> m_traits_ptr;
  };

  template<>
  ParsingDictionary<SpeciesAttribute::BasicTraits>  make_parsing_dictionary<SpeciesAttribute::BasicTraits>();

  inline
  jsonParser &to_json(SpeciesAttribute const &_attr, jsonParser &json) {
    return _attr.to_json(json);
  }
}
#endif
