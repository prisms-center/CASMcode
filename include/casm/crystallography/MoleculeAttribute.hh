#ifndef MOLECULEATTRIBUTE_HH
#define MOLECULEATTRIBUTE_HH

#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/misc/ParsingDictionary.hh"

namespace CASM {
  class jsonParser;
  class MasterSymGroup;
  class MoleculeAttribute;
  class SymOp;

  namespace MoleculeAttribute_impl {
    class BasicTraits {
    public:

      static std::string class_desc() {
        return "Molecule Attribute";
      }

      /// \brief allow destruction through base pointer
      virtual ~BasicTraits() {}

      /// \brief Name of attribute. Used globally to uniquely identify the attribute
      virtual std::string name() const = 0;

      /// \brief Populate @param _in from JSON
      virtual void from_json(MoleculeAttribute &_in, jsonParser const &_json) const = 0;

      /// \brief Output @param _out to JSON
      virtual jsonParser &to_json(MoleculeAttribute const &_out, jsonParser &_json) const = 0;

      /// \brief Generate a symmetry representation for this attribute

      /// if attribute is (M x 1) vector, the symrep will be (M x M)
      virtual SymGroupRepID generate_symrep(MasterSymGroup const &_group) const = 0;

      /// \brief non-virtual method to obtain copy through BasicTraits pointer
      std::unique_ptr<BasicTraits> clone() const {
        return std::unique_ptr<BasicTraits>(_clone());
      }
    private:
      virtual BasicTraits *_clone() const = 0;
    };


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /// \brief  Parsing dictionary for obtaining the correct MoleculeAttribute given a name
    using TraitsDictionary = ParsingDictionary<BasicTraits>;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class MoleculeAttribute {
  public:
    using BasicTraits = MoleculeAttribute_impl::BasicTraits;
    using KeyType = std::string;

    BasicTraits const &traits(KeyType const &key);

    MoleculeAttribute(std::string const &_name) :
      m_name(_name) {
      //_load_traits();
    }

    MoleculeAttribute(std::string const &_name,
                      Eigen::Ref<const Eigen::VectorXd> const &_value,
                      SymGroupRepID _rep_ID) :
      m_name(_name),
      m_value(_value),
      m_rep_ID(_rep_ID) {
      //_load_traits();
    }

    std::string const &name() const {
      return m_name;
    }

    Eigen::VectorXd const &value() const {
      return m_value;
    }

    bool identical(MoleculeAttribute const &other, double _tol) const;

    MoleculeAttribute &apply_sym(SymOp const &op);

    jsonParser &to_json(jsonParser &json) const {
      return _traits().to_json(*this, json);
    }

    void from_json(const jsonParser &json) {
      _traits().from_json(*this, json);
      return;
    }

  private:
    BasicTraits &_traits() const {
      return *m_traits_ptr;
    }

    void _generate_symrep(MasterSymGroup const &_group);

    std::string m_name;
    Eigen::VectorXd m_value;
    SymGroupRepID m_rep_ID;
    mutable notstd::cloneable_ptr<BasicTraits> m_traits_ptr;
  };

  template<>
  ParsingDictionary<MoleculeAttribute::BasicTraits>  make_parsing_dictionary<MoleculeAttribute::BasicTraits>();

  inline
  jsonParser &to_json(MoleculeAttribute const &_attr, jsonParser &json) {
    return _attr.to_json(json);
  }
}
#endif
