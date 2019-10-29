#ifndef CASM_SpeciesAttribute
#define CASM_SpeciesAttribute

#include "casm/CASM_global_Eigen.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/misc/unique_cloneable_map.hh"
#include "casm/symmetry/SymGroupRepID.hh"
#include "casm/misc/ParsingDictionary.hh"

namespace CASM {
  class jsonParser;
  namespace xtal {
    class SymOp;
    class SpeciesAttribute;

    namespace SpeciesAttribute_impl {
      class BasicTraits {
      public:

        static std::string class_desc() {
          return "Species Attribute";
        }

        BasicTraits(std::string const &_name) : m_name(_name) {}

        /// \brief allow destruction through base pointer
        virtual ~BasicTraits() {}

        /// \brief Name of attribute. Used globally to uniquely identify the attribute
        virtual std::string name() const {
          return m_name;
        }

        virtual bool time_reversal_active() const {
          return false;
        }

        /// \brief Populate @param _in from JSON
        virtual void from_json(SpeciesAttribute &_in, jsonParser const &_json) const = 0;

        /// \brief Output @param _out to JSON
        virtual jsonParser &to_json(SpeciesAttribute const &_out, jsonParser &_json) const = 0;

        /// \brief Generate a symmetry representation for the supporting vector space
        /// If attribute is Mx1 vector, res<ulting matrix will be MxM
        virtual SpeciesAttribute copy_apply(SymOp const &_op, SpeciesAttribute const &_attr) const = 0;

        /// \brief non-virtual method to obtain copy through BasicTraits pointer
        std::unique_ptr<BasicTraits> clone() const {
          return std::unique_ptr<BasicTraits>(_clone());
        }
      private:
        virtual BasicTraits *_clone() const = 0;

        std::string m_name;
      };


      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      /// \brief  Parsing dictionary for obtaining the correct SpeciesAttribute given a name
      using TraitsDictionary = ParsingDictionary<BasicTraits>;

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

    inline
    jsonParser &to_json(SpeciesAttribute const &_attr, jsonParser &json) {
      return _attr.to_json(json);
    }
  }
  template<>
  ParsingDictionary<xtal::SpeciesAttribute::BasicTraits>  make_parsing_dictionary<xtal::SpeciesAttribute::BasicTraits>();
}
#endif
