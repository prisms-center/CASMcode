#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/DataFormatterTools.hh"
#include "casm/symmetry/Orbit.hh"

namespace CASM {


  namespace PrimPeriodicDiffTransOrbitIO {

    template<typename ValueType>
    using GenericDiffTransOrbitFormatter = GenericDatumFormatter<ValueType, PrimPeriodicDiffTransOrbit>;

    class Contains : public BooleanAttribute<PrimPeriodicDiffTransOrbit> {

    public:

      static const std::string Name;

      static const std::string Desc;

      Contains() : BooleanAttribute<PrimPeriodicDiffTransOrbit>(Name, Desc) { };

      Contains(std::vector<std::string> search_list):
        BooleanAttribute<PrimPeriodicDiffTransOrbit>(Name, Desc),
        m_search_list(search_list) {};

      // --- Required implementations -----------

      /// \brief Returns true if all of the species in m_search_list are in the prototype
      ///  of the orbit
      bool evaluate(const PrimPeriodicDiffTransOrbit &orbit) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Contains> clone() const {
        return std::unique_ptr<Contains>(this->_clone());
      }

      /// \brief Expects 'contains("Specie1","Specie2",...)'
      bool parse_args(const std::string &args) override;

      /// \brief Short header returns: 'is_diff_trans_endpoint_of(diff_trans_name)', etc.
      std::string short_header(const PrimPeriodicDiffTransOrbit &_tmplt) const override {
        std::string big_string = m_search_list[0];
        auto it = m_search_list.begin();
        it++;
        for(; it != m_search_list.end(); ++it) {
          big_string += ("," + *it);
        }
        return "contains(" + big_string + ")";
      }
    private:
      /// \brief Clone using copy constructor
      Contains *_clone() const override {
        return new Contains(*this);
      }
      mutable std::vector<std::string> m_search_list;

    };

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<Index> multiplicity();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<Index> cluster_size();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> min_length();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> max_length();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<std::string> species_list();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<std::string> diff_trans_orbitname();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> min_dist_to_path();

    PrimPeriodicDiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> activation_barrier();
  }

  template<>
  StringAttributeDictionary<PrimPeriodicDiffTransOrbit> make_string_dictionary<PrimPeriodicDiffTransOrbit>();

  template<>
  BooleanAttributeDictionary<PrimPeriodicDiffTransOrbit> make_boolean_dictionary<PrimPeriodicDiffTransOrbit>();

  template<>
  IntegerAttributeDictionary<PrimPeriodicDiffTransOrbit> make_integer_dictionary<PrimPeriodicDiffTransOrbit>();

  template<>
  ScalarAttributeDictionary<PrimPeriodicDiffTransOrbit> make_scalar_dictionary<PrimPeriodicDiffTransOrbit>();

  template<>
  VectorXdAttributeDictionary<PrimPeriodicDiffTransOrbit> make_vectorxd_dictionary<PrimPeriodicDiffTransOrbit>();

  template<>
  VectorXiAttributeDictionary<PrimPeriodicDiffTransOrbit> make_vectorxi_dictionary<PrimPeriodicDiffTransOrbit>();


}
