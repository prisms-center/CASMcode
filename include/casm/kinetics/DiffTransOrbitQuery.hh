#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/Orbit_impl.hh"


namespace CASM {


	namespace DiffTransOrbitIO{

		template<typename ValueType>
		using GenericDiffTransOrbitFormatter = GenericDatumFormatter<ValueType,Kinetics::PrimPeriodicDiffTransOrbit>;

		class Contains : public BooleanAttribute<Kinetics::PrimPeriodicDiffTransOrbit> {
		
		public:

			static const std::string Name;

      static const std::string Desc;

      Contains() : BooleanAttribute<Kinetics::PrimPeriodicDiffTransOrbit>(Name,Desc){};

      Contains(std::vector<std::string> search_list):
      	BooleanAttribute<Kinetics::PrimPeriodicDiffTransOrbit>(Name,Desc),
      	m_search_list(search_list){};

      // --- Required implementations -----------

      /// \brief Returns true if all of the species in m_search_list are in the prototype
      ///  of the orbit
      bool evaluate(const Kinetics::PrimPeriodicDiffTransOrbit &orbit) const override;

      /// \brief Clone using copy constructor
      std::unique_ptr<Contains> clone() const {
        return std::unique_ptr<Contains>(this->_clone());
      }

      /// \brief Expects 'contains("Specie1","Specie2",...)'
      bool parse_args(const std::string &args) override;

    private:
    	/// \brief Clone using copy constructor
      Contains *_clone() const override {
        return new Contains(*this);
      }
    	mutable std::vector<std::string> m_search_list;

		};

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<Index> multiplicity();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<Index> cluster_size();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> min_length();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> max_length();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<std::string> species_list();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<std::string> difftransname();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> min_dist_to_path();

		DiffTransOrbitIO::GenericDiffTransOrbitFormatter<double> activation_barrier();
	}

	template<>
  StringAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_string_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>();

  template<>
  BooleanAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_boolean_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>();

  template<>
  IntegerAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_integer_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>();

  template<>
  ScalarAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_scalar_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>();

  template<>
  VectorXdAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_vectorxd_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>();
}