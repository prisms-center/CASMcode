#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/kinetics/DiffTransOrbitQuery.hh"

namespace CASM {

	namespace DiffTransOrbitIO {

		// Contains implementations

		const std::string Contains::Name = "contains";

		const std::string Contains::Desc = 
		"Checks to see if any of the listed species are in the prototype"
		"of the orbit"
		"Requires an argument which is a string of the form:"
		"contains(species1,species2,...)";

		/// \brief Returns true if all of the species in m_search_list are in the prototype
    ///  of the orbit
    bool Contains::evaluate(const Kinetics::PrimPeriodicDiffTransOrbit &orbit) const {
    	auto speciemap = orbit.prototype().specie_count();
    	for (auto it=this->m_search_list.begin();it!= this->m_search_list.end();it++){
    		for (auto it2 = speciemap.begin(); it2!= speciemap.end();it2++){
	    		if (it2->first.name==*it && it2->second==0){
	    			return false;
	    		}
    		}
    	}
    	return true;
    };

    /// \brief Expects 'contains("Specie1","Specie2",...)'
    bool Contains::parse_args(const std::string &args){
    	std::vector<std::string> splt_vec;
      boost::split(splt_vec, args, boost::is_any_of(","), boost::token_compress_on);
      m_search_list = splt_vec;
    	return true;
    };

    //Generic Formatter functions

    GenericDiffTransOrbitFormatter<Index> multiplicity() {
      return GenericDiffTransOrbitFormatter<Index>("multiplicity",
                                           "Symmetric multiplicity of the diffusion transformation, excluding translational equivalents.",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)->Index {
        return orbit.size();
      });
    }

    GenericDiffTransOrbitFormatter<Index> cluster_size() {
      return GenericDiffTransOrbitFormatter<Index>("cluster_size",
                                           "Size of the cluster corresponding to the diffusion transformation.",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)->Index {
        return orbit.prototype().size();
      });
    }

    GenericDiffTransOrbitFormatter<double> min_length() {
      return GenericDiffTransOrbitFormatter<double>("min_length",
                                           "Returns the minimum length between two sites in the cluser of the diffusion transformation",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> double {
        return orbit.prototype().min_length();
      });
    }

    GenericDiffTransOrbitFormatter<double> max_length() {
      return GenericDiffTransOrbitFormatter<double>("max_length",
                                           "Returns the maximum length between two sites in the cluser of the diffusion transformation",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> double {
        return orbit.prototype().max_length();
      });
    }

    GenericDiffTransOrbitFormatter<double> min_dist_to_path() {
      return GenericDiffTransOrbitFormatter<double>("min_dist_to_path",
                                           "Returns the distance of the closest atom to the diffusion transformation",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> double {
        return min_dist_to_path(orbit.prototype());
      });
    }

    GenericDiffTransOrbitFormatter<std::string> species_list() {
      return GenericDiffTransOrbitFormatter<std::string>("species_list",
                                           "Returns the species that move in this diffusion transformation",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> std::string {

      	auto speciemap = orbit.prototype().specie_count();
      	std::string result = "";
      	for (auto it=speciemap.begin(); it!=speciemap.end();++it){
      		for (int i= 0; i<it->second;i++){
      			result += (it->first.name + " ");
      		}
      	}
        return result;
      });
    }

    GenericDiffTransOrbitFormatter<std::string> difftransname() {
      return GenericDiffTransOrbitFormatter<std::string>("difftransname",
                                           "Returns the name of the diffusion transformation",
      [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> std::string {
        return orbit.prototype().name();
      });
    }

    GenericDiffTransOrbitFormatter<double> min_activation_barrier() {
    return GenericDiffTransOrbitFormatter<double>("min_activation_barrier",
                                         "Returns the locally cluster expanded activation barrier associated with this diffusion transformation",
    [](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> double {
    	//TEMPORARY FUNCTION UNTIL ACTIVATION BARRIER CAN BE CALCULATED
      return min_dist_to_path(orbit.prototype());
    },[](const Kinetics::PrimPeriodicDiffTransOrbit & orbit)-> bool {
      return false;
    });
  	}
	}

    

	template<>
  StringAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_string_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>() {

    using namespace DiffTransOrbitIO;
    StringAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      species_list(),
      difftransname()
    );

    return dict;
  }

	template<>
  BooleanAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_boolean_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>() {

    using namespace DiffTransOrbitIO;
    BooleanAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      Contains()
    );

    return dict;
  }

	template<>
  IntegerAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_integer_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>() {

    using namespace DiffTransOrbitIO;
    IntegerAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      multiplicity(),
      cluster_size()
    );

    return dict;
  }

	template<>
  ScalarAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> make_scalar_dictionary<Kinetics::PrimPeriodicDiffTransOrbit>() {

    using namespace DiffTransOrbitIO;
    ScalarAttributeDictionary<Kinetics::PrimPeriodicDiffTransOrbit> dict;

    dict.insert(
      min_length(),
      max_length(),
      min_dist_to_path()
    );

    return dict;
  }      
}