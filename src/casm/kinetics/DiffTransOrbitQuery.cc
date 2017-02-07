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

		/// \brief Returns if all of the species in m_search_list are in the prototype
    ///  of the orbit
    bool Contains::evaluate(const Kinetics::PrimPeriodicDiffTransOrbit &orbit) const {
    	return false;
    };

    /// \brief Expects 'contains("Specie1","Specie2",...)'
    bool Contains::parse_args(const std::string &args){
    	return false;
    };

	}
}