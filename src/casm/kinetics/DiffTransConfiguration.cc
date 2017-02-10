#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {

	namespace Kinetics {


		DiffTransConfiguration::DiffTransConfiguration(const Configuration &_from_config,const DiffusionTransformation &_diff_trans) :
    		m_diff_trans(_diff_trans), m_from_config(_from_config){
    			
    		} 
	}
}