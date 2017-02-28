#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {

	namespace Kinetics {


		DiffTransConfiguration::DiffTransConfiguration(const Configuration &_from_config,const DiffusionTransformation &_diff_trans) :
    		m_diff_trans(_diff_trans), m_from_config(_from_config){			
    }

    /// \brief sort DiffTransConfiguration in place
    DiffTransConfiguration &DiffTransConfiguration::sort(){
    	Configuration to = to_config();
    	if (to < m_from_config){
    		m_from_config = to;
    		m_diff_trans.reverse();
    	}
    	return *this;
    }

    /// \brief Returns a sorted version of this DiffTransConfiguration 
    DiffTransConfiguration DiffTransConfiguration::sorted() const{
      DiffTransConfiguration tmp {*this};
      return tmp.sort();
    }

    bool DiffTransConfiguration::is_sorted() const{
    	Configuration to = to_config();
    	return to < m_from_config;
    }

    DiffTransConfiguration &DiffTransConfiguration::canonical_form(){
    	SymOp op = m_from_config.to_canonical().sym_op();
    	m_diff_trans.apply_sym(op);
    	m_from_config.canonical_form();
    	return *this;
    }

   	DiffTransConfiguration &DiffTransConfiguration::canonical_form(){
			DiffTransConfiguration tmp {*this};
    	return tmp.canonical_form();
    }

	}
}