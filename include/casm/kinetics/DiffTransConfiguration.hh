#ifndef CASM_DiffTransConfiguration
#define CASM_DiffTransConfiguration

#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/Orbit_impl.hh"

namespace CASM {
  
  namespace Kinetics {
	  class DiffTransConfiguration {

	  public:

	  	/// \brief Constructor
	    DiffTransConfiguration(const Configuration &_from_config, const DiffusionTransformation &_diff_trans);

	    /// \brief Returns the initial configuration
	    const Configuration &from_config() const{
	    	return m_from_config;
	    }

	    /// \brief Returns the final configuration
	    const Configuration &to_config() const{
	    	return m_diff_trans.apply_to(m_from_config);
	    }

	    /// \brief Returns the diffusion transformation that is occurring
	    const DiffusionTransformation &diff_trans() const{
	    	return m_diff_trans;
	    }

	    /// \brief Compare DiffTransConfiguration
	    ///
	    /// - Comparison is made using the sorted forms
	    bool operator<(const DiffTransConfiguration &B) const {
	      return this->sorted()._lt(B.sorted());
	    }

	    /// ToDo:
	    ///   generate 'to_config' from 'from_config' and 'diff_trans'
	    ///   if 'to_config' is lt 'from_config', from_config = to_config

	    /// \brief sort DiffTransConfiguration in place
	    DiffTransConfiguration &sort();

	    /// \brief Returns a sorted version of this DiffTransConfiguration 
	    DiffTransConfiguration sorted() const;

	    /// \brief Returns true if the DiffTransConfiguration is sorted 
	    bool is_sorted() const;

	    const DiffTransConfiguration canonical_form() const;

	    bool is_canonical() const{
	    	return m_from_config.is_canonical();
	    }

	  private:

	    bool _lt(const DiffTransConfiguration &B) const {
	      // compare diff_trans
	      if(this->diff_trans() < B.diff_trans()) {
	        return true;
	      }
	      else if(B.diff_trans() < this->diff_trans()) {
	        return false;
	      }

	      // if diff_trans are equal, compare 'from_config'
	      if(B.from_config().is_equivalent(this->from_config())) {
	        // if equal, the DiffTransConfigurations must be equal, so return false
	        return false;
	      }
	      // else return (this->from_config() < B.from_config())
	      return this->from_config() < B.from_config();
	    }


	    mutable Configuration m_from_config;

	    // not necessary to store, could be determined by applying m_diff_trans
	    //Configuration m_to_config;

	    mutable DiffusionTransformation m_diff_trans;
	  };
	}
}

#endif
