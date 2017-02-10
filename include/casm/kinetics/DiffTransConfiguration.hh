#ifndef CASM_DiffTransConfiguration
#define CASM_DiffTransConfiguration

#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/Orbit_impl.hh"

namespace CASM {
  
  namespace Kinetics {
	  class DiffTransConfiguration {

	  public:

	    DiffTransConfiguration(const Configuration &_from_config, const DiffusionTransformation &_diff_trans);

	    const Configuration &from_config() const;

	    const Configuration &to_config() const;

	    const DiffusionTransformation &diff_trans() const;

	    /// \brief Compare DiffTransConfiguration
	    ///
	    /// - Comparison is made using the sorted forms
	    bool operator<(const DiffTransConfiguration &B) const {
	      return this->sorted()._lt(B.sorted());
	    }

	    /// ToDo:
	    ///   generate 'to_config' from 'from_config' and 'diff_trans'
	    ///   if 'to_config' is lt 'from_config', from_config = to_config
	    DiffTransConfiguration &sort();

	    DiffTransConfiguration sorted() const;

	    bool is_sorted() const;


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

	    Configuration m_from_config;

	    // not necessary to store, could be determined by applying m_diff_trans
	    //Configuration m_to_config;

	    DiffusionTransformation m_diff_trans;
	  };
	}
}

#endif
