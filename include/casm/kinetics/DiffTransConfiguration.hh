#ifndef CASM_DiffTransConfiguration
#define CASM_DiffTransConfiguration

#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {
  
  namespace Kinetics {
	  class DiffTransConfiguration : public Comparisons<DiffTransConfiguration> {

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

	    /// \brief sort DiffTransConfiguration in place
	    DiffTransConfiguration &sort();

	    /// \brief Returns a sorted version of this DiffTransConfiguration 
	    DiffTransConfiguration sorted() const;

	    /// \brief Returns true if the DiffTransConfiguration is sorted 
	    bool is_sorted() const;

	    /// \brief Converts this DiffTransConfiguration to canonical form
	    DiffTransConfiguration &canonical_form();

	    /// \brief Returns a DiffTransConfiguration that is the canonical form of this
	    DiffTransConfiguration canonical_equiv() const;

	    /// \brief States if this DiffTransConfiguration is in canonical form
	    bool is_canonical() const{
	    	return m_from_config.is_canonical();
	    }

	    /// \brief applies the symmetry op corresponding to the PermuteIterator to the 
	    /// DiffTransConfiguration
	    DiffTransConfiguration &apply_sym_impl(const PermuteIterator &it);

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

	  DiffTransConfiguration &apply_sym(DiffTransConfiguration &diff_trans_config, const PermuteIterator &it);

		DiffTransConfiguration copy_apply_sym(const DiffTransConfiguration &diff_trans_config, const PermuteIterator &it);
	}
}

#endif
