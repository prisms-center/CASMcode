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
	    const Configuration to_config() const{
	    	Configuration tmp {m_from_config};
	    	return m_diff_trans.apply_to(tmp);
	    }

	    /// \brief Returns the diffusion transformation that is occurring
	    const DiffusionTransformation &diff_trans() const{
	    	return m_diff_trans;
	    }

	    /// \brief Compare DiffTransConfiguration
	    /// Compares m_diff_trans first then
	    /// m_from_config if m_diff_trans are equal
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

	    /// \brief Returns a DiffTransConfiguration that is the canonical form of this
	    DiffTransConfiguration canonical_form() const; 

	    /// \brief Returns a PermuteIterator that takes this to canonical form
	    PermuteIterator to_canonical() const;

	    /// \brief Returns a PermuteIterator that takes the canonical form of this to this
	    PermuteIterator from_canonical() const{
	    	return to_canonical().inverse();
	    }

	    /// \brief States if this DiffTransConfiguration is in canonical form
	    bool is_canonical() const;

	    /// \brief applies the symmetry op corresponding to the PermuteIterator to the 
	    /// DiffTransConfiguration in place
	    DiffTransConfiguration &apply_sym(const PermuteIterator &it);

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
	      return this->from_config() < B.from_config();
	    }

	    Configuration m_from_config;

	    // not necessary to store, could be determined by applying m_diff_trans
	    //Configuration m_to_config;

	    DiffusionTransformation m_diff_trans;
	  };

	  DiffTransConfiguration &apply_sym(DiffTransConfiguration &diff_trans_config, const PermuteIterator &it);

		DiffTransConfiguration copy_apply_sym(const DiffTransConfiguration &diff_trans_config, const PermuteIterator &it);
	}
}

#endif
