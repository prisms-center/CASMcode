#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/kinetics/DiffTransConfiguration.hh"

namespace CASM {

	namespace Kinetics {


		DiffTransConfiguration::DiffTransConfiguration(const Configuration &_from_config,
      const DiffusionTransformation &_diff_trans) :
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
    	return m_from_config < to;
    }

    PermuteIterator DiffTransConfiguration::to_canonical() const{
    	// check which supercell factor group operations 
    	// when applied to m_diff_trans results in the greatest 
    	// DiffusionTransformation
    	DiffusionTransformation greatest {m_diff_trans};
    	std::vector<PermuteIterator> checklist;
    	ScelPeriodicDiffTransSymCompare symcompare(m_from_config.supercell().prim_grid(),
        m_from_config.supercell().crystallography_tol());
    	for (auto it = m_from_config.supercell().permute_begin(); 
        it != m_from_config.supercell().permute_begin(); it++){
        //
        DiffusionTransformation tmp = symcompare.prepare(copy_apply(it.sym_op(),m_diff_trans));

    		if (tmp == greatest){
    			checklist.push_back(it);
    		}
    		else if (tmp > greatest){
    			checklist.clear();
    			greatest = tmp;
    			checklist.push_back(it);
    		}
    	}
    	// of these operations check which one maximizes
    	// the result of applying to m_from_config
    	Configuration max_config {m_from_config};
    	PermuteIterator canon_op_it;
    	for (auto it = checklist.begin(); it != checklist.end();it++){
    		Configuration tmp = copy_apply(*it,m_from_config);
    		DiffTransConfiguration dtc_tmp(tmp,greatest);
    		dtc_tmp.sort();
    		if (it == checklist.begin() || dtc_tmp.m_from_config > max_config){
    			max_config = tmp;
    			canon_op_it = *it;
    		}
    	}
    	// return the operation that transforms this to canonical form
    	return canon_op_it;
    }

   	DiffTransConfiguration DiffTransConfiguration::canonical_form() const{
    	return copy_apply(this->to_canonical(),*this);
    }

    DiffTransConfiguration &DiffTransConfiguration::apply_sym(const PermuteIterator &it){
    	m_from_config = apply(it,m_from_config);
    	ScelPeriodicDiffTransSymCompare symcompare(m_from_config.supercell().prim_grid(),
        m_from_config.supercell().crystallography_tol());
    	m_diff_trans.apply_sym(it.sym_op());
    	m_diff_trans = symcompare.prepare(m_diff_trans);

    	return *this;
    }

   
    
	}
}