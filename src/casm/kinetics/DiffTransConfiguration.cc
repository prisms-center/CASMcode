#include "casm/kinetics/DiffusionTransformation.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
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
    	return m_from_config < to;
    }

    DiffTransConfiguration &DiffTransConfiguration::canonical_form(){
    	this->sort();

    	// check which supercell factor group operations 
    	// when applied to m_diff_trans results in the greatest 
    	// DiffusionTransformation
    	DiffusionTransformation greatest {m_diff_trans};
    	std::vector<PermuteIterator> checklist;
    	for (auto it = m_from_config.supercell().permute_begin(); it != m_from_config.supercell().permute_begin(); it++){
    		DiffusionTransformation tmp = copy_apply(it.sym_op(),m_diff_trans);
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
    		if (tmp > max_config){
    			max_config = tmp;
    			canon_op_it = *it;
    		}
    	}
    	// apply the operation to both m_diff_trans and m_from_config
    	SymOp op = canon_op_it.sym_op();
    	m_diff_trans.apply_sym(op);
    	apply(canon_op_it,m_from_config);
    	return *this;
    }

   	DiffTransConfiguration DiffTransConfiguration::canonical_equiv() const{
			DiffTransConfiguration tmp {*this};
    	return tmp.canonical_form();
    }

    DiffTransConfiguration &DiffTransConfiguration::apply_sym(const PermuteIterator &it){
    	m_from_config = apply(it,m_from_config);
    	m_diff_trans.apply_sym(it.sym_op());
    	return *this;
    }

   
    
	}
}