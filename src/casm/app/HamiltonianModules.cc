#include "casm/app/HamiltonianModules.hh"
#include "casm/misc/ParsingDictionary.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/crystallography/MoleculeAttribute.hh"

namespace CASM {

  HamiltonianModules::HamiltonianModules(const ProjectSettings &_set) :
    m_dof_dict(make_parsing_dictionary<DoFDictionary::value_type>().clone()),
    m_mol_attribute_dict(make_parsing_dictionary<MolAttributeDictionary::value_type>().clone()) {

    // add DoF plugins
    load_dof_plugins(
      _set,
      std::inserter(*m_dof_dict, m_dof_dict->end()),
      std::inserter(m_dof_lib, m_dof_lib.end()));

    // add attribute plugins
    load_mol_attribute_plugins(
      _set,
      std::inserter(*m_mol_attribute_dict, m_mol_attribute_dict->end()),
      std::inserter(m_mol_attribute_lib, m_mol_attribute_lib.end()));


  }

  HamiltonianModules::~HamiltonianModules() {
    // order of deletion matters
    m_dof_dict->clear();
    m_dof_lib.clear();
    // order of deletion matters
    m_mol_attribute_dict->clear();
    m_mol_attribute_lib.clear();
  }

  //std::unique_ptr<HamiltonianModules > HamiltonianModules::clone() const {
  //return std::unique_ptr<HamiltonianModules >(new HamiltonianModules(*this));
  //}

  HamiltonianModules::DoFDictionary &HamiltonianModules::dof_dict() {
    return *m_dof_dict;
  }

  HamiltonianModules::DoFDictionary const &HamiltonianModules::dof_dict()const {
    return *m_dof_dict;
  }

  HamiltonianModules::MolAttributeDictionary &HamiltonianModules::mol_attribute_dict() {
    return *m_mol_attribute_dict;
  }

  HamiltonianModules::MolAttributeDictionary const &HamiltonianModules::mol_attribute_dict()const {
    return *m_mol_attribute_dict;
  }


}
