#include "casm/app/HamiltonianModules_impl.hh"
#include "casm/misc/ParsingDictionary.hh"
#include "casm/crystallography/SymRepBuilder.hh"
namespace CASM {

  template<>
  HamiltonianModules::AnisoValDictionary make_parsing_dictionary<AnisoValTraits>() {
    HamiltonianModules::AnisoValDictionary dict;

    dict.insert(
      AnisoValTraits::disp(),
      AnisoValTraits::energy(),
      AnisoValTraits::cost(),
      AnisoValTraits::coordinate(),
      AnisoValTraits::latvec(),
      AnisoValTraits::selective_dynamics(),
      AnisoValTraits::Cmagspin(),
      AnisoValTraits::Cunitmagspin(),
      AnisoValTraits::NCmagspin(),
      AnisoValTraits::NCunitmagspin(),
      AnisoValTraits::SOmagspin(),
      AnisoValTraits::SOunitmagspin(),
      AnisoValTraits::isometry(),
      AnisoValTraits::strain("B"),
      AnisoValTraits::strain("U"),
      AnisoValTraits::strain("EA"),
      AnisoValTraits::strain("GL"),
      AnisoValTraits::strain("H"),
      AnisoValTraits::force(),
      AnisoValTraits::d_orbital_occupation(),
      AnisoValTraits::d_orbital_occupation_spin_polarized());

    return dict;
  }

  template<>
  HamiltonianModules::SymRepBuilderDictionary make_parsing_dictionary<SymRepBuilderInterface>() {
    HamiltonianModules::SymRepBuilderDictionary dict;
    dict.insert(
      SymRepBuilder::Identity(),
      SymRepBuilder::Cartesian(),
      SymRepBuilder::TimeReversal(),
      SymRepBuilder::AngularMomentum(),
      SymRepBuilder::Rank2Tensor(),
      SymRepBuilder::dOrbitalOccupation(),
      SymRepBuilder::dOrbitalOccupationSpinPolarized());
    return dict;
  }

  HamiltonianModules::HamiltonianModules(ProjectSettings const *_set) :
    m_dof_dict(make_parsing_dictionary<DoFDictionary::value_type>().clone()),
    m_aniso_val_dict(make_parsing_dictionary<AnisoValDictionary::value_type>().clone()),
    m_symrep_builder_dict(make_parsing_dictionary<SymRepBuilderDictionary::value_type>().clone()) {
    //std::cout << "Inside HamiltonianModules Constructor\n";
    // add DoF plugins
    if(_set) {
      load_dof_plugins(
        *_set,
        std::inserter(*m_dof_dict, m_dof_dict->end()),
        std::inserter(m_dof_lib, m_dof_lib.end()));

      // add attribute plugins
      load_symrep_builder_plugins(
        *_set,
        std::inserter(*m_symrep_builder_dict, m_symrep_builder_dict->end()),
        std::inserter(m_symrep_builder_lib, m_symrep_builder_lib.end()));
    }


  }

  HamiltonianModules::~HamiltonianModules() {
    m_aniso_val_dict->clear();

    // order of deletion matters
    m_dof_dict->clear();
    m_dof_lib.clear();

    // order of deletion matters
    m_symrep_builder_dict->clear();
    m_symrep_builder_lib.clear();
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

  HamiltonianModules::AnisoValDictionary &HamiltonianModules::aniso_val_dict() {
    return *m_aniso_val_dict;
  }

  HamiltonianModules::AnisoValDictionary const &HamiltonianModules::aniso_val_dict()const {
    return *m_aniso_val_dict;
  }

  HamiltonianModules::SymRepBuilderDictionary &HamiltonianModules::symrep_builder_dict() {
    return *m_symrep_builder_dict;
  }

  HamiltonianModules::SymRepBuilderDictionary const &HamiltonianModules::symrep_builder_dict()const {
    return *m_symrep_builder_dict;
  }


}
